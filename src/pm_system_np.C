// Parallel Finite Element-General Geometry Ewald-like Method.

// Copyright (C) 2015-2016 Xujun Zhao, Jiyuan Li, Xikai Jiang

// This code is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.


// This code is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.


// You should have received a copy of the GNU General Public
// License along with this code; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


// std C++
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <time.h>
#include <math.h>

// libmesh head files
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmv_io.h"
#include "libmesh/vtk_io.h"

// local header files
#include "pm_toolbox.h"
#include "brownian_system.h"
#include "pm_system_np.h"
#include "pm_system_poisson.h"

namespace libMesh {
// ======================================================================================
PMSystemNP::PMSystemNP(EquationSystems  & es,
                       const std::string& name,
                       const unsigned int number)
  : PMLinearImplicitSystem(es, name, number),
  _solver_np(es)
{
  // NP equation assembly
  _assemble_np = (new AssembleNP(es, "NP"));
  analytical_solution = _assemble_np->get_analytical_solution();
  set_init_cd = true;
}

// ==================================================================================
PMSystemNP::~PMSystemNP()
{
  // Clear data
  this->clear();
}

// ==================================================================================
void PMSystemNP::clear()
{
  // delete the pointer
  if (_assemble_np) {
    delete _assemble_np;
  }
}

// =================================================================
void PMSystemNP::attach_ion_type(const int& _ion_id,
                                 const std::string& _ion_name)
{
  START_LOG("attach_ion_type()", "PMSystemNP");

  ion_id = _ion_id;
  ion_name = _ion_name;
  if (this->name() != (std::string("NP")+":"+ion_name))
  {
    PMToolBox::output_message("Error: system name is not the same as "
                              "NP:'ion_name'. Exiting...", this->comm());
    libmesh_error();
  }

  STOP_LOG("attach_ion_type()", "PMSystemNP");
}

// ==================================================================
Real PMSystemNP::get_dt()
{
  START_LOG("get_dt()", "PMSystemNP");
  // Nernst-Planck dt (unit=tc=6*pi*eta*Rb^3/kbT)
  // fixme: we need to calculate dt_np automatically to ensure numerical
  // stability
  Real dt = 0.1;
  PMToolBox::output_message("Warning:: need to implement get_dt for NP "
                            "system, currently dt is set to be 0.1 constant",
                            this->comm());

  STOP_LOG("get_dt()", "PMSystemNP");
  return dt;
}

// =============================================================
Number PMSystemNP::init_solution(const Point & p,
                                 const Parameters & parameters,
                                 const std::string &,
                                 const std::string &)
{
  // special initial conditions for our test system
  if (parameters.get<std::string>("simulation_name")=="np_validation_analytic")
  {
    // the position of instant source at t = 0.
    Point p0(0.);
    // num = (p - velocity*t - p0)^2
    Real num = 0.;
    for (int dim_i=0; dim_i < 3; dim_i++)
    {
      num += pow(p(dim_i) - p0(dim_i), 2.);
    }
    // den = D * (4. * t + 1)
    const Real& D = parameters.get<std::vector<Real>>("ion_diffusivity")[0];
    const Real den = D;
    // return exp(-num/den) / (4. * t + 1)^1.5
    return exp(-num / den);
  }
  // in general the initial condition is just 0. everywhere
  else return 0.;
}

// =============================================================
//Number PMSystemNP::init_solution_bc (const Point & p,
//                                     const Parameters & parameters,
//                                     const std::string &,
//                                     const std::string &)
//{
//  return 1.;
//}

// =============================================================
void PMSystemNP::init_cd(const Real& relax_t_final)
{
  START_LOG("init_cd()", "PMSystemNP");

  // get dt for NP system
  np_dt = this->get_dt();

  // project initial conditions on both bulk and boundary elements. By
  // default, the initial ion concentration is set to be 0 across the domain
  PMToolBox::output_message("==> Initialize ion concentration for NP system",
    this->comm());
  PMToolBox::output_message("----> 1. init uniform concentration",
    this->comm());
  this->project_solution(init_solution, libmesh_nullptr,
    this->get_equation_systems().parameters);

  // write initial condition to output
  MeshBase& mesh = this->get_mesh();
  const std::string init_cd_name = "init_cd.e";
#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO(mesh).write_equation_systems(init_cd_name,
    this->get_equation_systems());
#endif

  // project Dirichlet BC on initial solutions
  // fixme: the bounary_project_solution is not able to work in parallel
  //  because of memory issue, need to figure out why if we need to assign
  //  different initial conditions on boundaries versus in bulk
//  PMToolBox::output_message("----> Set initial solution for boundary "
//                            "elements.", this->comm());
//  for (unsigned int i = 0; i < boundary_id_dirichlet_np.size(); i++)
//  {
//    PMToolBox::output_message(std::string("** for boundary id") + std::to_string
//    (i), this->comm());
//    // create a boundary id set for this boundary
//    std::set <boundary_id_type> b{boundary_id_type(boundary_id_dirichlet_np[i])};
//    // we only have one variable 'c'
//    std::vector<unsigned int> variables{0};
//    // project boundary value defined in init_solution_bc () to the
//    // boundary nodes in initial solution
//    PMToolBox::output_message(std::string("--------> On Dirichlet Boundary: ")
//                              + std::to_string
//                                (boundary_id_dirichlet_np[i]),
//                              this->comm());
//    this->boundary_project_solution(b, variables, init_solution_bc,
//                                    libmesh_nullptr,
//                                    this->get_equation_systems().parameters);
//  }

  // relax NP system with Poisson System on
  PMToolBox::output_message("----> 2. relax ion concentration",
                            this->comm());
  unsigned int relax_step_id = 0;
  Real relax_time = 0.;
  while(relax_time < relax_t_final)
  {
    // solve the NP system with the coupling of diffusion & electrostatics, i
    // .e., turn the fluid off, if "module_poisson" is true
    if(this->get_equation_systems().parameters.get<bool>("module_poisson"))
    {
      // need to solve the poisson system first before solving NP system
      this->get_equation_systems().get_system<PMSystemPoisson>("Poisson")
        .solve("unused");
      // solve the coupled NP system
      this->solve("diffusion&electrostatics");
    }
    // solve the NP system alone, i.e., diffusion only if "module_poisson" if
    // false
    else
    {
      this->solve("diffusion");
    }
    relax_time += np_dt;
    relax_step_id += 1;
    PMToolBox::output_message(std::string("relax system time : ") +
      std::to_string(relax_time), this->comm());
    // write the solution during relaxation
#ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO exo(mesh);
    exo.append(true);
    exo.write_timestep(init_cd_name,
                       this->get_equation_systems(),
                       relax_step_id,
                       relax_time);
#endif
  }

  // avoid setting initial condition again
  set_init_cd = false;

  STOP_LOG("init_cd()", "PMSystemNP");
}

// ===========================================================
void PMSystemNP::assemble_matrix(const std::string& system_name,
                                 const std::string& option)
{
  libmesh_assert(this->matrix);
  libmesh_assert(this->matrix->initialized());

  START_LOG("assemble_matrix()", "PMSystemNP");

  // zero the matrix before assembling
  this->matrix->zero();
  // call the assemble function to assemble the matrix
  _assemble_np->assemble_global_K(system_name, option);
  // close the matrix after assembling
  this->matrix->close();

  STOP_LOG("assemble_matrix()", "PMSystemNP");
}

// =========================================================================
void PMSystemNP::assemble_rhs(const std::string& system_name,
                              const std::string& option)
{
  libmesh_assert(this->rhs);
  libmesh_assert(this->rhs->initialized());

  START_LOG("assemble_rhs()", "PMSystemNP");

  // zero the rhs vector
  this->rhs->zero();
  // call the assemble function to assemble the vector
  _assemble_np->assemble_global_F(system_name, option);
  // close rhs vector
  this->rhs->close();

  STOP_LOG("assemble_rhs()", "PMSystemNP");
}

// =========================================================================
void PMSystemNP::solve(const std::string& option)
{
  START_LOG("solve()", "PMSystemNP");

  // When _re_init is true, Assemble the global matrix and pc matrix
  // this should only happen at the first timestep
  if (_re_init)
  {
    // Set the solver type for the NP equation
    const SystemSolverType solver_type =
      this->get_equation_systems().parameters.get<SystemSolverType>(
        "solver_type_np");
    _solver_np.set_solver_type(solver_type);

    // Assemble the global matrix
    this->assemble_matrix(this->name(), option);

    // init the KSP solver
    _solver_np.init_ksp_solver(this->name());

    //set _re_init to false once K matrix is built
    _re_init = false;
  }

  this->assemble_rhs(this->name(), option);

  // solve the problem
  _solver_np.solve();

  STOP_LOG("solve()", "PMSystemNP");
}

// =========================================================================
void PMSystemNP::test_concentration_profile()
{
  START_LOG("test_concentration_profile()", "PMSystemPoisson");

  // fixme:

  STOP_LOG("test_concentration_profile()", "PMSystemPoisson");
}
} // end of namespace
