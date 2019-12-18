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
  init_cd_set = false;
  relaxed = false;
  o_precision = 6;
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
  Real dt = 0.001;
  PMToolBox::output_message("Warning:: need to implement get_dt for NP "
                            "system, currently dt is set to be 0.1 constant",
                            this->comm());

  STOP_LOG("get_dt()", "PMSystemNP");
  return dt;
}

// ===============================================================
void PMSystemNP::init_params()
{
  // initialize some parameters of this NP system
  Parameters& params = this->get_equation_systems().parameters;
  ion_diffusivity = (params.get<std::vector<Real>>("ion_diffusivity"))[ion_id];
  ion_valence = (params.get<std::vector<int>>("ion_valence"))[ion_id];
  np_coeff = params.get<Real>("bjerrum_length") * ion_diffusivity *
    ion_valence;
  boundary_id_dirichlet_np = params.get<std::vector<unsigned int>>
    ("boundary_id_dirichlet_np");
  boundary_value_dirichlet_np = (params.get<std::vector<std::vector<Real>>>
    ("boundary_value_dirichlet_np"))[ion_id];
  // get dt from equation system instead from get_dt() function of this
  // system because we need to collect the minimum dt of all systems
  dt = params.get<Real>("dt");
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
    for (int dim_i=0; dim_i<3; dim_i++)
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
void PMSystemNP::init_cd()
{
  START_LOG("init_cd()", "PMSystemNP");

  // initialize system parameters
  this->init_params();

  // set initial ion concentration
  this->project_solution(init_solution, libmesh_nullptr,
    this->get_equation_systems().parameters);

  // sync solution to current_local_solution in case we need it
  this->update();

  init_cd_set = true;

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
//    PMToolBox::output_message("start assemble matrix", this->comm());
    this->assemble_matrix(this->name(), option);

    // init the KSP solver
//    PMToolBox::output_message("start init ksp solver", this->comm());
    _solver_np.init_ksp_solver(this->name());

    //set _re_init to false once K matrix is built
    _re_init = false;
  }
//  PMToolBox::output_message("start assemble rhs", this->comm());
  this->assemble_rhs(this->name(), option);

  // solve the problem
//  PMToolBox::output_message("start solve Ax=b", this->comm());
  _solver_np.solve(this->name());

  std::ostringstream oss;
  oss<<"* Solved '"<<this->name()<<"' with option '"<<option
     <<"'. Max solution = "<<this->solution->max();
  PMToolBox::output_message(oss, this->comm());

  STOP_LOG("solve()", "PMSystemNP");
}

// =========================================================================
void PMSystemNP::write_out_solution()
{
  // get system dimension
  MeshBase& mesh = this->get_mesh();
  const std::size_t dim = mesh.mesh_dimension();
  // get current time
  const Real& real_time = this->get_equation_systems().parameters.get<Real>
    ("real_time");
  // Output electrical potential profiles along xyz directions, global + local
  // solutions
  const Point& box_min = _point_mesh->pm_periodic_boundary()->box_min();
  const Point& box_len = _point_mesh->pm_periodic_boundary()->box_length();
  // nn: number of spacial points to check
  const unsigned int ns = 200;
  // define output file format parameters
  std::ofstream outfile;
  o_precision = this->get_equation_systems().parameters.get<int>
    ("o_precision");
  // define output filename for x, y, z directions
  std::vector<std::string> filenames;
  filenames.push_back("x");
  filenames.push_back("y");
  if (dim==3) filenames.push_back("z");
  // loop over all [dim_i] directions and write the output
  for (int dim_i=0; dim_i<dim; dim_i++)
  {
    std::ostringstream oss;
    oss << "output_concentration_profile_" << filenames[dim_i] << ".csv";
    // check if file already exists
    std::ifstream infile(oss.str());

    // open file for appending if real_time > 0. and file exists
    if (real_time>0. and infile.good())
    {
      outfile.open(oss.str(), std::ios_base::app);
    }
    else
    {
      outfile.open(oss.str(), std::ios_base::out);
      if (this->comm().rank()==0){
        std::cout << "create output file " << filenames[dim_i]
          << " at real_time: " << real_time << std::endl;
        outfile << "real_time,"
                << filenames[dim_i] << ","
                << "c,"
                << "c_exact\n";
      }
    }

    // set output format
    outfile.setf(std::ios::fixed);
    outfile.precision(o_precision);
    // write out solution at this step
    const Real ds = box_len(dim_i) / Real(ns);
    // create helper vector
    std::vector<int> helper_vec(dim, 0);
    helper_vec[dim_i] = 1;
    // loop over all [ns] space points along this direction, the coordinate of
    // other two directions are set to be zero
    for (std::size_t i = 0; i < ns + 1; ++i)
    {
      // create the i-th spacial point
      Point pt((box_min(0) + Real(i) * ds) * helper_vec[0],
               (box_min(1) + Real(i) * ds) * helper_vec[1],
               (box_min(2) + Real(i) * ds) * helper_vec[2]);
      // global potential from FEM
      const unsigned int c_var = 0; // c_var = 0
      // get global phi from FEM solution, this is slow, but we tolerate it
      // only for test
      Real c = this->point_value(c_var, pt);
      // get exact solution for an unbounded domain from analytical solution
      const Real c_exact = analytical_solution->exact_solution_infinite_domain(
        pt, real_time, ion_diffusivity);
      // write the result to output file
      if (this->comm().rank() == 0)
        outfile << real_time << ","
                << pt(dim_i) << ","
                << c << ","
                << c_exact << "\n";
    } // end for i-loop
    outfile.close();
  } // end loop dim_i
}

// =========================================================================
void PMSystemNP::test_concentration_profile()
{
  START_LOG("test_concentration_profile()", "PMSystemNP");

  // initialize system parameters
  this->init_params();

  // we only test the NP system for n_step at a constant dt
  const unsigned int n_step = 10;
  dt = 0.01;

  PMToolBox::output_message("====>2. Test concentration profile for NP system:",
    this->comm());

  // Make sure system time is 0. now
  if (this->get_equation_systems().parameters.get<Real>("real_time")>0.)
  {
    PMToolBox::output_message("Error: current real time is not zero before "
                              "testing concentration profile. Exiting...",
                              this->comm());
    libmesh_error();
  }
  // solve NP system to get concentration at t = 0. + dt
  _re_init = true;


  PMToolBox::output_message("----> 1. init system solution as analytic "
                            "concentration at time 0 (a gaussian distribution)",
                            this->comm());
  this->project_solution(init_solution, libmesh_nullptr,
                         this->get_equation_systems().parameters);
  PMToolBox::output_message(std::string("**** ion concentration at time 0: ")
    , this->comm());
  PMToolBox::output_message(std::string("max concentration = ") +
                            std::to_string(this->solution->max()), this->comm());

  // write initial condition to output

  this->write_out_solution();

  MeshBase& mesh = this->get_mesh();
  const std::string init_cd_name = "np_validation_analytic.e";
#ifdef LIBMESH_HAVE_EXODUS_API
      ExodusII_IO(mesh).write_equation_systems(init_cd_name,
    this->get_equation_systems());
#endif

  PMToolBox::output_message("----> 2. integrate np system with diffusion only",
                            this->comm());
  PMToolBox::output_message(std::string("dt = ") + std::to_string
    (dt), this->comm());
  // initialize step id and relax time
  unsigned int step_id = 0;
  Real real_time = 0.;
  while(step_id < n_step)
  {
    // update real_time since boundary condition is condition at next step
    real_time += dt;
    this->get_equation_systems().parameters.set<Real>("real_time") = real_time;

    // solve the diffusion system
    this->solve("diffusion");

//    this->rhs->print_matlab(std::string("rhs_step_")+std::to_string(step_id)
//    +".mat");
//    this->matrix->print_matlab(std::string("matrix_step_")+std::to_string
//    (step_id)+".mat");
//    this->solution->print_matlab(std::string("sol_step_")+std::to_string
//    (step_id+1)+".mat");

    // write solution to output
    step_id += 1;
    this->write_out_solution();

    // screen print out
    PMToolBox::output_message(
      std::string("**** ion concentration at real time:")
      + std::to_string(this->get_equation_systems().parameters.get<Real>
        ("real_time")), this->comm());
    PMToolBox::output_message(std::string("max concentration = ")
                              + std::to_string(this->solution->max()), this->comm());
    // write the solution during relaxation
#ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO exo(mesh);
exo.append(true);
exo.write_timestep(init_cd_name,
                   this->get_equation_systems(),
                   step_id,
                   real_time);
#endif
  }
  // reset real_time to zero after relax
  PMToolBox::output_message("----> test done.", this->comm());

  STOP_LOG("test_concentration_profile()", "PMSystemNP");
}
} // end of namespace
