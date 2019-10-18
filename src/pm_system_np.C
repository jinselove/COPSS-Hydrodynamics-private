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

namespace libMesh {
// ======================================================================================
PMSystemNP::PMSystemNP(EquationSystems  & es,
                       const std::string& name,
                       const unsigned int number)
  : PMLinearImplicitSystem(es, name, number),
  _solver_np(es)
{
  if (name != "NP")
  {
    PMToolBox::output_message("Error: system name in PMSystemNP initializer "
                              "is not 'NP'. Exiting...", this->comm());
    libmesh_error();
  }
  // Stokes equation assembly
  _assemble_np = (new AssembleNP(es, "NP"));
  analytical_solution = _assemble_np->get_analytical_solution();
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

// ==================================================================================
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

// ==================================================================================
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
    this->assemble_matrix("NP", option);
    // init the KSP solver
    _solver_np.init_ksp_solver();
    //set _re_init to false once K matrix is built
    _re_init = false;
  }

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
