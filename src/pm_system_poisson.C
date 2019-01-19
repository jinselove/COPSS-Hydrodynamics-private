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
#include "ggem_system.h"
#include "brownian_system.h"
#include "analytical_solution.h"
#include "pm_system_poisson.h"


namespace libMesh
{

// ======================================================================================
PMSystemPoisson::PMSystemPoisson(EquationSystems& es,
                                 const std::string& name,
                                 const unsigned int number)
: PMLinearImplicitSystem (es, name, number),
  _solver_poisson(es)
{
  // Poisson equation assembly
  _assemble_poisson = ( new AssemblePoisson(es) );
}



// ==================================================================================
PMSystemPoisson::~PMSystemPoisson()
{
  // Clear data
  this->clear();
}



// ==================================================================================
void PMSystemPoisson::clear ()
{
  // delete the pointer
  if(_assemble_poisson) {
    delete _assemble_poisson;
  }
}



// ===========================================================
void PMSystemPoisson::assemble_matrix(const std::string& system_name,
                                      const std::string& option)
{
  libmesh_assert (this->matrix);
  libmesh_assert (this->matrix->initialized());

  START_LOG("assemble_matrix()", "PMSystemPoisson");

  // init the matrices: global stiffness and PC matrix (if required)
  this->matrix->zero();
  const bool user_defined_pc = this->get_equation_systems().parameters.get<bool>("user_defined_pc");
  if( user_defined_pc ) {
    std::cout << "--->test in PMSystemPoisson::assemble_matrix(): "
              << "Initialize the preconditioning matrix (all zeros) \n";
    this->get_matrix("Preconditioner").zero();
  }

  // Call assemble function to assemble matrix
  _assemble_poisson->assemble_global_K(system_name, option);

  // close the matrices
  this->matrix->close();  // close the matrix
  if( user_defined_pc ) this->get_matrix("Preconditioner").close();

  STOP_LOG("assemble_matrix()", "PMSystemPoisson");
}



// ==================================================================================
void PMSystemPoisson::assemble_rhs(const std::string& system_name,
                                  const std::string& option)
{
  libmesh_assert (this->rhs);
  libmesh_assert (this->rhs->initialized());

  START_LOG("assemble_rhs()", "PMSystemPoisson");

  // init, assemble, and close the rhs vector
  this->rhs->zero();
  _assemble_poisson->assemble_global_F(system_name, option);
  this->rhs->close();

  STOP_LOG("assemble_rhs()", "PMSystemPoisson");
}



// ==================================================================================
void PMSystemPoisson::solve(const std::string& option,
                            const bool& re_init)
{
  START_LOG("solve()", "PMSystemPoisson");

  //Real t1, t2;
  //std::string msg = "---> solve Poisson";
  //PMToolBox::output_message(msg, this->comm());

  // Assemble the global matrix and pc matrix, and record the CPU wall time.
  if(re_init)
  {
    //t1 = MPI_Wtime();

    // Set the solver type for the Poisson equation
    const SystemSolverType solver_type  = this->get_equation_systems().parameters.get<SystemSolverType> ("solver_type");
    _solver_poisson.set_solver_type(solver_type);

    // Assemble the global matrix, and init the KSP solver
    this->assemble_matrix("Poisson",option);
    _solver_poisson.init_ksp_solver();

    //t2 = MPI_Wtime();
    //std::cout << "For Poisson equation, time used to assemble the global matrix and reinit KSP is " <<t2-t1<<" s\n\n";
  }

  // assemble the rhs vector, and record the CPU wall time.
  //t1 = MPI_Wtime();
  this->assemble_rhs ("Poisson",option);
  //t2 = MPI_Wtime();
  //std::cout << "For Poisson equation, time used to assemble the right-hand-side vector is " <<t2-t1<<" s\n";

  // solve the problem
  _solver_poisson.solve();

  STOP_LOG("solve()", "PMSystemPoisson");
}

} // end of namespace
