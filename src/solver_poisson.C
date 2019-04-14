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



#include "solver_poisson.h"

// make sure libMesh has been complied with PETSc
#ifdef LIBMESH_HAVE_PETSC 

// include PETSc KSP solver
EXTERN_C_FOR_PETSC_BEGIN
#  include <petscsys.h>
#  include <petscksp.h>
#  include <petscpc.h>
#  include <petscis.h>
#  include <petscviewer.h>
EXTERN_C_FOR_PETSC_END

// libmesh headers
#include "libmesh/libmesh_config.h"
#include "libmesh/dof_map.h"
#include "libmesh/system.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_linear_solver.h"

// local header
#include "pm_system_poisson.h"

using namespace libMesh;


// =======================================================================================
SolverPoisson::SolverPoisson(EquationSystems& es)
: Solver(es)
{
}



// =======================================================================================
SolverPoisson::SolverPoisson(EquationSystems& es,
                             const SystemSolverType solver_type)
: Solver(es, solver_type)
{
}



// =======================================================================================
SolverPoisson::~SolverPoisson()
{
  // destroy the PETSc context if it is created in this function.
  // NOTE: if init_ksp_solver() is not called, We don't destroy!
  if(_is_init)
  {
    int ierr;
    ierr = KSPDestroy(&_ksp);   CHKERRABORT(this->comm().get(), ierr);
  }
}



// ==================================================================================
void SolverPoisson::init_ksp_solver()
{
  START_LOG("init_ksp_solver()", "SolverPoisson");

  // Control parameters of KSP solver
  unsigned int max_iter = _equation_systems.parameters.get<unsigned int>("linear solver maximum iterations");
  Real linear_sol_rtol  = _equation_systems.parameters.get<Real> ("linear solver rtol");
  Real linear_sol_atol  = _equation_systems.parameters.get<Real> ("linear solver atol");
  _rtol    = static_cast<PetscReal>(linear_sol_rtol);
  _atol    = static_cast<PetscReal>(linear_sol_atol);
  _max_it  = static_cast<PetscInt>(max_iter);
 
  // Get a reference to the system Matrix
  PMSystemPoisson & system = _equation_systems.get_system<PMSystemPoisson> ("Poisson");
  PetscMatrix<Number>* matrix = cast_ptr<PetscMatrix<Number>*>( system.matrix );
  //this->petsc_view_matrix( matrix->mat() );
 
  // set up KSP in PETSc, resuse the Preconditioner
  //PetscPrintf(this->comm().get(), "--->test in SolverPoisson::init_ksp_solver(): ");
  //PetscPrintf(this->comm().get(), "Start to initialize the PETSc KSP solver...\n"); 
  int ierr;
  ierr = KSPCreate(this->comm().get(), &_ksp);               CHKERRABORT(this->comm().get(),ierr);
  ierr = KSPSetOperators(_ksp,matrix->mat(),matrix->mat() ); CHKERRABORT(this->comm().get(),ierr);
  ierr = KSPSetReusePreconditioner(_ksp,PETSC_TRUE);         CHKERRABORT(this->comm().get(),ierr);
  ierr = KSPSetTolerances (_ksp, _rtol, _atol, PETSC_DEFAULT, _max_it);
                                                             CHKERRABORT(this->comm().get(), ierr);
  ierr = KSPSetFromOptions(_ksp);                            CHKERRABORT(this->comm().get(), ierr);
  //ierr = KSPSetUp(_ksp); CHKERRABORT(this->comm().get(), ierr); // cause error
 
  // start from non-zero inital guess can reduce the total iteration steps of convergence?
  // actually NOT at all, and sometimes even worse !
  //ierr = KSPSetInitialGuessNonzero(_ksp,PETSC_TRUE); CHKERRABORT(this->comm().get(),ierr);
 
  // label the ksp solver is initialized
  _is_init = true;
 
  //PetscPrintf(this->comm().get(), "--->test in SolverPoisson::init_ksp_solver(): ");
  //PetscPrintf(this->comm().get(), "the PETSc KSP solver is initialized \n");
  STOP_LOG("init_ksp_solver()", "SolverPoisson");
}



// =======================================================================================
void SolverPoisson::solve()
{
  START_LOG("solve()", "SolverPoisson");

  int ierr;
  PetscInt      its     = 0;
  PetscReal final_resid = 0.;
  //Real t1, t2;
  //t1 = MPI_Wtime();
  //PetscPrintf(this->comm().get(), "--->test in SolverPoisson::solve(): ");
  //PetscPrintf(this->comm().get(), "Start the KSP solve... \n");
 
  // Get a reference to the Particle-Mesh linear implicit system object,
  // and the assembled matrix and the rhs vector, solution
  PMSystemPoisson & system = _equation_systems.get_system<PMSystemPoisson> ("Poisson");
  PetscVector<Number>*      rhs   = cast_ptr<PetscVector<Number>*>( system.rhs );
  NumericVector<Number>& sol_in   = *(system.solution);
  PetscVector<Number>* solution   = cast_ptr<PetscVector<Number>*>( &sol_in );
 
  // Look at the matrix for debug purpose
  //PetscMatrix<Number>* matrix     = cast_ptr<PetscMatrix<Number>*>( system.matrix );
  //this->petsc_view_matrix( matrix->mat() );
  //this->petsc_view_vector( rhs->vec() );
 
  // KSP solve
  ierr = KSPSolve (_ksp,rhs->vec(),solution->vec() ); CHKERRABORT(this->comm().get(), ierr);
  ierr = KSPGetIterationNumber (_ksp, &its);          CHKERRABORT(this->comm().get(), ierr);
  ierr = KSPGetResidualNorm (_ksp, &final_resid);     CHKERRABORT(this->comm().get(), ierr);

  // output the convergence infomation
  if(std::abs(final_resid)>_atol) {
    PetscPrintf(this->comm().get(), "   Linear solver does NOT converged after %d iteration,",its);
  }
  //else {
    //PetscPrintf(this->comm().get(), "   Linear solver converged after %d iteration,",its);
    //PetscPrintf(this->comm().get(), " and the residual norm is %E.\n\n",final_resid);
  //}
 
  // Update the system after the solve
  system.update();
  //t2 = MPI_Wtime();
  this->comm().barrier();
  //PetscPrintf(this->comm().get(),"   Time used to solve the linear equation Ax=b is %f\n",t2-t1);
  //std::cout << "\nTime used to solve the linear equation Ax=b is " <<t2-t1<<" s\n\n";
 
  STOP_LOG("solve()", "SolverPoisson");
}

#endif
