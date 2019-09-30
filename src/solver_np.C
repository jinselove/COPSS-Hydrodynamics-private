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


#include "solver_np.h"

// make sure libMesh has been complied with PETSc
#ifdef LIBMESH_HAVE_PETSC

// include PETSc KSP solver
EXTERN_C_FOR_PETSC_BEGIN
# include <petscsys.h>
# include <petscksp.h>
# include <petscpc.h>
# include <petscis.h>
# include <petscviewer.h>
EXTERN_C_FOR_PETSC_END

// libmesh headers
# include "libmesh/libmesh_config.h"
# include "libmesh/dof_map.h"
# include "libmesh/system.h"
# include "libmesh/linear_implicit_system.h"
# include "libmesh/petsc_matrix.h"
# include "libmesh/petsc_vector.h"
# include "libmesh/petsc_linear_solver.h"

// local header
# include "pm_system_np.h"

using namespace libMesh;


// =======================================================================================
SolverNP::SolverNP(EquationSystems& es)
  : Solver(es)
{}

// =======================================================================================
SolverNP::SolverNP(EquationSystems      & es,
                   const SystemSolverType solver_type)
  : Solver(es, solver_type)
{}

// =======================================================================================
SolverNP::~SolverNP()
{

}

// ==================================================================================
void SolverNP::init_ksp_solver()
{
  START_LOG("init_ksp_solver()", "SolverNP");


  STOP_LOG("init_ksp_solver()", "SolverNP");
}

// =======================================================================================
void SolverNP::solve()
{
  START_LOG("solve()", "SolverNP");


  STOP_LOG("solve()", "SolverNP");
}

#endif // ifdef LIBMESH_HAVE_PETSC
