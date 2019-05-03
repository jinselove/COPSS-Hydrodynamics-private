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


#include "solver.h"

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
# include "pm_system_stokes.h"

using namespace libMesh;


// =======================================================================================
Solver::Solver(EquationSystems& es)
  : ParallelObject(es),
  _equation_systems(es),
  _solver_type(field_split),
  _is_init(false)
{}

// =======================================================================================
Solver::Solver(EquationSystems      & es,
               const SystemSolverType solver_type)
  : ParallelObject(es),
  _equation_systems(es),
  _solver_type(solver_type),
  _is_init(false)
{}

// =======================================================================================
Solver::~Solver()
{}

// =======================================================================================
void Solver::set_solver_type(const SystemSolverType solver_type)
{
  _solver_type = solver_type;
}

// =======================================================================================
void Solver::petsc_view_is(IS is_p) const
{
  int ierr;
  PetscViewer is_viewer;

  ierr = PetscPrintf(this->comm().get(), "View IS info:\n"); CHKERRABORT(
    this->comm().get(),
    ierr);
  ierr = PetscViewerCreate(this->comm().get(), &is_viewer); CHKERRABORT(
    this->comm().get(),
    ierr);
  ierr = PetscViewerSetType(is_viewer, PETSCVIEWERASCII);   CHKERRABORT(
    this->comm().get(),
    ierr);
  ierr = ISView(is_p, is_viewer);                           CHKERRABORT(
    this->comm().get(),
    ierr);
  ierr = PetscViewerDestroy(&is_viewer);                    CHKERRABORT(
    this->comm().get(),
    ierr);
}

// =======================================================================================
void Solver::petsc_view_vector(Vec vector) const
{
  int ierr;
  PetscViewer vec_viewer;

  ierr = PetscPrintf(this->comm().get(), "View Vec info: \n"); CHKERRABORT(
    this->comm().get(),
    ierr);
  ierr = PetscViewerCreate(this->comm().get(), &vec_viewer); CHKERRABORT(
    this->comm().get(),
    ierr);
  ierr = PetscViewerSetType(vec_viewer, PETSCVIEWERASCII);   CHKERRABORT(
    this->comm().get(),
    ierr);
  ierr = VecView(vector, vec_viewer);                        CHKERRABORT(
    this->comm().get(),
    ierr);
  ierr = PetscViewerDestroy(&vec_viewer);                    CHKERRABORT(
    this->comm().get(),
    ierr);
}

/// =======================================================================================
// ///
void Solver::petsc_view_matrix(Mat matix) const
{
  int ierr;
  PetscViewer mat_viewer;

  ierr = PetscPrintf(this->comm().get(), "View Mat info: \n"); CHKERRABORT(
    this->comm().get(),
    ierr);
  ierr = PetscViewerCreate(this->comm().get(), &mat_viewer); CHKERRABORT(
    this->comm().get(),
    ierr);
  ierr = PetscViewerSetType(mat_viewer, PETSCVIEWERASCII);   CHKERRABORT(
    this->comm().get(),
    ierr);
  ierr = MatView(matix, mat_viewer);                         CHKERRABORT(
    this->comm().get(),
    ierr);
  ierr = PetscViewerDestroy(&mat_viewer);                    CHKERRABORT(
    this->comm().get(),
    ierr);
}

#endif // ifdef LIBMESH_HAVE_PETSC
