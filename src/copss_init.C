// Copyright (C) 2015-2016 Jiyuan Li, Xikai Jiang

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

#include "copss_init.h"

// PETSc
#ifdef LIBMESH_HAVE_PETSC
# include "petscsys.h"
#endif // ifdef LIBMESH_HAVE_PETSC

#ifdef LIBMESH_HAVE_OPENMP
# include <omp.h>
#endif // ifdef LIBMESH_HAVE_OPENMP

namespace libMesh {
CopssInit::CopssInit(int argc, char **argv)
  : LibMeshInit(argc, argv)
{
#ifdef LIBMESH_HAVE_PETSC
  PetscPopSignalHandler(); // get rid of Petsc error handler
#endif // ifdef LIBMESH_HAVE_PETSC

  // // Set the number of OpenMP threads to the same as the number of threads
  // libMesh is going to use
  // #ifdef LIBMESH_HAVE_OPENMP
  //   omp_set_num_threads(libMesh::n_threads());
  // #endif
}
} // end namespace
