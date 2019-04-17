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


#pragma once

#include "libmesh/petsc_macro.h" // ahead of LIBMESH_HAVE_PETSC
#include "libmesh/libmesh_config.h"

// make sure libMesh has been complied with PETSc
#ifdef LIBMESH_HAVE_PETSC

// C++ includes
# include <iostream>

// #include <memory>

// Local includes libmesh headers
# include "libmesh/equation_systems.h"
# include "libmesh/parallel_object.h"
# include "libmesh/reference_counted_object.h"

// include PETSc KSP solver
EXTERN_C_FOR_PETSC_BEGIN
# include <petscksp.h>
EXTERN_C_FOR_PETSC_END

using libMesh::EquationSystems;
using libMesh::ReferenceCountedObject;
using libMesh::ParallelObject;
using libMesh::Real;


/*
 * This class provides basic template for solving different
 * kinds of partial differential equations.
 */

enum SystemSolverType
{
  superLU_dist,
  field_split,
  user_define
};


class Solver : public ReferenceCountedObject<Solver>,
               public ParallelObject {
public:

  /*
   * Constructor
   */
  Solver(EquationSystems& es);


  /*
   * Constructor
   */
  Solver(EquationSystems      & es,
         const SystemSolverType solver_type);


  /*
   * Destructor
   */
  ~Solver();


  /*
   * Init the KSP solver:
   * The system matrix needs to be assembled before calling this init function!
   */
  virtual void init_ksp_solver() = 0;


  /*
   * Solve the equation system Ax = b
   */
  virtual void solve() = 0;


  /*
   * Return if the ksp solver is initialized or not.
   */
  const bool is_ksp_initialized() const {
    return _is_init;
  }

  /*
   * Set the solver type
   */
  void set_solver_type(const SystemSolverType solver_type);


  /*
   * Petsc View
   */
  void petsc_view_is(IS is_p) const;
  void petsc_view_vector(Vec vector) const;
  void petsc_view_matrix(Mat matix) const;

protected:

  // EquationSystems
  EquationSystems& _equation_systems;

  // the type of the solver used for Stokes
  SystemSolverType _solver_type;

  // Label if the system is initialized.
  // If not, the destructor cannot destroy PETSc objects.
  bool _is_init;
};

#endif // #ifdef LIBMESH_HAVE_PETSC
