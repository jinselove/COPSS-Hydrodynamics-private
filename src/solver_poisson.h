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

// Local includes
#include "solver.h"


/*
 * This class solves Poisson equation using either an iterative
 * solver with a preconditioner or direct solver (for example, superLU).
 */
class SolverPoisson : public Solver
{
public:

  /*
   * Constructor
   */
  SolverPoisson(EquationSystems& es);


  /*
   * Constructor
   */
  SolverPoisson(EquationSystems& es,
                const SystemSolverType solver_type);

 
  /*
   * Destructor
   */
  ~SolverPoisson();

 
  /*
   * Init the KSP solver:
   * The system matrix needs to be assembled before calling this init function! 
   */
  void init_ksp_solver();
 

  /*
   * solve the equation system Ax = b
   */
  void solve();

 
private:
 
  // solver relative tolerance
  PetscReal _rtol;
 
  // solver absolute tolerance
  PetscReal _atol;
 
  // (iterative) solver maximum iteration
  PetscInt _max_it;

  // KSP sover
  KSP _ksp;
};
