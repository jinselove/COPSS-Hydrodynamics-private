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
 * This class solves Stokes equation using either an iterative
 * solver with a preconditioner or direct solver (for example, superLU).
 * When building the preconditioner, it performs schur complement
 * reduction type solve for the saddle point problems arised from
 * mixed finite element method.
 */
class SolverStokes : public Solver {
public:

  /*
   * Constructor
   */
  SolverStokes(EquationSystems& es);


  /*
   * Constructor
   */
  SolverStokes(EquationSystems      & es,
               const SystemSolverType solver_type);


  /*
   * Destructor
   */
  ~SolverStokes();


  /*
   * Init the KSP solver:
   * The system matrix needs to be assembled before calling this init function!
   */
  void init_ksp_solver(const std::string& system_name);


  /*
   * solve the equation system Ax = b
   */
  void solve();


  /*
   * Build IS for u/p for PC field split
   *    -called by: solve()
   */
  void build_is(IS *is_v,
                IS *is_p);


  /*
   * Set up the PC for the schur complement reduction algorithm
   *    -called by solve()
   */
  void setup_schur_pc(KSP        ksp,
                      IS         is_v,
                      IS         is_p,
                      Mat       *pmat,
                      const bool userPC,
                      const bool userKSP);


  /*                                ^
   * Schur complement approximation S used as a PC for S*y1=y2, for example:
   *     S1 = A11 - A10 diag(A00)^(-1) A01;
   *     Sa = 1/v * Mp (pressure matrix, v is kinematic viscosity)
   *    -called by solve()
   */
  void setup_approx_schur_matrix(IS   is_v,
                                 IS   is_p,
                                 Mat *pmat);

private:

  // solver relative tolerance
  PetscReal _rtol;

  // solver absolute tolerance
  PetscReal _atol;

  // (iterative) solver maximum iteration
  PetscInt _max_it;

  // KSP sover
  KSP _ksp;

  // IS pointers for velocity and pressure, respectively
  IS _is_v;
  IS _is_p;

  // preconditioning matrix for Schur Complement
  Mat _schur_pmat;
};
