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

// Local Includes -----------------------------------
#include "assemble_nernst_planck.h"
#include "solver_np.h"
#include "pm_linear_implicit_system.h"

namespace libMesh {
/*
 * The PMSystemNP is designed to solve the Nernst-Planck
 * equation
 */

class PMSystemNP : public PMLinearImplicitSystem {
public:

  /**
   * Constructor.
   */
  PMSystemNP(EquationSystems  & es,
             const std::string& name,
             const unsigned int number); // number of systems


  /**
   * Destructor.
   */
  virtual ~PMSystemNP();


  /**
   * The type of system.
   */
  typedef PMSystemNP sys_type;


  /**
   * @returns a clever pointer to the system.
   */
  sys_type& system() {
    return *this;
  }

  /**
   * Clear all the data structures associated with the system.
   */
  void clear();


  /**
   * Init ion id for this NP system and validate if system name is correctly set
   */
  void attach_ion_type(const int& _ion_id, const std::string& _ion_name);

  /**
   * Assemble the system matrix.
   * option ==
   */
  void assemble_matrix(const std::string& system_name,
                       const std::string& option) override;


  /**
   * Assemble the system rhs.
   */
  void assemble_rhs(const std::string& system_name,
                    const std::string& option) override;


  /**
   * Solve the system.
   * option = ...
   * re_init = true => re-assemble the matrix and reinit the KSP solver.
   */
  void solve(const std::string& option) override;


  /**
   * Initialize parameters for this NP system
   * This function should be called in solve() function when re_init is true
   */
  void init_params();


  /**
   * Add the local solution to the global solution
   * This function does not apply to NP system
   * but need to re-implement it here to avoid compilation error
   */
  void add_local_solution() override {};


  /**
   * Compute the L2-error in an unbounded domain
   * This function does not apply to NP system
   */
  void test_l2_norm(bool& neighbor_list_update_flag) override {};


  /**
   * update system solution for output equation systems
   * fixme:implement this function after NP system is built
   */
  void update_solution_for_output(const std::string& solution_name = "total")
    override {};


  /**
   * Test the concentration profile for a preset test systems
   * this function is debug and validation purpose
   */
   void test_concentration_profile();


  /**
  * Return the NPSolver
  */
  SolverNP& solver_np() {
    return _solver_np;
  }

  // time stepping for NP system
  Real dt_np;

  // ion id
  int ion_id;

  // ion name
  std::string ion_name;

  // ion diffusivity
  Real ion_diffusivity;

  // ion valence
  int ion_valence;

  // Dirichlet Boundary id for this ion
  std::vector<unsigned int> boundary_id_dirichlet_np;

  // Dirichlet Boundary value for this ion
  std::vector<Real> boundary_value_dirichlet_np;

  // equilibrium tolerance for ion concentration
  Real equil_tol;

private:

  // NP solver
  SolverNP _solver_np;

  // Assemble NP system
  AssembleNP *_assemble_np = nullptr;

  // Get a pointer to AnalyticalSolutionNP
  AnalyticalSolutionNP *analytical_solution = nullptr;



};
} // end namespace libMesh
