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
   * override get_dt function in parent class since this system needs more
   * conditions on dt to make the system stable
   */
  Real get_dt();


  /**
   * Init ion id for this NP system and validate if system name is correctly set
   */
  void attach_ion_type(const int& _ion_id, const std::string& _ion_name);


  /*
   * Set up the initial condition for NP system
   * For general systems, we initialize the NP system solution to be 0
   * everywhere
   * fixme: potentially we can reinitialize the solution to be a finite
   * number specified in the input file
   */

  void init_cd();


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
   * update system solution to total
   * the FEM solution of NP is already the total, thus we don't do anything
   */
  void update_solution_to_total() override{};

  /**
   * resume solution_to_global
   * no need to do anything in NP
   */
  void resume_solution_to_global() override {};


  /**
   * write out solution to csv file
   */
   void write_out_solution();

  /**
   * Test the concentration profile for a preset test systems
   * this function is debug and validation purpose
   * Currently test the NP system solution for first 10 steps at dt = 0.01
   */
   void test_concentration_profile();


  /**
  * Return the NPSolver
  */
  SolverNP& solver_np() {
    return _solver_np;
  }


  /**
   * Initial solution function prototype.  This gives the exact
   * solution as a function of space and time.  In this case the
   * initial condition will be taken as the exact solution at time 0
   */

  static Number init_solution (const Point & p,
                              const Parameters & parameters,
                              const std::string &,
                              const std::string &);

  /**
   * The initial solution at time 0 on the boundaries will be set by this
   * function
   */
  //fixme: implement this function to assign different initial conditions for
  //elements on boundaries versus in bulk
//  static Number init_solution_bc (const Point & p,
//                                 const Parameters & parameters,
//                                 const std::string &,
//                                 const std::string &);

  // ion id
  int ion_id;

  // np system dt
  Real dt;

  // ion name
  std::string ion_name;

  // ion diffusivity, unit = [1]
  Real ion_diffusivity;

  // ion valence, unit = [1]
  int ion_valence;

  // Dirichlet Boundary id for this ion
  std::vector<unsigned int> boundary_id_dirichlet_np;

  // Dirichlet Boundary value for this ion
  std::vector<Real> boundary_value_dirichlet_np;

  // if initial condition is set
  bool init_cd_set;

  // if NP system is relaxed without fluid
  bool relaxed;

  // NP coefficient = bjerrum_length * ion_diffusivity * ion_valence, where
  // bjerrum length = e^2/(4*pi*epsilon_0*epsilon*KB*T*bead_radius). Unit = [1]
  Real np_coeff;


private:

  // NP solver
  SolverNP _solver_np;

  // Assemble NP system
  AssembleNP *_assemble_np = nullptr;

  // Get a pointer to AnalyticalSolutionNP
  AnalyticalSolutionNP *analytical_solution = nullptr;

  // output precision (defined in input file, default is 6)
  int o_precision;
};
} // end namespace libMesh
