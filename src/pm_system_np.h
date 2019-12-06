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


  /**
   * Set up the initial condition for NP system
   * By default, we need to initialize the NP system solution at time 0
   * When input parameter relax_t_final > 0., we will relax the NP system
   * until relax_t_final by time step dt with poisson system
   * on but stokes system off. The purpose of doing such relaxation is to
   * make sure ion cloud accommodate with charged particles so that the
   * simulation does not takes forever in the initial stages. Rule of thumb
   * of relax_t_final is ~2 * diffusion time of ions over the particle radius
   */
   // fixme: current implementation treats boundary elements and bulk
   //  elements the same due to a technical issue. Need to fix the issue if
   //  thinking differently.
  void init_cd(const Real& relax_t_final = 0.);


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
   * update system solution for output equation systems
   * There is only 'total' solution for NP system, no need to update
   */
  void update_solution_before_output(const std::string& solution_name = "total")
    override
    {
      // we don't need solution_name input, but we have to put it in
      // function
      // definition since it's overriding the virtual function in parent class
      (void) solution_name;
    };

  /**
   * override the resume_solution_after_output defined in
   * PMLinearImplicitSystem class.
   * We don't need to do anything in this resume function for NP system since
   * we didn't do anything in the update_solution_before_output function
   */
  void resume_solution_after_output() override {};


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

  // ion diffusivity
  Real ion_diffusivity;

  // ion valence
  int ion_valence;

  // Dirichlet Boundary id for this ion
  std::vector<unsigned int> boundary_id_dirichlet_np;

  // Dirichlet Boundary value for this ion
  std::vector<Real> boundary_value_dirichlet_np;

  // reinitialize initial condition
  bool set_init_cd;

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
