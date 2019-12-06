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
#include "assemble_system.h"
#include "analytical_solution_np.h"

/*! \brief This class provides the basic components
 * for assembling the matrix and vector when solving
 * Nernst-Planck equations.
 *
 * For details of the implementation, refer to ...
 *
 */

class AssembleNP : public AssembleSystem {
public:

  /*! \brief Constructor

     @param[in,out] es EquationSystem
   */
  AssembleNP(EquationSystems& es,
             const std::string& name);


  /*! \brief Destructor

   */
  ~AssembleNP();


  /*! \brief Assemble the Global Matrix K

     \param[in] system_name Name of the system (should be "NP")
     \param[in] Option options of assembling the system
     \param[out] Ke Add element matrix to system

   */
  void assemble_global_K(const std::string& system_name,
                         const std::string& option) override;


  /*! \brief Assemble the Global force vector F

     @param[in] system_name Name of the system (should be "NP")
     @param[in] option Options of assembling the system
     @param[out] Fe Add rhs vector to system.
   */
  void assemble_global_F(const std::string& system_name,
                         const std::string& option) override;


  /*! \brief select sides on the boundary for all elements
   *
   */
  void select_boundary_side(const Elem *elem,
                            const std::string& system_name) override;

  /*! \brief Pointer to analytical_solution
  */
  AnalyticalSolutionNP* get_analytical_solution() {
    return analytical_solution;
  }

private:
  // Boundary sides that DirichletBCs are applied.
  // for each tuple, the first element is the side id, the second element is the
  // associated boundary id of this side, the third element is the corresponding
  // Dirichlet BC values associated with this side. For NP system, each
  // boundary value is the concentration of the ion associated with this NP
  // system
  std::vector<std::vector<std::tuple<dof_id_type,
                                     dof_id_type,
                                     Real
                                     >>> _boundary_sides_dirichlet_np;

  // For regular NP simulations, Dirichlet Boundary value is constant. Thus
  // we store penalty value of boundary surface nodes in a vector. This
  // vector does not get reinitialized unless _reinit_node_penalty is True.
  // If boundary value changes with time, then we need to reinitialize this
  // vector at every step; Notice that each node modifies part of this vector
  // and the information never gets communicated, however, since we will only
  // access to the local information, we don't need to sync the whole vector
  // over multiple cpus
  std::vector<std::vector<Real>> _node_penalty;

  // initialize reinit_node_penalty as true. It will be reset to false if
  // necessary
  bool _reinit_node_penalty = true;

  // Get a reference to AnalyticalSolutionNP
  AnalyticalSolutionNP *analytical_solution = nullptr;

  // A big number for penalty
  const Real penalty = 1.e10;

};
