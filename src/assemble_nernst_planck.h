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
  AssembleNP(EquationSystems& es, const std::string& name);


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


  /*! \brief Assemble function for the right-hand-side in NP equation.

      This calculates each element's contribution to the right-hand-side vector.
   */
  void compute_element_rhs(const Elem  *elem,
                           const unsigned int&  n_dofs,
                           const std::vector<Real>& JxW,
                           const std::vector<std::vector<Real>>& c,
                           const std::vector<Point> q_xyz,
                           const std::vector<std::size_t>n_list,
                           const bool& pf_flag,
                           const std::string& option,
                           DenseVector<Number>& Fe);


  /*! \brief select sides on the boundary for all elements
   *
   */
  void select_boundary_side(const Elem *elem) override;

  /*! \brief Apply BCs by penalty method.

   */
  void apply_bc_by_penalty(const Elem          *elem,
                           const std::string  & matrix_or_vector,
                           DenseMatrix<Number>& Ke,
                           DenseVector<Number>& Fe,
                           const std::string  & option) override;

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
  // boundary value is an vector that contains the ion concentration of all
  // ion species at this boundary.
  std::vector<std::vector<std::tuple<unsigned int,
                                     unsigned int,
                                     std::vector<Real>
                                     >>> _boundary_sides_dirichlet_poisson;

  // Get a reference to AnalyticalSolutionNP
  AnalyticalSolutionNP *analytical_solution = nullptr;
};
