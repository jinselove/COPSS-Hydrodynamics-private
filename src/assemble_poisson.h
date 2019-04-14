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
#include "ggem_poisson.h"
#include "analytical_solution_poisson.h"

/*! \brief This class provides the basic components
 * for assembling the matrix and vector for solving
 * Poisson equation with point charges suspended in
 * low-Reynolds number fluids.
 *
 * For details of the implementation, refer to
 * J. Chem. Phys. 142, 014108 (2015)
 * J. Chem. Theory Comput. 2018, 14, 4901-4913
 *
 * The method uses Ewald-split to seperate total
 * electrical potential in to local and global
 * components. This class is for solving global
 * electrical potential using finite element method.
 *
 */

class AssemblePoisson : public AssembleSystem
{
public:
  /*! \brief Constructor

  @param[in,out] es EquationSystem
  */
  AssemblePoisson(EquationSystems& es,
                  const std::string& name);


  /*! \brief Destructor

  */
  ~AssemblePoisson();


  /*! \brief Assemble the Global Matrix K

    \param[in] system_name Name of the system (should be "Poisson")
    \param[out] Ke Add element matrix to global matrix K

  */
  void assemble_global_K(const std::string& system_name,
                         const std::string& option) override;


  /*! \brief Assemble the Global force vector F

    @param[in] system_name Name of the system (should be "Poisson")
    @param[out] Fe Add rhs vector to global vector F
  */
  void assemble_global_F(const std::string& system_name,
                         const std::string& option) override;


  /*! \brief Assemble function on each element for the right-hand-side in Poisson equation.

      This calculates each element's contribution to the right-hand-side vector.
  */
  void compute_element_rhs(const Elem* elem,
                           const unsigned int n_u_dofs,
                           FEBase& fe_v,
                           const std::vector<std::size_t> n_list,
                           const bool& pf_flag,
                           const std::string& option,
                           DenseVector<Number>& Fe) override;


  /*! \brief Select sides on Dirichlet and Neumann boundaries for all elements
  *
  */
  void select_boundary_side(const Elem* elem) override;


  /*! \brief Apply Dirichlet BC by penalty method to impose electrical potential on relevant boundaries.

  */
  void apply_bc_by_penalty(const Elem* elem,
                           const std::string& matrix_or_vector,
                           DenseMatrix<Number>& Ke,
                           DenseVector<Number>& Fe,
                           const std::string& option) override;


  /*! \brief Apply Neumann BC to impose surface charge density on relevant boundaries.

  */
  void apply_bc_neumann(const Elem* elem,
                        FEBase& fe_phi,
                        FEBase& fe_face,
                        DenseVector<Number>& Fe);


  /*! \brief Initialize ggem_poisson for local field calcualtions

  */
  void init_ggem_poisson(const std::string& system_name);


  /*! \brief Pointer to analytical_solution
  */
  AnalyticalSolutionPoisson* get_analytical_solution() {return analytical_solution;};
  

  /*! \brief Pointer to ggem_poisson for local field calculations

  */
  GGEMPoisson* get_ggem_poisson() {return ggem_poisson;};



private:

  // Boundary sides that Dirichlet and Neumann BCs are applied.
  std::vector<std::vector<unsigned int> > _boundary_sides_dirichlet_poisson, _boundary_sides_neumann_poisson;

  // Get a reference to GGEMPoisson
  GGEMPoisson* ggem_poisson;
  
  // Get a reference to AnalyticalSolutionPoisson
  AnalyticalSolutionPoisson* analytical_solution;
 
};
