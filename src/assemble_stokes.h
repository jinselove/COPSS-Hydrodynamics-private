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
#include "analytical_solution_stokes.h"
#include "ggem_stokes.h"

/*! \brief This class provides the basic components
 * for assembling the matrix and vector when solving
 * Stokes equations.
 *
 * For details of the numerical discretization, refer to
 * The finite element method in heat transfer and fluid dynamics (3rd ed)
 * J.N. Reddy and D.K. Gartling. 2010, CRC Press
 */

class AssembleStokes : public AssembleSystem {
public:

  /*! \brief Constructor

     @param[in,out] es EquationSystem
   */
  AssembleStokes(EquationSystems  & es,
                 const std::string& name);


  /*! \brief Destructor

   */
  ~AssembleStokes();


  /*! \brief Assemble the Global Matrix K

     \param[in] system_name Name of the system (should be "Stokes")
     \param[in] Option options of assembling the system ("disturbed" or
        "undisturbed")
     \param[out] Ke Add element matrix to system

   */
  void assemble_global_K(const std::string& system_name,
                         const std::string& option) override;


  /*! \brief Assemble the Global force vector F

     @param[in] system_name Name of the system (should be "Stokes")
     @param[in] option Options of assembling the system ("disturbed" or
        "undisturbed")
     @param[out] Fe Add rhs vector to system.
   */
  void assemble_global_F(const std::string& system_name,
                         const std::string& option,
                         const bool is_brownian = false);


  /*! \brief Assemble the element matrix K_IJ

      Reinit and compute the element matrix K_ij, which will be added into K
      matrix after calling assemble_global_K(). Size of this submatrix is
      n_u_dofs * n_u_dofs = n_v_dofs * n_v_dofs = n_w_dofs * n_w_dofs
   */
  void assemble_element_KIJ(const std::vector<Real>                      & JxW,
                            const std::vector<std::vector<RealGradient> >& dphi,
                            const unsigned int                             n_u_dofs,
                            const unsigned int                             I,
                            const unsigned int                             J,
                            DenseMatrix<Number>                          & Kij);


  /*! \brief Assemble the element matrices Q_I, i.e., kup, kvp, kwp, for
     pressure

      These element matrices will be added to Ke after calling
         assemble_global_K().
      This function only exists for Stokes equation.
   */
  void assemble_element_QI(const std::vector<Real>                      & JxW,
                           const std::vector<std::vector<RealGradient> >& dphi,
                           const std::vector<std::vector<Real> >        & psi,
                           const unsigned int                             n_v_dofs,
                           const unsigned int                             n_p_dofs,
                           const unsigned int                             I,
                           DenseMatrix<Number>                          & Qi);


  /*! \brief Assemble the element mass matrix M for preconditioning matrix.

      This function only exists for Stokes equation.
   */
  void assemble_element_MIJ(const std::vector<Real>              & JxW,
                            const std::vector<std::vector<Real> >& phi,
                            const unsigned int                     n_v_dofs,
                            DenseMatrix<Number>                  & Mij);


  /*! \brief Assemble function for the right-hand-side in Stokes equation.

      This calculates each element's contribution to the right-hand-side vector.
   */
  void compute_element_rhs(const Elem                   *elem,
                           const unsigned int&            n_u_dofs,
                           FEBase                      & fe_v,
                           const std::vector<dof_id_type>& n_list,
                           const bool                  & pf_flag,
                           const std::string           & option,
                           DenseVector<Number>         & Fe);


  /*! \brief select sides on the boundary for all elements
   *
   */
  void select_boundary_side(const Elem *elem,
                            const std::string& system_name) override;


  /*! \brief Apply BCs by penalty method.

   */
  void apply_bc_by_penalty(const Elem          *elem,
                           const std::string  & matrix_or_vector,
                           DenseMatrix<Number>& Ke,
                           DenseVector<Number>& Fe,
                           const std::string  & option);


  /*! \brief Define the pressure jump at the inlet and outlet of the channel

   */
  Real boundary_pressure_jump(const std::string& which_side) const;


  /*! \brief Define the pressure jump at the inlet and outlet of the channel

   */
  Real boundary_traction(const std::string& which_side) const;


  /*! \brief Pointer to ggem_stokes
   */
  void init_ggem_stokes(const std::string& system_name);


  /*! \brief Pointer to analytical_solution
   */
  AnalyticalSolutionStokes* get_analytical_solution() {
    return analytical_solution;
  }

  /*! \brief Pointer to ggem_stokes
   */
  GGEMStokes* get_ggem_stokes() {
    return ggem_stokes;
  }

private:

  // ! Get a reference to AnalyticalSolution
  AnalyticalSolutionStokes *analytical_solution = nullptr;

  // ! Get a reference to GGEMStokes
  GGEMStokes *ggem_stokes = nullptr;

  // vector stores dof sizes for all elems
  std::vector<unsigned int> _n_dofs;
  
  // dof indices
  std::vector<std::vector<dof_id_type> > _dof_indices;
  std::vector<unsigned int> _n_u_dofs;
  std::vector<unsigned int> _n_p_dofs;
  std::vector<unsigned int> _n_uvw_dofs;
  std::vector<std::vector<dof_id_type>> _dof_indices_u;
  std::vector<std::vector<dof_id_type>> _dof_indices_p;
};
