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

// C++ includes
#include <stdio.h>
#include <iostream>
#include <cstring>
#include <utility>
#include <vector>
#include <map>

// Libmesh includes
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/point.h"
#include "pm_linear_implicit_system.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;


/*! \brief This class provides the basic components
 * for assembling the matrix and vector when solving
 * Navier-Stokes equations.
 * 
 * For details of the numerical discretization, refer to
 * The finite element method in heat transfer and fluid dynamics (3rd ed)
 * J.N. Reddy and D.K. Gartling. 2010, CRC Press
 */

class AssembleNS : public ReferenceCountedObject<AssembleNS>
{
public:
  /*! \brief Constructor
  
  @param[in,out] es EquationSystem
  */
  AssembleNS(EquationSystems& es);
  
  
  /*! \brief Destructor
  
  */
  ~AssembleNS();


  /*! \brief Assemble the Global Matrix K
  
    \param[in] system_name Name of the system (should be "Stokes")
    \param[in] Option options of assembling the system ("disturbed" or "undisturbed")
    \param[out] Ke Add element matrix to system

  */
  void assemble_global_K(const std::string& system_name,
                         const std::string& option);
  
  
  /*! \brief Assemble the Global force vector F
  
    @param[in] system_name Name of the system (should be "Stokes")
    @param[in] option Options of assembling the system ("disturbed" or "undisturbed")
    @param[out] Fe Add rhs vector to system.
  */
  void assemble_global_F(const std::string& system_name,
                         const std::string& option);
  
  
  /*! \brief Assemble the element matrix K_IJ
      
      Reinit and compute the element matrix K_ij, which will be added into K
      matrix after calling assemble_global_K(). Size of this submatrix is
      n_u_dofs * n_u_dofs = n_v_dofs * n_v_dofs = n_w_dofs * n_w_dofs
  */
  void assemble_element_KIJ(const std::vector<Real>& JxW,
                            const std::vector<std::vector<RealGradient> >& dphi,
                            const Real&     mu,
                            const unsigned int n_u_dofs,
                            const unsigned int I,
                            const unsigned int J,
                            DenseMatrix<Number>& Kij);

  
  /*! \brief Assemble the element matrices Q_I, i.e., kup, kvp, kwp
  
      These element matrices will be added to Ke after calling assemble_global_K()
  */
  void assemble_element_QI(const std::vector<Real>& JxW,
                           const std::vector<std::vector<RealGradient> >& dphi,
                           const std::vector<std::vector<Real> >& psi,
                           const unsigned int n_v_dofs,
                           const unsigned int n_p_dofs,
                           const unsigned int I,
                           DenseMatrix<Number>& Qi);
  
  
  /*! \brief Assemble the element mass matrix M.
  
  */
  void assemble_element_MIJ(const std::vector<Real>& JxW,
                            const std::vector<std::vector<Real> >& phi,
                            const unsigned int n_v_dofs,
                            DenseMatrix<Number>& Mij);
  
  
  /*! \brief  assemble function for slit channel
  
  */
  void compute_element_rhs(const Elem*     elem,
                           const unsigned int n_u_dofs,
                           FEBase& fe_v,
                           const std::vector<std::size_t> n_list,
                           const bool& pf_flag,
                           const std::string& option,
                           const Real& alpha,
                           DenseVector<Number>& Fe);
  
  /*! \brief Apply BCs by penalty method.
  
  */
  void apply_bc_by_penalty(const Elem* elem,
                           const std::string& matrix_or_vector,
                           DenseMatrix<Number>& Ke,
                           DenseVector<Number>& Fe,
                           const std::string& option);
    
  
  /*! \brief Penalize element matrix or vector with a large number.

  */
  void penalize_elem_matrix_vector(DenseMatrix<Number>& Ke,
                                   DenseVector<Number>& Fe,
                                   const std::string & matrix_or_vector,
                                   const unsigned int& var_number,     // variable number
                                   const unsigned int& local_node_id,  // local node id to be penalized
                                   const unsigned int& n_nodes_elem,   // vel-node number of elem!
                                   const Real& penalty,
                                   const Real& value);
  
  
  /*! \brief Define the pressure jump at the inlet and outlet of the channel
  
  */
  Real boundary_pressure_jump(const std::string& which_side) const;
  
  
  /*! \brief Define the pressure jump at the inlet and outlet of the channel
  
  */
  Real boundary_traction(const std::string& which_side) const;
  
  
private:
  // Equation systems
  EquationSystems& _eqn_sys;
 
  // system dimension
  unsigned int _dim;

  // mesh
  MeshBase& _mesh;
  // the particle-mesh (linear implicit) system

  // Boundary ids - build_square(): (2D);  build_cube(): (3D)
  const std::vector<boundary_id_type> _boundary_id_3D = {4,2,1,3,0,5}; // left_id(x-)=4, right_id(x+)=2, bottom_id=1(y-), top_id=3(y+), back_id=0(z-), top_id=5(z+)

  const std::vector<std::string> _boundary_name_3D = {"left", "right", "bottom", "top", "back", "front"};
};

