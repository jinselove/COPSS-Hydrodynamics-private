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


/*! \brief This base class provides the template for
 * assembling the matrix and right-hand-side vector
 * when solving any partial differential equations
 * using finite element method.
 */

class AssembleSystem : public ReferenceCountedObject<AssembleSystem>
{
public:
  /*! \brief Constructor
 
  @param[in,out] es EquationSystem
  */
  AssembleSystem(EquationSystems& es);
 
 
  /*! \brief Destructor
 
  */
  ~AssembleSystem();


  /*! \brief Assemble the Global Matrix K
 
    \param[in] system_name Name of the system (could be "Stokes", "Poisson", or "NP")
    \param[in] Option options of assembling the system ("disturbed" or "undisturbed") for Stokes equation
    \param[out] Ke Add element matrix to system

  */
  virtual void assemble_global_K(const std::string& system_name,
                                 const std::string& option){}; //{} defines virtual function
 
 
  /*! \brief Assemble the Global force vector F
 
    @param[in] system_name Name of the system (could be "Stokes", "Poisson", or "NP")
    @param[in] option Options of assembling the system ("disturbed" or "undisturbed") for Stokes equation
    @param[out] Fe Add rhs vector to system.
  */
  virtual void assemble_global_F(const std::string& system_name,
                                 const std::string& option){};
 
 
  /*! \brief Assemble the element matrix K_IJ
 
      Reinit and compute the element matrix K_ij, which will be added into K
      matrix after calling assemble_global_K(). For Stokes equation, size of
      this submatrix is n_u_dofs * n_u_dofs = n_v_dofs * n_v_dofs = n_w_dofs * n_w_dofs
  */
  virtual void assemble_element_KIJ(const std::vector<Real>& JxW,
                                    const std::vector<std::vector<RealGradient> >& dphi,
                                    const unsigned int n_u_dofs,
                                    const unsigned int I,
                                    const unsigned int J,
                                    DenseMatrix<Number>& Kij){};
 
 
  /*! \brief  Assemble function for calculating each element's contribution to
   *          the right-hand-side vector. It is used in assemble_global_F()
 
  */
  virtual void compute_element_rhs(const Elem* elem,
                                   const unsigned int n_u_dofs,
                                   FEBase& fe_v,
                                   const std::vector<std::size_t> n_list,
                                   const bool& pf_flag,
                                   const std::string& option,
                                   const Real& alpha,
                                   DenseVector<Number>& Fe){};

 
  /*! \brief Apply Boundary Conditions by penalty method.
 
  */
  virtual void apply_bc_by_penalty(const Elem* elem,
                                   const std::string& matrix_or_vector,
                                   DenseMatrix<Number>& Ke,
                                   DenseVector<Number>& Fe,
                                   const std::string& option){};

 
  /*! \brief Assemble int_force matrix for every element,
   * this includes Gaussian quadrature weights multiplied by shape functions.
   * The product is calculated once and is stored in _int_force.

  */
  void assemble_int_force(const Elem*     elem,
                          const unsigned int n_u_dofs,
                          FEBase& fe_v);


  /*! \brief select sides on the boundary for all elements
   *
   */
  void select_boundary_side(const Elem* elem);

 
  /*! \brief Universal function to penalize element matrix or vector with a large number.

  */
  void penalize_elem_matrix_vector(DenseMatrix<Number>& Ke,
                                   DenseVector<Number>& Fe,
                                   const std::string & matrix_or_vector,
                                   const unsigned int& var_number,     // variable number
                                   const unsigned int& local_node_id,  // local node id to be penalized
                                   const unsigned int& n_nodes_elem,   // vel-node number of elem!
                                   const Real& penalty,
                                   const Real& value);
 
 
 
protected:
  // Equation systems
  EquationSystems& _eqn_sys;
 
  // system dimension
  unsigned int _dim;

  // Fluid mesh
  MeshBase& _mesh;

  // Boundary ids
  // defaults in cubic geometry generated from libMesh
  const std::vector<boundary_id_type> _boundary_id_3D = {4,2,1,3,0,5};
  const std::vector<std::string> _boundary_name_3D = {"left", "right", "bottom", "top", "back", "front"};

  // int_force matrix
  // this matrix stores the product of JxW[qp] * phi[k][qp]
  // size = num_elem * (n_u_dofs * n_quad_points) 
  std::vector<std::vector<Real>> _int_force;

  // vector stores q_xyz size
  // size = num_elem * q_xyz.size()
  std::vector<std::vector<Point>> _q_xyz;

  // vector stores dof sizes for all elems
  std::vector<unsigned int> _n_dofs;
  std::vector<unsigned int> _n_u_dofs;
  std::vector<unsigned int> _n_p_dofs;
  std::vector<unsigned int> _n_uvw_dofs;

  // dof indices
  std::vector<std::vector<dof_id_type> > _dof_indices;
  std::vector<std::vector<dof_id_type> > _dof_indices_u;
  std::vector<std::vector<dof_id_type> > _dof_indices_p;

  // sides on the boundary
  std::vector<std::vector<unsigned int> > _boundary_sides;

};
