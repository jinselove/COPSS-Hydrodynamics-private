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

class AssembleSystem : public ReferenceCountedObject<AssembleSystem> {
public:

  /*! \brief Constructor

     @param[in,out] es EquationSystem
   */
  AssembleSystem(EquationSystems& es);


  /*! \brief Destructor

   */
  ~AssembleSystem();


  /*! \brief Assemble the Global Matrix K

     \param[in] system_name Name of the system (could be "Stokes", "Poisson", or
        "NP")
     \param[in] Option options of assembling the system ("disturbed" or
        "undisturbed") for Stokes equation
     \param[out] Ke Add element matrix to system

   */
  virtual void assemble_global_K(const std::string& system_name,
                                 const std::string& option) = 0;


  /*! \brief Assemble the Global force vector F

     @param[in] system_name Name of the system (could be "Stokes", "Poisson", or
        "NP")
     @param[in] option Options of assembling the system ("disturbed" or
        "undisturbed") for Stokes equation
     @param[out] Fe Add rhs vector to system.
   */
//  virtual void assemble_global_F(const std::string& system_name,
//                                 const std::string& option) = 0;


  /*! \brief select sides on the boundary for all elements
   *
   */
  virtual void select_boundary_side(const Elem *elem,
                                    const std::string& system_name) = 0;


  /*! \brief Universal function to penalize element matrix or vector with a
     large number.

   */
  void penalize_elem_matrix_vector(DenseMatrix<Number>& Ke,
                                   DenseVector<Number>& Fe,
                                   const std::string  & matrix_or_vector,
                                   const unsigned int & var_number,    // variable
                                                                       // number
                                   const unsigned int & local_node_id, // local
                                                                       // node
                                                                       // id to
                                                                       // be
                                                                       // penalized
                                   const unsigned int & n_nodes_elem,  // vel-node
                                                                       // number
                                                                       // of
                                                                       // elem!
                                   const Real         & penalty,
                                   const Real         & value);

protected:

  // Equation systems
  EquationSystems& _eqn_sys;

  // System dimension
  unsigned int _dim;

  // Fluid mesh
  MeshBase& _mesh;

  // Boundary ids
  const std::vector<boundary_id_type> _boundary_id_3D = { 4, 6, 3, 5, 2, 1};
  
  const std::vector<std::string> _boundary_name_3D    =
  { "x_negative", "x_positive", "y_negative", "y_positive", "z_negative", 
    "z_positive" };

  // Sides on the boundary for different SubEquationSystems. Boundaries applied
  // with
  // Dirichlet BC using penalty method are different for different systems.
  // For example, in Poisson system, we don't use pressure inlet/outlet
  // condition
  // to decide whether this side has to be included in the Dirichlet BC.
  std::vector<std::vector<unsigned int> > _boundary_sides;

  // vector stores quadrature points for all elements
  // size = num_elem * q_xyz.size()
  std::vector<std::vector<Point> > _q_xyz;

  // vector store shape function
  std::vector<std::vector<std::vector<Real>>> _phi;

  // vector stores the gradient of the shape function
  std::vector<std::vector<std::vector<RealGradient>>> _dphi;

  // Jacobian of each element
  std::vector<std::vector<Real>> _JxW;

  // number of dofs for each element
  std::vector<unsigned int> _n_dofs;

  // dof indices of each element.
  std::vector<std::vector<dof_id_type> > _dof_indices;

};
