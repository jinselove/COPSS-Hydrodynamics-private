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


// C++ includes
#include <algorithm>
#include <math.h>

// Libmesh includes
#include "libmesh/exodusII_io.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/auto_ptr.h"

// For systems of equations the DenseSubMatrix and DenseSubVector provide
// convenient ways
// for assembling the element matrix and vector on a component-by-component
// basis.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/mesh.h"

// User defined header includes
#include "pm_toolbox.h"
#include "assemble_system.h"


// ==================================================================================
AssembleSystem::AssembleSystem(EquationSystems& es)
  : _eqn_sys(es),
  _mesh(es.get_mesh())
{
  _dim = es.get_mesh().mesh_dimension();
  _boundary_sides.resize(1); // intialize _boundary_sides matrix
}

// ==================================================================================
AssembleSystem::~AssembleSystem()
{
  // do nothing
}

// ==================================================================================


// =======================================================================================
void AssembleSystem::penalize_elem_matrix_vector(DenseMatrix<Number>& Ke,
                                                 DenseVector<Number>& Fe,
                                                 const std::string  & matrix_or_vector,
                                                 const unsigned int & var_number,
                                                 const unsigned int & local_node_id,
                                                 const unsigned int & n_nodes_elem,
// total vel-node #!
                                                 const Real         & penalty,
                                                 const Real         & value)
{
  START_LOG("penalize_elem_matrix_vector()", "AssembleSystem");

  // Penalize Matrix or vector.
  // n_nodes_elem: number of nodes in an element
  const unsigned int n = local_node_id;

  if (matrix_or_vector == "matrix")
  {
    Ke(var_number * n_nodes_elem + n, var_number * n_nodes_elem + n) += penalty;
  }
  else if (matrix_or_vector == "vector") // Penalize Right-hand-side.
  {
    Fe(var_number * n_nodes_elem + n) += penalty * value;
  }
  else if (matrix_or_vector == "both")
  {
    Ke(var_number * n_nodes_elem + n, var_number * n_nodes_elem + n) += penalty;
    Fe(var_number * n_nodes_elem + n) += penalty * value;
  }
  else
  {
    libmesh_assert("*** error in AssembleSystem::penalize_elem_matrix_vector():");
    libmesh_assert(               "*** ---> invalid argument: matrix_or_vector!");
    libmesh_error();
  }

  STOP_LOG("penalize_elem_matrix_vector()", "AssembleSystem");
}
