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
#include "assemble_nernst_planck.h"
#include "pm_system_np.h"

// ==================================================================================
AssembleNP::AssembleNP(EquationSystems& es,
                       const std::string& name)
  : AssembleSystem(es)
{
  libmesh_assert_equal_to(name, "NP");
  analytical_solution = new AnalyticalSolutionNP(name);
  // Initialize boundary sides for BC
  _boundary_sides_dirichlet_poisson.resize(_mesh.n_elem());
}

// ==================================================================================
AssembleNP::~AssembleNP()
{
  delete analytical_solution; analytical_solution = nullptr;
}

// ==================================================================================
void AssembleNP::assemble_global_K(const std::string& system_name,
                                   const std::string& option)
{
  START_LOG("assemble_global_K()", "AssembleNP");


  STOP_LOG("assemble_global_K()", "AssembleNP");
}

// ==================================================================================
void AssembleNP::assemble_global_F(const std::string& system_name,
                                   const std::string& option)
{
  START_LOG("assemble_global_F()", "AssembleNP");

  STOP_LOG("assemble_global_F()", "AssembleNP");
}

// ==================================================================================
void AssembleNP::compute_element_rhs(const Elem  *elem,
                                     const unsigned int&  n_dofs,
                                     const std::vector<Real>& JxW,
                                     const std::vector<std::vector<Real>>& c,
                                     const std::vector<Point> q_xyz,
                                     const std::vector<std::size_t>n_list,
                                     const bool& pf_flag,
                                     const std::string& option,
                                     DenseVector<Number>& Fe)
{
  START_LOG("compute_element_rhs()", "AssembleNP"); // libMesh log


  STOP_LOG("compute_element_rhs()", "AssembleNP");
}


// ==================================================================================
void AssembleNP::select_boundary_side(const Elem *elem)
{
  START_LOG("select_boundary_side()", "AssembleNP");

  STOP_LOG("select_boundary_side()", "AssembleNP");
}

// =======================================================================================
void AssembleNP::apply_bc_by_penalty(const Elem          *elem,
                                     const std::string  & matrix_or_vector,
                                     DenseMatrix<Number>& Ke,
                                     DenseVector<Number>& Fe,
                                     const std::string  & option)
{
  START_LOG("apply_bc_by_penalty()", "AssembleNP"); // libMesh log

  STOP_LOG("apply_bc_by_penalty()", "AssembleNP");
} // end of function apply_bc_by_penalty()
