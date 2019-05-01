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

// ==================================================================================
AssembleNP::AssembleNP(EquationSystems& es)
  : AssembleSystem(es)
{
  // do nothing
}

// ==================================================================================
AssembleNP::~AssembleNP()
{
  // do nothing
}

// ==================================================================================
void AssembleNP::assemble_global_K(const std::string& system_name,
                                   const std::string& option)
{
  START_LOG("assemble_global_K()", "AssembleNP");

  /*! It is a good idea to make sure we are assembling the proper system.*/
  libmesh_assert_equal_to(system_name, "NP");

  // -------------------------------------------------------------------------------------------
  //  if (_pm_system.comm().rank()==0){
  //  printf("assemble_matrix_K(): The global matrix K has been assembled
  // ...\n");
  // }
  // -------------------------------------------------------------------------------------------
  return;

  STOP_LOG("assemble_global_K()", "AssembleNP");
}

// ==================================================================================
void AssembleNP::assemble_global_F(const std::string& system_name,
                                   const std::string& option)
{
  START_LOG("assemble_global_F()", "AssembleNP");

  // PerfLog perf_log("assemble_global_F");
  // perf_log.push("preparation");
  // It is a good idea to make sure we are assembling the proper system.
  libmesh_assert_equal_to(system_name, "NP");

  STOP_LOG("assemble_global_F()", "AssembleNP");

  // ---------------------------------------------------------------------------------------------
  // if (_pm_system.comm().rank()==0){
  //  printf("assemble_global_F(): The global RHS vector has been assembled
  // ...\n");
  // }
  // ---------------------------------------------------------------------------------------------
}

// ==================================================================================
void AssembleNP::compute_element_rhs(const Elem                   *elem,
                                     const unsigned int            n_u_dofs,
                                     FEBase                      & fe_v,
                                     const std::vector<std::size_t>n_list,
                                     const bool                  & pf_flag,
                                     const std::string           & option,
                                     DenseVector<Number>         & Fe)
{
  START_LOG("compute_element_rhs()", "AssembleNP"); // libMesh log

  // libmesh_assert_equal_to(system_name, "NP");

  STOP_LOG("compute_element_rhs()", "AssembleNP");
}

// ==================================================================================
void AssembleNP::assemble_element_KIJ(
  const std::vector<Real>                      & JxW,
  const std::vector<std::vector<RealGradient> >& dphi,
  const unsigned int                             n_u_dofs,
  const unsigned int                             i,
  const unsigned int                             j,
  DenseMatrix<Number>                          & Kij)
{
  START_LOG("assemble_element_KIJ()", "AssembleNP"); // libMesh log

  STOP_LOG("assemble_element_KIJ()", "AssembleNP");  // libMesh log
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
