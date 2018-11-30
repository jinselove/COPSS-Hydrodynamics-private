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

// For systems of equations the DenseSubMatrix and DenseSubVector provide convenient ways
// for assembling the element matrix and vector on a component-by-component basis.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/mesh.h"

// User defined header includes
#include "analytical_solution.h"
#include "ggem_system.h"
#include "pm_toolbox.h"


// ==================================================================================
AssembleSystem::AssembleSystem(EquationSystems& es)
: _eqn_sys(es),
_mesh(es.get_mesh())
{
  _dim = es.get_mesh().mesh_dimension();
  _int_force.resize(1); // initialize _int_force matrix
  _boundary_sides.resize(1); // intialize _boundary_sides matrix
}



// ==================================================================================
AssembleSystem::~AssembleSystem()
{
  // do nothing
}



// ==================================================================================
void AssembleSystem::assemble_int_force(const Elem* elem,
                                        const unsigned int n_u_dofs,
                                        FEBase& fe_v)
{
  START_LOG("assemble_int_force()", "AssembleSystem");
  const std::vector<Real>& JxW                = fe_v.get_JxW();
  const std::vector<std::vector<Real> >& phi  = fe_v.get_phi();
  const std::vector<Point>& q_xyz             = fe_v.get_xyz(); // xyz coords of quad pts
  //fe_v.reinit(elem);
  const std::size_t elem_id = elem->id();
  _int_force[elem_id].resize(n_u_dofs*q_xyz.size(),0.); // resize this row
  _q_xyz[elem_id] = q_xyz;
  for (unsigned int k=0; k<n_u_dofs; ++k){
    // loop over all gauss quadrature points
    for(unsigned int qp=0; qp<q_xyz.size(); qp++){
      _int_force[elem_id][k*q_xyz.size()+qp] = JxW[qp]*phi[k][qp];
    //_int_force[elem_id][k] += JxW[qp]*phi[k][qp];
    //printf("_int_force[elem_id=%d][k*q_xyz.size()+qp=%d] = JxW[qp=%d]*phi[k=%d][qp=%d] = %f *%f = %f\n", elem_id, k*q_xyz.size()+qp, qp,k,qp, JxW[qp],phi[k][qp], _int_force[elem_id][k*q_xyz.size()+qp]);
    }
  }
  STOP_LOG("assemble_int_force()", "AssembleSystem");
  return;
}



// ==================================================================================
void AssembleSystem::select_boundary_side(const Elem* elem)
{
  START_LOG("select_boundary_side()", "AssembleSystem");

  // Get a reference to the Particle-Mesh linear implicit system object.
  PMLinearImplicitSystem & pm_system = _eqn_sys.get_system<PMLinearImplicitSystem> ("Stokes");
  const std::vector<bool>& periodicity = pm_system.point_mesh()->pm_periodic_boundary()->periodic_direction();
  const std::vector<bool>& inlet_direction = pm_system.point_mesh()->pm_periodic_boundary()->inlet_direction();
  const std::size_t elem_id = elem->id();

  // The following loops over the sides of the element. If the element has NO
  // neighbors on a side then that side MUST live on a boundary of the domain.
  for (unsigned int s=0; s<elem->n_sides(); s++)
  {
    bool apply_bc = false;
    if (elem->neighbor(s) == NULL)
    {
      apply_bc = true; 
      // if this side is on the periodic side or pressure side, don't apply penalty
      for (int i = 0; i < _dim; i++)
      {
        if(periodicity[i] or inlet_direction[i]){
          if (_mesh.get_boundary_info().has_boundary_id(elem, s, _boundary_id_3D[2*i+0]) or
              _mesh.get_boundary_info().has_boundary_id(elem, s, _boundary_id_3D[2*i+1]))
          { apply_bc = false; }
        }
      }   
    }
    if (apply_bc) _boundary_sides[elem_id].push_back(s);
  }
  STOP_LOG("select_boundary_side()", "AssembleSystem");
}



// =======================================================================================
void AssembleSystem::penalize_elem_matrix_vector(DenseMatrix<Number>& Ke,
                                                 DenseVector<Number>& Fe,
                                                 const std::string & matrix_or_vector,
                                                 const unsigned int& var_number,
                                                 const unsigned int& local_node_id,
                                                 const unsigned int& n_nodes_elem, // total vel-node #!
                                                 const Real& penalty,
                                                 const Real& value)
{
  START_LOG("penalize_elem_matrix_vector()", "AssembleSystem");
  
  // Penalize Matrix or vector.
  // n_nodes_elem: number of nodes in an element
  const unsigned int n = local_node_id;
  if (matrix_or_vector=="matrix")
  {
    Ke(var_number*n_nodes_elem+n,var_number*n_nodes_elem+n) += penalty;
  }
  else if (matrix_or_vector=="vector")  // Penalize Right-hand-side.
  {
    Fe(var_number*n_nodes_elem + n) += penalty*value;
  }
  else if (matrix_or_vector=="both")
  {
    Ke(var_number*n_nodes_elem+n,var_number*n_nodes_elem+n) += penalty;
    Fe(var_number*n_nodes_elem + n) += penalty*value;
  }
  else
  {
    libmesh_assert("*** error in AssembleSystem::penalize_elem_matrix_vector():");
    libmesh_assert("*** ---> invalid argument: matrix_or_vector!");
    libmesh_error();
  }
  
  STOP_LOG("penalize_elem_matrix_vector()", "AssembleSystem");
}
