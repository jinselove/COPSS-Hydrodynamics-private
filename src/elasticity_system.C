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


// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <string.h>
#include <math.h>


// Basic include file needed for the mesh functionality.
#include "libmesh/mesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/linear_implicit_system.h"

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


// local include
#include "pm_toolbox.h"
#include "elasticity_system.h"


namespace libMesh {
// ======================================================================================
ElasticitySystem::ElasticitySystem(EquationSystems  & es,
                                   const std::string& name,
                                   const unsigned int number)
  : Parent(es, name, number),
  _matrix_assembled(false)
{
  // do nothing
  // this->init_ksp_solver();
}

// ==================================================================================
ElasticitySystem::~ElasticitySystem()
{
  // Clear data
  this->clear();
}

// ==================================================================================
void ElasticitySystem::clear()
{
  // clear the parent data
  Parent::clear();

  //  this->destroy_ksp_solver(); // cause errors
}

// ==================================================================================
void ElasticitySystem::init_ksp_solver()
{
  KSPCreate(this->comm().get(), &_ksp);
}

// ==================================================================================
void ElasticitySystem::destroy_ksp_solver()
{
  if (_ksp) KSPDestroy(&_ksp);
}

// ==================================================================================
void ElasticitySystem::build_nodal_force_gravity(const std::vector<Real>& f)
{
  START_LOG("build_nodal_force_gravity()", "ElasticitySystem");

  // Get a constant reference to the mesh object and mesh dimension.
  const MeshBase& mesh   = this->get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  this->rhs->zero();

  /* ------------- validation: dim = 2 for surface element -------------- */

  //  std::ostringstream  oss;
  //  oss << "--->test: ElasticitySystem::build_nodal_force_gravity() mesh dim =
  // "<<dim;
  //  PMToolBox::output_message(oss.str(), this->comm()); oss.str("");
  // oss.clear();

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  // Numeric ids corresponding to each variable in the system
  const unsigned int u_var = this->variable_number("U"); // u_var = 0
  const unsigned int v_var = this->variable_number("V"); // v_var = 1
  unsigned int w_var       = 0;

  if (dim == 3) w_var = this->variable_number("W");      // w_var = 2


  // Build a FEBase object of the specified type
  // Note: typedef FEGenericBase<Real> FEBase in libMesh;
  FEType fe_vel_type = this->variable_type(u_var);
  UniquePtr<FEBase> fe_vel(FEBase::build(dim, fe_vel_type));

  // fe_vel->print_info(libMesh::out);

  // Gauss Quadrature rule and tell the finite element objects
  QGauss qrule(dim, fe_vel_type.default_quadrature_order()); // 3^dim pts
  fe_vel->attach_quadrature_rule(&qrule);

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>& JxW               = fe_vel->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();

  // const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();

  // A reference to the DofMap object for this system.
  const DofMap& dof_map = this->get_dof_map();
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u, dof_indices_v, dof_indices_w;

  // element nodal vector
  DenseVector<Number>    Fe;
  DenseSubVector<Number> Fu(Fe),    Fv(Fe),    Fw(Fe);


  // Now we will loop over all the elements in the mesh that live
  // on the local processor, and compute the element rhs Fe.
  MeshBase::const_element_iterator el =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently working on.
    const Elem *elem = *el;

    // Get the degree of freedom indices for the current element.
    dof_map.dof_indices(elem, dof_indices);
    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_v, v_var);

    const unsigned int n_dofs   = dof_indices.size();
    const unsigned int n_u_dofs = dof_indices_u.size();
    const unsigned int n_v_dofs = dof_indices_v.size();
    unsigned int n_w_dofs       = 0;

    if (dim == 3)
    {
      dof_map.dof_indices(elem, dof_indices_w, w_var);
      n_w_dofs = dof_indices_w.size();
    }

    // NOTE: here JxW and dphi and other element quantities are not computed.
    // This will be done in the elem loop after fe->reinit()
    fe_vel->reinit(elem);

    // elem->print_info(libMesh::out); libMesh::out<<"\n\n";
    // fe_vel->print_info(libMesh::out); libMesh::out<<"\n";

    // Zero the element right-hand-side vector before summing them.
    // DenseSubVector.reposition () member takes the (row_offset, row_size)
    Fe.resize(n_dofs); // Ke.resize (n_dofs, n_dofs);
    Fu.reposition(u_var * n_u_dofs, n_u_dofs);
    Fv.reposition(v_var * n_u_dofs, n_v_dofs);

    if (dim == 3) Fw.reposition(w_var * n_u_dofs, n_w_dofs);

    // Loop over each Gauss pt
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      // Assemble the subvector
      for (unsigned int i = 0; i < n_u_dofs; i++)
      {
        Fu(i) += JxW[qp] * phi[i][qp] * f[0];
        Fv(i) += JxW[qp] * phi[i][qp] * f[1];

        if (dim == 3) Fw(i) += JxW[qp] * phi[i][qp] * f[2];
      }
    }


    // If this assembly program were to be used on an adaptive mesh,
    // we would have to apply any hanging node constraint equations.
    dof_map.constrain_element_vector(Fe, dof_indices);

    // Add the element matrix and rhs vector to the global system.
    this->rhs->add_vector(Fe, dof_indices);
  } // end for el-loop
  this->rhs->close();


  /* Note the force vector computed above is not ordered
   * according to node #, but according to the degree of freedom.
   * We hence need to reorder the vector and distribute it to all processors.
   */
  this->rhs->localize(_nodal_force);

  // test view vector
  this->solution = this->rhs->clone();

  //  this->solution->print_global();

  STOP_LOG("build_nodal_force_gravity()", "ElasticitySystem");
}

// ==================================================================================
std::vector<Real>ElasticitySystem::mesh_size() const
{
  START_LOG("mesh_size()", "ElasticitySystem");

  Real hmin = 0.0, hmax = 0.0;
  std::size_t count = 0;

  // Now we will loop over all the elements in the mesh that live
  // on the local processor, and compute the element size.
  const MeshBase& mesh                = this->get_mesh();
  MeshBase::const_element_iterator el =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently working on.
    const Elem *elem = *el;
    Real hmin_t      = elem->hmin();
    Real hmax_t      = elem->hmax();

    if (count == 0)
    {
      hmin = hmin_t;
      hmax = hmax_t;
    }
    else
    {
      hmin = std::min(hmin, hmin_t);
      hmax = std::min(hmax, hmax_t);
    }
    count++;
  } // end for el-loop

  this->comm().min(hmin);
  this->comm().max(hmax);
  std::vector<Real> min_max(2);
  min_max[0] = hmin;
  min_max[1] = hmax;

  STOP_LOG("mesh_size()", "ElasticitySystem");
  return min_max;
}
} // end of namespace libMesh
