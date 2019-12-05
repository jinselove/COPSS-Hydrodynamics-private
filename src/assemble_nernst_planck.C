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
  _boundary_sides_dirichlet_np.resize(_mesh.n_elem());
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

  (void) option;

  // Get a Reference to the LinearImplicitSystem we are solving
  PMSystemNP& np_system = _eqn_sys.get_system<PMSystemNP> (system_name);

  // Get a const reference to the Finite element type for the first (and only
  // for NP system) variable in the system
  FEType fe_type = np_system.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe (FEBase::build(_dim, fe_type));
  UniquePtr<FEBase> fe_face (FEBase::build(_dim, fe_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  QGauss qrule (_dim, fe_type.default_quadrature_order());
  QGauss qface (_dim-1, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);
  fe_face->attach_quadrature_rule(&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.  We will start
  // with the element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW      = fe->get_JxW();
  const std::vector<Real> & JxW_face = fe_face->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();
  const std::vector<std::vector<Real> > & psi = fe_face->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // The XY locations of the quadrature points used for face integration
  const std::vector<Point> & qface_points = fe_face->get_xyz();

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  // in future examples.
  const DofMap & dof_map = np_system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  // we will not update Fe in assemble_matrix() function, but we need it when
  // we call apply_bc_by_penalty() method to modify Ke for BC
//  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.
  MeshBase::const_element_iterator el = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh
    .active_local_elements_end();

  for(; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem * elem = *el;

    // Fill _boundary_sides_dirichlet_np
//    this->select_boundary_side(elem, system_name);

    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices (elem, dof_indices);

    // Compute the element-specific data for the current
    // element.  This involves computing the location of the
    // quadrature points (q_point) and the shape functions
    // (phi, dphi) for the current element.
    fe->reinit (elem);

    // Zero the element matrix and right-hand side before
    // summing them.  We use the resize member here because
    // the number of degrees of freedom might have changed from
    // the last element.  Note that this will be the case if the
    // element type is different (i.e. the last element was a
    // triangle, now we are on a quadrilateral).
    Ke.resize (dof_indices.size(),
               dof_indices.size());
//    Fe.resize (dof_indices.size());

    // Now we will build the element matrix.
    // Constructing the RHS requires the solution and its
    // gradient from the previous timestep.  This myst be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    // loop over Gauss quadrature points
    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {
      // loop over rows of Ke
      for (unsigned int i=0; i<phi.size(); i++)
      {
        // loop over cols of Ke
        for (unsigned int j=0; j<phi.size(); j++)
        {
          Ke(i,j) += JxW[qp] * (
                                  // Mass term
                                  phi[i][qp]*phi[j][qp]
                                  // Diffusion term when using semi implicit
                                  // Euler for Diffusion
                                  + 0.5 * np_system.dt * np_system
                                  .ion_diffusivity * dphi[i][qp] *
                                  dphi[j][qp]
                               );
        } // end loop over cols of Ke
      } // end loop over rows of Ke
    } // end loop for qp

    // apply bc by penalty on K matrix
    {
      // The penalty value.
      const Real penalty = 1.e10;

      // The following loops over the sides of the element.
      // If the element has no neighbor on a side then that
      // side MUST live on a boundary of the domain.
      for (unsigned int s=0; s<elem->n_sides(); s++)
        if (elem->neighbor_ptr(s) == libmesh_nullptr)
        {
          fe_face->reinit(elem, s);

          for (unsigned int qp=0; qp<qface.n_points(); qp++)
          {
            // Matrix contribution
            for (unsigned int i=0; i<psi.size(); i++)
              for (unsigned int j=0; j<psi.size(); j++)
                Ke(i,j) += penalty*JxW_face[qp]*psi[i][qp]*psi[j][qp];
          }
        }
    }
    // If this assembly program were to be used on an adaptive mesh,
    // we would have to apply any hanging node constraint equations.
    dof_map.constrain_element_matrix(Ke, dof_indices);

    // Add the element matrix to the global system.
    np_system.matrix->add_matrix(Ke, dof_indices);
  } // end for loop over local elements

  STOP_LOG("assemble_global_K()", "AssembleNP");
}

// ==================================================================================
void AssembleNP::assemble_global_F(const std::string& system_name,
                                   const std::string& option)
{
  START_LOG("assemble_global_F()", "AssembleNP");

  // Get a Reference to the LinearImplicitSystem we are solving
  PMSystemNP& np_system = _eqn_sys.get_system<PMSystemNP> (system_name);

  // Get a const reference to the Finite element type for the first (and only
  // for NP system) variable in the system
  FEType fe_type = np_system.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe (FEBase::build(_dim, fe_type));
  UniquePtr<FEBase> fe_face (FEBase::build(_dim, fe_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  QGauss qrule (_dim, fe_type.default_quadrature_order());
  QGauss qface (_dim-1, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);
  fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.  We will start
  // with the element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW      = fe->get_JxW();
  const std::vector<Real> & JxW_face = fe_face->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();
  const std::vector<std::vector<Real> > & psi = fe_face->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // The XY locations of the quadrature points used for face integration
  const std::vector<Point> & qface_points = fe_face->get_xyz();

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  // in future examples.
  const DofMap & dof_map = np_system.get_dof_map();

  // we will not update Fe in assemble_matrix() function, but we need it when
  // we call apply_bc_by_penalty() method to modify Ke for BC
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.

  MeshBase::const_element_iterator el = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh
    .active_local_elements_end();

  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem * elem = *el;

    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices (elem, dof_indices);

    // Compute the element-specific data for the current
    // element.  This involves computing the location of the
    // quadrature points (q_point) and the shape functions
    // (phi, dphi) for the current element.
    fe->reinit (elem);

    // Zero the element matrix and right-hand side before
    // summing them.  We use the resize member here because
    // the number of degrees of freedom might have changed from
    // the last element.  Note that this will be the case if the
    // element type is different (i.e. the last element was a
    // triangle, now we are on a quadrilateral).
    Fe.resize (dof_indices.size());

    // Now we will build the element matrix and right-hand-side.
    // Constructing the RHS requires the solution and its
    // gradient from the previous timestep.  This myst be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {
      // Values to hold the old system solution & its gradient at this qp point
      Number c_old = 0.;
      Gradient gradient_c_old;

      if (option == "diffusion")
      {
        // compute the old solution & its gradient at this qp point
        for (unsigned int l = 0; l < phi.size(); l++)
        {
          // get old solution at this qp point
          // notice that we are not using Libmesh::TransientSystem, thus the
          // current_solution of our system is actually solution from previous
          // time step
          c_old += phi[l][qp] * (*np_system.current_local_solution)(dof_indices[l]);
//          c_old += phi[l][qp] * 1.;

          // get solution gradient at this qp_point
          gradient_c_old.add_scaled(dphi[l][qp], (*np_system.current_local_solution)(dof_indices[l]));
//          gradient_c_old.add_scaled(dphi[l][qp], 1.);
        }

        // Now compute the element rhs
        for (unsigned int i = 0; i < phi.size(); i++)
        {
//          std::cout << "JxW[qp]="<<JxW[qp]<<"; c_old="<<c_old
//            <<"; phi[i][qp]="<< phi[i][qp] << "; np_system.dt" << np_system.dt
//            <<"; np_system.ion_diffusivity="<<np_system.ion_diffusivity
//            <<"; gradient_c_old="<<gradient_c_old
//            <<"; dphi[i][qp]="<<dphi[i][qp]<<std::endl;
          Fe(i) += JxW[qp] * (
            // Mass term
            c_old * phi[i][qp]
            // diffusion term when using semi-implicit Euler for diffusion
            - 0.5 * np_system.dt * np_system.ion_diffusivity * gradient_c_old
            * dphi[i][qp]);
//          std::cout<<"Fe(i)="<<Fe(i)<<std::endl;
        } // end loop i < phi.size()
      }
      else if (option == "diffusion&convection")
      {
        // fixme
      }
      else if (option == "diffusion&electrostatics")
      {
        // fixme
      }
      else if (option == "diffusion&electrostatics&convection")
      {
        // fixme
      }
      else
      {
        std::cout << "Error: unsupported option to assembel global F for NP "
                     "system. Exiting..." << std::endl;
        libmesh_error();
      }
    } // end loop qp<q_rule.n_points()

    // apply bc by penalty
    {
      // The penalty value.
      const Real penalty = 1.e10;

      // The following loops over the sides of the element.
      // If the element has no neighbor on a side then that
      // side MUST live on a boundary of the domain.
      for (unsigned int s=0; s<elem->n_sides(); s++)
      {
        if (elem->neighbor_ptr(s) == libmesh_nullptr)
        {
          fe_face->reinit(elem, s);

          for (unsigned int qp=0; qp<qface.n_points(); qp++)
          {
            const Real value =
              analytical_solution->exact_solution_infinite_domain
                (qface_points[qp],
                 _eqn_sys.parameters.get<Real>("real_time"),
                 _eqn_sys.parameters.get<std::vector<Real>>("ion_diffusivity")[0]);

            // RHS contribution
            for (unsigned int i=0; i<psi.size(); i++)
              Fe(i) += penalty*JxW_face[qp]*value*psi[i][qp];
          } // end loop qp<qface.n_points()
        } // end if elem->neighbor_ptr(s) == libmesh_nullptr
      } // end loop s<elem->n_sides()
    }// end apply penalty
//      this->apply_bc_by_penalty(elem, "vector", Ke, Fe, option);

    // If this assembly program were to be used on an adaptive mesh,
    // we would have to apply any hanging node constraint equations.
    // perf_log.push("finish");
    dof_map.constrain_element_vector(Fe, dof_indices);

    // Add the element matrix and rhs vector to the global system.
    // PMToolBox::zero_filter_dense_vector(Fe, 1e-10);
    // PMToolBox::output_dense_vector(Fe);
//    std::cout<<"----> elem id="<<elem->id()<<std::endl;
//    for (int ii=0; ii<dof_indices.size(); ii++)
//    {
//      std::cout<<"dof = " << dof_indices[ii] << "; Fe = "<< Fe(ii) <<std::endl;
//    }
    np_system.rhs->add_vector(Fe, dof_indices);
  } // end loop over local elements


  STOP_LOG("assemble_global_F()", "AssembleNP");
}

// ==================================================================================
void AssembleNP::select_boundary_side(const Elem *elem,
                                      const std::string& system_name)
{
  START_LOG("select_boundary_side()", "AssembleNP");
  // Get a reference to the Particle-Mesh System.
  PMSystemNP & pm_system = _eqn_sys.get_system<PMSystemNP>(system_name);
  const std::size_t elem_id   = elem->id();
  _boundary_sides_dirichlet_np[elem_id].resize(0);

  // The following loops over the sides of the element. If the element has NO
  // neighbors on a side then that side MUST live on a boundary of the domain.
//  std::cout<<"elem_id="<<elem_id<<"; n_sides="<<elem->n_sides()<<std::endl;
  for (unsigned int s = 0; s < elem->n_sides(); s++)
  {
    // If this side is on the boundary
    if (elem->neighbor(s) == NULL)
    {
//      std::cout << "--> s="<<s<<"; has no neighbor"<<std::endl;
      // Check if the s-th face of this element is associated with any Dirichlet boundary id
      // if any, stop looking
      for (unsigned int i = 0; i < pm_system.boundary_id_dirichlet_np.size(); i++)
      {
        if (_mesh.get_boundary_info().has_boundary_id(elem, s,
                                                      pm_system.boundary_id_dirichlet_np[i]))
        {
          _boundary_sides_dirichlet_np[elem_id].push_back(
            std::make_tuple(s,
                            pm_system.boundary_id_dirichlet_np[i],
                            pm_system.boundary_value_dirichlet_np[i]));
          break;
        }
      }
    }
  }

  STOP_LOG("select_boundary_side()", "AssembleNP");
}

// =======================================================================================
void AssembleNP::apply_bc_by_penalty(const Elem          *elem,
                                     const std::string  & matrix_or_vector,
                                     DenseMatrix<Number>& Ke,
                                     DenseVector<Number>& Fe,
                                     const std::string  & option)
{
  START_LOG("apply_bc_by_penalty()", "AssembleNP");

  // we don't need option for NP system
  (void) option;

  // set penalty value to be a large number
  const Real penalty = 1.e10;

  // get the element id
  const std::size_t elem_id = elem->id();

  // Loop through sides in this element that sits on system's boundary
  for (unsigned int s=0; s<_boundary_sides_dirichlet_np[elem_id].size(); s++)
  {
    // extract the tuple for this side of this element, use referece to avoid data copying
    // get<0>(t) : side id
    // get<1>(t) : the dirichlet boundary id associated with this side
    // get<2>(t) : the dirichlet value associated with this side
    const auto &t = _boundary_sides_dirichlet_np[elem_id][s];
    // get the Dirichlet Boundary Value on this side
    Real boundary_value = std::get<2>(t);

    // Build the full-order side element for Dirichlet BC at the walls.
    UniquePtr<Elem> side(elem->build_side(std::get<0>(t)));

    // loop over all nodes on this side element
    for (unsigned int nn=0; nn<side->n_nodes(); nn++)
    {
      // update boundary value for this node only for our validation system
      const Point& ptx = side->point(nn);
      if (_eqn_sys.parameters.get<std::string>("simulation_name") ==
        "np_validation_analytic")
      {
        boundary_value = analytical_solution->exact_solution_infinite_domain
          (ptx, _eqn_sys.parameters.get<Real>("real_time"), _eqn_sys.parameters
          .get<std::vector<Real>>("ion_diffusivity")[0]);
      }

      // Find the node on the element matching this node on the side.
      // That defined where in the element matrix the BC will be applied.
      const unsigned int& local_node_id = elem->local_node(side->node_id(nn));
      // Penalize phi at the current node, var_number is 0 for 'phi'
      std::cout<<"elem_id="<<elem_id
        <<"; id of sides on dirichlet boundaries=" << std::get<0>(t)
        <<"; dirichlet boundary id associated with this side="<<std::get<1>(t)
        <<"; boundary_value="<<boundary_value << std::endl;
      this->penalize_elem_matrix_vector(Ke,
                                        Fe,
                                        matrix_or_vector,
                                        0,
                                        local_node_id,
                                        elem->n_nodes(),
                                        penalty,
                                        boundary_value);
    }
  }


  STOP_LOG("apply_bc_by_penalty()", "AssembleNP");
} // end of function apply_bc_by_penalty()
