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
#include "ggem_poisson.h"
#include "pm_toolbox.h"
#include "assemble_poisson.h"
#include "pm_system_poisson.h"

// ==================================================================================
AssemblePoisson::AssemblePoisson(EquationSystems  & es,
                                 const std::string& name)
  : AssembleSystem(es)
{
  libmesh_assert_equal_to(name, "Poisson");
  analytical_solution = new AnalyticalSolutionPoisson(name);
  ggem_poisson        = new GGEMPoisson();
  // Initialize boundary sides for Dirichlet and Neumann BC
  _boundary_sides_dirichlet_poisson.resize(_mesh.n_elem());
  _boundary_sides_neumann_poisson.resize(_mesh.n_elem());
}

// ==================================================================================
AssemblePoisson::~AssemblePoisson()
{
  delete analytical_solution; analytical_solution = nullptr;
  delete ggem_poisson; ggem_poisson = nullptr;
}


// ==================================================================================
void AssemblePoisson::assemble_global_K(const std::string& system_name,
                                        const std::string& option)
{
  START_LOG("assemble_global_K()", "AssemblePoisson");

 // Get a reference to the LinearImplicitSystem we are solving
  PMSystemPoisson& _pm_system = _eqn_sys.get_system<PMSystemPoisson>(system_name);
  
  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers. 
  const DofMap& dof_map = _pm_system.get_dof_map();

  // Numeric id that corresponds to variable (potential) in the system
  const unsigned int phi_var = _pm_system.variable_number("phi");

  // Get Finite Element type for "phi"
  FEType fe_phi_type = _pm_system.variable_type(phi_var);

  // Build a Finite Element object of the specified type for the potential
  // variable
  UniquePtr<FEBase> fe_phi(FEBase::build(_dim, fe_phi_type));

  // Gauss quadrature rule for numerical integration
  QGauss qrule(_dim, fe_phi_type.default_quadrature_order());

  // Tell finite element objects to use the quadrature rule
  fe_phi->attach_quadrature_rule(&qrule);

  // Element Jacobian*quadrature weight at each integration point
  const std::vector<Real>& JxW = fe_phi->get_JxW();

  // Element shape function gradients for potential variable evaluated at
  // quadrature points
  const std::vector<std::vector<RealGradient> >& dphi = fe_phi->get_dphi();

  // Element shape function for potential variable evaluated at quadrature
  // points
  const std::vector<std::vector<Real> >& phi = fe_phi->get_phi();

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Element matrix contribution
  DenseMatrix<Number> Ke;

  // Element vector contribution
  DenseVector<Number> Fe;

  // attach PointMesh to AnalyticalSolution (only do this once for
  // "ggem_validation")
  if ((_eqn_sys.parameters.get<std::string>("simulation_name") ==
       "ggem_validation_poisson") &&
      !analytical_solution->get_point_mesh())
  {
    analytical_solution->attach_point_mesh(_pm_system.point_mesh());
  }

  // Initialize (set) parameters in ggem_poisson
  this->init_ggem_poisson(system_name);

  // Loop over all the elements in the mesh that live on the local processor.
  // We will compute the element matrix Ke.
  MeshBase::const_element_iterator el = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
          _mesh.active_local_elements_end();
  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently working on.
    const Elem *elem = *el;
    
    // Fill _boundary_sides_neumann_poisson and _boundary_sides_dirichlet_poisson 
    this->select_boundary_side(elem, system_name);
    
    // Get the degree of freedom indices for the current element. These define
    // where in the global matrix and right-hand-side this element will contribute to.
    dof_map.dof_indices(elem, dof_indices);

    // Cache the number of degrees of freedom on this element, for
    // use as a loop bound later.  We use cast_int to explicitly
    // convert from size() (which may be 64-bit) to unsigned int
    // (which may be 32-bit but which is definitely enough to count
    // *local* degrees of freedom.
    const unsigned int n_dofs = cast_int<unsigned int> (dof_indices.size());

    // Compute the element-specific data for the current element. This involves
    // computing the location of the quadrature points (q_point) and the 
    // shape functions (phi, dphi) for the current element.
    fe_phi->reinit(elem);

    // Zero the element matrix before summing them.  We use the resize member here because
    // the number of degrees of freedom might have changed from
    // the last element.  Note that this will be the case if the
    // element type is different (i.e. the last element was a
    // triangle, now we are on a quadrilateral).
    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);

    // Now loop over the quadrature points. This handles the numeric
    // integration.
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      // Add the matrix contribution
      for (std::size_t i = 0; i < phi.size(); i++)
        for (std::size_t j = 0; j < phi.size(); j++)
        { 
          Ke(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
        }
    }

    // apply BCs by penalty method
    this->apply_bc_by_penalty(elem, "matrix", Ke, Fe, option);

    // If this assembly program were to be used on an adaptive mesh,
    // we would have to apply any hanging node constraint equations.
    dof_map.constrain_element_matrix(Ke, dof_indices);

    // Add the element matrix to the global system.
    // GeomTools::zero_filter_dense_matrix(Ke, 1e-10);
    // PMToolBox::output_dense_matrix(Ke);
    _pm_system.matrix->add_matrix(Ke, dof_indices);
  }

  STOP_LOG("assemble_global_K()", "AssemblePoisson");
  return;
}

// ==================================================================================
void AssemblePoisson::assemble_global_F(const std::string& system_name,
                                        const std::string& option)
{
  START_LOG("assemble_global_F()", "AssemblePoisson");

  // Make sure we are assembling the proper system.
  libmesh_assert_equal_to(system_name, "Poisson");

  // Get a reference to the LinearImplicitSystem we are solving
  PMSystemPoisson& _pm_system = _eqn_sys.get_system<PMSystemPoisson>(system_name);

  if (_eqn_sys.parameters.get<bool>("module_np")){
    // if np_systems vectors are not filled with NP systems pointers, we do it
    if (np_systems.size()==0)
    {
      // loop over all systems and push the references of all NP systems and
      // their dof maps to np_systems and np_dof_maps vector respectively
      unsigned int n_sys = _eqn_sys.n_systems();
      for (unsigned int s_id=0; s_id<n_sys; s_id++){
        const std::string& s_name = _eqn_sys.get_system(s_id).name();
        if (s_name.rfind("NP:", 0)==0){
          // np system
          np_systems.push_back(&(_eqn_sys.get_system<PMSystemNP>(s_name)));
          // np system dof map
          np_dof_maps.push_back(&(_eqn_sys.get_system<PMSystemNP>(s_name)
            .get_dof_map()));
        }
      }
    }
  }

  // A reference to the DofMap object for this system.
  const DofMap& dof_map = _pm_system.get_dof_map();

  // Numeric id corresponding to potential variable
  const unsigned int phi_var = _pm_system.variable_number("phi");

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a std::unique_ptr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.  
  FEType fe_phi_type = _pm_system.variable_type(phi_var);
  UniquePtr<FEBase> fe_phi(FEBase::build(_dim, fe_phi_type));

  // Define Gauss quadrature rule for numerical integration.
  QGauss qrule(_dim, fe_phi_type.default_quadrature_order());
  fe_phi->attach_quadrature_rule(&qrule);

  // Declare a special finite element object for
  // boundary integration.
  UniquePtr<FEBase> fe_face(FEBase::build(_dim, fe_phi_type));
  // Boundary integration requires one quadrature rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  QGauss qface(_dim - 1, fe_phi_type.default_quadrature_order());
  fe_face->attach_quadrature_rule(&qface);

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>& JxW               = fe_phi->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe_phi->get_phi();
  const std::vector<Point>& q_xyz            = fe_phi->get_xyz();
  const std::vector<std::vector<RealGradient>>& dphi = fe_phi->get_dphi(); 

  // Define data structures to contain the element matrix Ke and vector Fe
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;
  DenseVector<Number> sigma_e;
  DenseVector<Number> sigma_e_global;
  DenseVector<Number> sigma_e_local;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
 
 // Now we will loop over all the elements in the mesh that live
  // on the local processor, and compute the element vector Fe.
  MeshBase::const_element_iterator el =
    _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    _mesh.active_local_elements_end();
  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently working on.
    const Elem *elem           = *el;
    const unsigned int elem_id = elem->id();
    // std::cout<<"-------------------elem_id = "<<elem_id <<"------------------\n";
    // Get the degree of freedom indices for the current element. These define
    // where in the global matrix and right-hand-side this element will contribute to.
    dof_map.dof_indices(elem, dof_indices);
    
    // Cache the number of degrees of freedom on this element, for
    // use as a loop bound later.  We use cast_int to explicitly
    // convert from size() (which may be 64-bit) to unsigned int
    // (which may be 32-bit but which is definitely enough to count
    // *local* degrees of freedom.
    const unsigned int n_dofs = cast_int<unsigned int> (dof_indices.size());
    
    // Compute the element-specific data for the current element. This involves
    // computing the location of the quadrature points (q_point) and the 
    // shape functions (phi, dphi) for the current element.
    fe_phi->reinit(elem);
   
    // update Fe size
    Fe.resize(n_dofs);
    sigma_e.resize(n_dofs);
    sigma_e_local.resize(n_dofs);
    sigma_e_global.resize(n_dofs);
    

    // if elem_neighbor_list is pre-built, we can access it directly
    const std::vector<dof_id_type>& n_list =
      _pm_system.point_mesh()->get_elem_point_neighbor_list(elem->id());

    // if this elem has no neighboring particle we turn the pc_flag to 'false'
    bool pc_flag = true;                      // a flag for the point charge
    if (n_list.size() == 0) pc_flag = false;  // No point force because no point
                                              // list.

    // Now compute Fe caused by the regularized point charge.
    // FIXME: need to implement space charge density when considering
    // Nernst-Planck solver
    this->compute_element_rhs(elem, n_dofs, JxW, phi, q_xyz, n_list, pc_flag,
                              option, Fe);

    // Impose the Dirichlet BC (electrical potential) via the penalty method.
    this->apply_bc_by_penalty(elem, "vector", Ke, Fe, option);
    // PMToolBox::output_dense_vector(Fe);

    // Impose Neumann BC (surface charge density) to the right hand side
    this->apply_bc_neumann(elem, *fe_phi, *fe_face, Fe, sigma_e, sigma_e_global, sigma_e_local);
    // PMToolBox::output_dense_vector(Fe);

    // std::cout << "sigma_e:\n";
    // PMToolBox::output_dense_vector(sigma_e);
    // std::cout << "sigma_e_global:\n";
    // PMToolBox::output_dense_vector(sigma_e_global);
    // std::cout << "sigma_e_local:\n";
    // PMToolBox::output_dense_vector(sigma_e_local);
    // If this assembly program were to be used on an adaptive mesh,
    // we would have to apply any hanging node constraint equations.
    dof_map.constrain_element_vector(Fe, dof_indices);

    // PMToolBox::output_dense_vector(Fe);

    // _n_dofs[elem_id] = _dof_indices[elem_id].size();
   // Add the element rhs vector to the global system.
    // PMToolBox::zero_filter_dense_vector(Fe, 1e-10);
    // PMToolBox::output_dense_vector(Fe);
    _pm_system.rhs->add_vector(Fe, dof_indices);
  } // end for elem-loop
  STOP_LOG("assemble_global_F()", "AssemblePoisson");
}


// ==================================================================================
void AssemblePoisson::compute_element_rhs(const Elem *elem,
                                          const unsigned int& n_dofs,
                                          const std::vector<Real>& JxW,
                                          const std::vector<std::vector<Real>>& phi,
                                          const std::vector<Point>& q_xyz,
                                          const std::vector<dof_id_type>&n_list,
                                          const bool& pc_flag,
                                          const std::string& option,
                                          DenseVector<Number>& Fe)
{
  START_LOG("compute_element_rhs()", "AssemblePoisson");

  // we don't need option for Poisson system
  (void) option;

  PMSystemPoisson& _pm_system = _eqn_sys.get_system<PMSystemPoisson>("Poisson");

  // Get a reference to the PointMesh, PMPeriodicBoundary, and PointParticles
  PointMesh<3> *_point_mesh                 = _pm_system.point_mesh();
  PMPeriodicBoundary *_pm_periodic_boundary = _point_mesh->pm_periodic_boundary();
  std::vector<PointParticle *> _particles   = _point_mesh->particles();

  // 1. Add the regularized point charge, first examine if this element has
  // neighboring point charge sources
  if (pc_flag)
  {
    // Number of particles in the neighbor list
    const std::size_t n_pts = n_list.size();

    // Initialize variables for position, charge, distance, etc.
    Point np_pos(0.);
    Real charge_val = 0., np_charge = 0., pi_4 = 4. * libMesh::pi;
    Point r;
    unsigned int qp_size = q_xyz.size();

    // Assemble int_force to avoid duplicate calculations
    // Fixme: avoid calculating int_force at every time step
    std::vector<Real> int_force(n_dofs * qp_size, 0.);
    for (unsigned int l=0; l<n_dofs; l++)
    {
      for (unsigned int qp=0; qp<qp_size; qp++)
      {
        int_force[l * qp_size + qp] = JxW[qp] * phi[l][qp];
      }
    }

    // get dof indices for this element in the NP systems. Assuming all NP
    // systems have the same dof_indices for this element
    std::vector<dof_id_type> np_dof_indices;
    if (np_dof_maps.size()>0)
      np_dof_maps[0]->dof_indices(elem, np_dof_indices);

    // loop over all qp points to compute rhs
    for (unsigned int qp=0; qp<qp_size; qp++)
    {
      // calculate the global charge density at this qp
      Real rho_global = 0.;

      // loop over all neighbor points of this qp point to accumulate rho_global
      for (unsigned int np=0; np<n_pts; np++)
      {
        // get a constant pointer to this neighbor particle
        const Point& pt = _particles[n_list[np]]->point();
        const Real& charge = _particles[n_list[np]]->charge();
        // calculate the contribution of this neighbor particle to global
        // charge density of this qp point
        rho_global += ggem_poisson->smoothed_charge_exp
          (_pm_periodic_boundary->point_vector(q_xyz[qp], pt)) * charge;
      }

      // calculate the contribution of ion cloud to rho_global if NP systems
      // exists. If there are no np_systems, np_systems.size() is 0, so this
      // loop will be skipped
      Real rho_global_ion = 0.;
      for (short int s_id=0; s_id<np_systems.size(); s_id++)
      {
        // interpolate ion concentration at this qp point using the ion
        // concentration at all nodes of this element
        Real tmp = 0.;
        for (unsigned int l=0; l<n_dofs; l++)
        {
          tmp += (phi[l][qp] *
            np_systems[s_id]->current_local_solution->operator()
            (np_dof_indices[l]));
        }
        // scale tmp by the coefficients of this np system
        tmp *= (_eqn_sys.parameters.get<Real>("coeff_ion_charge_density") *
          np_systems[s_id]->ion_valence);

        // add the contribution of this np system to global charge density
        rho_global_ion += tmp;
      }
//      std::cout<<"rho_global_ion = "<<rho_global_ion<<std::endl;
      rho_global += rho_global_ion;

      // scale rho_global by 4*pi, this comes from the dimensionless process
      rho_global *= pi_4;

      // add the contribution of this qp point to the rhs vector
      for (unsigned int l=0; l<n_dofs; l++)
      {
        Fe(l) += int_force[l * qp_size + qp] * rho_global;
      }
    }
  }       // end if( pc_flag )

  STOP_LOG("compute_element_rhs()", "AssemblePoisson");
}


// ==================================================================================
void AssemblePoisson::select_boundary_side(const Elem *elem,
                                           const std::string& system_name)
{
  START_LOG("select_boundary_side()", "AssemblePoisson");

  // Get a reference to the Particle-Mesh System.
  PMSystemPoisson & pm_system = _eqn_sys.get_system<PMSystemPoisson>(system_name);
  const std::size_t elem_id   = elem->id();
  _boundary_sides_dirichlet_poisson[elem_id].resize(0);
  _boundary_sides_neumann_poisson[elem_id].resize(0);

  // Get boundary ids for each BC
  const std::vector<unsigned int> boundary_id_dirichlet_poisson =
    _eqn_sys.parameters.get<std::vector<unsigned int> >(
      "boundary_id_dirichlet_poisson");
  const std::vector<unsigned int> boundary_id_neumann_poisson =
    _eqn_sys.parameters.get<std::vector<unsigned int> >(
      "boundary_id_neumann_poisson");
  const std::vector<Real> boundary_value_dirichlet_poisson =
    _eqn_sys.parameters.get<std::vector<Real> >(
      "boundary_value_dirichlet_poisson");
  const std::vector<Real> boundary_value_neumann_poisson =
    _eqn_sys.parameters.get<std::vector<Real> >(
      "boundary_value_neumann_poisson");

  // The following loops over the sides of the element. If the element has NO
  // neighbors on a side then that side MUST live on a boundary of the domain.
  for (unsigned int s = 0; s < elem->n_sides(); s++)
  {
    // bool apply_bc_dirichlet = false, apply_bc_neumann = false;

    // If this side is on the boundary
    if (elem->neighbor(s) == NULL)
    {
      // Check if the s-th face of this element is associated with any Dirichlet boundary id
      // if any, stop looking
      for (unsigned int i = 0; i < boundary_id_dirichlet_poisson.size(); i++) 
      {
        if (_mesh.get_boundary_info().has_boundary_id(elem, s,
                                                      boundary_id_dirichlet_poisson[i])) 
        {
          // apply_bc_dirichlet = true;
          _boundary_sides_dirichlet_poisson[elem_id].push_back(
            std::make_tuple(s, 
              boundary_id_dirichlet_poisson[i], 
              boundary_value_dirichlet_poisson[i]));
          break;
        }
      }

      // Check if the s-th face of this element is associated with any Neumann boundary id
      // if any, stop looking
      for (unsigned int i = 0; i < boundary_id_neumann_poisson.size(); i++) 
      {
        if (_mesh.get_boundary_info().has_boundary_id(elem, s,
                                                      boundary_id_neumann_poisson[
                                                        i])) {
          // apply_bc_neumann = true;
          _boundary_sides_neumann_poisson[elem_id].push_back(
            std::make_tuple(s, 
              boundary_id_neumann_poisson[i], 
              boundary_value_neumann_poisson[i]));
          break;
        }
      }
    
    }
  }

  STOP_LOG("select_boundary_side()", "AssemblePoisson");
}

// =======================================================================================
void AssemblePoisson::apply_bc_by_penalty(const Elem          *elem,
                                          const std::string  & matrix_or_vector,
                                          DenseMatrix<Number>& Ke,
                                          DenseVector<Number>& Fe,
                                          const std::string  & option)
{
  START_LOG("apply_bc_by_penalty()", "AssemblePoisson");

  // we don't need option for Poisson system
  (void) option;

  const std::size_t  elem_id = elem->id();
  
  // If the dirichlet boundary sides of this elem is none, return
  if (_boundary_sides_dirichlet_poisson[elem_id].size() == 0)
  {
    return;
  }
 
  // Get a reference to the Particle-Mesh System
  PMSystemPoisson& pm_system = _eqn_sys.get_system<PMSystemPoisson>("Poisson");
  
  // The penalty value
  const Real penalty         = 1E6;

  // Loop through sides in this element that sits on system's boundary
  for (unsigned int s = 0; s < _boundary_sides_dirichlet_poisson[elem_id].size(); s++)
  {
    // extract the tuple for this side of this element, use referece to avoid data copying
    // get<0>(t) : side id
    // get<1>(t) : the dirichlet boundary id associated with this side
    // get<2>(t) : the dirichlet potential associated with this side
    const auto &t = _boundary_sides_dirichlet_poisson[elem_id][s];
    // Build the full-order side element for "potential" Dirichlet BC at the walls.
    // UniquePtr<Elem> side(elem->build_side(
                           // _boundary_sides_dirichlet_poisson[elem_id][s]));
    UniquePtr<Elem> side(elem->build_side(std::get<0>(t)));

    // Total electrical potential on this side
    Real phi_total = std::get<2>(t);

    // Loop through nodes on this side element
    for (unsigned int nn = 0; nn < side->n_nodes(); nn++)
    {
      // Calculate local part of the electrical potential from GGEM
      const Point ptx       = side->point(nn); // Coordinate of the node
      const Real  phi_local = pm_system.local_potential_field(ptx,
        "regularized", elem->id());

      // When we run ggem_validation_poisson, Dirichelet boundary condition is
      // applied on all
      // boundaries, the electrical potential on boundaries due to point charges
      // are evaluated
      // using free-space Green's function.
      if (_eqn_sys.parameters.get<std::string>("simulation_name") ==
          "ggem_validation_poisson")
      {
        phi_total = analytical_solution->exact_solution_infinite_domain(
          *ggem_poisson,
          ptx);
      }

      // Because of Ewald split: global_potential = total_potential -
      // ggem_local_potential
      Real phi_global = phi_total - phi_local;
      if ((phi_local>0)||(phi_local<0))
        std::cout<<"vrai"<<"\n";

      // Find the node on the element matching this node on the side.
      // That defined where in the element matrix the BC will be applied.
      const unsigned int& local_node_id = elem->local_node(side->node_id(nn));
      // Penalize phi at the current node, var_number is 0 for 'phi'
      this->penalize_elem_matrix_vector(Ke,
                                        Fe,
                                        matrix_or_vector,
                                        0,
                                        local_node_id,
                                        elem->n_nodes(),
                                        penalty,
                                        phi_global);
    }     // end for nn-loop
  }       // end for s-loop

  STOP_LOG("apply_bc_by_penalty()", "AssemblePoisson");
}

// =======================================================================================
void AssemblePoisson::apply_bc_neumann(const Elem          *elem,
                                       FEBase             & fe_phi,
                                       FEBase             & fe_face,
                                       DenseVector<Number>& Fe,
                                       DenseVector<Number>& sigma_e,
                                       DenseVector<Number>& sigma_e_global,
                                       DenseVector<Number>& sigma_e_local
                                      )
{
  START_LOG("apply_bc_neumann()", "AssemblePoisson");
  
  const std::size_t  elem_id = elem->id();
  
  // If the Neumann boundaries of this element is None, return None
  if (_boundary_sides_neumann_poisson[elem->id()].size() == 0)
  {
    return;
  }
  /*std::cout<<"Warning: needs more testing/validation on the Neumann BC of "
             "Poisson system. Exiting..." <<std::endl;*/
//  libmesh_error();

  // Get a reference to the Particle-Mesh System
  PMSystemPoisson& pm_system = _eqn_sys.get_system<PMSystemPoisson>("Poisson");
  
  // Integration terms for side element
  const std::vector<std::vector<Real>>& phi_face = fe_face.get_phi();
  const std::vector<std::vector<RealGradient>> dphi_face = fe_face.get_dphi();
  const std::vector<Real>& JxW_face               = fe_face.get_JxW();
  const std::vector<Point>& qface_point = fe_face.get_xyz();
  const std::vector<Point>& face_normals          = fe_face.get_normals();

//      std::ofstream outfile;
//  outfile.open("debug.csv", std::ios_base::out);
//  outfile<<"face normal x,face normal y,face normal z,Potential gradient x,Potential gradient y,Potential gradient z,sigma_total,dphi_local/dn,phi_local\n";


  // Extract the shape function derivatives to be evaluated at the nodes
  // const std::vector<std::vector<RealGradient> >& dphi = fe_phi.get_dphi();
  // Extract the element node coordinates in the reference frame
  // std::vector<Point> nodes;
  // fe_phi.get_refspace_nodes(elem->type(), nodes);
  // Evaluate the shape functions derivatives at the nodes
  // fe_phi.reinit(elem, &nodes);

  
//  fe_phi.reinit(elem);
  // Pre-calculate 'Phi_local' and 'dphi_local' on all nodes of this element
  // so that we can directly use them later when we loop over nodes on the 
  // boundary surfaces
//   std::vector<Real> phi_local(elem->n_nodes());
  // std::vector<RealGradient> dphi_local(elem->n_nodes());
  // // Loop through nodes on this element
  // for (unsigned int nn = 0; nn < elem->n_nodes(); nn++)
  // {
  //   // Coordinate of the node
  //   const Point ptx = elem->point(nn);
  //   // Evaluate local electrical potential (phi) on each node
  //   phi_local[nn] = pm_system.local_potential_field(elem, ptx, "regularized");
  //   // phi_local[nn] = pm_system.local_potential_field(ptx, "regularized");
  // }
  // // Gradient of electrical potential on all nodes in this element
  // for (unsigned int nn = 0; nn < elem->n_nodes(); nn++)
  // {
  //   for (unsigned int i = 0; i < phi_local.size(); i++)
  //   {
  //     dphi_local[nn] += phi_local[i] * dphi[i][nn];
  //   }
  // }
  // Loop through sides in this element that sits on system's boundary
  // std::cout << "number of boundary faces of this elem = " << _boundary_sides_neumann_poisson[elem_id].size() <<"\n";
  
  for (unsigned int s = 0; s < _boundary_sides_neumann_poisson[elem_id].size(); s++)
  {
    // Build the full-order side element for "potential" Dirichlet BC at the
    // walls.
    
    // get<0>(t) : side id
    // get<1>(t) : the neumann boundary id associated with this side
    // get<2>(t) : the neumann potential associated with this side
    const auto &t = _boundary_sides_neumann_poisson[elem_id][s];
    const unsigned int& side_id = std::get<0>(t);
    const unsigned int& boundary_id = std::get<1>(t);
    const Real& sigma_total = std::get<2>(t);
    // std::cout << "side id = " << side_id << "\n";
    // std::cout << "boundary id of this side = " << boundary_id << "\n";
    // std::cout << "sigma_total of this side = " << sigma_total <<"\n";

    // reinit fe_face on this side 
    fe_face.reinit(elem, side_id);
    // std::cout << "number of quadrature points of this side = " << qface_point.size() <<"\n";
    // std::cout << "size of JxW of this side = " << JxW_face.size() <<"\n";
    // std::cout << "size of phi of this side (evaludate on quadrature points) = " << phi_face.size() <<"(of "<< phi_face[0].size()<<")"<<"\n";
    // std::cout << "size of dphi of this side (evaludate on quadrature points) = " << dphi_face.size() <<"\n";
    Real val = 0.;
      std::pair<Real, Point> local_potential_old;
    for (unsigned int qp=0; qp<qface_point.size(); qp++)
    {
      // calculate the neumann boundary value at the quadrature point


        pm_system.local_potential_field(qface_point[qp], "regularized",
                                        "grad", local_potential_old, elem_id);

        const Real dphi_local=local_potential_old.second(0)*face_normals[0](0)+
                local_potential_old.second(1)*face_normals[0](1)+
                local_potential_old.second(2)*face_normals[0](2);


        val = sigma_total*(4*libMesh::pi)-dphi_local;

      for (unsigned int i=0; i<phi_face.size(); i++)
      {
        Fe(i) += JxW_face[qp] * val * phi_face[i][qp];
      }
    }

    // // build this side
    // // UniquePtr<Elem> side(elem->build_side(std::get<0>(t)));
    // 
    // // dphi / dn (surface charge density) on this side
    // 
    // // Calcualte Jacobian*Weights and shape functions on this side
    // fe_face.reinit(elem, s);
    // std::cout << "number of nodes of this side = " << side->n_nodes() <<"\n";
    // std::cout << "number of quadrature points of this side = " << q_xyz_face.size() <<"\n";
    // std::cout << "size of JxW of this side = " << JxW_face.size() <<"\n";
    // std::cout << "size of phi of this side (evaludate on quadrature points) = " << phi_face.size() <<"(of "<< phi_face[0].size()<<")"<<"\n";
    // std::cout << "size of dphi of this side (evaludate on quadrature points) = " << dphi_face.size() <<"\n";
    // 
    // // Calculate sigma_local, i.e., (dphi / dn)_local on side nodes
    // std::vector<Real> sigma_local(side->n_nodes(), 0.);
    // std::vector<Real> sigma_global(side->n_nodes(), 0.);
    // // std::cout << "sigma_local (of this face) size = " << sigma_local.size() << "\n";
    // // std::cout << "sigma_global (of this face) size = " << sigma_global.size() << "\n";
    // // Loop through nodes on this side element
    // // std::cout << "local ids of nodes on the side are: ";
    // std::cout << "size of normals of this side = " << face_normals.size() << "\n";
    // std::cout<< "Normals of the face:\n";
    // for (int i=0; i<face_normals.size(); i++) 
    // {
    //   face_normals[i].print();
    //   std::cout<<"\n";
    // }
    
    // for (unsigned int sn = 0; sn < side->n_nodes(); sn++)
    // {
    //   // // Find the node on the element matching this node on the side.
    //   // // That defined where in the element matrix the BC will be applied.
    //   // for (unsigned int nn = 0; nn < elem->n_nodes(); nn++)
    //   // {
    //   //   if (elem->node(nn) == side->node(sn))
    //   //   {
    //   //     // Use the normal vector at a quadrature point, since
    //   //     // normally the side element has zero curvature
    //   //     sigma_local[sn] = dphi_local[nn](0) * face_normals[0](0)
    //   //                       + dphi_local[nn](1) * face_normals[0](1)
    //   //                       + dphi_local[nn](2) * face_normals[0](2);
    //   // 
    //   //   }
    //   // }
    // 
    //   // get the local node id of the elem using the global node id
    //   const unsigned int& nn = elem->local_node(side->node_id(sn));
    //   sigma_local[sn] = dphi_local[nn](0) * face_normals[0](0)
    //                     + dphi_local[nn](1) * face_normals[0](1)
    //                     + dphi_local[nn](2) * face_normals[0](2);
    //   // Calculate (dphi / dn)_global
    //   sigma_global[sn] = sigma_total - sigma_local[sn];
    //   // std::cout << "local_node_id = " << nn << "; sigma_total = " << sigma_total << "; sigma_global = " << sigma_global[sn] << "; sigma_local = "<<sigma_local[sn] << "\n";
    // }
    // std::cout << "\n";

    // Integrate on this side
    // FIXME: this can be optimized by store JxW_face * phi_face for face
    // element
    // in a vector similar to that in _int_force
    // for (unsigned int i=0; i<phi_face.size(); i++)
    // {
    //   for (unsigned int qp=0; qp<JxW_face.size(); qp++)
    //   {
    //     Fe(i) += JxW_face[qp] * phi_face[i][qp] * sigma_global[i];
    //   }
    //   // sigma_e(i) = sigma_total;
    //   // sigma_e_global(i) = sigma_global[i];
    //   // sigma_e_local(i) = sigma_local[i];
    // }
  } // end for s-loop
//    outfile.close();
  STOP_LOG("apply_bc_neumann()", "AssemblePoisson");
}

// ============================================================================================
void AssemblePoisson::init_ggem_poisson(const std::string& system_name)
{
  START_LOG("init_ggem_poisson()", "AssemblePoisson");

  //libmesh_assert_equal_to(system, "Poisson");

  PMSystemPoisson& _pm_system =
    _eqn_sys.get_system<PMSystemPoisson>(system_name);
  PointMesh<3>   *_point_mesh = _pm_system.point_mesh();
  const PointType point_type  = _point_mesh->particles()[0]->point_type();
  ggem_poisson->set_alpha(_eqn_sys.parameters.get<Real>("alpha"));
  ggem_poisson->set_br0(_eqn_sys.parameters.get<Real>("br0"));
  ggem_poisson->set_point_type(point_type);
  ggem_poisson->set_ksi();
  ggem_poisson->set_phi0(_eqn_sys.parameters.get<Real>("phi0"));

  STOP_LOG("init_ggem_poisson()", "AssemblePoisson");
}
