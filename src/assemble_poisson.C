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
#include "ggem_poisson.h"
#include "pm_toolbox.h"
#include "assemble_poisson.h"
#include "pm_system_poisson.h"

// ==================================================================================
AssemblePoisson::AssemblePoisson(EquationSystems& es,
                                 const std::string& name)
: AssembleSystem(es)
{
  if (name != "Poisson") libmesh_error();
  analytical_solution = new AnalyticalSolutionPoisson(name);
  ggem_poisson = new GGEMPoisson();
}



// ==================================================================================
AssemblePoisson::~AssemblePoisson()
{
  if (ggem_poisson){
    delete ggem_poisson;
  }
}



// ==================================================================================
void AssemblePoisson::assemble_global_K(const std::string& system_name,
                                        const std::string& option)
{
  START_LOG ("assemble_global_K()", "AssemblePoisson");

  // Make sure we are assembling the proper system
  libmesh_assert_equal_to (system_name, "Poisson");

  const unsigned int n_mesh_elem = _mesh.n_elem();

  PMSystemPoisson& _pm_system = _eqn_sys.get_system<PMSystemPoisson>(system_name);

  // Numeric id that corresponds to variable (potential) in the system
  const unsigned int phi_var = _pm_system.variable_number ("phi");

  // Get Finite Element type for "phi"
  FEType fe_phi_type = _pm_system.variable_type(phi_var);

  // Build a Finite Element object of the specified type for the potential variable
  UniquePtr<FEBase> fe_phi (FEBase::build(_dim, fe_phi_type) );

  // Gauss quadrature rule for numerical integration
  QGauss qrule (_dim, SECOND);

  // Tell finite element objects to use the quadrature rule
  fe_phi->attach_quadrature_rule (&qrule);

  // Element Jacobian*quadrature weight at each integration point
  const std::vector<Real>& JxW = fe_phi->get_JxW();

  // Element shape function gradients for potential variable evaluated at quadrature points
  const std::vector<std::vector<RealGradient>>& dphi = fe_phi->get_dphi();

  // Element shape function for potential variable evaluated at quadrature points
  const std::vector<std::vector<Real>>& phi = fe_phi->get_phi();

  // Reference to the DofMap object for this system
  const DofMap & dof_map = _pm_system.get_dof_map();
  std::vector<dof_id_type> dof_indices;

  // Element matrix contribution
  DenseMatrix<Number> Ke;

  // Element vector contribution
  DenseVector<Number> Fe;

  // Build _boundary_sides_dirichlet_poisson and
  //       _boundary_sides_neumann_poisson vectors at beginning of simulation
  if(_boundary_sides_dirichlet_poisson.size() == 0 and _boundary_sides_neumann_poisson.size() == 0 ){
    _boundary_sides_dirichlet_poisson.resize(n_mesh_elem);
    _boundary_sides_neumann_poisson.resize(n_mesh_elem);
    MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = _mesh.active_local_elements_end();
    for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently working on.
      const Elem* elem = *el;
      this->select_boundary_side(elem);
      //printf("finished select_boundary_side\n");
    }
  }

  // attach PointMesh to AnalyticalSolution (only do this once for "ggem_validation")
  if (_eqn_sys.parameters.get<std::string> ("simulation_name") == "ggem_validation_poisson" &&
      !analytical_solution -> get_point_mesh())
  {
      analytical_solution -> attach_point_mesh(_pm_system.point_mesh());
  }

  // Initialize (set) parameters in ggem_poisson
  this->init_ggem_poisson(system_name);

  // Loop over all the elements in the mesh that live on the local processor.
  // We will compute the element matrix Ke.
  MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_local_elements_end();
  for ( ; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently working on.
    const Elem* elem = *el;

    // Get the degree of freedom indices for the current element. These define where
    // in the global matrix and right-hand-side this element will contribute to.
    dof_map.dof_indices (elem, dof_indices);

    const unsigned int n_dofs = dof_indices.size();

    // Zero the element matrix Ke before summing them. We use the resize member here because
    // the number of degrees of freedom might have changed from the last element. Note that
    // this will be the case if the element type is different (i.e. the last element was a
    // triangle, now we are on a quadrilateral.)
    Ke.resize (n_dofs, n_dofs);

    // Compute the element-specific data for the current element. This involves computing
    // the location of the quadrature points (q_point) and the shape functions (phi, dphi)
    // for the current element.
    fe_phi->reinit (elem);

    // Now loop over the quadrature points. This handles the numeric integration.
    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {
      // Add the matrix contribution
      for (std::size_t i=0; i<phi.size(); i++)
        for (std::size_t j=0; j<phi.size(); j++)
          Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
    }

    // apply BCs by penalty method
    this->apply_bc_by_penalty(elem, "matrix", Ke, Fe, option);

    // If this assembly program were to be used on an adaptive mesh,
    // we would have to apply any hanging node constraint equations.
    dof_map.constrain_element_matrix (Ke, dof_indices);

    // Add the element matrix to the global system.
    //GeomTools::zero_filter_dense_matrix(Ke, 1e-10);
    //PMToolBox::output_dense_matrix(Ke);
    _pm_system.matrix->add_matrix (Ke, dof_indices);
  }
  //if (_pm_system.comm().rank()==0){
  // printf("assemble_matrix_K(): The global matrix K has been assembled ...\n");
  //}

  return;
  STOP_LOG ("assemble_global_K()", "AssemblePoisson");
}



// ==================================================================================
void AssemblePoisson::assemble_global_F(const std::string& system_name,
                                        const std::string& option)
{
  START_LOG ("assemble_global_F()", "AssemblePoisson");

  // Make sure we are assembling the proper system.
  libmesh_assert_equal_to (system_name, "Poisson");

  PMSystemPoisson& _pm_system = _eqn_sys.get_system<PMSystemPoisson> (system_name);

  // Numeric id corresponding to potential variable
  const unsigned int phi_var = _pm_system.variable_number ("phi");

  // Get FE type for "phi"
  FEType fe_phi_type = _pm_system.variable_type(phi_var);
  UniquePtr<FEBase> fe_phi (FEBase::build(_dim, fe_phi_type));
  // Define Gauss quadrature rule for numerical integration.
  QGauss qrule (_dim, SECOND);
  fe_phi->attach_quadrature_rule (&qrule);

  // build the face element for boundary integration
  UniquePtr<FEBase> fe_face (FEBase::build(_dim, fe_phi_type));
  QGauss qface(_dim-1, SECOND);
  fe_face->attach_quadrature_rule (&qface);

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>& JxW              = fe_phi->get_JxW();
  const std::vector<std::vector<Real>>& phi = fe_phi->get_phi();
  const std::vector<Point>& q_xyz           = fe_phi->get_xyz(); // xyz coords of quad pts

  // A reference to the DofMap object for this system.
  const DofMap & dof_map = _pm_system.get_dof_map();

  // Define data structures to contain the element matrix Ke and vector Fe
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // Build _int_force vector at the beginning of Poisson solver
  // This can work for both Poisson and Stokes equations, since there is only one
  // forcing term on the right-hand-side of the equations.
  if(_int_force.size() == 1){
    if(_pm_system.comm().rank()==0){
      printf("\nassemble_int_force() for Poisson solver at the beginning of simulation\n\n");
    }

    const unsigned int n_mesh_elem = _mesh.n_elem();
    _int_force.resize(n_mesh_elem);
    _q_xyz.resize(n_mesh_elem);
    _n_dofs.resize(n_mesh_elem);
    _dof_indices.resize(n_mesh_elem);

    // Now we will loop over all the elements in the mesh that live on the local processor.
    MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = _mesh.active_local_elements_end();
    for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently working on.
      const Elem* elem = *el;
      const unsigned int elem_id = elem->id();

      // Get the degree of freedom indices for the current element.
      dof_map.dof_indices (elem, _dof_indices[elem_id]);

      _n_dofs[elem_id] = _dof_indices[elem_id].size();
      Fe.resize (_n_dofs[elem_id]);

      // NOTE: here JxW and dphi and other element quantities are not computed up to now,
      // and these will be done in the elem loop after fe->reinit()
      fe_phi->reinit (elem);

      //qrule.print_info();
      this->assemble_int_force(elem, _n_dofs[elem_id], *fe_phi);
      // printf("finished assemble_int_force\n");
    }
  }

  // Now we will loop over all the elements in the mesh that live
  // on the local processor, and compute the element vector Fe.
  MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_local_elements_end();
  for ( ; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently working on.
    const Elem* elem = *el;
    const unsigned int elem_id = elem->id();

    // Get the degree of freedom indices for the current element.
    // FIXME:Why do we update dof_indices again?
    //dof_map.dof_indices (elem, _dof_indices[elem_id]);

    // const unsigned int n_dofs = dof_indices.size();
    Fe.resize (_n_dofs[elem_id]);

    // NOTE: here JxW and dphi and other element quantities are not computed up to now,
    // and these will be done in the elem loop after fe->reinit()
    //fe_vel->reinit (elem);
    //qrule.print_info();

    // if elem_neighbor_list is pre-built, we can access it directly
    const std::vector<std::size_t>& n_list = _pm_system.point_mesh()->elem_neighbor_list(elem);

    // if this elem has no neighboring particle we turn the pc_flag to 'false'
    bool pc_flag = true;                  // a flag for the point charge
    if(n_list.size()==0) pc_flag = false; // No point force because no point list.

    // Now compute Fe caused by the regularized point charge.
    // FIXME: need to implement space charge density when considering Nernst-Planck solver
    this->compute_element_rhs(elem, _n_dofs[elem_id], *fe_phi, n_list, pc_flag, option, Fe);

    // Impose the Dirichlet BC (electrical potential) via the penalty method.
    this->apply_bc_by_penalty(elem, "vector", Ke, Fe, option);

    // Impose Neumann BC (surface charge density) to the right hand side 
    this->apply_bc_neumann(elem, *fe_phi, *fe_face, Fe);

    // If this assembly program were to be used on an adaptive mesh,
    // we would have to apply any hanging node constraint equations.
    dof_map.constrain_element_vector (Fe, _dof_indices[elem_id]);

    // Add the element rhs vector to the global system.
    //PMToolBox::zero_filter_dense_vector(Fe, 1e-10);
    //PMToolBox::output_dense_vector(Fe);
    _pm_system.rhs->add_vector (Fe, _dof_indices[elem_id]);
  }// end for elem-loop

  STOP_LOG ("assemble_global_F()", "AssemblePoisson");

  //if (_pm_system.comm().rank()==0){
  //  printf("assemble_global_F(): The global RHS vector has been assembled ...\n");
  //}
  return;
}



// ==================================================================================
void AssemblePoisson::compute_element_rhs(const Elem* elem,
                                          const unsigned int n_u_dofs,
                                          FEBase& fe_v,
                                          const std::vector<std::size_t> n_list,
                                          const bool& pc_flag,
                                          const std::string& option,
                                          DenseVector<Number>& Fe)
{
  START_LOG("compute_element_rhs()", "AssemblePoisson");

  libmesh_assert_equal_to (system_name, "Poisson");

  PMSystemPoisson& _pm_system = _eqn_sys.get_system<PMSystemPoisson> ("Poisson");

  // Get a reference to the PointMesh, PMPeriodicBoundary, and PointParticles
  PointMesh<3>* _point_mesh = _pm_system.point_mesh();
  PMPeriodicBoundary* _pm_periodic_boundary = _point_mesh->pm_periodic_boundary();
  std::vector<PointParticle*> _particles = _point_mesh->particles();

  // The element Jacobian * quadrature weight at each quad pt(high order Qgauss).
  //const std::vector<Real>& JxW                = fe_v.get_JxW();
  //const std::vector<std::vector<Real> >& phi  = fe_v.get_phi();
  //const std::vector<Point>& q_xyz             = fe_v.get_xyz(); // xyz coords of quad pts
  //printf("q_xyz size = %d\n", q_xyz.size());
  //fe_v.reinit(elem);

  const unsigned int elem_id = elem->id();
  const std::vector<Point>& q_xyz = _q_xyz[elem_id]; // xyz coords of quad pts

  // 1. Add the regularized point charge, first examine if this element has
  // neighboring point charge sources
  if( pc_flag )
  {
    // Number of particles in the neighbor list
    const std::size_t n_pts = n_list.size();

    // Initialize variables for position, charge, distance, etc.
    Point np_pos(0.);
    Real r=0., charge_val=0., np_charge=0.;
    unsigned int qp_size = q_xyz.size();
    //printf("qp_size = %d\n", q_xyz.size());

    // Now we will build the element RHS using gauss quadrature integration.
    // first loop over all neighboring particles near this element
    for(unsigned int np = 0; np<n_pts; ++np){
      // Charge on this bead
      np_charge = _particles[n_list[np]]->charge();

      // Get the location of this bead
      np_pos = _particles[n_list[np]]->point();

      for(unsigned int qp=0; qp<qp_size; qp++){
        // Distance from quadrature point to the charge point
        r = _pm_periodic_boundary->point_distance(q_xyz[qp], np_pos);

        // Evaluate the value of regularized gaussian charge at this quadrature point
        charge_val = ggem_poisson->smoothed_charge_exp(r) * np_charge;
//std::cout<<"np_charge="<<np_charge<<", r="<<r<<", smoothed_function="<<ggem_poisson->smoothed_charge_exp(r)<<", charge_val="<<charge_val<<std::endl;
        // FIXME:Need to add nodal space charge density from ion concentration fields
        // this will need to access PMSystemNP, and approximate ion concentration on
        // quadrature points?
        //Real space_charge_density = 0.;

       	for(unsigned int k=0; k<n_u_dofs; ++k){
          //Fe(k) += JxW[qp]*phi[k][qp]*charge_val;
          Fe(k) += _int_force[elem_id][k*qp_size+qp]*charge_val;
       	} // end loop over nodes (dofs)
      } // end loop over quadrature points
    } // end loop over beads
  } // end if( pc_flag )

  STOP_LOG("compute_element_rhs()", "AssemblePoisson");
  return;
}



// ==================================================================================
void AssemblePoisson::select_boundary_side(const Elem* elem)
{
  START_LOG("select_boundary_side()", "AssemblePoisson");

  // Get a reference to the Particle-Mesh System.
  PMSystemPoisson& pm_system = _eqn_sys.get_system<PMSystemPoisson> ("Poisson");
  const std::size_t elem_id = elem->id();

  // Get boundary ids for each BC
  const std::vector<unsigned int> boundary_id_dirichlet_poisson = _eqn_sys.parameters.get<std::vector<unsigned int>> ("boundary_id_dirichlet_poisson");
  const std::vector<unsigned int> boundary_id_neumann_poisson   = _eqn_sys.parameters.get<std::vector<unsigned int>> ("boundary_id_neumann_poisson");

  // The following loops over the sides of the element. If the element has NO
  // neighbors on a side then that side MUST live on a boundary of the domain.
  for(unsigned int s=0; s<elem->n_sides(); s++)
  {
    bool apply_bc_dirichlet = false, apply_bc_neumann = false;

    // If this side is on the boundary
    if(elem->neighbor(s) == NULL)
    {
      // Find this side's corresponding Dirichlet BoundaryID. If found, stop looping all BoundaryIDs 
      for(unsigned int i=0; i<boundary_id_dirichlet_poisson.size(); i++){
        if(_mesh.get_boundary_info().has_boundary_id(elem, s, boundary_id_dirichlet_poisson[i])){
          apply_bc_dirichlet = true;
          break;
        }
      }

      // Find this side's corresponding Neumann BoundaryID. If found, stop looping all BoundaryIDs 
      for(unsigned int i=0; i<boundary_id_neumann_poisson.size(); i++){
        if(_mesh.get_boundary_info().has_boundary_id(elem, s, boundary_id_neumann_poisson[i])){
          apply_bc_neumann = true;
          break;
        }
      }
    }

    // If BC is applied, store this side_id for this element
    if(apply_bc_dirichlet) _boundary_sides_dirichlet_poisson[elem_id].push_back(s);
    if(apply_bc_neumann)   _boundary_sides_neumann_poisson[elem_id].push_back(s);
  }

  STOP_LOG("select_boundary_side()", "AssemblePoisson");
}



// =======================================================================================
void AssemblePoisson::apply_bc_by_penalty(const Elem* elem,
                                          const std::string& matrix_or_vector,
                                          DenseMatrix<Number>& Ke,
                                          DenseVector<Number>& Fe,
                                          const std::string& option)
{
  START_LOG("apply_bc_by_penalty()", "AssemblePoisson");

  // Get a reference to the Particle-Mesh System
  PMSystemPoisson& pm_system = _eqn_sys.get_system<PMSystemPoisson> ("Poisson");

  const Real penalty = 1E6; //The penalty value
  const unsigned int n_nodes = elem->n_nodes();
  const std::size_t elem_id = elem->id();

  // Dirichlet boundaries
  const std::vector<unsigned int> boundary_id_dirichlet_poisson = _eqn_sys.parameters.get<std::vector<unsigned int>> ("boundary_id_dirichlet_poisson");
  const std::vector<Real> boundary_value_dirichlet_poisson = _eqn_sys.parameters.get<std::vector<Real>> ("boundary_value_dirichlet_poisson");

  // Loop through sides in this element that sits on system's boundary
  for (unsigned int s=0; s<_boundary_sides_dirichlet_poisson[elem_id].size(); s++)
  {
    // Build the full-order side element for "potential" Dirichlet BC at the walls.
    UniquePtr<Elem> side (elem->build_side(_boundary_sides_dirichlet_poisson[elem_id][s]));

    // Total electrical potential on this side
    Real phi_total = 0.;
    for(unsigned int i=0; i<boundary_value_dirichlet_poisson.size(); i++){
      if(_mesh.get_boundary_info().has_boundary_id(elem, _boundary_sides_dirichlet_poisson[elem_id][s] , boundary_id_dirichlet_poisson[i])){
        phi_total = boundary_value_dirichlet_poisson[i];
        break;
      }
    }

    // Loop through nodes on this side element
    for (unsigned int nn=0; nn<side->n_nodes(); nn++)
    {
      // Calculate local part of the electrical potential from GGEM
      const Point ptx = side->point(nn); // Coordinate of the node
      const Real phi_local = pm_system.local_potential_field(elem, ptx, "regularized");
      //const Real phi_local = pm_system.local_potential_field(ptx, "regularized");

      // When we run ggem_validation_poisson, Dirichelet boundary condition is applied on all
      // boundaries, the electrical potential on boundaries due to point charges are evaluated
      // using free-space Green's function.
      if(_eqn_sys.parameters.get<std::string> ("simulation_name") == "ggem_validation_poisson")
      {
        phi_total = analytical_solution -> exact_solution_infinite_domain(*ggem_poisson, ptx);
      }

      // Because of Ewald split: global_potential = total_potential - ggem_local_potential
      Real phi_global = phi_total - phi_local;

      // Find the node on the element matching this node on the side.
      // That defined where in the element matrix the BC will be applied.
      for (unsigned int n=0; n<elem->n_nodes(); n++)
      {
        if (elem->node(n) == side->node(nn))
        {
          // Penalize phi at the current node, var_number is 0 for 'phi'
          this->penalize_elem_matrix_vector(Ke,Fe,matrix_or_vector,0,n,n_nodes,penalty,phi_global);
        } // end if (elem->node(n) == side->node(nn))
      } // enf for n-loop
    } // end for nn-loop
  } // end for s-loop

  STOP_LOG("apply_bc_by_penalty()", "AssemblePoisson");
}



// =======================================================================================
void AssemblePoisson::apply_bc_neumann(const Elem* elem,
                                       FEBase& fe_phi,
                                       FEBase& fe_face,
                                       DenseVector<Number>& Fe)
{
  START_LOG("apply_bc_neumann()", "AssemblePoisson");

  // Get a reference to the Particle-Mesh System
  PMSystemPoisson& pm_system = _eqn_sys.get_system<PMSystemPoisson> ("Poisson");

  const unsigned int n_nodes = elem->n_nodes();
  const std::size_t elem_id = elem->id();

  // Neumann boundaries
  const std::vector<unsigned int> boundary_id_neumann_poisson = _eqn_sys.parameters.get<std::vector<unsigned int>> ("boundary_id_neumann_poisson");
  const std::vector<Real> boundary_value_neumann_poisson = _eqn_sys.parameters.get<std::vector<Real>> ("boundary_value_neumann_poisson");

  // Integration terms for side element
  const std::vector<Real>               & JxW_face = fe_face.get_JxW();
  const std::vector<std::vector<Real> > & phi_face = fe_face.get_phi();
  const std::vector<Point>          & face_normals = fe_face.get_normals();

  // Extract the shape function derivatives to be evaluated at the nodes
  const std::vector<std::vector<RealGradient> > & dphi = fe_phi.get_dphi();
  // Extract the element node coordinates in the reference frame
  std::vector<Point> nodes;
  fe_phi.get_refspace_nodes(elem->type(), nodes);
  // Evaluate the shape functions derivatives at the nodes
  fe_phi.reinit(elem, &nodes);

  // Local electrical potential and its gradient on all nodes in this element
  std::vector<Real> phi_local;
  phi_local.resize(elem->n_nodes());
  std::vector<RealGradient> dphi_local;
  dphi_local.resize(elem->n_nodes());
  // Loop through nodes on this element
  for (unsigned int nn=0; nn<elem->n_nodes(); nn++)
  {
    // Coordinate of the node
    const Point ptx = elem->point(nn);
    // Evaluate local electrical potential (phi) on each node
    phi_local[nn] = pm_system.local_potential_field(elem, ptx, "regularized");
    //phi_local[nn] = pm_system.local_potential_field(ptx, "regularized");
  }

  // Gradient of electrical potential on all nodes in this element
  for (unsigned int nn=0; nn<elem->n_nodes(); nn++)
    for (unsigned int i=0; i<phi_local.size(); i++)
      dphi_local[nn] += phi_local[i] * dphi[i][nn];

  // Loop through sides in this element that sits on system's boundary
  for (unsigned int s=0; s<_boundary_sides_neumann_poisson[elem_id].size(); s++)
  {
    // Build the full-order side element for "potential" Dirichlet BC at the walls.
    UniquePtr<Elem> side (elem->build_side(_boundary_sides_neumann_poisson[elem_id][s]));

    // dphi / dn (surface charge density) on this side
    Real sigma_total = 0.;
    for(unsigned int i=0; i<boundary_value_neumann_poisson.size(); i++){
      if(_mesh.get_boundary_info().has_boundary_id(elem, _boundary_sides_neumann_poisson[elem_id][s], boundary_id_neumann_poisson[i])){
        sigma_total = boundary_value_neumann_poisson[i];
        break;
      }
    }

    // Calcualte Jacobian*Weights and shape functions on this side
    fe_face.reinit(elem, s);

    // Calculate (dphi / dn)_local on side nodes
    std::vector<Real> sigma_local;
    sigma_local.resize(side->n_nodes());
    // Loop through nodes on this side element
    for (unsigned int nn=0; nn<side->n_nodes(); nn++)
    {
      // Find the node on the element matching this node on the side.
      // That defined where in the element matrix the BC will be applied.
      for (unsigned int n=0; n<elem->n_nodes(); n++)
      {
        if (elem->node(n) == side->node(nn))
        {
          // Use the normal vector at a quadrature point, since
          // normally the side element has zero curvature
          sigma_local[nn] = dphi_local[n](0) * face_normals[0](0)
                          + dphi_local[n](1) * face_normals[0](1)
                          + dphi_local[n](2) * face_normals[0](2);
        }
      }
    }

    // Calculate (dphi / dn)_global
    std::vector<Real> sigma_global;
    sigma_global.resize(side->n_nodes());
    for (unsigned int nn=0; nn<side->n_nodes(); nn++)
      sigma_global[nn] = sigma_total - sigma_local[nn];

    // Integrate on this side
    // FIXME: this can be optimized by store JxW_face * phi_face for face element
    // in a vector similar to that in _int_force
    for (unsigned int qp=0; qp<JxW_face.size(); qp++)
      for (unsigned int i=0; i<phi_face.size(); i++)
        Fe(i) += JxW_face[qp] * phi_face[i][qp] * sigma_global[i];
  }//end for s-loop

  STOP_LOG("apply_bc_neumann()", "AssemblePoisson");
}



// ============================================================================================
void AssemblePoisson::init_ggem_poisson(const std::string& system_name)
{
  START_LOG("init_ggem_poisson()", "AssemblePoisson");

  libmesh_assert_equal_to(system, "Poisson");

  PMSystemPoisson& _pm_system = _eqn_sys.get_system<PMSystemPoisson> (system_name);
  PointMesh<3>* _point_mesh = _pm_system.point_mesh();
  const PointType point_type = _point_mesh->particles()[0]->point_type();
  ggem_poisson -> set_alpha(_eqn_sys.parameters.get<Real> ("alpha"));
  ggem_poisson -> set_br0(_eqn_sys.parameters.get<Real> ("br0"));
  ggem_poisson -> set_point_type(point_type);
  ggem_poisson -> set_ksi();
  ggem_poisson -> set_phi0(_eqn_sys.parameters.get<Real> ("phi0"));

  STOP_LOG("init_ggem_poisson()", "AssemblePoisson");
}
