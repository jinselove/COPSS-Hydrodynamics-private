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
#include "assemble_navier_stokes.h"
#include "pm_linear_implicit_system.h"

// ==================================================================================
AssembleNS::AssembleNS(EquationSystems& es)
: _eqn_sys(es),
_mesh(es.get_mesh())

{
  // do nothing
  _dim  = es.get_mesh().mesh_dimension();
  
}



// ==================================================================================
AssembleNS::~AssembleNS()
{
  // do nothing
}







// ==================================================================================
void AssembleNS::assemble_global_K(const std::string& system_name,
                                   const std::string& option)
{
  /*! It is a good idea to make sure we are assembling the proper system.
  
  */
  libmesh_assert_equal_to (system_name, "Stokes");
//  const MeshBase& _mesh = _eqn_sys.get_mesh();
  PMLinearImplicitSystem& _pm_system = _eqn_sys.get_system<PMLinearImplicitSystem> (system_name);
  // Numeric ids corresponding to each variable in the system
  const unsigned int  u_var = _pm_system.variable_number ("u");      // u_var = 0
  const unsigned int  v_var = _pm_system.variable_number ("v");      // v_var = 1
  unsigned int w_var = 0;
  if(_dim==3)         w_var = _pm_system.variable_number ("w");      // w_var = 2 if dim==3
  const unsigned int  p_var = _pm_system.variable_number ("p");      // p_var = 2(dim=2); = 3(dim=3)
  
  // Get the Finite element type for "u".
  // Note "u" is the same as the type for "v".
  FEType fe_vel_type  = _pm_system.variable_type(u_var);
  
  // Get the Finite element type for "p"
  FEType fe_pres_type = _pm_system.variable_type(p_var);
  
  // Build a Finite element object of the specified type for the velocity variables
  UniquePtr<FEBase> fe_vel  (FEBase::build(_dim, fe_vel_type) );

  // Build a Finite element object of the specified type for the pressure variables
  UniquePtr<FEBase> fe_pres (FEBase::build(_dim, fe_pres_type));
   
  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  QGauss qrule (_dim, fe_vel_type.default_quadrature_order());   //3^dim pts
  
  // Tell the finite element objects to use our quadrature rule.
  fe_vel ->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);
  
  // Here we define some references to cell-specific data that will be used to assemble
  // linear system.

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>& JxW                        = fe_vel->get_JxW();
  
  // The element shape function gradients for the velocity variables evaluated at the 
  // quadrature points. 
  const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();
  
  // 
  const std::vector<std::vector<Real> >& phi          = fe_vel->get_phi();
  
  // The element shape functions for the pressure variable evaluated at the quadrature
  //points.
  const std::vector<std::vector<Real> >& psi          = fe_pres->get_phi();
  
  
  // A reference to the DofMap object for this system.
  const DofMap & dof_map = _pm_system.get_dof_map();
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u, dof_indices_p;
  
  // Element matrix contribution 
  DenseMatrix<Number> Ke; 

  // Input of compute_element_rhs() 
  DenseVector<Number> Fe; 
  
  // Element pressure mass matrix
  DenseMatrix<Number> Mp;  
  
  // retrieve system parameters
  const std::string schur_pc_type =  _eqn_sys.parameters.get<std::string> ("schur_pc_type");
  const bool user_defined_pc      =  _eqn_sys.parameters.get<bool>        ("user_defined_pc");
  const Real mu                   =  _eqn_sys.parameters.get<Real>        ("viscosity_0");
  
  // -------------------------------------------------------------------------------------------
  //if(system.comm().rank()==0){
  //  printf("assemble_global_K(): Start to assemble the global matrix K ...\n");
  //}
  // -------------------------------------------------------------------------------------------
  
  
  // Now we will loop over all the elements in the mesh that live
  // on the local processor. We will compute the element matrix Ke. In case users
  // later modify the program to include refinement, we will be safe and will only
  // consider the active elements; hence we use a variant of the active_elem_iterator
  MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_local_elements_end();
  

  for ( ; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently working on.
    const Elem* elem = *el;
    
    // Get the degree of freedom indices for the current element. Tnese define where
    // in the global matrix and right-hand-side this element will contribute to.
    dof_map.dof_indices (elem, dof_indices);
    dof_map.dof_indices (elem, dof_indices_u, u_var);
    dof_map.dof_indices (elem, dof_indices_p, p_var);
    
    const unsigned int n_dofs   = dof_indices.size();
    const unsigned int n_u_dofs = dof_indices_u.size();
    const unsigned int n_p_dofs = dof_indices_p.size();
    const unsigned int n_uvw_dofs = n_u_dofs*_dim;
    
    // Zero the element matrix Ke before summing them. We use the resize member here because
    // the number of degrees of freedom might have changed from the last element. Note that
    // this will be the case if the element type is different (i.e. the last element was a 
    // triangle, now we are on a quadrilateral.)
    Ke.resize (n_dofs, n_dofs);
    
    // Compute the element-specific data for the current element. This involves computing
    // the location of the quadrature points (q_point) and the shape functions (phi, dphi)
    // for the current element.
    fe_vel->reinit (elem);
    fe_pres->reinit(elem);
    
    DenseMatrix<Number> Ktt(n_u_dofs,n_u_dofs);

    for (unsigned int i=0; i<_dim; ++i)
    {
      DenseMatrix<Number> Kij;
      this->assemble_element_KIJ(JxW,dphi,mu,n_u_dofs,i,i,Kij);
      Ktt += Kij;
    }
    
    // Compute KIJ and add it to Ke
    for (unsigned int i=0; i<_dim; i++)
    {
      for (unsigned int j=0; j<_dim; j++)
      {
        // -------------------------------------------------------------
        DenseMatrix<Number> Kij;
        this->assemble_element_KIJ(JxW,dphi,mu,n_u_dofs,i,j,Kij);
        
        // add Kij to Ke
        for(unsigned int k=0; k<n_u_dofs; ++k)
          for(unsigned int l=0; l<n_u_dofs; ++l)
            Ke(n_u_dofs*i+k,n_u_dofs*j+l) += Kij(k,l);
        // -------------------------------------------------------------
        
        if(i==j)
        {
          for(unsigned int k=0; k<n_u_dofs; ++k)
            for(unsigned int l=0; l<n_u_dofs; ++l)
              Ke(n_u_dofs*i+k,n_u_dofs*j+l) += Ktt(k,l);
        }
      } // end for J-loop
      
      // Compute Q matrix (Kup, Kvp, Kwp), and add it to Ke
      DenseMatrix<Number> Qi;
      this->assemble_element_QI(JxW,dphi,psi,n_u_dofs,n_p_dofs,i,Qi);
      for(unsigned int k=0; k<n_u_dofs; ++k)
      {
        for(unsigned int l=0; l<n_p_dofs; ++l)
        {
          Ke(n_u_dofs*i+k,n_uvw_dofs+l) += -Qi(k,l);
          Ke(n_uvw_dofs+l,n_u_dofs*i+k) += -Qi(k,l);
        }
      }
      
    } // end for I-loop
    
    
//    if(_eqn_sys.comm().rank()==0){
//      printf("----->TEST: element id = %u\n", elem->id());
//      PMToolBox::output_dense_matrix(Ktt);
//    }
    
    
    // apply BCs by penalty method
    this->apply_bc_by_penalty(elem, "matrix", Ke, Fe, option);
    
    
    // If this assembly program were to be used on an adaptive mesh,
    // we would have to apply any hanging node constraint equations.
    dof_map.constrain_element_matrix (Ke, dof_indices);
    
    
    // Add the element matrix and rhs vector to the global system.
    //GeomTools::zero_filter_dense_matrix(Ke, 1e-10);
    //PMToolBox::output_dense_matrix(Ke);
    _pm_system.matrix->add_matrix (Ke, dof_indices);
    
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     * The following will modify the elem pc matrix, then add it to the global PCMatrix
     * It has to be done after the penalty, because Ke and Fe have been changed!
     * FIXME: how can we modify Mp according to the change of Ke due to penalty?
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    if( user_defined_pc )
    {
      // pressure mass matrix (negative) Mp
      if ( schur_pc_type=="SMp" || schur_pc_type=="SMp_lump" )
      {
        DenseMatrix<Number> Mij;
        this->assemble_element_MIJ(JxW, psi, n_p_dofs, Mij);
        
        Mp.resize(n_p_dofs,n_p_dofs);
        for (unsigned int i=0; i<n_p_dofs; i++)
          for (unsigned int j=0; j<n_p_dofs; j++)
            Mp(i,j) += -Mij(i,j)/mu;
      } // end if ( schur_pc_type=="SMp" )
      
      // Mp pc
      if(schur_pc_type=="SMp")
        _pm_system.get_matrix("Preconditioner").add_matrix (Mp, dof_indices_p);
      // -------------------------------------------------------------------------------
      
      //
      if(schur_pc_type=="SMp_lump")
      {
        for (unsigned int i=0; i<n_p_dofs; i++)
          for (unsigned int j=0; j<n_p_dofs; j++)
            if (i != j) { Mp(i,i) += Mp(i,j);   Mp(i,j)  = 0.0;  }
        
        _pm_system.get_matrix("Preconditioner").add_matrix (Mp, dof_indices_p);
      }
      // -------------------------------------------------------------------------------
    } // end if( user_defined_pc )
    
  } // end of elem-loop
  
  
  // That's it.
  // -------------------------------------------------------------------------------------------
  //  if (_pm_system.comm().rank()==0){
  //  printf("assemble_matrix_K(): The global matrix K has been assembled ...\n");
  // }
  // -------------------------------------------------------------------------------------------
  return;
  
}




// ==================================================================================
void AssembleNS::assemble_global_F(const std::string& system_name,
                                   const std::string& option)
{
  // It is a good idea to make sure we are assembling the proper system.
  libmesh_assert_equal_to (system_name, "Stokes");
//  const MeshBase& _mesh = _eqn_sys.get_mesh();
  PMLinearImplicitSystem& _pm_system = _eqn_sys.get_system<PMLinearImplicitSystem> (system_name);
  
  // Numeric ids corresponding to each variable in the system
  const unsigned int  u_var = _pm_system.variable_number ("u");      // u_var = 0
  const unsigned int  p_var = _pm_system.variable_number ("p");      // p_var = 2(dim=2); = 3(dim=3)
  
  // Get the FE type for "u" and "p".  Note "u" is the same as the type for "v".
  FEType fe_vel_type  = _pm_system.variable_type(u_var);
  UniquePtr<FEBase> fe_vel  (FEBase::build(_dim, fe_vel_type) );
  
  // Define Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  QGauss qrule (_dim, fe_vel_type.default_quadrature_order());   //3^dim pts
  //QGauss qrule (_dim, SIXTH);   // SIXTH
  fe_vel->attach_quadrature_rule (&qrule);
  
  // build the face element for boundary traction
  UniquePtr<FEBase> fe_face (FEBase::build(_dim, fe_vel_type));
  QGauss qface(_dim-1, fe_vel_type.default_quadrature_order());
  fe_face->attach_quadrature_rule (&qface);
  
  
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>& JxW                = fe_vel->get_JxW();
  const std::vector<std::vector<Real> >& phi  = fe_vel->get_phi();
  const std::vector<Point>& q_xyz             = fe_vel->get_xyz(); // xyz coords of quad pts
  
  // A reference to the DofMap object for this system.
  const DofMap & dof_map = _pm_system.get_dof_map();
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u, dof_indices_p;
  
  
  // Define data structures to contain the element matrix Ke and vector Fe
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;
  
  
  // retrieve system parameters
  const Real alpha   =  _eqn_sys.parameters.get<Real> ("alpha");
  
  // -------------------------------------------------------------------------------------------
  //if(system.comm().rank()==0){
  //  printf("\nassemble_global_F(): Start to assemble the global force vector F ...\n");
  //}
  // -------------------------------------------------------------------------------------------

  
  // Now we will loop over all the elements in the mesh that live
  // on the local processor, and compute the element matrix Ke.
  MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_local_elements_end();
  for ( ; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently working on.
    const Elem* elem = *el;
    
    // Get the degree of freedom indices for the current element.
    dof_map.dof_indices (elem, dof_indices);
    dof_map.dof_indices (elem, dof_indices_u, u_var);
    dof_map.dof_indices (elem, dof_indices_p, p_var);
    
    const unsigned int n_dofs   = dof_indices.size();
    const unsigned int n_u_dofs = dof_indices_u.size();
    const unsigned int n_p_dofs = dof_indices_p.size();
    const unsigned int n_uvw_dofs = n_u_dofs*_dim;
    Fe.resize (n_dofs);
    
    // NOTE: here JxW and dphi and other element quantities are not computed up to now,
    // and these will be done in the elem loop after fe->reinit()
    fe_vel->reinit (elem);
    //qrule.print_info();
    
    // if elem_neighbor_list is pre-built, we can access it directly
    const std::vector<std::size_t> n_list = _pm_system.point_mesh()->elem_neighbor_list(elem);
    
    
    // Now compute Fe caused by the regularized point force and boundary traction.
    // if this elem has no neighboring particle or only the undisturbed field
    // is required, we turn the pf_flag to 'false'!
    bool pf_flag = true;                          // a flag for the point force
    if (n_list.size()==0)       pf_flag = false;  // No point force because no point list.
    if(option == "undisturbed") pf_flag = false;  // No point force
    
    this->compute_element_rhs(elem, n_u_dofs, *fe_vel, n_list,
                              pf_flag, option, alpha, Fe);
    
    
    // imposed the Dirichlet BC at no-slip walls & pressure jump at the inlet/outlet
    // via the penalty method.
    this->apply_bc_by_penalty(elem, "vector", Ke, Fe, option);
    
    
    // If this assembly program were to be used on an adaptive mesh,
    // we would have to apply any hanging node constraint equations.
    dof_map.constrain_element_vector (Fe, dof_indices);
    
    
    // Add the element matrix and rhs vector to the global system.
    //PMToolBox::zero_filter_dense_vector(Fe, 1e-10);
    //PMToolBox::output_dense_vector(Fe);
    _pm_system.rhs->add_vector (Fe, dof_indices);
  } // end for elem-loop
  
  
  // ---------------------------------------------------------------------------------------------
  //if (_pm_system.comm().rank()==0){
  //  printf("assemble_global_F(): The global RHS vector has been assembled ...\n");
  // }
  // ---------------------------------------------------------------------------------------------
  return;
}


// ==================================================================================
void AssembleNS::compute_element_rhs(const Elem*     elem,
                                     const unsigned int n_u_dofs,
                                     FEBase& fe_v,
                                     const std::vector<std::size_t> n_list,
                                     const bool& pf_flag,
                                     const std::string& option,
                                     const Real& alpha,
                                     DenseVector<Number>& Fe)
{
  START_LOG("compute_element_rhs()", "AssembleNS");  // libMesh log
  
  libmesh_assert_equal_to (system_name, "Stokes");
//  const MeshBase& _mesh = _eqn_sys.get_mesh();
  PMLinearImplicitSystem& _pm_system = _eqn_sys.get_system<PMLinearImplicitSystem> ("Stokes");
 /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Get a reference to the Particle-Mesh linear implicit system.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PointMesh<3>* _point_mesh = _pm_system.point_mesh();
  PMPeriodicBoundary* _pm_periodic_boundary = _point_mesh -> pm_periodic_boundary();
  const std::vector<bool>& _inlet_direction = _pm_periodic_boundary -> inlet_direction();
  const std::vector<Real>& _inlet_pressure = _pm_periodic_boundary->inlet_pressure();
  std::vector<PointParticle*> _particles = _point_mesh->particles();
 // The element Jacobian * quadrature weight at each quad pt(high order Qgauss).
  const std::vector<Real>& JxW                = fe_v.get_JxW();
  const std::vector<std::vector<Real> >& phi  = fe_v.get_phi();
  const std::vector<Point>& q_xyz             = fe_v.get_xyz(); // xyz coords of quad pts
  fe_v.reinit(elem);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   1. add the regularized point force for "disturbed" flow field!
   first examine if this element is "close to" the point force sources
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if( pf_flag && (option == "disturbed") )
  {
    // get the location of points in the neighbor list
    GGEMSystem ggem_sys;
    const std::size_t n_pts = n_list.size();
    // Now we will build the element RHS using high order gauss quadrature.
    // first loop over all neighboring particles near this element
    Point np_pos(0.);
    std::vector<Real> np_force(_dim, 0.);
    Real r = 0. , force_val = 0.;
    std::vector<Real> fvalues(_dim,0.);
    for(unsigned int np=0; np<n_pts; ++np)
    {
      // force vector on each particle(note they are different!)
     // std::vector<Real> pf(_dim,0.0);    // point force vector of the current particle
     np_force = _particles[n_list[np]]->particle_force();
     np_pos = _particles[n_list[np]]->point();  

    // Next loop over gauss points
      for (unsigned int qp=0; qp<q_xyz.size(); qp++)
      {
        // distance from Gaussian point to the force point
        r = _pm_periodic_boundary->point_distance(q_xyz[qp], np_pos);
        
        // evaluate the value of gauss force at this quad pt.
        force_val = ggem_sys.smoothed_force_exp(r, alpha);
        
        // force magnitudes on each particle(note they are different!)
        for (unsigned int j=0; j<_dim; ++j) {
          fvalues[j] = np_force[j]*force_val;
        }
        
        // compute the element rhs vector
        for (unsigned int j=0; j<_dim; ++j){
          for(unsigned int k=0; k<n_u_dofs; ++k){
            Fe(j*n_u_dofs + k) += JxW[qp]*phi[k][qp]*fvalues[j];
          } // end for k-loop
        } // end for j-loop
        
      } // end for qp-loop
    } // end for np-loop
    
  } // end if( pf_flag )

  if( option == "undisturbed" )
  {
    FEType fe_vel_type = fe_v.get_fe_type();
    UniquePtr<FEBase> fe_face (FEBase::build(_dim, fe_vel_type));
    QGauss qface(_dim-1, fe_vel_type.default_quadrature_order());
    fe_face->attach_quadrature_rule (&qface);
     
    const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
    const std::vector<Real>&                JxW_face = fe_face->get_JxW();
    const std::vector<Point>&             q_xyz_face = fe_face->get_xyz(); // xyz of quad pts

            // The following loops over the sides of the element. If the element has NO
            // neighbors on a side then that side MUST live on a boundary of the domain.
            for (unsigned int s=0; s<elem->n_sides(); s++)
            {
              if (elem->neighbor(s) == NULL)
              {
                for (int i = 0; i < _dim; i++)
                {
                  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                  // This operation applies the traction jump on the inlet and outlet boundaries.
                  // when s is on either face in direction i and inlet_direction[i]==true
                  //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                  if((_mesh.get_boundary_info().has_boundary_id(elem,s,_boundary_id_3D[2*i]) or
		      _mesh.get_boundary_info().has_boundary_id(elem,s,_boundary_id_3D[2*i+1])) and
		      _inlet_direction[i])
		  {
                    fe_face->reinit(elem, s);  
                    // traction: t = sigma*n.
                    for (unsigned int qp=0; qp<JxW_face.size(); qp++)
                    {
                        // traction jump at the inlet and outlet
                        for (unsigned int j=0; j<n_u_dofs; j++)
                        {
                              Fe(i*n_u_dofs + j) += JxW_face[qp]* _inlet_pressure[i] * phi_face[j][qp];
                        }
                    } // end for qp-loop
                   //PMToolBox::zero_filter_dense_vector(Fe,1E-10);
                   //PMToolBox::output_dense_vector(Fe);
                  }
                }
              } // end if(elem->neighbor(s) == NULL)     
            } // end for s-loop
  } // end if( option == "undisturbed" )

  STOP_LOG("compute_element_rhs()", "AssembleNS");
  return;
}


// ==================================================================================
void AssembleNS::assemble_element_KIJ(const std::vector<Real>& JxW,
                                      const std::vector<std::vector<RealGradient> >& dphi,
                                      const Real&     mu,
                                      const unsigned int n_u_dofs,
                                      const unsigned int i,
                                      const unsigned int j,
                                      DenseMatrix<Number>& Kij)
{
  START_LOG("assemble_element_KIJ()", "AssembleNS");  // libMesh log
  
  Kij.resize(n_u_dofs, n_u_dofs);
  
  // reinit and compute the element matrix K_ij
  for (unsigned int qp=0; qp<JxW.size(); qp++)
  {
    // Assemble the u-velocity row uu coupling
    for (unsigned int m=0; m<n_u_dofs; m++)
      for (unsigned int n=0; n<n_u_dofs; n++)
        Kij(m,n) += JxW[qp]*(dphi[m][qp](i)*dphi[n][qp](j))*mu; // with viscosity
  }
  
  STOP_LOG("assemble_element_KIJ()", "AssembleNS");  // libMesh log
  return;
}



// ==================================================================================
void AssembleNS::assemble_element_QI(const std::vector<Real>& JxW,
                                     const std::vector<std::vector<RealGradient> >& dphi,
                                     const std::vector<std::vector<Real> >& psi,
                                     const unsigned int n_v_dofs,
                                     const unsigned int n_p_dofs,
                                     const unsigned int I,
                                     DenseMatrix<Number>& Qi)
{
  START_LOG("assemble_element_QI()", "AssembleNS");  // libMesh log
  
  // reinit and compute the element matrix K_ij
  Qi.resize(n_v_dofs,n_p_dofs);
  for (unsigned int qp=0; qp<JxW.size(); qp++)
  {
    // Assemble the u-velocity row uu coupling
    for (unsigned int m=0; m<n_v_dofs; m++)
      for (unsigned int n=0; n<n_p_dofs; n++)
        Qi(m,n) += JxW[qp]*(dphi[m][qp](I)*psi[n][qp]);
  }
  
  STOP_LOG("assemble_element_QI()", "AssembleNS");  // libMesh log
}



// ==================================================================================
void AssembleNS::assemble_element_MIJ(const std::vector<Real>& JxW,
                                      const std::vector<std::vector<Real> >& phi,
                                      const unsigned int n_v_dofs,
                                      DenseMatrix<Number>& Mij)
{
  START_LOG("assemble_element_MIJ()", "AssembleNS");  // libMesh log
  
  // reinit and compute the element matrix K_ij
  Mij.resize(n_v_dofs, n_v_dofs);
  for (unsigned int qp=0; qp<JxW.size(); qp++)
  {
    // Assemble the u-velocity row uu coupling
    for (unsigned int m=0; m<n_v_dofs; m++)
      for (unsigned int n=0; n<n_v_dofs; n++)
        Mij(m,n) += JxW[qp]*(phi[m][qp]*phi[n][qp]); // with density
  }
  
  STOP_LOG("assemble_element_MIJ()", "AssembleNS");  // libMesh log
}

/// ======================================================================================= ///
void AssembleNS::apply_bc_by_penalty(const Elem* elem,
                                     const std::string& matrix_or_vector,
                                     DenseMatrix<Number>& Ke,
                                     DenseVector<Number>& Fe,
                                     const std::string& option)
{
  START_LOG("apply_bc_by_penalty()", "AssembleNS");  // libMesh log
  
  // Get a reference to the Particle-Mesh linear implicit system object.
  PMLinearImplicitSystem & pm_system = _eqn_sys.get_system<PMLinearImplicitSystem> ("Stokes");
//  const MeshBase& mesh   = pm_system.get_mesh();
  const Real penalty     = 1E6;    // The penalty value.
  const Real tol         = 1E-6;   // The tolerence value.
  const unsigned int n_nodes = elem->n_nodes();
  const std::vector<bool>& periodicity = pm_system.point_mesh()->pm_periodic_boundary()->periodic_direction();
  const std::vector<bool>& inlet_direction = pm_system.point_mesh()->pm_periodic_boundary()->inlet_direction();
  AnalyticalSolution analytical_solution(pm_system);
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   The following loops over the sides of the element. If the element has NO
   neighbors on a side then that side MUST live on a boundary of the domain.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  for (unsigned int s=0; s<elem->n_sides(); s++)
  {
    
    if (elem->neighbor(s) == NULL)
    { 
      bool _continue = false;
      for (int i = 0; i < _dim; i++)
      {
        if(periodicity[i] or inlet_direction[i]){
          if (_mesh.get_boundary_info().has_boundary_id(elem, s, _boundary_id_3D[2*i+0]) or
              _mesh.get_boundary_info().has_boundary_id(elem, s, _boundary_id_3D[2*i+1]))
          { _continue = true; }
        }
      }   
      // For any other periodic boundaries, which is not a wall, we do nothing
      if(_continue) continue;
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       // if this is neither the inlet nor the outlet, it MUST be no-slip walls
       // the following part -2- will impose the no-slip BC by penalty method.
       
       // -2.- build the full-order side element for "vel" Dirichlet BC at the walls.
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
      UniquePtr<Elem> side (elem->build_side(s));
      for (unsigned int ns=0; ns<side->n_nodes(); ns++)
      {
        // Otherwise if this is the wall of the channel, set u/v/w = 0
        std::vector<Real> uvw(_dim,0.0);
        
        // In GGEM, vel bc at the no-slip boundaries is u_bc = -u_local
        // In the case of "undisturbed" flow field, set u_bc = 0.
        // Note this only influence the rhs vector
        if(option=="disturbed")
        {
          const Point ptx = side->point(ns);
          const std::vector<Real> u_local = pm_system.local_velocity_fluid(ptx,"regularized");
          // ---------------- setup for validation test 01 -------------
          if(_eqn_sys.parameters.get<std::string> ("test_name") == "ggem_validation")
          {
            const std::vector<Real> u_boundary = analytical_solution.exact_solution_infinite_domain(ptx);
            for(unsigned int k=0;k<_dim;++k) uvw[k] = u_boundary[k]-u_local[k];
          }
          // ---------------- setup for other test -------------
          else
          {
            for(unsigned int k=0;k<_dim;++k) uvw[k] = -u_local[k];
          }
        }
        
        // Find the node on the element matching this node on the side.
        // That defined where in the element matrix the BC will be applied.
        for (unsigned int n=0; n<elem->n_nodes(); n++)
        {
          if (elem->node(n) == side->node(ns))
          {
            // Penalize u, v and w at the current node.
            for(unsigned int k=0;k<_dim;++k){
              this->penalize_elem_matrix_vector(Ke,Fe,matrix_or_vector,k,n,n_nodes,penalty,uvw[k]);
            } // end for k-loop
          } // end if (elem->node(n) == side->node(ns))
        } // enf for n-loop
        
        
      } // end for ns-loop
      
    } // end if (elem->neighbor(side) == NULL)
    
  } // end for s-loop
  
  STOP_LOG("apply_bc_by_penalty()", "AssembleNS");
} // end of function apply_bc_by_penalty()



/// ======================================================================================= ///
void AssembleNS::penalize_elem_matrix_vector(DenseMatrix<Number>& Ke,
                                             DenseVector<Number>& Fe,
                                             const std::string & matrix_or_vector,
                                             const unsigned int& var_number,
                                             const unsigned int& local_node_id,
                                             const unsigned int& n_nodes_elem, // total vel-node #!
                                             const Real& penalty,
                                             const Real& value)
{
  START_LOG("penalize_elem_matrix_vector()", "AssembleNS");
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Penalize Matrix or vector.
   n_nodes_elem: number of nodes in an element
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
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
    libmesh_assert("*** error in AssembleNS::penalize_elem_matrix_vector():");
    libmesh_assert("*** ---> invalid argument: matrix_or_vector!");
    libmesh_error();
  }
  
  STOP_LOG("penalize_elem_matrix_vector()", "AssembleNS");
}



// ============================================================================================
Real AssembleNS::boundary_pressure_jump(const std::string& which_side) const
{
  // set the pressure on the left boundary as a time dependent function
  if(which_side=="left")
  {
    return +0.00;
    //return 10.0*std::sin(20.0*pi*t);
  }
  // set the pressure on the right boundary as a constant
  else if(which_side=="right")
  {
    return +0.00;
  }
  else
  {
    std::cout <<"********* error in AssembleNS::boundary_pressure_jump(): *********"<<std::endl
              <<"********** the side can only be left or right boundary !**********"<<std::endl
              <<"******************************************************************"<<std::endl;
    return 1E100;
  } // end if-else
  
}



// ============================================================================================
Real AssembleNS::boundary_traction(const std::string& which_side) const
{
  // set the traction on the left boundary as a time dependent function
  if(which_side=="left")
  {
    return +1;
    //return 10.0*std::sin(20.0*pi*t);
  }
  // set the traction on the right boundary as a constant
  else if(which_side=="right")
  {
    return +1;
  }
  else
  {
    std::cout <<"************ error in AssembleNS::boundary_traction(): ***********"<<std::endl
              <<"********** the side can only be left or right boundary !**********"<<std::endl
              <<"******************************************************************"<<std::endl;
    return 1E100;
  } // end if-else
  
}






/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 -0.- If periodic BC is along all three directions, we need to introduce a
 constraint to avoid rigid motion
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
//      unsigned int periodic_count = 0;
//      if( x_periodic ) periodic_count++;
//      if( y_periodic ) periodic_count++;
//      if( z_periodic ) periodic_count++;
//      if(periodic_count==3)
//      {
//        // build side element and loop over its each node
//        UniquePtr<Elem> side(elem->build_side(s)); // side elem with the same order as the bulk elem
//        for (unsigned int ns=0; ns<side->n_nodes(); ns++)
//        {
//          // Check if this is the corner node
//          const Point& Pcorner = side->point(ns);
//          const Point  Pab = Pcorner - ref_point;
//          if(Pab.size()<tol)  // There may be no node at the ref_point!
//          {
//            // Because we don't know the local node # in this elem
//            // So we find it by the following loop
//            for (unsigned int n=0; n<elem->n_nodes(); n++) {
//              if (elem->node(n) == side->node(ns)) {
//                if(x_periodic)
//                  this->penalize_elem_matrix_vector(Ke,Fe,matrix_or_vector,0,n,n_nodes,penalty,0.0);
//                if(y_periodic)
//                  this->penalize_elem_matrix_vector(Ke,Fe,matrix_or_vector,1,n,n_nodes,penalty,0.0);
//                if(z_periodic)
//                  this->penalize_elem_matrix_vector(Ke,Fe,matrix_or_vector,2,n,n_nodes,penalty,0.0);
//              }
//            } // end for n-loop
//          } // end if (Pab.size()<tol)
//        } // end for ns-loop
//      } // end if(x_periodic && y_periodic && z_periodic)
//#include "libmesh/unique_ptr.hpp"
