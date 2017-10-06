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




#include <iomanip>
#include <stdio.h>
#include <iostream>

#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/boundary_info.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/auto_ptr.h"
//#include "libmesh/unique_ptr.hpp"

#include "pm_toolbox.h"
#include "rigid_particle.h"

  
namespace libMesh
{ 

// ======================================================================
RigidParticle::RigidParticle(const dof_id_type& particle_id,
                             const int& particle_type,
                             const Point& pt,
                             const Point& mag,
                             const Point& rot,
                             const Real& charge,
                             const Real& epsilon_in,
                             const Point& sedimentation_body_force_density,
                             const Parallel::Communicator &comm_in)
: ParallelObject(comm_in),
  _particle_id(particle_id),
  _particle_type(particle_type),
  _center0(pt),
  _radius(0.),
  _mag(mag),
  _rot(rot),
  _charge(charge), _epsilon_in(epsilon_in),
  _sedimentation_body_force_density(sedimentation_body_force_density),
  _volume0(0.), _area0(0.),
  _processor_id(-1),
  _mesh(comm_in),
  _mesh_spring_network(NULL)
{
}
 

// ======================================================================
RigidParticle::~RigidParticle()
{
  // do nothing
}


//=======================================================================
void RigidParticle::init()
{
  START_LOG("init()", "RigidParticle");
  this->compute_volume();
  this->compute_area();
  this->compute_centroid();
  _area0 = _area;
  _volume0 = _volume;
  _centroid0 = _centroid;
  STOP_LOG("init()", "RigidParticle");
}
  
// ======================================================================
void RigidParticle::extract_surface_mesh(const std::string& vmesh,
                                         const std::string& smesh) const
{
  START_LOG("extract_surface_mesh()", "RigidParticle");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   read the bulk mesh
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  SerialMesh bulk_mesh(this->comm());
  SerialMesh boundary_mesh(this->comm());
  bulk_mesh.read(vmesh);
  
  // print out bulk mesh info.
//  std::string msg = "--->test in RigidParticle::extract_surface_mesh(): print bulk mesh info:";
//  PMToolBox::output_message(msg, this->comm());
//  bulk_mesh.print_info();
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   sync the bulk and boundary mesh
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  bulk_mesh.get_boundary_info().sync(boundary_mesh);
  
  // print out surface mesh info:
//  msg = "--->test in RigidParticle::extract_surface_mesh(): print surface mesh info:";
//  PMToolBox::output_message(msg, this->comm());
//  boundary_mesh.print_info();
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write out the surface mesh
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  boundary_mesh.write(smesh);
  
  STOP_LOG("extract_surface_mesh()", "RigidParticle");
}
  
  
  
  
// ======================================================================
void RigidParticle::read_mesh(const std::string& filename,
                              const std::string& mesh_type)
{
  START_LOG("read_mesh()", "RigidParticle");
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   read the mesh file
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  _mesh_type = mesh_type;
  if(_particle_type == 0){
    this->read_mesh_sphere(filename);
  }
  else if(_particle_type == 1){
    this->read_mesh_cylinder(filename);
  }
  else{
    std::cout <<"Warning: unsupported rigid particle type: " <<_particle_type << "; only 1(sphere), 2(cylinder) is supported for now" <<std::endl; 
  }
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   * init mesh characteristics
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  this->init();
  STOP_LOG("read_mesh()", "RigidParticle");
}
  
  
  
// ======================================================================
void RigidParticle::read_mesh_sphere(const std::string& filename)
{
  START_LOG("read_mesh_sphere()", "RigidParticle");
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   read the mesh data file
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  _mesh.read(filename);   //_mesh.all_first_order();
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   change the coordinates of the surface mesh according to its radius
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  _radius = 1. * _mag(0);
  _dim    = _mesh.spatial_dimension();
  _mesh_dim = _mesh.mesh_dimension();
  MeshBase::node_iterator       nd     = _mesh.active_nodes_begin();
  const MeshBase::node_iterator end_nd = _mesh.active_nodes_end();
  for ( ; nd != end_nd; ++nd)
  {
    // Store a pointer to the element we are currently working on.
    Node* node = *nd;
    for(unsigned int i=0; i<_dim; ++i){
      (*node)(i) =  ((*node)(i)*_radius + _center0(i) );
    }
  } // end for
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Compute the volume of the sphere.
   FIXME: the volume for sphere can be computed analytically!
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  _particle_shape = "sphere";
  STOP_LOG("read_mesh_sphere()", "RigidParticle");
}
  

  
// ======================================================================
void RigidParticle::read_mesh_cylinder(const std::string& filename)
{
  START_LOG("read_mesh_cylinder()", "RigidParticle");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Check the existance of the file
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  const bool file_exist = PMToolBox::file_exist(filename);
  if( !file_exist )
  {
    printf("***error in RigidParticle::read_mesh_cylinder(): file does NOT exist!");
    libmesh_error();
  }
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   read the data file
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  _mesh.read(filename);
  _dim    = _mesh.spatial_dimension();
  _mesh_dim = _mesh.mesh_dimension();
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Magnify, rotate and shift
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  PMToolBox::magnify_serial_mesh(_mesh, _mag);
  PMToolBox::rotate_serial_mesh(_mesh, _rot);
  std::vector<Real> shift_dist(3);
  for(unsigned int i=0; i<3; ++i){
    shift_dist[i] =  _center0(i);
  }
  PMToolBox::shift_serial_mesh(_mesh, shift_dist);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Compute the volume of the cylinder
   FIXME: the volume for cylinder can be computed analytically!
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  _particle_shape = "cylinder";
  STOP_LOG("read_mesh_cylinder()", "RigidParticle");
}
  
  
  
// ======================================================================
void RigidParticle::write_mesh(const std::string& filename)
{
  START_LOG("write_mesh()", "RigidParticle");
  
  // write the mesh data file
  _mesh.write(filename);
  
  STOP_LOG("write_mesh()", "RigidParticle");
}

  
  
// ======================================================================
void RigidParticle::update_mesh(const std::vector<Point>& nodal_vec)
{
  START_LOG("update_mesh()", "RigidParticle");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   We update the position of nodes on the surface mesh. Algorithm I:
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  MeshBase::node_iterator       nd     = _mesh.active_nodes_begin();
  const MeshBase::node_iterator end_nd = _mesh.active_nodes_end();
  for ( ; nd != end_nd; ++nd)
  {
    
    // Store a pointer to the current node
    Node* node = *nd;
    const dof_id_type n_id = node->id();
    for(unsigned int i=0; i<_dim; ++i) {
      (*node)(i) = nodal_vec[n_id](i);
    }
  } // end for
  // update particle centroid
  this->compute_centroid();
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Algorithm II:
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
//  const std::size_t n_points = _mesh.n_nodes();
//  for (std::size_t i=0; i<n_points; ++i)
//  {
//    Node& node = _mesh.node(i);
//    for(unsigned int j=0; j<3; ++j)
//      node(j) = nodal_vec[i](j);
//  }
  
  STOP_LOG("update_mesh()", "RigidParticle");
}

// ======================================================================
void RigidParticle::translate_mesh(const Point& translationDist)
{
  START_LOG("translate_mesh()", "RigidParticle");
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   We apply a translation move to the position of nodes of particle surface mesh 
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  MeshBase::node_iterator       nd     = _mesh.active_nodes_begin();
  const MeshBase::node_iterator end_nd = _mesh.active_nodes_end();
  for ( ; nd != end_nd; ++nd)
  {
    // Store a pointer to the current node
    Node* node = *nd;
    for(unsigned int i=0; i<_dim; ++i) {
      (*node)(i) = translationDist(i);
    }
  } // end for
  _centroid += translationDist;
  
  STOP_LOG("translate_mesh()", "RigidParticle");
}
  


// ======================================================================
void RigidParticle::extract_nodes(std::vector<Point>& node_xyz)
{
  START_LOG("extract_nodes()", "RigidParticle");
  
  // init the vector
  node_xyz.resize(_mesh.n_nodes());
  
  // extract the nodal coordinates
  MeshBase::node_iterator       nd     = _mesh.active_nodes_begin();
  const MeshBase::node_iterator end_nd = _mesh.active_nodes_end();
  for ( ; nd != end_nd; ++nd)
  {
    // Store a pointer to the current node
    Node* node = *nd;
    for(unsigned int i=0; i<_dim; ++i)
    node_xyz[node->id()](i) =  (*node)(i) ;
  } // end for
  
  STOP_LOG("extract_nodes()", "RigidParticle");
}

  
  
// ======================================================================
SerialMesh& RigidParticle::mesh()
{
  return _mesh;
}
  
  
  
// ======================================================================
std::vector<Real> RigidParticle::mesh_size() const
{
  return PMToolBox::mesh_size(_mesh);
}

  
// ======================================================================
const Point& RigidParticle::mesh_point(const std::size_t i) const
{
  return _mesh.point(i);
}



// ======================================================================
std::size_t RigidParticle::num_mesh_nodes() const
{
  return _mesh.n_nodes();
}



// ======================================================================
std::size_t RigidParticle::num_mesh_elem() const
{
  return _mesh.n_elem();
}


// ======================================================================
bool RigidParticle::on_the_periodic_boundary()
{
  START_LOG("on_the_periodic_boundary()", "RigidParticle");
  // init flag
  bool flag = false;
  PMPeriodicBoundary* pbc_tmp = _mesh_spring_network->periodic_boundary();

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   First we check if there is periodic boundary condition
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  bool pbc  = false;
  for(std::size_t i=0; i<_dim; ++i)
  {  
    if( pbc_tmp->periodic_direction(i) )
    {
      pbc = true;  
      break;  // return pbc = true and break the for loop!
    }  
  }  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   If there is a periodic boundary, we will exam if this particle cuts the pb.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  if( pbc )
  {
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Loop over each element (operated on all the processors)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    MeshBase::const_element_iterator       el     = _mesh.active_elements_begin();
    const MeshBase::const_element_iterator end_el = _mesh.active_elements_end();
    for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently working on.
      const Elem* elem = *el;

      // check if this element is on the periodic boundary
      std::vector<bool> elem_flag = pbc_tmp->image_elem(elem);
      
      // If this element is on the pb, label the flag to be true
      for(std::size_t i=0; i<elem_flag.size(); ++i){
        if(elem_flag[i]){
          flag = true;  // modify flag value
          break;        // break i-loop
        }
      }
      
      // If flag==true, break the el-loop
      if(flag) break;
    } // end for el-loop
  } // end if(pbc)  
  STOP_LOG("on_the_periodic_boundary()", "RigidParticle");
  _on_pb = flag;
  return flag;
}


// ======================================================================
void RigidParticle::rebuild_periodic_mesh()
{
  START_LOG("rebuild_periodic_mesh()", "RigidParticle");
  
  //---- test : output the original mesh
  //ExodusII_IO(_mesh).write("test_rebuild_mesh0.e");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over all the nodes (operated on all the processors)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  MeshBase::node_iterator       nd     = _mesh.active_nodes_begin();
  const MeshBase::node_iterator end_nd = _mesh.active_nodes_end();
  for ( ; nd != end_nd; ++nd)
  {
    // Store a pointer to the node we are currently working on.
    Node* node = *nd;
    
    // Loop in each direction, and find the periodic boundary
    PMPeriodicBoundary* pbc = _mesh_spring_network->periodic_boundary();
    for(std::size_t i=0; i<_dim; ++i)
    {
      if( pbc->periodic_direction(i) )
      {
        const Real dist = (*node)(i) - pbc->box_min()(i);
        const Real Lc = pbc->box_length()(i);
        if( dist <= Lc/2. ){
          (*node)(i) += Lc;
        }
      } // end if
    } // end for i-loop
    
  } // end for nd-loop
  STOP_LOG("rebuild_periodic_mesh()", "RigidParticle");
}
  
  
  
// ======================================================================
void RigidParticle::restore_periodic_mesh()
{
  START_LOG("restore_periodic_mesh()", "RigidParticle");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over all the nodes (operated on all the processors)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  MeshBase::node_iterator       nd     = _mesh.active_nodes_begin();
  const MeshBase::node_iterator end_nd = _mesh.active_nodes_end();
  for ( ; nd != end_nd; ++nd)
  {
    // Store a pointer to the node we are currently working on.
    Node* node = *nd;
    
    // Loop in each direction, and find the periodic boundary
    PMPeriodicBoundary* pbc = _mesh_spring_network->periodic_boundary();
    pbc->correct_position(*node);
    
  } // end for nd-loop
  STOP_LOG("restore_periodic_mesh()", "RigidParticle");
}

// // ======================================================================
// void RigidParticle::add_particle_force(const std::vector<Real>& pforce)
// {
//   for(std::size_t i=0; i<pforce.size(); ++i)
//   {
//     _force[i] += pforce[i];
//   }
// }

void RigidParticle::add_body_force_density(Point& body_force_density) 
{
  START_LOG("RigidParticle::add_body_force_density()", "RigidParticle");
  _body_force_density += body_force_density;
  STOP_LOG("RigidParticle::add_body_force_density()", "RigidParticle");
}

void RigidParticle::add_sedimentation_body_force_density()
{
  START_LOG("RigidParticle::add_sedimentation_body_force_density()", "RigidParticle");
  _body_force_density += _sedimentation_body_force_density;
  STOP_LOG("RigidParticle::add_sedimentation_body_force_density()", "RigidParticle");
}
  
// ======================================================================
void RigidParticle::build_nodal_force(std::vector<Point>& nf)
{
  START_LOG("build_nodal_force()", "RigidParticle");
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   We will build the nodal force vector through the EquationSystems.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  _surface_force_density = _body_force_density * _volume0 / _area0;
  _centroid_force = _body_force_density * _volume0;
  EquationSystems equation_systems (_mesh);
  ExplicitSystem& system = equation_systems.add_system<ExplicitSystem>("Elasticity");

  unsigned int u_var = 0, v_var = 0, w_var = 0;
  u_var = system.add_variable ("U", FIRST);
  v_var = system.add_variable ("V", FIRST);
  if(_dim==3)  w_var  = system.add_variable ("W", FIRST);
  equation_systems.init();
  system.rhs->zero();
  // equation_systems.print_info();
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Build a FEBase object and qrule (typedef FEGenericBase<Real> FEBase)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  FEType fe_vel_type  = system.variable_type(u_var);
  UniquePtr<FEBase> fe_vel (FEBase::build(_mesh_dim, fe_vel_type) );
  
  // Gauss Quadrature rule and tell the finite element objects
  QGauss qrule (_mesh_dim, fe_vel_type.default_quadrature_order());   //3^dim pts
  fe_vel->attach_quadrature_rule (&qrule);
  
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>& JxW                        = fe_vel->get_JxW();
  const std::vector<std::vector<Real> >& phi          = fe_vel->get_phi();
  
  // A reference to the DofMap object for this system.
  const DofMap & dof_map = system.get_dof_map();
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u, dof_indices_v, dof_indices_w;
  
  // element nodal vector
  DenseVector<Number>     Fe;
  DenseSubVector<Number>  Fu(Fe),    Fv(Fe),    Fw(Fe);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   loop over all the elements and compute the element rhs vector Fe.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_local_elements_end();
  for ( ; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently working on.
    const Elem* elem = *el;
    
    // Get the degree of freedom indices for the current element.
    dof_map.dof_indices (elem, dof_indices);
    dof_map.dof_indices (elem, dof_indices_u, u_var);
    dof_map.dof_indices (elem, dof_indices_v, v_var);
    
    const unsigned int n_dofs   = dof_indices.size();
    const unsigned int n_u_dofs = dof_indices_u.size();
    const unsigned int n_v_dofs = dof_indices_v.size();
    unsigned int n_w_dofs = 0;
    if(_dim==3)
    {
      dof_map.dof_indices (elem, dof_indices_w, w_var);
      n_w_dofs = dof_indices_w.size();
    }
    
    // NOTE: JxW, phi/dphi and other element quantities are not computed until reinit().
    fe_vel->reinit (elem);
    
    // Zero the element right-hand-side vector before summing them.
    // DenseSubVector.reposition () member takes the (row_offset, row_size)
    Fe.resize (n_dofs);   // Ke.resize (n_dofs, n_dofs);
    Fu.reposition (u_var*n_u_dofs, n_u_dofs);
    Fv.reposition (v_var*n_u_dofs, n_v_dofs);
    if(_dim==3)    Fw.reposition (w_var*n_u_dofs, n_w_dofs);
    
    // Loop over each Gauss pt
    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {
      // Assemble the subvector
      for (unsigned int i=0; i<n_u_dofs; i++)
      {
        Fu(i) += JxW[qp]*phi[i][qp]*_surface_force_density(0);
        Fv(i) += JxW[qp]*phi[i][qp]*_surface_force_density(1);
        if(_dim==3) Fw(i) += JxW[qp]*phi[i][qp]*_surface_force_density(2);
      }
    }
    
    // apply any hanging node constraint equations.
    dof_map.constrain_element_vector (Fe, dof_indices);
    
    // add the element vector
    system.rhs->add_vector (Fe, dof_indices);
    
  } // end for el-loop
  system.rhs->close();
//  system.update();
//  ExodusII_IO(_mesh).write_equation_systems("test_es0.e",equation_systems);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   localize the vector
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  std::vector<Real> local_f;
  system.rhs->localize(local_f);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   convert it to the output format
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  nf.resize(_mesh.n_nodes());
  MeshBase::node_iterator       nd     = _mesh.active_nodes_begin();
  const MeshBase::node_iterator end_nd = _mesh.active_nodes_end();
  for ( ; nd != end_nd; ++nd)
  {
    // Store a pointer to the element we are currently working on.
    Node* node = *nd;
    const dof_id_type node_id = node->id();
  
    // get the dof numbers at this node
    for(unsigned int k=0; k<_dim; ++k){ // sys#; var(k) = 0, 1, 2; component=0
      const dof_id_type dof_num = node->dof_number(0,k,0);
      nf[node_id](k) = local_f[dof_num];
    }
  }
  _nf = nf;
  STOP_LOG("build_nodal_force()", "RigidParticle");
}


  
// ======================================================================
void RigidParticle::compute_volume()
{
  START_LOG("compute_volume()", "RigidParticle");  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Define FE and quadrature rules
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  FEType fe_type;         // default: first order and Lagrange element
  UniquePtr<FEBase>  fe_base (FEBase::build(_mesh_dim, fe_type) );
  QGauss qrule (_mesh_dim, fe_type.default_quadrature_order());
  fe_base->attach_quadrature_rule (&qrule);  
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>&  JxW        = fe_base->get_JxW();
  const std::vector<Point>& q_xyz      = fe_base->get_xyz(); // xyz of quad pts 
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over each element and compute the volume(area)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  Real val = 0.0;
  MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_local_elements_end();
  for ( ; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently working on.
    const Elem* elem = *el;
    
    // If this particle has volume mesh, just sum up the volumes of all elem
    if (_mesh_dim==3){
      val += elem->volume();
    }
    
    // If this particle has surface mesh, compute \int (x*n)
    else if(_mesh_dim==2)
    {
      fe_base->reinit(elem);
      
      // Get the unit normal value at the center of the element
      const Point center(-0.,-0.) ;// != elem->centroid();// ref coordinate
      const Point n = this->elem_surface_normal(*elem,center);
      
      Real elem_val = 0.;
      for(std::size_t i=0; i<JxW.size(); ++i){
        const Real n_dot_x = q_xyz[i](0)*n(0) + q_xyz[i](1)*n(1) + q_xyz[i](2)*n(2);
        elem_val += n_dot_x*JxW[i];
      }
      val += elem_val;
    }
    else {
      libmesh_error();
    } // end if-else
    
  } // end for el-loop
  
  // sum the values on all the processors
  this->comm().sum(val);
  _volume = val/3.; 

  STOP_LOG("compute_volume()", "RigidParticle");
}
  
  
  
  
// ======================================================================
void RigidParticle::compute_area()
{
  START_LOG("compute_area()", "RigidParticle");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over each element and compute the volume(area)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  Real val = 0.0;
  MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_local_elements_end();
  for ( ; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently working on.
    const Elem* elem = *el;
    
    // If this particle has surface mesh, just sum up the area of all elem
    if (_mesh_dim==2)
    {
      val += elem->volume();
    }
    
    // If this particle has volume mesh, only compute the area of its external side
    else if(_mesh_dim==3)
    {
      // loop over each side, and find the surface side.
      for (unsigned int s=0; s<elem->n_sides(); s++)
      {
        if (elem->neighbor(s) == NULL){
          UniquePtr<Elem>  elem_side( elem->build_side(s) );
          val += elem_side->volume();
        }
      } // end for s-loop
    }
    else
    {
      libmesh_error();
    }
    
  } // end for el-loop
  
  // sum the values on all the processors
  this->comm().sum(val);
  
  _area = val;  
  STOP_LOG("compute_area()", "RigidParticle");
}

  
  
// ======================================================================
void RigidParticle::compute_centroid()
{
  START_LOG("compute_centroid()", "RigidParticle");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over each element and compute the element volume(area) and center
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_local_elements_end();
  Real  val_ax = 0., val_ay = 0., val_az = 0.;
  Real  val_b = 0.;
  for ( ; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently working on.
    const Elem* elem = *el;
    const Real  elem_v = elem->volume();
    const Point elem_c = elem->centroid();
    val_ax += elem_c(0)*elem_v;
    val_ay += elem_c(1)*elem_v;
    val_az += elem_c(2)*elem_v;
    val_b += elem_v;
  } // end for el-loop
  
  // Sum across processors
  this->comm().sum(val_ax);
  this->comm().sum(val_ay);
  this->comm().sum(val_az);
  this->comm().sum(val_b);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Compute the centroid
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  _centroid(0) = val_ax/val_b;
  _centroid(1) = val_ay/val_b;
  _centroid(2) = val_az/val_b;
  
  STOP_LOG("compute_centroid()", "RigidParticle");
}
  
  
  
  
// ======================================================================
Point RigidParticle::elem_surface_normal(const Elem& elem,
                                         const Point& pt0) const
{
  START_LOG("elem_surface_normal()", "RigidParticle");
  
  // - - - - - - - - - - - - - - - - - - - - - - - - -
//  const Point& pt0 = elem->point(0);
//  const Point& pt1 = elem->point(1);
//  const Point& pt2 = elem->point(2);
//  const Point v1 = pt0 - pt1;
//  const Point v2 = pt2 - pt1;
//  elem_n(0) = v1(1)*v2(2) - v1(2)*v2(1);
//  elem_n(1) = v1(2)*v2(0) - v1(0)*v2(2);
//  elem_n(2) = v1(0)*v2(1) - v1(1)*v2(0);
//  elem_n /= -elem_n.size(); // normalize
  // - - - - - - - - - - - - - - - - - - - - - - - - -
  
  Point elem_n;
  const unsigned int     dim = 2;
  const unsigned int n_nodes = elem.n_nodes();
  FEType fe_type;   // default: first order and Lagrange element
  FE<dim,LAGRANGE>  fe_base(fe_type);  // 2D FE<Dim,FEFamily>
  
  // Compute derivatives of shape functions, and nodal coordinates
  std::vector<Point> node_xyz(n_nodes), dNdx(n_nodes);
  for (unsigned int i=0; i<n_nodes; ++i)
  {
    node_xyz[i] = elem.point(i);
    for (unsigned int j=0; j<dim; ++j){ // dN/dr and dN/ds
      dNdx[i](j) = fe_base.shape_deriv(&elem,fe_type.order,i,j,pt0);
    }
  }
  
  // dx/dr = sum(i=0:N)dN(i)*x(i)
  Point dxdr, dxds;
  for (unsigned int i=0; i<n_nodes; ++i){
    for (unsigned int j=0; j<3; ++j){
      dxdr(j) += dNdx[i](0)*node_xyz[i](j);
      dxds(j) += dNdx[i](1)*node_xyz[i](j);
    }
  }
  
  // n = dxdr X dxds
  elem_n(0) = dxdr(1)*dxds(2) - dxdr(2)*dxds(1);
  elem_n(1) = dxdr(2)*dxds(0) - dxdr(0)*dxds(2);
  elem_n(2) = dxdr(0)*dxds(1) - dxdr(1)*dxds(0);
  elem_n /= elem_n.norm(); // normalize
  
  
  STOP_LOG("elem_surface_normal()", "RigidParticle");
  return elem_n;
}
  
  
  
// ======================================================================
std::vector<Point> RigidParticle::compute_surface_normal()
{
  START_LOG("compute_surface_normal()", "RigidParticle");
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   check the input argument. Only allows surface mesh!
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  if(_mesh_type != "surface_mesh")
  {
    printf("--->error in RigidParticle::compute_surface_normal()\n");
    printf("    ONLY surface mesh is supported in this function!\n");
    libmesh_error();
  }
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   init the vector
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  std::vector<Point> s_normal( _mesh.n_nodes() );
  std::vector<std::size_t> node_count(_mesh.n_nodes());
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over each element (operated on all the processors)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  MeshBase::const_element_iterator       el     = _mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_elements_end();
  for ( ; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently working on.
    const Elem* elem = *el;
    const unsigned int n_nodes = elem->n_nodes();
    
    // compute the normal of this element
    const Point center(0.,0.);  // != elem->centroid();// ref coordinate
    const Point n = this->elem_surface_normal(*elem,center);
    
    // add the value to each node
    for (unsigned int i=0; i<n_nodes; ++i) {
      s_normal[elem->node(i)] += n;
      node_count[elem->node(i)]++;
    }
    
  } // end for el-loop
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Average nodal values
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  for(std::size_t i=0; i<_mesh.n_nodes(); ++i)
  {
    s_normal[i] /= Real(node_count[i]);
//    printf("--->TEST: pt = (%f,%f,%f), n = (%f,%f,%f), node_count = %lu\n",
//           _mesh.point(i)(0),_mesh.point(i)(1),_mesh.point(i)(2),
//           s_normal[i](0),s_normal[i](1),s_normal[i](2), node_count[i]);
  }
  
  STOP_LOG("compute_surface_normal()", "RigidParticle");
  return s_normal;
}
  
  
  
// ======================================================================
void RigidParticle::volume_conservation()
{
  START_LOG("volume_conservation()", "RigidParticle");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Check if this particle is on the periodic boundary. If so, move the nodes
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  const bool on_pb = this->on_the_periodic_boundary();
  if( on_pb ) {
    this->rebuild_periodic_mesh();
    this->comm().barrier();
  }
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Compute nodal position increment
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  std::vector<Point> n_vec = this->compute_surface_normal();
  this->compute_volume();
  this->compute_area();
  const Real aa = -3.*(_volume - _volume0) / _area; // 3 or 6?
  if(this->comm().rank()==0){
    printf("--->TEST in RigidParticle::volume_conservation(): volume0 = %f, volume = %f, area = %f, aa = %f\n",
           _volume0,_volume, _area0, aa);
  }
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Form the nodal position vector, and update the mesh coordinates.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  std::vector<Point> node_xyz;
  this->extract_nodes(node_xyz);
  for (std::size_t i=0; i<node_xyz.size(); ++i) {
//    printf("         %lu-th node  x0 = (%f,%f,%f), dx =  (%f,%f,%f)\n",
//           i,node_xyz[i](0),node_xyz[i](1),node_xyz[i](2),
//           aa*nodal_vec[i](0),aa*nodal_vec[i](1),aa*nodal_vec[i](2));
    
    node_xyz[i] += aa*n_vec[i];
  }
  this->update_mesh(node_xyz);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if this particle is on the periodic boundary, restore(move back) the nodes
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  if( on_pb ) this->restore_periodic_mesh();
  
  STOP_LOG("volume_conservation()", "RigidParticle");
}


// ======================================================================  
void RigidParticle::debug_body_force() const
{
  START_LOG("RigidParticle::debug_body_force()", "RigidParticle");
  std::cout << "============debug_body_force()================= " <<std::endl;
  std::cout << "==> 1) body force density = ";
  for (int i=0; i<_dim; i++){
    std::cout <<_body_force_density(i) <<"; ";
  }
  std::cout << "\n total force = _body_force_density * _volume0 = ";
  for (int i=0; i<_dim; i++){
    std::cout <<_centroid_force(i) << "; ";
  }
  std::cout << std::endl;

  std::cout << "==> 2) surface force density = _surface_force_density * _area0 = ";
  for (int i=0; i<_dim; i++){
    std::cout <<_surface_force_density(i) <<"; ";
  }
  std::cout << "\n total force = ";
  for (int i=0; i<_dim; i++){
    std::cout <<_surface_force_density(i) * _area0 << "; ";
  }
  std::cout << std::endl;


  Point total_force;
  for (int i =0; i<_nf.size(); i++) total_force += _nf[i];
  std::cout << "==> 3) total force by adding up node force = sum (_node_force) = ";
  for(int i=0; i<_dim; i++){
    std::cout << total_force(i) << ";";
  }    
  std::cout << std::endl;

  STOP_LOG("RigidParticle::debug_body_force()", "RigidParticle");
}



// ======================================================================
void RigidParticle::print_info() const
{
  /* --------------------------------------------------------------------------
   * Scheme 1: using printf. Then every process will print out info on its own.
   * --------------------------------------------------------------------------*/
      const std::vector<Real> hs = this->mesh_size();
      std::cout << "------>" << _particle_id <<"-th particle: " << std::endl
           << "                  particle_type = " << _particle_type << std::endl
           << "                  particle shape = "<<_particle_shape <<std::endl
           << "                  charge = " <<_charge <<std::endl
           << "                  dielectric constant = " <<_epsilon_in <<std::endl 
           << "                  center0 = "<< _center0(0) <<"," <<_center0(1) <<"," <<_center0(2) <<std::endl
           << "                  centroid0 = ("<< _centroid0(0) <<","<< _centroid0(1) <<","<< _centroid0(2) <<")" <<std::endl
           << "                  magnitude = " <<_mag(0) <<"," <<_mag(1)<<"," <<_mag(2) <<std::endl
           << "                  rotation = "<<_rot(0) <<"," <<_rot(1) <<"," <<_rot(2) <<std::endl
           << "                  area0   = "<< _area0 << std::endl
           << "                  volume0 = "<< _volume0 << std::endl
           << "                  hmin = "<< hs[0] <<", hmax = "<< hs[1] << std::endl
           << "                  number of elements = " << this->num_mesh_elem() << std::endl
           << "                  number of nodes = " << this->num_mesh_nodes() << std::endl;  
  // // output the neighbor list
  // if (print_neighbor_list)
  // {
  //   if( _neighbor_list.size()>0 )
  //   {
  //     printf("              the neighbor list includes: \n");
  //     for (std::size_t i=0; i<_neighbor_list.size(); ++i )
  //       printf("              ---particle %lu,   distance = %f\n",
  //              _neighbor_list[i].first, _neighbor_list[i].second);
  //   }
  //   else
  //     printf("              there are no neighbors around this particle!\n");
  // }
  // printf("\n");
  
} // the end of print_info()
  
} // end of namespace
  
