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


#pragma once

#include <stdio.h>


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/parallel_object.h"
#include "libmesh/point.h"
#include "libmesh/id_types.h"
#include "libmesh/mesh.h"
#include "libmesh/serial_mesh.h"
#include "mesh_spring_network.h"


namespace libMesh {
/*
 * Shape of rigid particles
 */
enum ParticleShape {
  SPHERE    = 0,
  CYLINDER  = 1,
  ELLIPSOID = 2,
  CUBE      = 3,
  GENERAL   = 100
};


/*
 * This class only defines the rigid particle with general geometries
 */


class RigidParticle : public ReferenceCountedObject<RigidParticle>,
                      public ParallelObject {
public:

  // Constructor for a spherical particle with radius, density, total charge,
  // and relative permittivity
  RigidParticle(const dof_id_type           & particle_id,
                const int                   & particle_type,
                const Point                 & pt,
                const Point                 & mag,
                const Point                 & rot,
                const Real                  & charge,
                const Real                  & epsilon_in,
                const Point                 & sedimentation_body_force_density,
                const Parallel::Communicator& comm_in);

  // ~ Destructor
  ~RigidParticle();


  /*
   * Initialize particle volume0, volume, area0, area, centroid0, centroid, etc.
   */
  void              init();

  /*
   * The counter that counts how many times a particle has crossed boundaries
   */
  std::vector<int>& counter() {
    return _counter;
  }

  /*
   * Return particle id, which is unique for a particle.
   * This id is set during the initialization,and not allowed to reset
   *afterwards.
   */
  dof_id_type id() const {
    return _particle_id;
  }

  /*
   * particle radius.
   * This parameter is only for spherical particles,
   * but meaningless for particles with general shapes!
   */
  Real radius() const {
    return _radius;
  }

  /*
   * free charge carried by the particle
   */
  const Real& charge() const  {
    return _charge;
  }

  /*
   * relative permittivity of the particle
   */
  const Real& epsilon_in() const {
    return _epsilon_in;
  }

  /*
   * The processor that the particle belongs to, which is the same as the
   * processor id of its hosting element;
   *
   * For a finite sized particle whose parts can be on different processors,
   * this becomes meaningless. So it is here set as -1 at this time
   */
  int processor_id() const {
    return _processor_id;
  }

  /*
   * (re-)set the processor id for the particle
   */
  void set_processor_id(const int pid) {
    _processor_id = pid;
  }

  /*
   * The mesh type of the particle
   */
  const std::string& mesh_type() const {
    return _mesh_type;
  }

  /*
   * Set the neighbor list of the particle
   * This is set by the member function in the class "ParticleMesh"
   */
  void set_neighbor_list(const std::vector<std::pair<std::size_t,
                                                     Real> >& nei_list)
  {
    _neighbor_list = nei_list;
  }

  /*
   * Return the neighbor list of the particle.
   * NOTE, this includes the particle ids and distance values to this particle.
   */
  std::vector<std::pair<std::size_t, Real> >neighbor_list() const
  {
    return _neighbor_list;
  }

  /*
   * Attach MeshSpringNetwork & return the pointer
   */
  void attach_mesh_spring_network(MeshSpringNetwork *msn)
  {
    _mesh_spring_network = msn;
  }

  MeshSpringNetwork* mesh_spring_network() const {
    return _mesh_spring_network;
  }

  /*
   * This extract the surface mesh from a volume mesh,
   * and write the surface mesh to the local file
   */
  void extract_surface_mesh(const std::string& vmesh,
                            const std::string& smesh) const;


  /*
   * Read the mesh data for the current particle without any modification,
   * Then compute the particle's volume according to the mesh.
   * It can be either a surface mesh of a volume mesh!
   */
  void read_mesh(const std::string& filename,
                 const std::string& mesh_type);


  /*
   * Read the mesh data for a spherical particle.
   * and modify the mesh size according to the radius value
   * and center position.
   */
  void read_mesh_sphere(const std::string& filename);


  /*
   * Read the mesh data for a cylindrical particle.
   * and modify the mesh size according to the radius value
   * height value, and center position.
   * (assume the original size of the particle is : r=0.5; h=1;
   * center = (0,0,0) and aligned along z-direction )
   *
   * Note this can also be used to read other particles, e.g. spheres
   */
  void read_mesh_cylinder(const std::string& filename);

  /*
   * Write the mesh data for the current particle.
   */
  void write_mesh(const std::string& filename);


  /*
   * Update the particle position, which is achieved
   * by updating the coordinates of each node on the mesh.
   */
  void update_mesh(const std::vector<Point>& nodal_pos,
                   const std::vector<Point>& nodal_vel);


  /*
   * Update the particle position, which is achieved
   * by updating the coordinates of each node on the mesh.
   */
  void update_mesh(const std::vector<Point>& nodal_pos);


  /*
   * Rebuild mesh when rigidparticle is forced to make a translation move
   *
   */
  void translate_mesh(const Point& translationDist);

  /*
   * Extract nodal coordinates of the particle mesh.
   * Note, this only extract the nodal values, but not change them!
   */
  void extract_nodes(std::vector<Point>& node_xyz);


  /*
   * Return the mesh associated with this particle
   */
  SerialMesh& mesh();


  /*
   * Return the mesh size (hmin/hmax) associated with this particle
   */
  std::vector<Real>mesh_size() const;


  /*
   * Return the coordinate of the i-th node (tracking point)
   */
  const Point& mesh_point(const std::size_t i) const;

  /*
   * Return node to center distance
   */
  const Real   node_center_dist(const std::size_t i) const;


  /*
   * total number of nodes of the mesh
   */
  std::size_t num_mesh_nodes() const;


  /*
   * total number of elements of the mesh
   */
  std::size_t num_mesh_elem() const;


  /*
   * Check if this particle is sitting on the periodic boundary
   */
  bool        on_the_periodic_boundary();

  /*
   * get the information that if a particle is on the periodic boundary, this
   *function is used in FixRigid::check_pbc_post_fix()
   */
  const bool& if_on_the_periodic_boundary() const {
    return _on_pb;
  }

  /*
   * When a particle is sitting on the periodic boundary, we need to
   * rebuild the periodic mesh by moving nodes from one side to the other.
   *
   * This is necessary to correctly compute the particle quantities:
   * volume/area/center/normal ...
   */
  void rebuild_periodic_mesh();


  /*
   * Restore the periodic mesh: moving nodes back to the right positions.
   * This is usually used together with rebuild_periodic_mesh()
   */
  void restore_periodic_mesh();


  /*
   * This function builds the nodal force vector from a constant force density
   *\f
   * f = (fx,fy,fz) => nf = (f1x,f1y,f1z; f2x,f2y,f2z; ... fnx,fny,fnz;)
   *
   * It is an area density for surface mesh, and volume density of volume mesh!
   * NOTE: we don't use "()" const for constructing es
   */
  void build_nodal_sedimentation_force();

  /*
   * This function adds forces to the center of mass of rigid particles
   * These forces will be distributed over the nodal points
   */

  // void add_particle_force(const std::vector<Real>& pforce);


  /*
   * compute the volume of the particle.
   * The algorithm is different for surface mesh and volume mesh.
   */
  void compute_volume();


  /*
   * compute the area of the particle.
   * The algorithm is different for surface mesh and volume mesh.
   */
  void compute_area();


  /*
   * compute the centroid of the particle.
   *       \int x*dV        sum Ci*Vi
   * xc =  ----------  =  -------------
   *        \int dV          sum Vi
   *
   * This function doesn't care what mesh type is used.
   */
  void compute_centroid();


  /*
   * Get volume0
   */
  const Real& get_volume0() const {
    return _volume0;
  }

  /*
   * Get volume
   */
  const Real& get_volume() const {
    return _volume;
  }

  /*
   * Get area0
   */
  const Real& get_area0() const {
    return _area0;
  }

  /*
   * Get area
   */
  const Real& get_area() const {
    return _area;
  }

  /*
   * Get center0
   */
  const Point& get_center0() const {
    return _center0;
  }

  /*
   * Get centroid
   */
  const Point& get_centroid0() {
    return _centroid0;
  }

  /*
   * Get centroid0
   */
  const Point& get_centroid() {
    return _centroid;
  }

  /*
   * Get particle shape
   */
  const std::string& get_particle_shape() {
    return _particle_shape;
  }

  // this function has not implemented
  // Point compute_centroid(const std::string& mesh_type = "surface_mesh"){};

  /*
   * Construct the unit surface normal at point pt0 for a surface element.
   * NOTE: this is only for surface mesh.
   */
  Point elem_surface_normal(const Elem & s_elem,
                            const Point& pt0) const;


  /*
   * Construct the unit surface normal
   * NOTE: this is only for surface mesh.
   */
  std::vector<Point>compute_surface_normal();


  /*
   * Correct the position of tracking points to conserve the volume!
   * NOTE: this is only for surface mesh.
   */
  void volume_conservation();


  /*
   * Print information of this particle
   */
  void print_info() const;

  /*
   * set node velocity
   */
  void set_node_velocity(std::vector<Point>& node_velocity) {
    _node_velocity = node_velocity;
  }

  /*! \brief force on surface nodes
   *
   */
  void zero_node_force();

  /*! \brief add force to surface nodes
   *
   */
  void add_node_force(std::vector<Point>& node_force);

  /*! \brief get force on surface nodes
   *
   */
  void get_node_force(std::vector<Point>& node_force) const {
    node_force = _node_force;
  }

  /*! \brief compute body velocity
   * each node of this rigid particle does not have member velocity
   * so we can only calculate centroid_velocity in point_mesh.C and assign it
   * to this rigid particle
   */
  void         compute_centroid_velocity();

  /*! \brief get centroid velocity
   *
   */
  const Point& centroid_velocity() {
    return _centroid_velocity;
  }

  /*! \brief compute body force
   *
   */
  void         compute_centroid_force();

  /*! get centroid force
   *
   */
  const Point& centroid_force() {
    return _centroid_force;
  }

private:

  // particle id
  dof_id_type _particle_id;

  // particle type
  int _particle_type;

  // number of surface nodes
  std::size_t _num_mesh_node;

  // number of surface element
  std::size_t _num_mesh_elem;

  // mesh size
  std::vector<Real>_mesh_size;

  // dim
  std::size_t _dim;

  // mesh dim
  std::size_t _mesh_dim;

  // how many time the particle has crossed the boundary
  std::vector<int>_counter;

  // radius of particle (for spherical particles only)
  Real _radius;

  // free charge of particle
  Real _charge;

  // relative permittivity of particle
  Real _epsilon_in;

  // the processor that the particle belongs to
  // *** Currently, this is not used for the finite size particles.
  int _processor_id;

  // neighbor particles around the present particle: particle id and distance
  // value.
  std::vector<std::pair<std::size_t, Real> >_neighbor_list;

  // mesh of the particle, which can be eigther surface mesh or volume mesh
  SerialMesh _mesh;

  // Type of the particle's mesh: surface_mesh or volume_mesh
  std::string _mesh_type;


  // Mesh - spring network
  MeshSpringNetwork *_mesh_spring_network = nullptr;

  // Initial values (if it is the perfect shape)
  Real _volume0;
  Real _area0;
  Point _center0;
  Point _centroid0;


  // discretized values
  Real _volume;
  Real _area;
  Point _centroid;

  // particle shape
  std::string _particle_shape;

  // particle rotation
  Point _rot;

  // particle magnitude
  Point _mag;

  // check if particle is on the periodic boundary
  bool _on_pb;

  // // body force density
  // Point _body_force_density;

  // // surface force density
  // Point _surface_force_density;
  // sedimentation force
  Point _sedimentation_force;

  // sedimentation body force density
  Point _body_sedimentation_force_density;

  // sedimentation surface force density
  Point _surface_sedimentation_force_density;

  // forces on each nodes
  std::vector<Point>_node_force;

  // velocity of each nodes
  std::vector<Point>_node_velocity;

  // velocity of the centroid (mean over all node velocity)
  Point _centroid_velocity;

  // force if the centroid (mean over all node force)
  Point _centroid_force;
}; // end of class defination
}  // end of namespace
