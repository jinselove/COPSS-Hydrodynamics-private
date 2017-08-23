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



namespace libMesh
{
  
  /*
   * Shape of rigid particles
   */
  enum ParticleShape {
    SPHERE    = 0,
    ELLIPSOID = 1,
    CYLINDER  = 2,
    CUBE      = 3,
    GENERAL   = 100
  };
  
  
  
  /*
   * This class only defines the rigid particle with general geometries
   */
  
  
class RigidParticle : public ReferenceCountedObject<RigidParticle>,
                      public ParallelObject
{
public:
  
  // Constructor for a spherical particle with radius and (relative) density
  RigidParticle(const Point pt,
                const dof_id_type particle_id,
                const Real r,
                const Real density,
                const Parallel::Communicator &comm_in);

  
  // A simple Constructor for general shaped particles without r (general geometry)
  RigidParticle(const Point pt,
                const dof_id_type particle_id,
                const Real density,
                const Parallel::Communicator &comm_in);
  
  
  // A simpler Constructor for general shaped particles without r and density.
  RigidParticle(const Point pt,                         // center
                const dof_id_type particle_id,          // id
                const Parallel::Communicator &comm_in); // comm
  
  
  // Constructor for general shaped particles from a given mesh.
  // NOTE: we don't use const SerialMesh& for constructing es
  RigidParticle(SerialMesh& pmesh,                // particle mesh
                const std::string& mesh_type,
                const dof_id_type particle_id,          // id
                const Parallel::Communicator &comm_in); // comm
  
  
  // Constructor for a spherical particle with radius, density, total charge, and relative permittivity
  RigidParticle(const Point pt,
                const dof_id_type point_id,
                const Real r,
                const Real density,
                const Real charge,
                const Real epsilon_in,
                const Parallel::Communicator &comm_in);
  
  // ~ Destructor
  ~RigidParticle();
  
  
  /*
   * The coordinate of the particle center which is a writeable reference
   * Note that this value will be changed when the particle moves
   */
  //Point& point()  {  return _center;  };
  Point& center() {  return _center;  };
  
  
  /*
   * Return particle id, which is unique for a particle.
   * This id is set during the initialization,and not allowed to reset afterwards.
   */
  dof_id_type id() const {  return _id; };
  
  
  /*
   * particle radius.
   * This parameter is only for spherical particles,
   * but meaningless for particles with general shapes!
   */
  Real radius() const { return _radius; }
  
  
  /*
   * density, which may be used for inertial flow.
   */
  Real density() const { return _density; }


  /*
   * free charge carried by the particle 
   */
  const Real& charge() const  { return _charge; }


  /*
   * relative permittivity of the particle
   */
  const Real& epsilon_in() const { return _epsilon_in; }
  
  
  /*
   * The processor that the particle belongs to, which is the same as the
   * processor id of its hosting element;
   *
   * For a finite sized particle whose parts can be on different processors,
   * this becomes meaningless. So it is here set as -1 at this time
   */
  int processor_id() const {  return _processor_id;  };
  
  
  /*
   * (re-)set the processor id for the particle
   */
  void set_processor_id(const int pid) { _processor_id = pid; };
  
  
  /*
   * The mesh type of the particle
   */
  const std::string& mesh_type() const { return _mesh_type; };
  
  
  
  /*
   * Set the neighbor list of the particle
   * This is set by the member function in the class "ParticleMesh"
   */
  void set_neighbor_list(const std::vector<std::pair<std::size_t,Real> >& nei_list)
  { _neighbor_list = nei_list;  };
  
  
  /*
   * Return the neighbor list of the particle.
   * NOTE, this includes the particle ids and distance values to this particle.
   */
  std::vector<std::pair<std::size_t,Real> > neighbor_list() const
  { return _neighbor_list;  };
  
  
  /*
   * Set the force and torque vectors on the particle
   * This is set by the member function in the class "ParticleMesh"
   */
  void set_particle_force(const std::vector<Real>& pforce){ _force = pforce; };
  
  
  /*
   * Return the force and torque vector on the particle
   */
  const std::vector<Real>& particle_force() const {  return _force;  };
  
  
  
  /*
   * Attach MeshSpringNetwork & return the pointer
   */
  void attach_mesh_spring_network(MeshSpringNetwork* msn)
  { _mesh_spring_network = msn;  }
  
  MeshSpringNetwork* mesh_spring_network() const {  return _mesh_spring_network;  }
  
  
  
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
  void read_mesh_sphere(const std::string& filename,
                        const std::string& mesh_type);
  
  
  /*
   * Read the mesh data for a cylindrical particle.
   * and modify the mesh size according to the radius value
   * height value, and center position.
   * (assume the original size of the particle is : r=0.5; h=1;
   * center = (0,0,0) and aligned along z-direction )
   *
   * Note this can also be used to read other particles, e.g. spheres
   */
  void read_mesh_cylinder(const std::string& filename,
                          const std::string& mesh_type,
                          const std::vector<Real>& mag_factor,
                          const std::vector<Real>& angles);
  
  
  /*
   * Write the mesh data for the current particle.
   */
  void write_mesh(const std::string& filename);

  
  /*
   * Update the particle position, which is achieved
   * by updating the coordinates of each node on the mesh.
   */
  void update_mesh(const std::vector<Point>& nodal_vec);
  
  
  /*
   * Extract nodal coordinates of the particle mesh.
   * Note, this only extract the nodal values, but not change them!
   */
  void extract_nodes(std::vector<Point>& node_xyz) ;
  

  /*
   * Return the mesh associated with this particle
   */
  SerialMesh& mesh();
  
  
  /*
   * Return the mesh size (hmin/hmax) associated with this particle
   */
  std::vector<Real> mesh_size() const;
  
  
  /*
   * Return the coordinate of the i-th node (tracking point)
   */
  const Point& mesh_point(const std::size_t i) const;
  
  
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
  bool on_the_periodic_boundary() const;
  
  
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
   * This function builds the nodal force vector from a constant force density \f
   * f = (fx,fy,fz) => nf = (f1x,f1y,f1z; f2x,f2y,f2z; ... fnx,fny,fnz;)
   *
   * It is an area density for surface mesh, and volume density of volume mesh!
   * NOTE: we don't use "()" const for constructing es
   */
  void build_nodal_force(const std::vector<Real>& f,
                         std::vector<Point>& nf);
  
  /*
   * This function adds forces to the center of mass of rigid particles
   * These forces will be distributed over the nodal points 
   */
  void add_particle_force(const std::vector<Real>& pforce);

  
  
  /*
   * compute the volume of the particle.
   * The algorithm is different for surface mesh and volume mesh.
   */
  Real compute_volume();
  
  
  /*
   * compute the area of the particle.
   * The algorithm is different for surface mesh and volume mesh.
   */
  Real compute_area();
  
  
  /*
   * compute the centroid of the particle.
   *       \int x*dV        sum Ci*Vi
   * xc =  ----------  =  -------------
   *        \int dV          sum Vi
   *
   * This function doesn't care what mesh type is used.
   */
  Point compute_centroid(); 
 

  // this function has not implemented
  //Point compute_centroid(const std::string& mesh_type = "surface_mesh"){};

   /*
   * Construct the unit surface normal at point pt0 for a surface element.
   * NOTE: this is only for surface mesh.
   */
  Point elem_surface_normal(const Elem& s_elem,
                            const Point& pt0) const;
  
  
  /*
   * Construct the unit surface normal
   * NOTE: this is only for surface mesh.
   */
  std::vector<Point> compute_surface_normal(const std::string& mesh_type);
  
  
  /*
   * Correct the position of tracking points to conserve the volume!
   * NOTE: this is only for surface mesh.
   */
  void volume_conservation(const std::string& mesh_type);
  
  
  /*
   * Print information of this particle
   */
  void print_info(const bool & print_neighbor_list = true) const;
  
  
  
private:
  
  // The coordinate of the particle center
  Point _center;
  
  // particle id
  dof_id_type _id;
  
  // radius of particle (for spherical particles only)
  Real _radius;
  
  // (relative) density of the particle rel_den = density_particle/density_fluid
  Real _density;

  // free charge of particle, Xikai
  Real _charge;

  // relative permittivity of particle, Xikai
  Real _epsilon_in;
  
  // Initial volume of the particle.
  // This is used to monitor the volume change of the particle, and correct the volume value.
  Real _volume0;
  
  // the processor that the particle belongs to
  // *** Currently, this is not used for the finite size particles.
  int _processor_id;
  
  // neighbor particles around the present particle: particle id and distance value.
  std::vector<std::pair<std::size_t,Real> > _neighbor_list;
  
  // the force vector excerted on this particle(non-hydrodynamic and non-Brownian)
  std::vector<Real> _force;
  
  // mesh of the particle, which can be eigther surface mesh or volume mesh
  SerialMesh _mesh;
  
  // Type of the particle's mesh: surface_mesh or volume_mesh
  std::string _mesh_type;
  
  
  // Mesh - spring network
  MeshSpringNetwork* _mesh_spring_network;

}; // end of class defination
  
  
  
}  // end of namespace




