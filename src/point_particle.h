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

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/point.h"
#include "libmesh/id_types.h"



//
namespace libMesh
{
  
  /*
   * Type of a point
   */
  enum PointType {
    POLYMER_BEAD      = 0,
    LAGRANGIAN_POINT  = 1,
    POINT_PARTICLE    = 2,
    NOT_DEFINED       = 100
  };
  
  
  /*
   * This basic class only defines the point-type particle,
   * which is only a point without radius.
   * It is also the key component to define other types of structures
   */
  
  
class PointParticle : public ReferenceCountedObject<PointParticle>
{
public:
  
  // Constructor
  PointParticle(const Point pt,
                const dof_id_type point_id);
  
  
  // Constructor
  PointParticle(const Point pt,
                const dof_id_type point_id,
                const PointType point_type);
  
  // Constructor
  PointParticle(const Point pt,
                const dof_id_type point_id,
                const PointType point_type,
                const std::vector<Real>& rot_vec);
  
  
  // Copy constructor
  PointParticle(const PointParticle& particle);
  
  
  // ~ Destructor
  ~PointParticle();
  
  
  /*
   * The coordinate of the particle center(point) which is set writeable!
   */
  Point& point()  {  return _center;  };
  Point& center() {  return _center;  };

 

  /*
   * The counter that counts how many times a point has crossed boundaries
   */
  std::vector<int> & counter() {  return _counter; };


  /*
   * Set values for counter
   */
  void set_pbc_counter(int nx, int ny, int nz) { _counter[0]=nx; _counter[1]=ny; _counter[2]=nz; };

 
  /*
   * Return particle id, which is an unique id for a given particle.
   * This id is set during the initialization,and NOT allowed to be changed.
   */
  dof_id_type id() const {  return _id; };
  
  
  /*
   * Get and set parent id.
   * A point may be part of another upper level object, for example,
   * it may be a bead of a chain, or a node of a solid particle, which
   * may also have its id.
   *
   * It is also possible that this point doesn't have a parent. In this case,
   * this parent_id is set -1 (by default)
   */
  int parent_id() const { return _parent_id; };
  void  set_parent_id(const int pid) { _parent_id = pid; };
  
  
  /*
   * Get the point type: 0 - polymer point; 1 - tracking point; ...
   * This can be either a point on immersed bodies or a bead on the polymer chains.
   *
   * Note that points with the same type may have different parent id. For example,
   * They can be bead type points on two different polymer chains, or tracking points
   * on two different immersed bodies, whose parent ids are different.
   */
  PointType point_type() const {   return _point_type; };
  void set_point_type(const PointType p_type) { _point_type = p_type; };
  
  
  /*
   * The processor that the particle belongs to, which is the same as the
   * processor id of its hosting element.
   */
  int processor_id() const {  return _processor_id;  };
  void set_processor_id(const int pid) { _processor_id = pid; };
  
  
  /*
   * Elem id that the particle resides in.
   * Note this id will change when a particle moves from one element to another.
   * Therefore, it needs to be reset in this situation.
   */
  dof_id_type elem_id() const { return _elem_id;  };
  void set_elem_id(const dof_id_type e_id) {  _elem_id = e_id;  };
  
  
  /*
   * Set the neighbor list of the particle
   * This is used by the member function in the class "ParticleMesh"
   */
  void set_neighbor_list(const std::vector<std::pair<std::size_t,Real> >& nei_list)
  {
    if (nei_list.size() > 0) {
      _neighbor_list = nei_list;
    }
  };
  
  
  /*
   * Return the neighbor list of the particle.
   * NOTE, this includes the particle ids and distance values to this particle.
   */
  std::vector<std::pair<std::size_t,Real> > neighbor_list() const
  { return _neighbor_list;  };
  

  /*
   * Reinit neighbor list distances
   */ 
  void reinit_neighbor_distance();
  
  /*
   * Set the force vector on the particle
   * This is set by the member function in the class "ParticleMesh"
   */
  void set_particle_force(const Point& pforce);
  

  /*
   * Add the force vector
   */
  void add_particle_force(const Point& pforce);
  
  
  /*
   * set the force vector equal to zeros
   */
  void zero_particle_force();

  /*!\ brief particle force
   *
   */  
  const Point& particle_force() {return _force;};

  /*
   * Set particle velocity
   */ 
  void set_particle_velocity(const Point& p_velocity) {_velocity = p_velocity;};

  /*
   * Return the velocity vector of the particle
   */
  const Point& particle_velocity() const {return _velocity; };


  
  /*
   * Return the orientation vector of the particle
   */
  const std::vector<Real> orientation() const {  return _orientation;  };
  void set_orientation(const std::vector<Real>& rot_vec);
  
  
  
  /*
   * Re-init the point particle quantities:
   *  clean neighbor_list
   * _processor_id = -1;
   * _elem_id      = -1
   * The following members are NOT changed during the reinit:
   * (1)_center; (2) force; (3)_counter; (4)_id; (5)_parent_id; (6)_point_type; (7)_orientation; 
   */
  void reinit_particle();


  /*
   * Print information of this particle
   */
  void print_info(const bool & print_neighbor_list = true) const;
  
  
  
private:
  
  // The coordinate of the particle center
  Point _center;

  // The velocity of the particle center
  Point _velocity;

  // the force vector excerted on this particle(non-hydrodynamic and non-Brownian)
  Point _force;
  

  // Count how many times the point has crossed a boundary
  // Used to unfold point's coordinate with periodic boundaries
  std::vector<int> _counter;
  
  // particle id
  dof_id_type _id;
  
  
  // type of the point
  PointType _point_type;
  
  
  // id of its parent that it belongs to
  // Its parents can be either an immersed body or a polymer chain.
  // Knowing the parent's property will facilitate the computation of force field.
  // -1 by default, which means it has no parent.
  int _parent_id;
  
  
  // the processor that the particle belongs to
  int _processor_id;
  
  
  // Elem id that the particle reside in
  dof_id_type _elem_id; //unsigned int _elem_id (uint32 by default!);
  
  
  // neighbor particles around the present particle: particle id and distance value.
  std::vector<std::pair<std::size_t,Real> > _neighbor_list;
  

  
  // Define the orientation of this point(bending, torque)
  // Use the definition of quaternions to describe the orientation
  // as specified in lammps(http://lammps.sandia.gov/doc/set.html)
  std::vector<Real> _orientation;
};
  
  
} // end of namespace
