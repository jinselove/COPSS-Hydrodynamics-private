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
#include <utility>

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/parallel_object.h"
#include "libmesh/point.h"
#include "libmesh/elem.h"
#include "libmesh/mesh.h"


#include "pm_periodic_boundary.h"

namespace libMesh {
/*
 * This class is used to build the spring network
 * for a given (surface or volume) mesh
 *
 * This is only for one finite-sized particle, because
 * particles usually use the same surface/volume mesh,
 * but they may have different centers and diameters (or shapes).
 * In these cases, the spring network can be reused for the purpose
 * of saving memory and operation time, and it is not necessary
 * to build all the springs for each particle.
 *
 */


class MeshSpringNetwork : public ReferenceCountedObject<MeshSpringNetwork>,
                          public ParallelObject {
public:

  // Type to store the neighboring ids and the equilibrium distance(for spring)
  // For each node, there are typically more than one neighboring nodes.
  typedef std::vector<std::pair<std::size_t, Real> > neighbor_connection;


  // Constructor
  MeshSpringNetwork(MeshBase          & surface_mesh,
                    PMPeriodicBoundary& pbc);


  // ~ Destructor
  ~MeshSpringNetwork();


  /*
   * Build the spring network for all the nodes of the mesh.
   * This also include the spring connection between the nodes and its center,
   * which constrains the movement along the radial directions.
   *
   * This function typically is called only once at the beginning. It
   *initializes:
   *  - _nodes_neighbors
   *  - _node_center_equilibrium_dist
   */
  void build_spring_network(const Point& center);


  /*
   * Return the i-th node neighbors
   */
  const neighbor_connection& nodes_neighbors(const std::size_t i) const
  {
    return _nodes_neighbors[i];
  }

  /*
   * Return the (i-th) node-center equilibrium distance
   */
  const Real& node_center_equilibrium_dist(const std::size_t i) const
  {
    return _node_center_equilibrium_dist[i];
  }

  /*
   * Return pointer to the Periodic boundary condition
   */
  PMPeriodicBoundary* periodic_boundary() {
    return _periodic_boundary;
  }

  /*
   * print the MeshSpringNetwork info
   */
  void print_info();

private:

  /*
   * Find the nodal neighbors.
   * This function requires nodes_to_elem_map, which can be pre-computed by
   * MeshTools::build_nodes_to_elem_map ( mesh, nodes_to_elem_map );
   */
  void find_nodal_neighbors(const Node                                   & node,
                            const std::vector<std::vector<const Elem *> >& nodes_to_elem_map,
                            std::vector<const Node *>                    & neighbors)
  const;


  // Mesh base of the particle
  MeshBase& _p_mesh;

  // Pointer to the Periodic boundary condition
  PMPeriodicBoundary *_periodic_boundary;

  // node neighbors: neighbor id and equilibrium distance
  std::vector<neighbor_connection>_nodes_neighbors;

  // node-center equilibrium distance
  std::vector<Real>_node_center_equilibrium_dist;

  // nodal force vector of _p_mesh
  //  std::vector<Point> _nodal_force;
};
}
