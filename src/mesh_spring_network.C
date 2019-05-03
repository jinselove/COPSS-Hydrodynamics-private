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


#include "libmesh/mesh_tools.h"


#include "mesh_spring_network.h"


namespace libMesh {
// ======================================================================
MeshSpringNetwork::MeshSpringNetwork(MeshBase          & p_mesh,
                                     PMPeriodicBoundary& pbc)
  : ParallelObject(p_mesh),
  _p_mesh(p_mesh),
  _periodic_boundary(&pbc)
{}

// ======================================================================
MeshSpringNetwork::~MeshSpringNetwork()
{
  // do nothing
}

// ======================================================================
void MeshSpringNetwork::build_spring_network(const Point& center)
{
  START_LOG("build_spring_network()", "MeshSpringNetwork");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     init the variables, and build the nodes to element map
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  const std::size_t n_nodes = _p_mesh.n_nodes();
  _nodes_neighbors.resize(n_nodes);
  _node_center_equilibrium_dist.resize(n_nodes);
  std::vector<std::vector<const Elem *> > nodes_to_elem_map;
  MeshTools::build_nodes_to_elem_map(_p_mesh, nodes_to_elem_map);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Loop over each node, and built connected springs between neighboring nodes
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  MeshBase::node_iterator nd           = _p_mesh.active_nodes_begin();
  const MeshBase::node_iterator end_nd = _p_mesh.active_nodes_end();

  for (; nd != end_nd; ++nd)
  {
    // node pointer -> node id, node point, node neighbors
    Node *node            = *nd;
    const std::size_t nid = node->id();         // node id
    const Point pt0       = _p_mesh.point(nid); // node coord
    std::vector<const Node *> neighbors;        // node's neighboring nodes
    this->find_nodal_neighbors(*node, nodes_to_elem_map, neighbors);

    // build connections between the neighboring nodes
    const std::size_t n_neighbors = neighbors.size();
    _nodes_neighbors[nid].resize(n_neighbors);

    for (std::size_t i = 0; i < n_neighbors; ++i)
    {
      const std::size_t neigh_id = neighbors[i]->id();
      const Point neigh_pt       = _p_mesh.point(neigh_id);
      const Real  pt_dist        = _periodic_boundary->point_distance(pt0,
                                                                      neigh_pt);
      std::pair<std::size_t, Real> pp = std::make_pair(neigh_id, pt_dist);
      _nodes_neighbors[nid][i] = pp;
    }

    // build connections between the nodes and the center of the mesh
    _node_center_equilibrium_dist[nid] = _periodic_boundary->point_distance(pt0,
                                                                            center);
  } // end for
  STOP_LOG("build_spring_network()", "MeshSpringNetwork");
}

// ======================================================================
void MeshSpringNetwork::print_info()
{
  START_LOG("print_info()", "MeshSpringNetwork");

  // Loop over all the nodes, and print their spring connections
  MeshBase::node_iterator nd           = _p_mesh.active_nodes_begin();
  const MeshBase::node_iterator end_nd = _p_mesh.active_nodes_end();

  for (; nd != end_nd; ++nd)
  {
    // node pointer -> node id, node point, node neighbors
    Node *node                    = *nd;
    const std::size_t nid         = node->id();
    const Point pt0               = _p_mesh.point(nid);
    const std::size_t n_neighbors = _nodes_neighbors[nid].size();

    printf("Node # %lu, (x,y,z) = (%f, %f, %f) has %lu neigbhor nodes:\n",
           nid, pt0(0), pt0(1), pt0(2), n_neighbors);

    for (std::size_t i = 0; i < n_neighbors; ++i)
    {
      const std::size_t neigh_id = _nodes_neighbors[nid][i].first;
      const Real pt_dist         = _nodes_neighbors[nid][i].second;
      printf("     Node %lu, Equilibrium length = %f\n", neigh_id, pt_dist);
    }
    printf("     and Equilibrium distance to the center is %f\n",
           _node_center_equilibrium_dist[nid]);
  } // end for

  STOP_LOG("print_info()",
           "MeshSpringNetwork");
}

/* **********************************************************************

                       Priviate functions:

********************************************************************** */


// ======================================================================
void MeshSpringNetwork::find_nodal_neighbors(
  const Node                                   & node,
  const std::vector<std::vector<const Elem *> >& nodes_to_elem_map,
  std::vector<const Node *>                    & neighbors)
const
{
  START_LOG("find_nodal_neighbors()", "MeshSpringNetwork");

  MeshTools::find_nodal_neighbors(_p_mesh, node, nodes_to_elem_map, neighbors);

  STOP_LOG("find_nodal_neighbors()", "MeshSpringNetwork");
}
} // end of namespace
