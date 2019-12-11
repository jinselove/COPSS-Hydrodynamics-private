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


// #include "libmesh/point.h"
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/perf_log.h"

// C++ Includes   -----------------------------------
#include <stdio.h>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <unordered_set>

#include "pm_toolbox.h"
#include "particle_mesh.h"
#include "point_mesh.h"


namespace libMesh {
// ======================================================================
template<unsigned int KDDim>
PointMesh<KDDim>::PointMesh(ParticleMesh<KDDim>& particle_mesh,
                            const Real         & search_radius_p,
                            const Real         & search_radius_e)
  : ParallelObject(particle_mesh),
  _mesh(particle_mesh.mesh()),
  _search_radius_p(search_radius_p),
  _search_radius_e(search_radius_e),
  _point_list_adaptor(_particles),
  _polymer_chain(NULL),
  _is_sorted(true),
  _periodic_boundary(NULL)
{
  // total # of tracking points on the surface of all the finite-sized particles
  _num_point_particles = particle_mesh.num_mesh_points();
  _num_rigid_particles = particle_mesh.num_particles();
  _particles.resize(_num_point_particles);
  _rigid_particles.resize(_num_rigid_particles);
  _rigid_particles = particle_mesh.particles();
  _velocity_magnitude.resize(_num_point_particles);

  // create point particle list
  std::size_t count = 0;

  for (std::size_t i = 0; i < _num_rigid_particles; ++i)
  {
    // extract the nodal coordinates of the current particle
    std::vector<Point> node_xyz;
    _rigid_particles[i]->extract_nodes(node_xyz);

    // Create the PointParticle to form the particle list
    const std::size_t n_nodes = node_xyz.size();

    for (std::size_t j = 0; j < n_nodes; ++j)
    {
      Point pt(node_xyz[j]);
      PointParticle *particle = new PointParticle(pt, count);
      particle->set_parent_id(int(i));
      particle->set_point_type(LAGRANGIAN_POINT); // 1 - Lagrangian point; 0 -
                                                  // Polymer bead;
      _particles[count] = particle;
      count++;
    } // enf for j-loop
  } // enf for i-loop
  // Add the periodic boundary conditions
  this->add_periodic_boundary(*particle_mesh.pm_periodic_boundary());
  std::string msg =
    "PointMesh has been constructed from the extracted nodal points of particle's mesh!";

  //  msg += " count = " + std::to_string(count);
  PMToolBox::output_message(msg, this->comm());
}

// ======================================================================
template<unsigned int KDDim>
PointMesh<KDDim>::PointMesh(MeshBase                 & mesh,
                            std::vector<PolymerChain>& polymer_chains,
                            const Real               & search_radius_p,
                            const Real               & search_radius_e)
  : ParallelObject(mesh),
  _mesh(mesh),
  _search_radius_p(search_radius_p),
  _search_radius_e(search_radius_e),
  _point_list_adaptor(_particles),
  _polymer_chain(NULL),
  _is_sorted(true),
  _periodic_boundary(NULL)
{
  // Total # of beads and chains
  const std::size_t n_chains = polymer_chains.size();

  _num_point_particles = 0;

  for (std::size_t i = 0; i < n_chains; ++i)
  {
    _num_point_particles += polymer_chains[i].n_beads();
  }
  _particles.resize(_num_point_particles);
  _velocity_magnitude.resize(_num_point_particles);

  // Construct the PointParticle vector from the polymer chains
  std::size_t count                         = 0;
  std::vector<PointParticle *>::iterator it = _particles.begin();

  for (std::size_t i = 0; i < n_chains; ++i)
  {
    _particles.insert(it + count,
                      polymer_chains[i].beads().begin(),
                      polymer_chains[i].beads().end());

    count += polymer_chains[i].n_beads();
  }
}

// ======================================================================
template<unsigned int KDDim>
PointMesh<KDDim>::PointMesh(MeshBase    & mesh,
                            PolymerChain& polymer_chain,
                            const Real  & search_radius_p,
                            const Real  & search_radius_e)
  : ParallelObject(mesh),
  _mesh(mesh),
  _search_radius_p(search_radius_p),
  _search_radius_e(search_radius_e),
  _point_list_adaptor(_particles),
  _polymer_chain(&polymer_chain),
  _is_sorted(true),
  _periodic_boundary(NULL)
{
  // Get the point particles
  _particles           = polymer_chain.beads();
  _num_point_particles = _particles.size();
  _velocity_magnitude.resize(_num_point_particles);
}

// ======================================================================
template<unsigned int KDDim>
PointMesh<KDDim>::~PointMesh()
{
  // delete the particle pointers
  for (std::size_t i = 0; i < _num_point_particles; ++i)
  {
    // Only delete the Lagrangian points, polymer beads will be
    // destructed in PolymerChain class
    const PointType point_type = _particles[i]->point_type();

    if ((point_type != POLYMER_BEAD) && _particles[i])
    {
      delete _particles[i];
    }
  }
  _particles.clear();
//  _elem_neighbor_list.clear();
//  _local_elem_neighbor_list.clear();
//  _elem_elem_neighbor_list.clear();
//  _elem_point_containing_list.clear();
//  _elem_point_neighbor_list.clear();
}

// // ======================================================================
template<unsigned int KDDim>
void PointMesh<KDDim>::construct_kd_tree()
{
#ifdef LIBMESH_HAVE_NANOFLANN

  START_LOG("construct_kd_tree()", "PointMesh<KDDim>");

  // Initialize underlying KD tree if this is not constructed.
  if (_kd_tree.get() == NULL)
    _kd_tree.reset(new kd_tree_t(KDDim,
                                 _point_list_adaptor,
                                 nanoflann::KDTreeSingleIndexAdaptorParams(
                                   10 /* max leaf */)));

  libmesh_assert(_kd_tree.get() != NULL);

  _kd_tree->buildIndex();

  STOP_LOG("construct_kd_tree()",
           "PointMesh<KDDim>");
#endif // ifdef LIBMESH_HAVE_NANOFLANN
}

// ======================================================================
template<unsigned int KDDim>
void PointMesh<KDDim>::clear_kd_tree()
{
#ifdef LIBMESH_HAVE_NANOFLANN

  if (_kd_tree.get()) // If exist, delete the KD Tree and start fresh
    _kd_tree.reset(NULL);
#endif // ifdef LIBMESH_HAVE_NANOFLANN
}


// ======================================================================
template<unsigned int KDDim>
void PointMesh<KDDim>::reinit_neighbor_vector()
{
  START_LOG("reinit_neighbor_vector()", "PointMesh<KDDim>");

//  PMToolBox::output_message(">>>>>>>>>>> Start reinit neighbor vector",
//    this->comm
//  ());
  for (std::size_t p_id = 0; p_id < _num_point_particles; p_id++) {
    const Point pti = _particles[p_id]->point();
    const std::vector<dof_id_type>& n_list = _particles[p_id]->neighbor_list();
    std::vector<Point> r_ij(n_list.size());

    for (std::size_t j = 0; j < n_list.size(); j++) {
      const dof_id_type& n_id = n_list[j];
      const Point ptj         = _particles[n_id]->point();
      r_ij[j] = _periodic_boundary->point_vector(pti, ptj);
    }
    _particles[p_id]->set_neighbor_vector(r_ij);
  }
//  PMToolBox::output_message("Done !!!", this->comm());
  STOP_LOG("reinit_neighbor_vector()", "PointMesh<KDDim>");
}

// ======================================================================
template<unsigned int KDDim>
void PointMesh<KDDim>::reinit(bool& neighbor_list_update_flag)
{
  START_LOG("reinit()", "PointMesh<KDDim>");
  // zero particle force at every time step
  for (std::size_t j = 0; j < _num_point_particles; ++j) {
    _particles[j]->zero_particle_force();
  }

  // when reinit_neighbor_list is set to be true, which is done every some
  // diffusion
  // steps instead of every time step, we need to reinit particle neighbor list
  if (neighbor_list_update_flag) {
    // --> 1. build elem-elem neighbor list. Notice that this neighbor list
    // should and should only be built once since element locations don't
    // change during the simulation
    if(_elem_elem_neighbor_list.size()==0)
      this->build_elem_elem_neighbor_list();
    // --> 2. build elem-point containing list. This will update the points
    // contained in each element, which is necessary to build other types of
    // neighbor list later. Notice that this function will also update the
    // element id of each point.
    this->build_elem_point_containing_list();
    // --> 3. build elem-point neighbor list. This will update the neighbor
    // points of each element, which is necessary to build other types of
    // neighbor list later. Notice that this function will utilize the
    // information in elem_elem_neighbor_list and elem_point_containing_list.
    this->build_elem_point_neighbor_list();
    // --> 4. build point_point neighbor list. This will update the neighbor
    // points of each point. Notice that this function will utilize the elem
    // id updated in step 2 and elem_point_neighbor_list.
    this->build_point_point_neighbor_list();

    // after reinit neighbor list, set the flag to false
    neighbor_list_update_flag = false;
  } // end if (reinit_neighbor_list)

  // update particle neighbor distance at each reinit step
  this->reinit_neighbor_vector();
  STOP_LOG("reinit()", "PointMesh<KDDim>");
}

// =====================================================================
template<unsigned int KDDim>
void PointMesh<KDDim>::build_elem_elem_neighbor_list()
{
  // let's be sure we properly initialize on every processor at once: this
  // function is only executed in debug mode
  parallel_object_only();
  START_LOG("build_elem_elem_neighbor_list()", "PointMesh<KDDim>");

//  std::ostringstream oss;
//  oss << "\n>>>>>>>>>>> Building elem-elem neighbor list";
//  PMToolBox::output_message(oss, this->comm());
  // get some reference for convenience
  const dof_id_type& n_elem = _mesh.n_elem();

  //Initialize the global container
  _elem_elem_neighbor_list.clear();

  // create a helper map to store local mapping, we need this local mapping
  // to fill _elem_elem_neighbor_list globally.
  std::map<dof_id_type, std::vector<dof_id_type>>
  _local_elem_elem_neighbor_list;

  // loop over each element in parallel to build _elem_elem_neighbor_list
  MeshBase::const_element_iterator el = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh
    .active_local_elements_end();

  for(; el!=end_el; ++el)
  {
    // store a pointer to the current element
    const Elem *elem = *el;

    // get some references for convinence
    const std::size_t& elem_id = elem->id();
    const Point& p0 = elem->centroid();

    // create a list to store all element neighbors of this element
    std::vector<dof_id_type> nb_list;

    // Start to find all element neighbors of this element; we do this
    // naively since we only need to do this once
    for (dof_id_type i=0; i<n_elem; i++)
    {
      // neighbor cannot be itself
      if (elem_id != i)
      {
        //get a reference to this element
        const Elem* elem_i = _mesh.elem_ptr(i);
        const Point& pi = elem_i->centroid();

        // get the distance between there two element
        Point pt_ij;
        if (_periodic_boundary)
          pt_ij = _periodic_boundary->point_vector(p0, pi);
        else
          pt_ij = pi - p0;

        // if distance is within element search radius, then i is a neighbor
        if (pt_ij.norm() < _search_radius_e)
        {
          nb_list.push_back(i);
        }
      }
    }

    // insert the mapping of elem_id and nb_list to
    // _local_elem_elem_neighbor_list
    _local_elem_elem_neighbor_list.insert(std::make_pair(elem_id, nb_list));

  } // end loop elem in parallel

  // Until now, each processor has a different version of
  // _local_elem_elem_neighbor_list that stores the mapping of local elements
  // to their neighbor elements. For example, if element 0 is on processor 0
  // and element 1 is processor 1, then the mapping of element 0 is stored on
  // _local_elem_elem_neighbor_list on processor 0 and the mapping of element
  // 1 is store in _local_elem_elem_neighbor_list on processor 1. The problem
  // is that if we want to get access to the mapping of element 0 on
  // processor 1, we will not able to find it. Thus we need a global mapping

  // if in serial, _local_elem_elem_neighbor_list is the same as
  // _elem_elem_neighbor_list
  if (this->comm().size()<2)
  {
    _elem_elem_neighbor_list.insert(_local_elem_elem_neighbor_list.begin(),
      _local_elem_elem_neighbor_list.end());
  }

  // start creating global mapping

  // 1. create buffer vectors to help creating global mapping
  const dof_id_type len0 = _local_elem_elem_neighbor_list.size();
  // a vector that stores all element ids, size = len0
  std::vector<dof_id_type> buffer_elem_id_list(len0);
  // a vector that stores all neighbor element ids of all elements
  std::vector<dof_id_type> buffer_elem_nb_list;
  // a vector that stores all neighbor list length of all elements, size = len0
  std::vector<dof_id_type> buffer_elem_nb_len(len0);

  dof_id_type k=0;
  std::map<const dof_id_type, std::vector<dof_id_type>>::const_iterator p;

  for(p=_local_elem_elem_neighbor_list.begin();
    p!=_local_elem_elem_neighbor_list.end(); ++p)
  {
    buffer_elem_nb_len[k] = p->second.size();
    buffer_elem_id_list[k] = p->first;
    k++;
    // append neighbor list of current element to the end of buffer_elem_nb
    buffer_elem_nb_list.insert(buffer_elem_nb_list.end(), p->second.begin(),
      p->second.end());
  }

  // (2) all gather the local data. Check this link to understand how
  // allgather works: https://libmesh.github.io/doxygen/classlibMesh_1_1Parallel_1_1Communicator.html#a10632832fb05b53d1e7e0882c21fd6f3
  this->comm().allgather(buffer_elem_id_list);
  this->comm().allgather(buffer_elem_nb_len);
  this->comm().allgather(buffer_elem_nb_list);

  // (3) unpacked the gathered buffer to _elem_elem_neighbor_list
  k = 0;
  for(dof_id_type i=0; i<buffer_elem_id_list.size(); ++i)
  {
    const dof_id_type& elem_id = buffer_elem_id_list[i];
    const dof_id_type& elem_nb_len = buffer_elem_nb_len[i];
    std::vector<dof_id_type> elem_nb_list(elem_nb_len);

    for(dof_id_type j=0; j<elem_nb_len; ++j)
    {
      elem_nb_list[j] = buffer_elem_nb_list[k];
      k++;
    }

    _elem_elem_neighbor_list.insert(std::make_pair(elem_id, elem_nb_list));
  }

  // print out debug information to screen
//  oss << "Done!!!";
//  PMToolBox::output_message(oss, this->comm());
//  for (dof_id_type elem_id=0; elem_id<n_elem; elem_id++)
//  {
//    oss << "elem id " << elem_id << " neighbor element list: ";
//    const std::vector<dof_id_type>& elem_nb_list =
//      this->get_elem_elem_neighbor_list(elem_id);
//    for (dof_id_type i=0; i<elem_nb_list.size(); i++)
//      oss << elem_nb_list[i] << ", ";
//    oss << "\n";
//  }
//  PMToolBox::output_message(oss, this->comm());

  STOP_LOG("build_elem_elem_neighbor_list()", "PointMesh<KDDim>");
}

// =====================================================================
template<unsigned int KDDim>
void PointMesh<KDDim>::build_elem_point_containing_list()
{
  START_LOG("build_elem_point_containing_list()", "PointMesh<KDDim>");

//  std::ostringstream oss;
//  oss << "\n>>>>>>>>>>> Building elem-point containing list";
//  PMToolBox::output_message(oss, this->comm());

  // clear the mapping
  _elem_point_containing_list.clear();

//  std::cout<<"_particles.size() = "<<_particles.size()<<std::endl;
  for (std::size_t p_id=0; p_id<_particles.size(); p_id++)
  {
    bool use_locator = false;

    // if current element has not been assigned, we use locator
    if (_particles[p_id]->elem_id()==-1)
      use_locator = true;
    // if the current element doesn't contain this point anymore, we use locator
    else if (!_mesh.elem(_particles[p_id]->elem_id())->contains_point
      (_particles[p_id]->point()))
      use_locator = true;

    // use locator. If found, reset point elem id
    if (use_locator)
    {
      const Elem* elem = _mesh.point_locator().operator()
        (_particles[p_id]->point());
      if (elem== nullptr)
      {
        std::ostringstream oss;
        oss << "Error: cannot find the element that current particle is "
               "within : "
            <<"particle id = " << p_id <<", "
            <<"particle location = (" << _particles[p_id]->point()(0) << ", "
            << _particles[p_id]->point()(1) << ", "
            << _particles[p_id]->point()(2) << "). Exiting ...";
        PMToolBox::output_message(oss, this->comm());
        libmesh_error();
      }
      else
        _particles[p_id]->set_elem_id(elem->id());
    }

    // get the new elem_id
    const dof_id_type& elem_id = _particles[p_id]->elem_id();

    // insert a new pair if this elem_id has not been mapped to any particles
    // yet
    if (_elem_point_containing_list.count(elem_id)==0)
    {
      _elem_point_containing_list.insert
      (std::make_pair(elem_id, std::vector<dof_id_type>()));
    }

    // update mapping
    _elem_point_containing_list[elem_id].push_back(p_id);
  }

//  oss << "Done!!!";
//  PMToolBox::output_message(oss, this->comm());
//  for (dof_id_type elem_id=0; elem_id<_mesh.n_elem(); elem_id++) {
////    if (_elem_point_containing_list.count(elem_id)){
//      oss << "Elem id " << elem_id << " contains points: ";
//      const std::vector <dof_id_type> &contained_points =
//        this->get_elem_point_containing_list(elem_id);
//      for (dof_id_type i = 0; i < contained_points.size(); i++)
//        oss << contained_points[i] << ", ";
//      oss << "\n";
////    }
//  }
//  PMToolBox::output_message(oss, this->comm());

  STOP_LOG("build_elem_point_containing_list()", "PointMesh<KDDim>");
}

// =====================================================================
template<unsigned int KDDim>
void PointMesh<KDDim>::build_elem_point_neighbor_list()
{
  // let's be sure we properly initialize on every processor at once: this
  // function is only executed in debug mode
  parallel_object_only();
  START_LOG("build_elem_point_neighbor_list", "PointMesh<KDDim>");

//  std::ostringstream oss;
//  oss << "\n>>>>>>>>>>> Building elem-point neighbor list (using information "
//         "in elem-elem neighbor list and elem-point containing list)";
//  PMToolBox::output_message(oss, this->comm());
  // get some reference for convenience
  const dof_id_type& n_elem = _mesh.n_elem();

  //Initialize the global container
  _elem_point_neighbor_list.clear();

  // create a helper map to store local mapping, we need this local mapping
  // to fill _elem_elem_neighbor_list globally.
  std::map<dof_id_type, std::vector<dof_id_type>>
    _local_elem_point_neighbor_list;

  // loop over each element in parallel to build _elem_elem_neighbor_list
  MeshBase::const_element_iterator el = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh
    .active_local_elements_end();

  for(; el!=end_el; ++el)
  {
    // store a pointer to the current element
    const Elem *elem = *el;

    // get some references for convenience
    const std::size_t& elem_id = elem->id();

    // get the reference to all neighbor elements of this element
    const std::vector<dof_id_type>& elem_nb_list =
      this->get_elem_elem_neighbor_list(elem_id);

    // create a list to store all neighbor particles of this element
    std::vector<dof_id_type> nb_point_list;

    // start collecting all neighbor particles of this element
    // -> 0. first include points within this element itself
    const std::vector<dof_id_type>& containing_points =
      this->get_elem_point_containing_list(elem_id);
    for(dof_id_type j=0; j<containing_points.size(); j++)
      nb_point_list.push_back(containing_points[j]);
    // -> 1. loop over all neighbor elements of this element
    for (dof_id_type i=0; i<elem_nb_list.size(); i++)
    {
      const dof_id_type& nb_elem_id = elem_nb_list[i];
      // get a reference to the points within this neighbor element
      const std::vector<dof_id_type>& containing_points =
        this->get_elem_point_containing_list(nb_elem_id);
      // -> 2. loop over all points contained in this neighbor element
      for (dof_id_type j=0; j<containing_points.size(); j++)
        nb_point_list.push_back(containing_points[j]);
    };

    // insert the mapping of elem_id and nb_list to
    // _local_elem_elem_neighbor_list
    _local_elem_point_neighbor_list.insert(std::make_pair(elem_id,
      nb_point_list));
  } // end loop elem in parallel

  // Until now, each processor has a different version of
  // _local_elem_point_neighbor_list that stores the mapping of local elements
  // to their neighbor points. For example, if element 0 is on processor 0
  // and element 1 is processor 1, then the mapping of element 0 is stored on
  // _local_elem_point_neighbor_list on processor 0 and the mapping of element
  // 1 is store in _local_elem_point_neighbor_list on processor 1. The problem
  // is that if we want to get access to the mapping of element 0 on
  // processor 1, we will not able to find it. Thus we need a global mapping

  // if in serial, _local_elem_point_neighbor_list is the same as
  // _elem_point_neighbor_list
  if (this->comm().size()<2)
  {
    _elem_point_neighbor_list.insert(_local_elem_point_neighbor_list.begin(),
                                    _local_elem_point_neighbor_list.end());
  }

  // start creating global mapping

  // 1. create buffer vectors to help creating global mapping
  const dof_id_type len0 = _local_elem_point_neighbor_list.size();
  // a vector that stores all element ids, size = len0
  std::vector<dof_id_type> buffer_elem_id_list(len0);
  // a vector that stores all neighbor element ids of all elements
  std::vector<dof_id_type> buffer_point_nb_list;
  // a vector that stores all neighbor list length of all elements, size = len0
  std::vector<dof_id_type> buffer_point_nb_len(len0);

  dof_id_type k=0;
  std::map<const dof_id_type, std::vector<dof_id_type>>::const_iterator p;

  for(p=_local_elem_point_neighbor_list.begin();
      p!=_local_elem_point_neighbor_list.end(); ++p)
  {
    buffer_point_nb_len[k] = p->second.size();
    buffer_elem_id_list[k] = p->first;
    k++;
    // append neighbor list of current element to the end of buffer_elem_nb
    buffer_point_nb_list.insert(buffer_point_nb_list.end(), p->second.begin(),
                               p->second.end());
  }

  // (2) all gather the local data. Check this link to understand how
  // allgather works: https://libmesh.github.io/doxygen/classlibMesh_1_1Parallel_1_1Communicator.html#a10632832fb05b53d1e7e0882c21fd6f3
  this->comm().allgather(buffer_elem_id_list);
  this->comm().allgather(buffer_point_nb_len);
  this->comm().allgather(buffer_point_nb_list);

  // (3) unpacked the gathered buffer to _elem_elem_neighbor_list
  k = 0;
  for(dof_id_type i=0; i<buffer_elem_id_list.size(); ++i)
  {
    const dof_id_type& elem_id = buffer_elem_id_list[i];
    const dof_id_type& point_nb_len = buffer_point_nb_len[i];
    std::vector<dof_id_type> point_nb_list(point_nb_len);

    for(dof_id_type j=0; j<point_nb_len; ++j)
    {
      point_nb_list[j] = buffer_point_nb_list[k];
      k++;
    }

    _elem_point_neighbor_list.insert(std::make_pair(elem_id, point_nb_list));
  }

  // print out debug information to screen
//  oss << "Done!!!";
//  PMToolBox::output_message(oss, this->comm());
//  for (dof_id_type elem_id=0; elem_id<n_elem; elem_id++)
//  {
//    oss << "elem id " << elem_id << ", neighbor point list: ";
//    const std::vector<dof_id_type>& point_nb_list =
//      this->get_elem_point_neighbor_list(elem_id);
//    for (dof_id_type i=0; i<point_nb_list.size(); i++)
//      oss << point_nb_list[i] << ", ";
//    oss << "\n";
//  }
//  PMToolBox::output_message(oss, this->comm());


  STOP_LOG("build_elem_point_neighbor_list", "PointMesh<KDDim>");
}

// =====================================================================
template<unsigned int KDDim>
void PointMesh<KDDim>::build_point_point_neighbor_list()
{
  START_LOG("build_point_point_neighbor_list()", "PointMesh<KDDim>");

//  std::ostringstream oss;
//  oss << "\n>>>>>>>>>>> Building point-point neighbor list (first get the "
//         "element id of this point and then find all neighbor points of the "
//         "element)";
//  PMToolBox::output_message(oss, this->comm());

  // clear the mapping
  _point_point_neighbor_list.clear();

  for (std::size_t p_id=0; p_id<_particles.size(); p_id++)
  {
    // get the elem_id of this point
    const dof_id_type& elem_id = _particles[p_id]->elem_id();

    // get a copy of all the neighbor points of this element
    std::vector<dof_id_type> point_nb_list =
      this->get_elem_point_neighbor_list(elem_id);
    // erase the point itself from the neighbor list
    point_nb_list.erase(std::remove(point_nb_list.begin(), point_nb_list.end
    (), p_id), point_nb_list.end());

    // insert the mapping
//    _point_point_neighbor_list.insert(std::make_pair(p_id, point_nb_list));
    // set particle neighbor list
    _particles[p_id]->set_neighbor_list(point_nb_list);
  }

//  oss << "Done!!!";
//  PMToolBox::output_message(oss, this->comm());
//  for (dof_id_type p_id=0; p_id<_particles.size(); p_id++) {
//    oss << "Point id " << p_id <<" (elem_id = " << _particles[p_id]->elem_id()
//    <<"), neighbor points: ";
//    const std::vector <dof_id_type> &neighbor_points =
//      _particles[p_id]->neighbor_list();
//    for (dof_id_type i = 0; i<neighbor_points.size(); i++)
//      oss << neighbor_points[i] << ", ";
//    oss << "\n";
//  }
//  PMToolBox::output_message(oss, this->comm());

  STOP_LOG("build_point_point_neighbor_list()", "PointMesh<KDDim>");
}

// ======================================================================
template<unsigned int KDDim>
void PointMesh<KDDim>::build_particle_neighbor_list(const Point& tgt,
                                                    const bool is_sorted,
                                                    std::vector<std::pair<std::size_t, Real> >& IndicesDists)
{
  #ifdef LIBMESH_HAVE_NANOFLANN
    START_LOG("build_particle_neighbor_list(point)", "PointMesh<KDDim>");

  // if the KD tree is not built, construct the KD tree first
  if (_kd_tree.get() == NULL) this->construct_kd_tree();

  /* Find all the neighbors to a query_point within a maximum radius.
  * The output: the first element is a point index and the second the
  *corresponding distance.
  *
  *  If searchParams.sorted==true, the output list is sorted by ascending
  *distances.
  *
  *  For a better performance, it is advisable to do a .reserve() on the vector
  *  if you have any wild guess about the number of expected matches.
  */
  Real query_pt[] = { tgt(0), tgt(1), tgt(2) };
  nanoflann::SearchParams params;
  params.sorted = is_sorted; // not sorted, but with disordered sequence
  const Real& r_l2 = _search_radius_p * _search_radius_p;
  _kd_tree->radiusSearch(&query_pt[0], r_l2, IndicesDists, params);

  // the distance is L2 form, which is the distance square, so we take sqrt()
  for (std::size_t j = 0; j < IndicesDists.size(); ++j) {
  IndicesDists[j].second = std::sqrt(IndicesDists[j].second);
  }

  /* ------------------------------------------------------------------------
  * if the periodic boundary condition is applied, we must find the neighbor
  * list around its image particles for computing the interaction forces.
  *
  * In the following implementation, we assume that the domain size MUST
  * be larger than 4X search_radius so that only one image in this direction
  * needs to be considered. This is typically reasonable in realistic
  *simulations!
  * ------------------------------------------------------------------------*/
  if (_periodic_boundary != NULL)
  {
    // loop over each direction to find its images
    // std::size_t  NImage = 0;
    // if(KDDim==2) NImage = 3;
    // if(KDDim==3) NImage = 7;
    std::size_t NImage = 3;

    for (std::size_t i = 0; i < NImage; ++i)
    {
      Point im_pt;
      const bool has_image = _periodic_boundary->get_image_point(tgt,
                                                                 _search_radius_p,
                                                                 i,
                                                                 im_pt);

      if (has_image)
      {
        Real query_pt_im[KDDim];

        for (std::size_t j = 0; j < KDDim; ++j) query_pt_im[j] = im_pt(j);

        // find the neighbor particles around the image point!
        // Note that the im_pt is outside the box, so the returned list
        // does not include the im_pt itself.
        std::vector<std::pair<std::size_t, Real> > IndicesDists_image;
        _kd_tree->radiusSearch(&query_pt_im[0], r_l2, IndicesDists_image, params);

        // Add these to the list
        for (std::size_t j = 0; j < IndicesDists_image.size(); ++j)
        {
          IndicesDists_image[j].second = std::sqrt(IndicesDists_image[j].second);
          IndicesDists.push_back(IndicesDists_image[j]);
        }
      }
    }
  }
  /* -----------------------------------------------------------------------*/
  STOP_LOG("build_particle_neighbor_list(point)", "PointMesh<KDDim>");
  #endif // ifdef LIBMESH_HAVE_NANOFLANN
}

// ======================================================================
template<unsigned int KDDim>
void PointMesh<KDDim>::update_particle_mesh(ParticleMesh<KDDim> *particle_mesh)
const
{
  START_LOG("update_particle_mesh()", "PointMesh<KDDim>");

  // Extract the point coordinates from the point_mesh object
  std::vector<Point> nodal_vec(_num_point_particles);
  std::vector<Point> nodal_force(_num_point_particles);

  for (std::size_t i = 0; i < _num_point_particles; ++i) {
    nodal_vec[i]   = _particles[i]->point();
    nodal_force[i] = _particles[i]->particle_velocity();
  }
  particle_mesh->update_mesh(nodal_vec, nodal_force);

  STOP_LOG("update_particle_mesh()", "PointMesh<KDDim>");
}

// ======================================================================
template<unsigned int KDDim>
Real PointMesh<KDDim>::search_radius(const std::string& p_e) const
{
  if (p_e == "p") return _search_radius_p;
  else if (p_e == "e") return _search_radius_e;
  else
  {
    printf("PointMesh::search_radius(): the input option must be either p or e!\n");
    libmesh_error();
  }
}

// ======================================================================
// =========================== print functions ==========================
// ======================================================================
template<unsigned int KDDim>
void PointMesh<KDDim>::print_point_info() const
{
  // const int o_precision = 9;
  // std::ofstream out_file;
  // out_file.precision(o_precision);
  // out_file.open("out.initial_force_field", std::ios_base::out);
  // out_file <<"bead_id pos_x pos_y pos_z force_x force_y force_z\n";
  // for (std::size_t i = 0; i<_num_point_particles; ++i){
  //   out_file <<i <<" ";
  //   for (int j=0; j<3; j++) out_file << _particles[i]->point()(j) <<" ";
  //   for (int j=0; j<3; j++) out_file << _particles[i]->particle_force()[j]
  // <<" ";
  //   out_file<<"\n";
  // }
  // out_file.close();
  // printf("==> done: output bead info after first system.reinit() to
  // 'out.initial_force_field'\n");
  printf(
         "======================= printing the point information: =======================\n\n");
  printf("There are totally %lu points (search radius = %f)\n",
         _particles.size(),
         _search_radius_p);

  for (std::size_t j = 0; j < _particles.size(); ++j)
    _particles[j]->print_info(this->comm());
  printf(
    "========================= end of the point information =========================\n\n");
}

// ======================================================================
//template<unsigned int KDDim>
//void PointMesh<KDDim>::print_elem_neighbor_list(std::ostream& out) const
//{
//  printf(
//    "======================= printing the element neighbor list: ========================\n\n");
//  std::map<const dof_id_type, std::vector<dof_id_type> >::const_iterator p;
//
//  for (p = _elem_point_neighbor_list.begin(); p != _elem_point_neighbor_list
//  .end(); ++p)
//  {
//    const std::size_t elem_id             = p->first;
//    const std::vector<dof_id_type> n_list = p->second;
//    const Elem *elem                      = _mesh.elem(elem_id);
//    const Point center_pt                 = elem->centroid();
//    const Real  hmax                      = elem->hmax();
//
//    // ----- Scheme 2: print out information using printf (from all the
//    // processors) -----
//    if (n_list.size() > 0)
//    {
//      printf(
//             "========== There are %lu neighbor points around the element %lu, output rank = %u :===\n",
//        n_list.size(),
//        elem_id,
//        this->comm().rank());
//      printf("element centroid = (%f, %f, %f), and hmax = %f \n",
//             center_pt(0), center_pt(1), center_pt(2), hmax);
//
//      if (!n_list.empty())
//        for (std::size_t i = 0; i < n_list.size(); ++i) // printf("%lu    ",
//                                                        // n_list[i]);
//          _particles[n_list[i]]->print_info(false);
//      else out << "There is no neighboring particle around this element! \n";
//      printf("\n");
//    } // end if
//  }   // end for p-loop
//  printf(
//    "======================= end of the element neighbor list ======================\n\n");
//}

// ======================================================================
template<unsigned int KDDim>
void PointMesh<KDDim>::set_bead_velocity(const std::vector<Real>& vel)
{
  START_LOG("PointMesh::set_bead_velocity()", "PointMesh");
  const int dim = 3;

  // assign velocity to 0-th to (n-1)-th particle
  Point p_velocity(0.);

  for (std::size_t p_id = 0; p_id < _num_point_particles; p_id++) {
    for (int i = 0; i < dim; i++) p_velocity(i) = vel[p_id * dim + i];
    _particles[p_id]->set_particle_velocity(p_velocity);
    _velocity_magnitude[p_id] = p_velocity.norm_sq();
  }
  auto minmax_velocity_magniture = std::minmax_element(
    _velocity_magnitude.begin(),
    _velocity_magnitude.end());
  _min_velocity_magnitude = std::sqrt(*minmax_velocity_magniture.first);
  _max_velocity_magnitude = std::sqrt(*minmax_velocity_magniture.second);
  STOP_LOG("PointMesh::set_bead_velocity()", "PointMesh");
}

// ======================================================================
template<unsigned int KDDim>
const void PointMesh<KDDim>::write_initial_surface_node_pos() const
{
  std::ostringstream oss;

  oss << "output_initial_bead_position.csv";
  std::ofstream out_file;

  if (this->comm().rank() == 0) {
    out_file.open(oss.str(), std::ios_base::out);

    // write out the csv file
    // POINT data
    out_file << "rigid_particle_id node_id x_coord y_coord z_coord\n";
    out_file.precision(6);

    for (int k = 0; k < _num_rigid_particles; k++) {
      RigidParticle *rigid_particle = _rigid_particles[k];

      for (std::size_t i = 0; i < rigid_particle->num_mesh_nodes(); ++i)
      {
        const Point node_pos = rigid_particle->mesh_point(i);
        out_file << k << " " << i << " ";

        // write position
        for (std::size_t j = 0; j < KDDim; ++j) {
          out_file << node_pos(j) << " ";
        }
        out_file << "\n";
      }
    }
    out_file.close();
  } // end if comm_in_rank == 0
}

// ------------------------------------------------------------
// Explicit Instantiations
template class PointMesh<1>;
template class PointMesh<2>;
template class PointMesh<3>;
} // end of namespace
