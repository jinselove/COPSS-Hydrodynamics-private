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

// libmesh Includes -----------------------------------
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/parallel_object.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/exodusII_io.h"
#ifdef LIBMESH_HAVE_NANOFLANN
# include "libmesh/nanoflann.hpp"
#endif // ifdef LIBMESH_HAVE_NANOFLANN


// Local Includes   -----------------------------------
#include "rigid_particle.h"
#include "pm_periodic_boundary.h"

// #include "point_mesh.h"

// C++ Includes   -------------------------------------
// #include <stdio.h>
#include <string>
#include <vector>


namespace libMesh {
// enum ParticleType {POINT, RIGID, DEFORMABLE};


// class RigidParticle;

template<unsigned int KDDim>
class PointMesh;


/*
 * This template class defines particle-particle and particle-mesh
 * mapping relations. It takes advantage of KD Tree library nanoflann
 * to build the neighbor list for each particle.
 *
 * Generally speaking,for general shaped particles, their interaction
 * is achieved by computing force values on each tracking point(mesh node).
 * Therefore, the neighbor list may not be used in this case. However, we
 * still provide this functionality for use in some special cases, for
 * example, rigid particles with electrostatic interaction (no polarization)
 */


template<unsigned int KDDim>
class ParticleMesh :  public ReferenceCountedObject<ParticleMesh<KDDim> >,
                      public ParallelObject {
protected:

#ifdef LIBMESH_HAVE_NANOFLANN

  /**
   * This 'embeded in' class adapts list of libMesh \p Particle types
   * for use in a nanoflann KD-Tree. For more on the basic idea, see
   * examples/pointcloud_adaptor_example.cpp in the nanoflann src.
   * This class is similar to the "PointCloud" in the example.
   */

  template<unsigned int PLDim>
  class PointListAdaptor {
private:

    /* use ref, which will point to the address of _particles in the parent
       class
     * when _particles is modified, this ref need not to be updated
     * when the values change
     */
    const std::vector<RigidParticle *>& _pts;

public:

    // Constructor of PointListAdaptor
    PointListAdaptor(const std::vector<RigidParticle *>& particles)
      : _pts(particles)
    {}

    /**
     * libMesh \p Point coordinate type and std::size_t type
     */
    typedef Real        coord_t;
    typedef std::size_t size_t;

    /**
     * Must return the number of data points
     */
    inline size_t kdtree_get_point_count() const {
      return _pts.size();
    }

    /**
     * Returns the distance between the vector "p1[0:size-1]"
     * and the data point with index "idx_p2" stored in the class
     */
    inline coord_t kdtree_distance(const coord_t *p1,
                                   const size_t   idx_p2,
                                   size_t         size) const
    {
      libmesh_assert_equal_to(size, PLDim);
      libmesh_assert_less(idx_p2, _pts.size());

      // retrieve the point data with index = idx_p2
      const Point& p2(_pts[idx_p2]->get_centroid());

      switch (size)
      {
      case 3:
      {
        const coord_t d0 = p1[0] - p2(0);
        const coord_t d1 = p1[1] - p2(1);
        const coord_t d2 = p1[2] - p2(2);
        return d0 * d0 + d1 * d1 + d2 * d2;
      }

      case 2:
      {
        const coord_t d0 = p1[0] - p2(0);
        const coord_t d1 = p1[1] - p2(1);
        return d0 * d0 + d1 * d1;
      }

      case 1:
      {
        const coord_t d0 = p1[0] - p2(0);
        return d0 * d0;
      }

      default:
        libmesh_error_msg("ERROR: unknown size " << size);
      }

      return -1.0;
    } // end of kdtree_distance()

    /**
     * Returns the dim'th component of the idx'th point in the class:
     * Since this is inlined and the "dim" argument is typically an immediate
     *value, the
     *  "if's" are actually solved at compile time.
     */
    inline coord_t kdtree_get_pt(const size_t idx, int dim) const
    {
      libmesh_assert_less(dim, (int)PLDim);
      libmesh_assert_less(idx, _pts.size());
      libmesh_assert_less(dim, 3);

      const Point& p(_pts[idx]->get_centroid());

      if (dim == 0) return p(0);

      if (dim == 1) return p(1);

      return p(2);
    } // end of kdtree_get_pt()

    /**
     * Optional bounding-box computation: return false to default to a standard
     *bbox computation loop.
     * Return true if the BBOX was already computed by the class and returned in
     *"bb" so it can be
     * avoided to redo it again. Look at bb.size() to find out the expected
     *dimensionality
     * (e.g. 2 or 3 for point clouds)
     */
    template<class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const {
      return false;
    }
  }; // end of class PointListAdaptor


  typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<Real,
                                                                           PointListAdaptor<KDDim> >,
                                              PointListAdaptor<KDDim>,
                                              KDDim>kd_tree_t;

  /*
   *  If a data member is declared mutable, then it is legal to assign a value
      to this data member from a const member function.
   */
  mutable UniquePtr<kd_tree_t>_kd_tree;

#endif // LIBMESH_HAVE_NANOFLANN


  /**
   * Build & initialize the KD tree through nanoFlann, if needed.
   */
  virtual void construct_kd_tree();

public:

  // Constructor
  // Note if this is used to construct the object,
  // users are responsible to set the search radius
  // using the member function set_search_radius().
  ParticleMesh(MeshBase& mesh);


  // Constructor
  ParticleMesh(MeshBase  & mesh,
               const Real& search_radius_p,
               const Real& search_radius_e);


  // Constructor
  ParticleMesh(MeshBase          & mesh,
               PMPeriodicBoundary& pmpb,
               const Real        & search_radius_p,
               const Real        & search_radius_e);

  // Destructor
  ~ParticleMesh();

  /**
   * Read particle coordinate data from file
   * mesh_type = "surface_mesh" or "volume_mesh"
   * defaulty, we turn the "Electrostatics" off
   */
  void read_particles_data(const std::string& filename,            // particle
                                                                   // xyz file
                           const std::string& particle_mesh_type); // mesh type
                                                                   // of the
                                                                   // particle);

  /**
   * Read particle data from restart mesh files
   * only surface mesh is allowed
   */
  void read_particles_data_restart(const std::string& filename,
                                   const std::string& particle_mesh_type);

  /**
   * Read particle data from output_surface_node_*.csv
   */
  void read_particles_data_restart(const std::string & filename,
                                   const unsigned int& o_step);

  /**
   * Read chromatin data (cylinder particle) from the local file.
   * mesh_type = "surface_mesh" or "volume_mesh"
   */

  // void read_chromatin_data(const std::string& filename,      // particle xyz
  // file
  //                          const std::string& vmesh_file,    // volume mesh
  // file name
  //                          const std::string& smesh_file,    // surface mesh
  // file name
  //                          const std::string& mesh_type);    // mesh type of
  // the particle


  /**
   * Generate random particle coordinate data within the bounding box.
   * with r = 1.0, den = 1.0;
   */
  void generate_random_particles(const std::size_t N,
                                 const Real        bbox_XA,
                                 const Real        bbox_XB,
                                 const Real        bbox_YA,
                                 const Real        bbox_YB,
                                 const Real        bbox_ZA,
                                 const Real        bbox_ZB);


  void generate_random_particles(const std::size_t N,
                                 const Point     & bbox_min,
                                 const Point     & bbox_max);


  /**
   * Return the information of particles
   */
  std::vector<RigidParticle *>particles() const {
    return _particles;
  }

  /**
   * Return the total number of particles
   */
  std::size_t num_particles() const {
    return _n_rigid_particles;
  }

  /**
   * Check if the neighbor list is constructed sortedly
   */
  bool is_sorted() const {
    return _is_sorted;
  }

  /** Should we move this and above function to the ParticleMeshSystem???
  ** Because we need dist-function to calculate the particle-wall force!
  * Reinit the particle-mesh system. This includes:
  * (1) build the particle-particle neighbor list according to search radius;
  * (2) compute the force vectors on particles;
  * (3) build the element-particle neighbor list according to search radius;
  * (4) set the elem_id and proc_id for particles
  * ---- not incorporated:
  *   particle-wall repulsive force,
  *   which will be rebuilt in PMLinearSystem::reinit_particle_mesh()
  */
  void reinit();


  /**
   * print out the particle information
   */
  void print_particle_info() const;


  /**
   * Set the values of search radii
   */
  void set_search_radius(const Real rp, const Real re)
  {
    _search_radius_p = rp;  _search_radius_e = re;
  }

  /**
   * Return the search radius
   */
  Real search_radius(const std::string& p_e) const;


  /**
   * Add periodic boundary condition
   */
  void add_periodic_boundary(PMPeriodicBoundary& _periodic_bdry)
  {
    _periodic_boundary = &_periodic_bdry;
  }

  /**
   * Retrun the pointer to the periodic boundary for use
   */
  PMPeriodicBoundary* pm_periodic_boundary()
  {
    return _periodic_boundary;
  }

  /*
   * Return the Mesh
   */
  MeshBase& mesh() {
    return _mesh;
  }

  /*
   * Total number of the tracking points on the mesh,
   * which is the sum of number of nodes on each particle's mesh
   */
  std::size_t num_mesh_points() const;


  /*
   * Update the mesh for each particle
   * This is achieved by updating each nodal coordinates.
   */
  void update_mesh(const std::vector<Point>& nodal_pos,
                   const std::vector<Point>& nodal_vel);

  // void update_mesh(const std::vector<Real> & nodal_vec);


  // /*
  //  * update the (input) point_mesh object according to the particle_mesh
  //  * particle_mesh doesn't change, but point_mesh is changed.
  //  *
  //  * This function updates the input class, but doesn't change its own data.
  //  */
  // void update_point_mesh(PointMesh<KDDim>* point_mesh) const;


  // /*
  //  * update the particle_mesh object according to the (input) point_mesh
  //  * The input point_mesh doesn't change, but particle_mesh is changed.
  //  */
  // void update_particle_mesh(const PointMesh<KDDim>* point_mesh);


  /*
   * zero particle force density
   */
  void zero_node_force();


  /*
   * Return the mesh size (hmin/hmax) associated with this particle
   */
  std::vector<Real>mesh_size() const;


  /*
   * Correct the position of tracking points to avoid volume change!
   * NOTE: this is only for surface mesh.
   */
  void volume_conservation();


  /*
   * compute center of mass at step 0: center0, for rigid particles
   */
  void initial_particle_center_of_mass(std::vector<Point>& center0) const;

  /*
   * Write particle
   */
  void write_particle(const unsigned int            & step_id,
                      const unsigned int            & o_step,
                      const Real                    & real_time,
                      const std::vector<std::string>& output_file,
                      unsigned int                    comm_in_rank) const;

  /*
   * write out step and real time
   */
  void write_time(const unsigned int& step_id,
                  const unsigned int& o_step,
                  const Real        & real_time,
                  unsigned int        comm_in_rank) const;


  /*
   * Write out particle trajectories && forces && velocity to csv files
   * "output_bead_o_step.csv"
   */
  void write_particle_trajectory(const unsigned int& o_step,
                                 unsigned int        comm_in_rank) const;

  /*
   * Write out surface nodes positions to csv files
   */
  void write_surface_node(const unsigned int& o_step,
                          unsigned int        comm_in_rank) const;


  /*
   * Write out the particle's mesh(either surface mesh or volume mesh)
   */
  void        write_particle_mesh(const unsigned int& o_step) const;

  /*
   * Write out restart mesh (only save the latest step)
   */
  void        write_particle_mesh_restart() const;

  /*
   * Return a stitched mesh associated with all particles for electrostatic
   *solver
   * and assign particle's id to subdomain_id
   */
  SerialMesh& stitched_mesh();


  /********************************************************************************
      Uncomment play around the section below to construct neighbor_list for
      rigidParticle-rigidParticle and rigidParticle-elem
  ********************************************************************************/

  /**
   * Return the information of element neighbor list
   * around the search radius _search_radius_e
   * Mapping between elem id and element neighbor list of particles
   */
  std::map<const std::size_t, std::vector<std::size_t> >elem_neighbor_list() const
  {
    return _elem_neighbor_list;
  }

  /**
   * Return the particle neighbor list around a given element
   * around the search radius _search_radius_e
   */
  const std::vector<std::size_t>elem_neighbor_list(const Elem *elem)
  {
    const std::size_t elem_id = elem->id();

    return _elem_neighbor_list[elem_id];

    // return _elem_neighbor_list.at(elem_id); // not work well
  }

  /**
   * Return the particle neighbor list around a given LOCAL element
   * around the search radius _search_radius_e
   */
  const std::vector<std::size_t>local_elem_neighbor_list(const Elem *elem)
  {
    const std::size_t elem_id = elem->id();

    return _local_elem_neighbor_list[elem_id];
  }

  /**
   * Clear the KD-Tree if need
   * Note, different from the construct_kd_tree, this is a public function.
   */
  virtual void clear_kd_tree();


  /**
   * Build the neighbor list of a particle within the search_radius_p using KD
   *tree.
   * This includes the query particle itself, if the query_pt is also in the
   *particle list.
   * For example, the coords of the query_pt is the same as one of the particles
   *in the list.
   *
   * IndicesDists output pair<particle_id, distance>
   * Note we need the particle distance to compute their interaction forces.
   */
  virtual void build_particle_neighbor_list(
    const Point                  & query_pt,
    const bool                     is_sorted,
    std::vector<std::pair<std::size_t,
                          Real> >& IndicesDists);


  /**
   * Build the particle-particle neighbor list of particles within the
   *search_radius_p
   * using KD tree(logN complexity). This includes the original particle itself.
   * In fact, the particle neighbor list has been built in "reinit()", so this
   *is NOT USED!
   */
  virtual void build_particle_neighbor_list();


  /**NOT ACTUALLY USED NOW, written only for the purpose of comparing with
     KD-Tree!
   * Build the neighbor list of particles within the search radius directly!
   * (N^2 complexity)  This includes the original particle itself.
   * NOTE that the particle list constructed from this function is not sorted!
   */
  virtual void build_particle_neighbor_list_naively();


  /**
   * Build the neighbor list of partilces for a given elem
   * Output the n_list,
   * which gives the particle ids around elem's centroid within search_radius_e.
   *
   * Note we don't need the distance between particles and elem centroid here,
   * because this list is only used to evaluate element-wise load vector due to
   * the neighboring particles.
   */
  virtual void build_elem_neighbor_list(const Elem               *elem,
                                        const bool                is_sorted,
                                        std::vector<std::size_t>& n_list);


  /** FIXME: Should we build this locally and globally(by allgather)???
   *
   * Build the particle neighbor list for each element within the search radius
   * this will generate the _elem_neighbor_list and  _elem_neighbor_list_vector
   * NOTE this takes advantage of KD Tree
   */
  virtual void build_elem_neighbor_list();

  /**
   * print out the element neighbor list information
   */
  void         print_elem_neighbor_list(std::ostream& out = libMesh::out) const;

private:

  // number of rigid particles
  std::size_t _n_rigid_particles;

  // number of rigid particle types
  std::size_t _n_rigid_particle_types;

  // number of rigid particle mesh types
  std::size_t _n_rigid_particle_mesh_types;

  // mass values of rigid particles
  std::vector<Real>_mass;

  // mesh files
  std::vector<std::string>_rigid_particle_mesh_files;

  // A vector that store the pointers to Particle
  std::vector<RigidParticle *>_particles;

  // Mesh base: this is the domain(fluid) mesh, not the particle's mesh
  MeshBase& _mesh;

  // dimension
  std::size_t _dim;

  // Search radius (around a particle)
  Real _search_radius_p;

  // Search radius (around an element)
  Real _search_radius_e;

  // The point list adapter
  // - interface to the nanoflann in order to construct KD-Tree
  PointListAdaptor<KDDim>_point_list_adaptor;

  // Is the neighbor list sorted? set true defaultly (particle/elem neighbor
  // list)
  bool _is_sorted;


  // Element neighbor list: mapping between Elem and particle id's around it
  // This is NOT necessary to be on all the processors in implementation,
  // but we do this here only for test purpose. For local build
  std::map<const std::size_t, std::vector<std::size_t> >_elem_neighbor_list;
  std::map<const std::size_t, std::vector<std::size_t> >_local_elem_neighbor_list;

  // Pointer to the Periodic boundary condition
  PMPeriodicBoundary *_periodic_boundary;

  // precision of output
  const int o_precision = 6;
};
}
