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



//#include "libmesh/point.h"
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



namespace libMesh
{
  
// ======================================================================
template <unsigned int KDDim>
PointMesh<KDDim>::PointMesh(ParticleMesh<KDDim>& particle_mesh,
                            const Real& search_radius_p,
                            const Real& search_radius_e)
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
  for (std::size_t i=0; i<_num_rigid_particles; ++i)
  {
    // extract the nodal coordinates of the current particle
    std::vector<Point> node_xyz;
    _rigid_particles[i]->extract_nodes(node_xyz);    
    // Create the PointParticle to form the particle list
    const std::size_t n_nodes = node_xyz.size();
    for (std::size_t j=0; j<n_nodes; ++j)
    {
      Point pt ( node_xyz[j] );
      PointParticle* particle = new PointParticle(pt,count);
      particle->set_parent_id( int(i) );
      particle->set_point_type(LAGRANGIAN_POINT); // 1 - Lagrangian point; 0 - Polymer bead;
      _particles[count] = particle;
      count++;
    } // enf for j-loop
  } // enf for i-loop
  // Add the periodic boundary conditions
  this->add_periodic_boundary( *particle_mesh.pm_periodic_boundary() );
  std::string msg = "PointMesh has been constructed from the extracted nodal points of particle's mesh!";
//  msg += " count = " + std::to_string(count);
  PMToolBox::output_message(msg,this->comm());
  return;
}
  
  
  
  
// ======================================================================
template <unsigned int KDDim>
PointMesh<KDDim>::PointMesh(MeshBase& mesh,
                            std::vector<PolymerChain>& polymer_chains,
                            const Real& search_radius_p,
                            const Real& search_radius_e)
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
  for (std::size_t i=0; i<n_chains; ++i)
  {
    _num_point_particles += polymer_chains[i].n_beads();
  }
  _particles.resize(_num_point_particles);
  _velocity_magnitude.resize(_num_point_particles);
  // Construct the PointParticle vector from the polymer chains
  std::size_t count = 0;
  std::vector<PointParticle*>::iterator it = _particles.begin();
  for (std::size_t i=0; i<n_chains; ++i)
  {
    _particles.insert(it+count,
                      polymer_chains[i].beads().begin(),
                      polymer_chains[i].beads().end());
    
    count += polymer_chains[i].n_beads();
  }

}



// ======================================================================
template <unsigned int KDDim>
PointMesh<KDDim>::PointMesh(MeshBase& mesh,
                            PolymerChain& polymer_chain,
                            const Real& search_radius_p,
                            const Real& search_radius_e)
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
  _particles = polymer_chain.beads();
  _num_point_particles = _particles.size();
  _velocity_magnitude.resize(_num_point_particles);
}

  

// ======================================================================
template <unsigned int KDDim>
PointMesh<KDDim>::~PointMesh()
{
  // delete the particle pointers
  for (std::size_t i=0; i<_num_point_particles; ++i)
  {
    // Only delete the Lagrangian points, polymer beads will be
    // destructed in PolymerChain class
    const PointType point_type = _particles[i]->point_type();
    if( point_type!=POLYMER_BEAD && _particles[i] )
    {
      delete _particles[i];
    }
  }
  _particles.clear();
  _elem_neighbor_list.clear();
  _local_elem_neighbor_list.clear();
}



// ======================================================================
// template <unsigned int KDDim>
// void PointMesh<KDDim>::read_points_data(const std::string& filename)
// {
//   START_LOG ("read_points_data()", "PointMesh<KDDim>");
  
//   std::cout <<"\n###point data filename = "<<filename <<std::endl;
//   std::ifstream infile;
//   infile.open (filename, std::ios_base::in);
//   if( !infile.good() )
//   {
//     printf("***warning: read_points_data() can NOT read the particle coordinate file!");
//     libmesh_error();
//   }
  
//   // init
//   const std::size_t dim = _mesh.mesh_dimension();
//   const PointType point_type = POINT_PARTICLE;  // polymer bead; Lagrangian point; Point particle;
//   Real x=0., y=0., z=0., r=0., den=0.;  // initialize particle coords and radius
//   infile >> num_point_particles;  // total number of particles
  
  
//   // read particle data
//   _particles.resize(num_point_particles);
//   for (std::size_t i=0; i<num_point_particles; ++i)
//   {
//     infile >> x >> y >> z >> r >> den;
//     r = 0.0;  den = 0.0;  // zeros for the point particle
//     if (KDDim==2 || dim==2) z = 0.0;
//     Point pt(x,y,z);
//     PointParticle* particle = new PointParticle(pt, i);
//     particle->set_point_type(point_type);  // 1 - tracking point; 0 - bead point;
//     _particles[i] = particle;
    
//     // --------------------- test ------------------------------------------
// //    if(this->comm().rank()==0 )
// //    {
// //      printf("read_particles_data: x = %f, y = %f, z = %f. ",x, y, z );
// //      printf("radius = %f, relative density = %f. ",r, den );
// //      printf("MPI_rank = %d\n",this->comm().rank() );
// //    }
//     // --------------------- test ------------------------------------------
//   } // end for i-loop
  
//   infile.close();
  
//   this->comm().barrier();
//   std::cout << "Reading point data from "<<filename<<" is completed!\n\n";
  
//   STOP_LOG ("read_points_data()", "PointMesh<KDDim>");
// }



// // ======================================================================
// template <unsigned int KDDim>
// void PointMesh<KDDim>::generate_random_points(const std::size_t N,
//                                               const Real bbox_XA, const Real bbox_XB,
//                                               const Real bbox_YA, const Real bbox_YB,
//                                               const Real bbox_ZA, const Real bbox_ZB)
// {
//   START_LOG ("generate_random_points()", "PointMesh<KDDim>");
  
//   // problem dimension and domain size
//   const std::size_t dim = _mesh.mesh_dimension();
//   const Real max_range_x = bbox_XB - bbox_XA;
//   const Real max_range_y = bbox_YB - bbox_YA;
//   const Real max_range_z = bbox_ZB - bbox_ZA;
//   const Real r = 1.0, den = 1.0;
  
//   // generate random particle coordinates inside the domain, and write out the file
//   // We only use rank=0 proccessor to avoid generating different random numbers on each processor.
//   if( this->comm().rank()==0 )
//   {
//     printf("---> test in generate_random_points: Generating %lu random points ...\n",N);
    
//     // write the particle coordinates into a file
//     std::string filename = "random_points_file.txt";
//     int o_width = 5, o_precision = 9;
    
//     std::ofstream outfile;
//     outfile.open(filename,std::ios_base::out);
//     outfile << N << "\n";
//     for (size_t i=0;i<N;i++)
//     {
//       // generate random coordinates
//       Real x = max_range_x * (std::rand() % 1000) / 1000 + bbox_XA;
//       Real y = max_range_y * (std::rand() % 1000) / 1000 + bbox_YA;
//       Real z = max_range_z * (std::rand() % 1000) / 1000 + bbox_ZA;
//       if (KDDim==2 || dim==2) z = 0.0;
      
//       outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
//       outfile.precision(o_precision);   outfile.width(o_width);
//       outfile << x << "  " << y << "  " << z << "  " << r << "  "<< den << "  \n";
//     } // end loop-i
    
//     // close the file
//     outfile.close();
//     printf("---> test in generate_random_points: random particle file is created!\n");
//   }
  
//   this->comm().barrier();
  
//   STOP_LOG ("generate_random_points()", "PointMesh<KDDim>");
// }



// ======================================================================
// template <unsigned int KDDim>
// void PointMesh<KDDim>::generate_random_points(const std::size_t N,
//                                               const Point& bbox_min,
//                                               const Point& bbox_max)
// {
//   this->generate_random_points(N,
//                                bbox_min(0),bbox_max(0),
//                                bbox_min(1),bbox_max(1),
//                                bbox_min(2),bbox_max(2) );
// }



// // ======================================================================
template <unsigned int KDDim>
void PointMesh<KDDim>::construct_kd_tree ()
{
#ifdef LIBMESH_HAVE_NANOFLANN
  
  START_LOG ("construct_kd_tree()", "PointMesh<KDDim>");
  
  // Initialize underlying KD tree if this is not constructed.
  if (_kd_tree.get() == NULL)
    _kd_tree.reset(new kd_tree_t(KDDim,
                                 _point_list_adaptor,
                                 nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */)));
  
  libmesh_assert (_kd_tree.get() != NULL);
  
  _kd_tree->buildIndex();
  
  STOP_LOG ("construct_kd_tree()", "PointMesh<KDDim>");
#endif
}



// ======================================================================
template <unsigned int KDDim>
void PointMesh<KDDim>::clear_kd_tree()
{
#ifdef LIBMESH_HAVE_NANOFLANN
  if (_kd_tree.get())   // If exist, delete the KD Tree and start fresh
    _kd_tree.reset (NULL);
#endif
}



// ======================================================================
template <unsigned int KDDim>
void PointMesh<KDDim>::build_particle_neighbor_list(const Point &tgt,
                                                    const bool is_sorted,
                                                    std::vector<std::pair<std::size_t,Real> >& IndicesDists)
{
#ifdef LIBMESH_HAVE_NANOFLANN
  START_LOG ("build_particle_neighbor_list(point)", "PointMesh<KDDim>");
  
  // if the KD tree is not built, construct the KD tree first
  if (_kd_tree.get() == NULL)
    this->construct_kd_tree();
  
  /* Find all the neighbors to a query_point within a maximum radius.
   * The output: the first element is a point index and the second the corresponding distance.
   *
   *  If searchParams.sorted==true, the output list is sorted by ascending distances.
   *
   *  For a better performance, it is advisable to do a .reserve() on the vector
   *  if you have any wild guess about the number of expected matches.
   */
  Real query_pt[] = { tgt(0), tgt(1), tgt(2) };
  nanoflann::SearchParams params;
  params.sorted = is_sorted;  // not sorted, but with disordered sequence
  const Real& r_l2 = _search_radius_p*_search_radius_p;
  _kd_tree->radiusSearch(&query_pt[0], r_l2, IndicesDists, params);
  
  // the distance is L2 form, which is the distance square, so we take sqrt()
  for (std::size_t j=0; j<IndicesDists.size(); ++j){
    IndicesDists[j].second = std::sqrt( IndicesDists[j].second );
  }
  
  /* ------------------------------------------------------------------------
   * if the periodic boundary condition is applied, we must find the neighbor
   * list around its image particles for computing the interaction forces.
   *
   * In the following implementation, we assume that the domain size MUST
   * be larger than 4X search_radius so that only one image in this direction
   * needs to be considered. This is typically reasonable in realistic simulations!
   * ------------------------------------------------------------------------*/
  if(_periodic_boundary != NULL)
  {
    // loop over each direction to find its images
    // std::size_t  NImage = 0;
    // if(KDDim==2) NImage = 3;
    // if(KDDim==3) NImage = 7;
    std::size_t NImage = 3;
    for (std::size_t i=0; i<NImage; ++i)
    {
      Point im_pt;
      const bool has_image = _periodic_boundary->get_image_point(tgt,_search_radius_p,i,im_pt);
      
      if(has_image)
      {
        Real query_pt_im[KDDim];
        for(std::size_t j=0; j<KDDim; ++j) query_pt_im[j] = im_pt(j);
        
        // find the neighbor particles around the image point!
        // Note that the im_pt is outside the box, so the returned list
        // does not include the im_pt itself.
        std::vector<std::pair<std::size_t,Real> > IndicesDists_image;
        _kd_tree->radiusSearch(&query_pt_im[0], r_l2, IndicesDists_image, params);
        
        // Add these to the list
        for (std::size_t j=0; j<IndicesDists_image.size(); ++j)
        {
          IndicesDists_image[j].second = std::sqrt(IndicesDists_image[j].second);
          IndicesDists.push_back( IndicesDists_image[j] );
        }
        // -------------------------- test --------------------------
//        if (this->comm().rank()==0 && IndicesDists_image.size()>0)
//        {
//          printf("--->test in build_particle_neighbor_list() i=%lu:\n",i);
//          printf("    Original point = (%f %f %f) \n", tgt(0), tgt(1), tgt(2) );
//          printf("    Image point    = (%f %f %f) has the neighbors:\n", im_pt(0),im_pt(1),im_pt(2));
//          for (std::size_t j=0; j<IndicesDists_image.size(); ++j)
//            printf("      particle id %lu, distance to the image point is %f\n",
//                   IndicesDists_image[j].first, IndicesDists_image[j].second);
//        }
        // -------------------------- test --------------------------
      } // end if (has_image)
    } // end for i-loop
  } // end if(_periodic_boundary != NULL)
  /* -----------------------------------------------------------------------*/ 
  STOP_LOG ("build_particle_neighbor_list(point)", "PointMesh<KDDim>");
#endif
}



// ======================================================================
template <unsigned int KDDim>
void PointMesh<KDDim>::build_particle_neighbor_list()
{
#ifdef LIBMESH_HAVE_NANOFLANN
  
  // do a radius search to build the neighbor list for ALL particles.
  for (std::size_t j=0; j<_num_point_particles; ++j)
  {
    // get the neighbor indices & distance values
    const Point &tgt( _particles[j]->point() );
    std::vector<std::pair<std::size_t,Real> > IndicesDists0, IndicesDists;
    this->build_particle_neighbor_list(tgt, _is_sorted, IndicesDists);
    /* IndicesDists above returns <particle_id, distance>! */
    
    // Exclude the current particle itself!
    const std::size_t pid = _particles[j]->id();
    for (std::size_t i=0; i<IndicesDists.size(); ++i){
      if ( IndicesDists[i].first!= pid){
        IndicesDists0.push_back( IndicesDists[i] );
      }
    }
    
    // set the neighbor list of the j-th particle
    _particles[j]->set_neighbor_list (IndicesDists0);
  } // end for j-loop over particles
  
#endif
}


// ======================================================================
template <unsigned int KDDim>
void PointMesh<KDDim>::build_particle_neighbor_list_naively()
{
  START_LOG ("build_particle_neighbor_list_naively()", "PointMesh<KDDim>");
  
  // do a radius search to build the neighbor list for ALL particles.
  for (std::size_t i=0; i<_num_point_particles; ++i)
  {
    // The info of the i-th particle
    const Point pt0( _particles[i]->point() );
    std::vector<std::pair<std::size_t,Real> > IndicesDists;
    
    // j-loop over the particles
    for (std::size_t j=0; j<_num_point_particles; ++j)
    {
      if (i!=j)
      {
        const Point ptj( _particles[j]->point() );
        Point pt_ij;
        if(_periodic_boundary){
          pt_ij = _periodic_boundary->point_vector(pt0, ptj);
        }
        else{
          pt_ij = ptj - pt0;
        }
        const Real dist = pt_ij.norm(); // the real distance
        
        if( dist < _search_radius_p ){
          IndicesDists.push_back ( std::make_pair(_particles[j]->id(), dist) );
        } // end if
        
      } // end if(i!=j)
    } // end for j-loop
    
    // if needed, sort the particle neighbor list.
    if (_is_sorted)
    {
      // ***sort the particle neighbor list (NOT implemented here)
      printf("*** warning: the particle neighbor list is not sorted in this function!");
    }
    
    // set the neighbor list of the j-th particle (the result is not sorted)
    _particles[i]->set_neighbor_list (IndicesDists);
  }
  
  STOP_LOG ("build_particle_neighbor_list_naively()", "PointMesh<KDDim>");
}



// ======================================================================
template <unsigned int KDDim>
void PointMesh<KDDim>:: build_elem_neighbor_list(const Elem* elem,
                                                 const bool is_sorted,
                                                 std::vector<std::size_t>& n_list)
{
  // This function must be called on every processor.
  parallel_object_only();
  
#ifdef LIBMESH_HAVE_NANOFLANN
  START_LOG ("build_elem_neighbor_list(elem)", "PointMesh<KDDim>");
  
  
  // If the KD tree is not built, construct the KD tree first
  if (_kd_tree.get() == NULL)
    this->construct_kd_tree();
  
  // find the controid and size of this elem
  const Point center_pt = elem->centroid();
  const Real  hmax      = elem->hmax();
  
  // Build the nearest particle neighbor list of this elem within search radius.
  // The search radius is set the max of half elem hmax and the given radius
  // FIXME: the radiusSearch use L2 adapter, so the input radius is r^2 ??
  const Real query_pt[] = { center_pt(0), center_pt(1), center_pt(2) };
  const Real r = std::max( _search_radius_e, hmax/2. );
  const Real r_l2 = r*r;
  nanoflann::SearchParams params;
  params.sorted = is_sorted;      // with sorted/unsorted sequence
  std::vector<std::pair<std::size_t,Real> > IndicesDists;
  _kd_tree->radiusSearch(&query_pt[0], r_l2, IndicesDists, params);
  
  
  /* ------------------------------------------------------------------------
   * if the periodic boundary condition is applied, we must find the neighbor
   * list around its image particles for computing the interaction forces.
   *
   * In the following implementation, we assume that the domain size MUST
   * be larger than 4X search_radius so that only one image in this direction
   * needs to be considered. This is typically reasonable in realistic simulations!
   * ------------------------------------------------------------------------*/
  if(_periodic_boundary != NULL)
  {
    // loop over each direction to find its images
    // std::size_t  NImage = 0;
    // if(KDDim==2) NImage = 3;
    // if(KDDim==3) NImage = 7;
    std::size_t NImage = 3;
    for (std::size_t i=0; i<NImage; ++i)
    {
      Point im_pt;
      const bool has_image = _periodic_boundary->get_image_point(center_pt,_search_radius_e,i,im_pt);
      
      if(has_image)
      {
        Real query_pt_im[KDDim];
        for(std::size_t j=0; j<KDDim; ++j) query_pt_im[j] = im_pt(j);
        
        // find the neighbor particles around the image point!
        // Note that the im_pt is outside the box, so the returned list
        // does not include the im_pt itself.
        std::vector<std::pair<std::size_t,Real> > IndicesDists_image;
        _kd_tree->radiusSearch(&query_pt_im[0], r_l2, IndicesDists_image, params);
        
        // Add these to the list
        for (std::size_t j=0; j<IndicesDists_image.size(); ++j){
          IndicesDists.push_back( IndicesDists_image[j] );
        } // end for j-loop
        
      } // end if
      
    } // end for i-loop
  } // if(_periodic_boundary != NULL)
  /* -----------------------------------------------------------------------*/
  
  
  /* ------------------------------------------------------------------------
   We don't need to store the distance values in this elem-particle list!
   ------------------------------------------------------------------------*/
  const std::size_t np = IndicesDists.size();
  n_list.resize (np);
  for (std::size_t j=0; j<np; ++j){
    n_list[j] = IndicesDists[j].first;
  }
  
  STOP_LOG ("build_elem_neighbor_list(elem)", "PointMesh<KDDim>");
#endif
}



// ======================================================================
template <unsigned int KDDim>
void PointMesh<KDDim>::build_elem_neighbor_list()
{
  START_LOG ("build_elem_neighbor_list()", "PointMesh<KDDim>");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   init the container:
   the vector maps elem_id to the element neighbor list
   the local list stores the mapvector on each local processor
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  _elem_neighbor_list.clear();
  _local_elem_neighbor_list.clear();
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   prepare the send list for particle and element pairs
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  std::vector<dof_id_type> particle_id_send_list_vec, element_id_send_list_vec;
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   loop over each element to build its local neighbor list
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  std::size_t call_count = 0;
//  MeshBase::const_element_iterator       el     = _mesh.active_elements_begin();
//  const MeshBase::const_element_iterator end_el = _mesh.active_elements_end();
  MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_local_elements_end();
  for ( ; el != end_el; ++el)
  {
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Store a pointer to the element we are currently working on.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    const Elem* elem          = *el;
    const std::size_t elem_id = elem->id();
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Find the elem neighbor list of particles around the elem's centroid
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    std::vector<std::size_t> n_list;
    this->build_elem_neighbor_list(elem, _is_sorted, n_list);
    call_count++;
    
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Construct the local elem-neighbor list mapvector, and this will be
     distributed to all the processors through MPI_Allgather()
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    _local_elem_neighbor_list.insert (std::make_pair(elem_id,n_list) );
    
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Print out the element-particle neighbor list (for test purpose only)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//    if (n_list.size()>0 ) //if (pid==16 || pid==17 )
//    {
//      printf("--->test: in build_elem_neighbor_list()\n");
//      printf("--->elem id = %lu has %lu neighboring particles:\n", elem_id, n_list.size() );
//      printf("  neighbor list: \n");
//      for (std::size_t j=0; j<n_list.size(); ++j) printf("    %lu  ", n_list[j]);
//      printf("\n");
//    } // end if
    
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Find the particles sitting in this element. Note that a particle can sit in more
     than one elements if it is on the facets or corners of an element!
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//    printf("--->TEST: in PointMesh<KDDim>::build_elem_neighbor_list():\n");
//    printf("         Elem id = %lu, and volume is %f \n",elem_id, elem->volume() );
//    printf("         Element Nodes are : \n");
//    for(std::size_t k=0; k<elem->n_nodes(); ++k)
//    {
//      printf("            node %lu = (%f,%f,%f)\n", k,
//                        elem->point(k)(0), elem->point(k)(1), elem->point(k)(2));
//    }
    
    
    for (std::size_t j=0; j<n_list.size(); ++j)
    {
      const std::size_t pid = n_list[j];
      const Point pt        = _particles[pid]->point();
//      printf("         pt[%lu] = (%f,%f,%f)",  j, pt(0),pt(1),pt(2));
      const bool  inside    = elem->contains_point(pt);
//      if (inside) printf(": inside \n");
//      else printf(": outside\n");
      
      if (inside ){
        particle_id_send_list_vec.push_back(pid);
        element_id_send_list_vec. push_back(elem_id);
        //printf("*****pid = %lu, eid = %lu, processor_id = %i\n",pid,elem_id,elem->processor_id() );
        //printf("*****point xyz = (%f %f %f)"\n\n",pt(0),pt(1),pt(2)); //elem->print_info();
      } // end if
      
    } // end for j-loop
    
  } // end for elem-loop
  
//  this->comm().barrier();
//  printf("--->test in build_elem_neighbor_list() call_count = %lu on rank = %d\n",
//         call_count,this->comm().rank() );
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Allgather the particle_send_list and element_send_list
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  this->comm().allgather(particle_id_send_list_vec);
  this->comm().allgather(element_id_send_list_vec);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set the element id for each particle.
   This operation is on all the processors
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for(std::size_t i=0; i<particle_id_send_list_vec.size(); ++i)
  {
    const std::size_t pid = particle_id_send_list_vec[i];
    const std::size_t eid =  element_id_send_list_vec[i];
    _particles[pid]->set_elem_id(eid);
  }
  particle_id_send_list_vec.clear();  // clear the space
  element_id_send_list_vec.clear();
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   If this elem neighbor list is constructed serially, NO communication required.
   This function will return. Otherwise Allgather the local data!
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (this->comm().size() < 2 )
  {
    _elem_neighbor_list.insert(_local_elem_neighbor_list.begin(),
                               _local_elem_neighbor_list.end() );
    STOP_LOG ("build_elem_neighbor_list()", "PointMesh<KDDim>");
    return;
  }
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   If the above loop is performed over local elem on multiple processors,
   we need to allgather the local mapvector to all the processes.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const std::size_t len0 = _local_elem_neighbor_list.size();
  std::vector<std::size_t> buffer_elem_id_list(len0);
  std::vector<std::size_t> buffer_neighbor_list_len(len0);
  std::vector<std::size_t> buffer_elem_neighbor_list;
  //printf("---- length of elem neighbor list on process %d is %lu\n",this->comm().rank(),len0);
  
  // (1) pack the local_elem_neighbor_list, which will be used as "send buffer"
  std::size_t k = 0;
  std::map<const std::size_t, std::vector<std::size_t> >::const_iterator p;
  for (p=_local_elem_neighbor_list.begin(); p!=_local_elem_neighbor_list.end(); ++p)
  {
    buffer_neighbor_list_len[k] = p->second.size();
    buffer_elem_id_list[k]      = p->first;           k++;
    buffer_elem_neighbor_list.insert (buffer_elem_neighbor_list.end(),
                                      p->second.begin(), p->second.end());
  }
  
  
  // (2) allgather the local data
  this->comm().allgather(buffer_elem_id_list);
  this->comm().allgather(buffer_neighbor_list_len);
  this->comm().allgather(buffer_elem_neighbor_list);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   test output on other processes
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//  if (this->comm().rank() == 1)
//  {
//    printf("============ size of buffer_elem_id_list = %lu: printed by rank %d\n",
//           buffer_elem_id_list.size(), this->comm().rank() );
//    for (std::size_t i=0; i<buffer_elem_id_list.size(); ++i)
//      printf("buffer_elem_id_list[%lu] = %lu \n", i, buffer_elem_id_list[i]);
//    printf("\n");
//
//    printf("============ size of buffer_elem_neighbor_list = %lu: printed by rank %d\n",
//           buffer_elem_neighbor_list.size(), this->comm().rank() );
//    for (std::size_t i=0; i<buffer_elem_neighbor_list.size(); ++i)
//      printf("buffer_elem_neighbor_list[%lu] = %lu \n", i, buffer_elem_neighbor_list[i]);
//    printf("\n");
//  }
  
  
  // (3) unpack the buffer data on all the processes
  k = 0;
  for (std::size_t i=0; i<buffer_elem_id_list.size(); ++i)
  {
    const std::size_t elem_id = buffer_elem_id_list[i];
    const std::size_t ilen    = buffer_neighbor_list_len[i];
    std::vector<std::size_t> n_list(ilen);
    for (std::size_t j=0; j<ilen; ++j)
    {
      n_list[j] = buffer_elem_neighbor_list[k];
      k++;
    } // end for j-loop
    
    // make pair and reconstruct _elem_neighbor_list on all processors.
    _elem_neighbor_list.insert ( std::make_pair(elem_id,n_list) );
  } // end for i-loop
  
  
  STOP_LOG ("build_elem_neighbor_list()", "PointMesh<KDDim>");
}



// ======================================================================
template <unsigned int KDDim>
void PointMesh<KDDim>::reinit(const bool& with_hi,
                              bool& neighbor_list_update_flag)
{
  START_LOG ("reinit()", "PointMesh<KDDim>");

  // zero particle force at every time step
  for (std::size_t j=0; j<_num_point_particles; ++j){
    _particles[j]->zero_particle_force();
  }

  // when reinit_neighbor_list is set to be true, which is done every some diffusion
  // steps instead of every time step, we need to reinit particle neighbor list
  if(neighbor_list_update_flag){
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     if the KD tree is already built, we need to re-construct the KD tree in reinit()
     This is necessary at each time step when particles move!
     FIXME: is this necessary when the total number of particles does not change?
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    if ( _kd_tree.get() ) this->clear_kd_tree();
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Then we re-construct the kd-tree after clearing it!
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    this->construct_kd_tree();  
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Construct the particle-particle, elem-particle neighbor list and particle force.
     Loop over ALL particles to calculate their neighbor lists and forces.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    
    for (std::size_t j=0; j<_num_point_particles; ++j)
    {
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Re-init the particle j 
       This includes: set processor_id = -1
                      set elem_id = -1
                      clean neighbor_list
       (note that, this reinit_particle function does not include
       zero_particle_force())
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      _particles[j]->reinit_particle();
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       get the neighbor indices & distance values
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      const Point &tgt( _particles[j]->point() );

      std::vector<std::pair<std::size_t,Real> > IndicesDists0, IndicesDists;

      this->build_particle_neighbor_list(tgt, _is_sorted, IndicesDists);
      /* IndicesDists above returns <particle_id, distance>! */

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Exclude the current particle itself, so that the particle's neighbor list
       does NOT contain itself.
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      const std::size_t pid = _particles[j]->id();

      for (std::size_t i=0; i<IndicesDists.size(); ++i){
        if ( IndicesDists[i].first!= pid){
          IndicesDists0.push_back( IndicesDists[i] );
        }
      }

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       set the neighbor list of the j-th particle
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      _particles[j]->set_neighbor_list (IndicesDists0);
      
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       FIXME: Check if the j-th point particle is out of domain
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    } // end for j-loop over particles
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Construct the element-particle neighbor list and particle-element id map
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    if(with_hi){
      this->build_elem_neighbor_list();
    }    
    // after reinit neighbor list, set the flag to false
    neighbor_list_update_flag = false;
  } // end if (reinit_neighbor_list)


  STOP_LOG ("reinit()", "PointMesh<KDDim>");
}


  
// ======================================================================
template <unsigned int KDDim>
void PointMesh<KDDim>::update_particle_mesh(ParticleMesh<KDDim>* particle_mesh) const
{
  START_LOG("update_particle_mesh()", "PointMesh<KDDim>");
  
  // Extract the point coordinates from the point_mesh object
  std::vector<Point> nodal_vec(_num_point_particles);
  for(std::size_t i=0; i<_num_point_particles; ++i){
    nodal_vec[i] = _particles[i]->point();
  }
  particle_mesh->update_mesh(nodal_vec);
  
  STOP_LOG("update_particle_mesh()", "PointMesh<KDDim>");
}
  
  
  
// ======================================================================
template <unsigned int KDDim>
void PointMesh<KDDim>::update_point_mesh(const ParticleMesh<KDDim>* particle_mesh)
{
  START_LOG("update_point_mesh()", "PointMesh<KDDim>");
  
  // Loop over each Particle
  const std::size_t n_rigid_particles = particle_mesh->num_particles();
  std::size_t start_id = 0;
  for(std::size_t i=0; i<n_rigid_particles; ++i)
  {
    // Extract the point(node) coordiantes from particle_mesh
    std::vector<Point> node_xyz;
    particle_mesh->particles()[i]->extract_nodes(node_xyz);
    
    // Assign the values to the point_mesh
    const std::size_t n_points = node_xyz.size();
    for(std::size_t j=0; j<n_points; ++j)
    {
      _particles[start_id+j]->point() = node_xyz[j];
    }
    
    // Update the start id
    start_id += n_points;
  }
  
  STOP_LOG("update_point_mesh()", "PointMesh<KDDim>");
}
  
  
  

// ======================================================================
template <unsigned int KDDim>
Real PointMesh<KDDim>::search_radius(const std::string & p_e) const
{
  if (p_e=="p") return _search_radius_p;
  else if (p_e=="e") return _search_radius_e;
  else
  {
    printf("PointMesh::search_radius(): the input option must be either p or e!\n");
    libmesh_error();
  }
}





// ======================================================================
// =========================== print functions ==========================
// ======================================================================
template <unsigned int KDDim>
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
  //   for (int j=0; j<3; j++) out_file << _particles[i]->particle_force()[j] <<" ";        
  //   out_file<<"\n";
  // }
  // out_file.close();
  // printf("==> done: output bead info after first system.reinit() to 'out.initial_force_field'\n");
  printf("======================= printing the point information: =======================\n\n");
  printf("There are totally %lu points (search radius = %f)\n",_particles.size(),_search_radius_p);
  for (std::size_t j=0; j<_particles.size(); ++j)
    _particles[j]->print_info();
 printf("========================= end of the point information =========================\n\n");

}



// ======================================================================
template <unsigned int KDDim>
void PointMesh<KDDim>::print_elem_neighbor_list(std::ostream &out) const
{
  printf("======================= printing the element neighbor list: ========================\n\n");
  std::map<const std::size_t, std::vector<std::size_t> >::const_iterator p;
  for (p = _elem_neighbor_list.begin(); p != _elem_neighbor_list.end(); ++p)
  {
    const std::size_t elem_id             = p->first;
    const std::vector<std::size_t> n_list = p->second;
    const Elem* elem            = _mesh.elem(elem_id);
    const Point center_pt       = elem->centroid();
    const Real  hmax            = elem->hmax();
    
    // ----- Scheme 2: print out information using printf (from all the processors) -----
    if (n_list.size()>0)
    {
      printf("========== There are %lu neighbor points around the element %lu, output rank = %u :===\n",
             n_list.size(), elem_id, this->comm().rank() );
      printf("element centroid = (%f, %f, %f), and hmax = %f \n",
             center_pt(0), center_pt(1), center_pt(2), hmax);
      
      if (!n_list.empty())
        for (std::size_t i=0; i<n_list.size(); ++i)  //printf("%lu    ", n_list[i]);
          _particles[ n_list[i] ]->print_info(false);
      else
        out << "There is no neighboring particle around this element! \n";
      printf("\n");
    } // end if
    
  } // end for p-loop
  printf("======================= end of the element neighbor list ======================\n\n");
}


// ======================================================================
template <unsigned int KDDim>
void PointMesh<KDDim>::set_bead_velocity(const std::vector<Real>& vel)
{ 
  START_LOG("PointMesh::set_bead_velocity()", "PointMesh");
  const int dim = 3;
  // assign velocity to 0-th to (n-1)-th particle 
  Point p_velocity(0.);
  for(std::size_t p_id=0; p_id < _num_point_particles; p_id++){
    for (int i=0; i<dim; i++) p_velocity(i) = vel[p_id*dim+i];
    _particles[p_id]->set_particle_velocity(p_velocity);
    _velocity_magnitude[p_id] = p_velocity.norm_sq();
  }
  auto minmax_velocity_magniture = std::minmax_element(_velocity_magnitude.begin(), _velocity_magnitude.end());
  _min_velocity_magnitude = std::sqrt(*minmax_velocity_magniture.first);
  _max_velocity_magnitude = std::sqrt(*minmax_velocity_magniture.second);

  // calculate rigid particle centroid velocity if there are rigid particles
  // create point particle list
  std::size_t point_id = 0;
  for (std::size_t i=0; i<_rigid_particles.size(); ++i)
  {
    const std::size_t n_nodes = _rigid_particles[i]->num_mesh_nodes();
    std::vector<Point> node_velocity(n_nodes);
    for (std::size_t j=0; j<n_nodes; ++j)
    {
      node_velocity[j]=_particles[point_id+j]->particle_velocity();
    } // enf for j-loop
    _rigid_particles[i]->set_node_velocity(node_velocity);
    point_id += n_nodes;
  } // enf for i-loop
  STOP_LOG("PointMesh::set_bead_velocity()", "PointMesh");
}

// ------------------------------------------------------------
// Explicit Instantiations
template class PointMesh<1>;
template class PointMesh<2>;
template class PointMesh<3>;
  
} // end of namespace
