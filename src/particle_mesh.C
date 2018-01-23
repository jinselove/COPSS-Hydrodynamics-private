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



#include "libmesh/libmesh_logging.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"

// C++ Includes
#include <stdio.h>
#include <cstdlib>
#include <iomanip>  // providing parametric manipulators: setw
#include <fstream>
#include <vector>

#include "pm_toolbox.h"
#include "rigid_particle.h"
#include "point_mesh.h"
#include "particle_mesh.h"



namespace libMesh
{
  
  
  
// ======================================================================
template <unsigned int KDDim>
ParticleMesh<KDDim>::ParticleMesh(MeshBase& mesh)
  : ParallelObject(mesh),
    _mesh(mesh),
    _point_list_adaptor(_particles),
    _is_sorted(true),
    _periodic_boundary(NULL)
{
  _dim = _mesh.mesh_dimension();
}


  
// ======================================================================
template <unsigned int KDDim>
ParticleMesh<KDDim>::ParticleMesh(MeshBase& mesh,
                                  const Real& search_radius_p,
                                  const Real& search_radius_e)
: ParallelObject(mesh),
  _mesh(mesh),
  _search_radius_p(search_radius_p),
  _search_radius_e(search_radius_e),
  _point_list_adaptor(_particles),
  _is_sorted(true),
  _periodic_boundary(NULL)
{
  _dim = _mesh.mesh_dimension();
}

  
  
// ======================================================================
template <unsigned int KDDim>
ParticleMesh<KDDim>::ParticleMesh(MeshBase& mesh,
                                  PMPeriodicBoundary& pmpb,
                                  const Real& search_radius_p,
                                  const Real& search_radius_e)
  : ParallelObject(mesh),
  _mesh(mesh),
  _search_radius_p(search_radius_p),
  _search_radius_e(search_radius_e),
  _point_list_adaptor(_particles),
  _is_sorted(true),
  _periodic_boundary(&pmpb)
{
  // do nothing
}
  
  
  
// ======================================================================
template <unsigned int KDDim>
ParticleMesh<KDDim>::~ParticleMesh()
{
  // delete the particle pointers
  for (std::size_t i=0; i<_particles.size(); ++i)
    delete _particles[i];
  
  _particles.clear();
  _elem_neighbor_list.clear();
  _local_elem_neighbor_list.clear();
}

// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::read_particles_data(const std::string& filename,
                                              const std::string& particle_mesh_type)
{
  std::cout << std::endl << "###particle coordinate filename = " << filename << std::endl;
  // Check the existence of the particle input file
  std::ifstream infile;
  infile.open (filename, std::ios_base::in);
  if( !infile.good() )
  {
    std::cout << "***error in read_particles_data(): particle coordinate file does NOT exist!" << std::endl;
    libmesh_error();
  }

  // init variables
  const std::size_t dim = _mesh.mesh_dimension(); // fluid mesh dimension
  unsigned int p_id;
  unsigned  p_type;
  Point p_center;
  Point mag;
  Point rot;
  Real charge=0., epsilon_in=0.;  // parameters for Electrostatics
  Point sedimentation_body_force_density;
  std::size_t n1;

  // Read file line by line
  std::string line_str, str_tmpt;
  std::getline(infile, line_str); // 0. Header line
  infile >> _n_rigid_particles >> line_str; // # number of particles
  infile >> _n_rigid_particle_types >> line_str >> str_tmpt; // # number of particle type ( == number of mesh types == number of mesh files)

  // init
  _mass.resize(_n_rigid_particle_types);
  _particles.resize(_n_rigid_particles);
  _rigid_particle_mesh_files.resize(_n_rigid_particle_types);

  // Read the particle mass
  while (std::getline(infile, line_str)){
    if(line_str=="Masses") break;
  }
  std::getline(infile, line_str); // skip the empty line
  for(std::size_t i=0; i<_n_rigid_particle_types; ++i)
  {
    infile >> n1 >> _mass[i];
  }

  // Read mesh files
  while (std::getline(infile, line_str)){
    if(line_str=="Particle Meshes") break;
  }
  std::getline(infile, line_str); // skip the empty line
  for(std::size_t i=0; i<_n_rigid_particle_types; ++i)
  {
    infile >> n1 >> _rigid_particle_mesh_files[i];
  }
  // Read the rigid particles
  while (std::getline(infile, line_str)){
    if(line_str=="Particles") break;
  }
  std::getline(infile, line_str); // skip the empty line
  for (std::size_t i=0; i<_n_rigid_particles; ++i)
  {   
    infile >> p_id
           >> p_type 
           >> p_center(0) >> p_center(1) >> p_center(2) 
           >> mag(0) >> mag(1) >> mag(2) 
           >> rot(0) >> rot(1) >> rot(2) 
           >> charge 
           >> epsilon_in 
           >> sedimentation_body_force_density(0) >> sedimentation_body_force_density(1) >>sedimentation_body_force_density(2);
    RigidParticle* particle = new RigidParticle(p_id-1, p_type-1, p_center, mag, rot, charge, epsilon_in, sedimentation_body_force_density, this->comm());
    // When particles have different shapes, read mesh from different mesh files!
      //if( (!smesh_exist[p_type]) || i==0 )  particle->extract_surface_mesh( vmesh_file[p_type], smesh_file[p_type] );
    if(particle_mesh_type == "surface_mesh" or particle_mesh_type == "volume_mesh"){
      const std::string& particle_mesh_filename = _rigid_particle_mesh_files[p_type-1];
      bool particle_mesh_exist = PMToolBox::file_exist(particle_mesh_filename);
      if( !particle_mesh_exist )
      {
        std::cout << "***error in read_particles_data(): particle mesh file '" << particle_mesh_filename << "' does NOT exist!" << std::endl;
        libmesh_error();
      }
      else{
        particle->read_mesh(particle_mesh_filename, particle_mesh_type); 
      } 
    } // end if
    else
    {
      std::cout << "***error in read_particles_data(): invalid mesh type!" << std::endl;
      libmesh_error();
    } // end else
    // Assignment
    _particles[i] = particle; 
   // --------------------- test ------------------------------------------
   // if(this->comm().rank()==0 )
   // {
   //   printf("ParticleMesh::read_particles_data: x = %f, y = %f, z = %f. \n",x, y, z );
   //   printf("        radius = %f, relative density = %f. \n",r, den );
   //   printf("        mgx = %f, mgy = %f, mgz = %f. \n",mgx, mgy, mgz);
   //   printf("        rotation angle = (%f, %f, %f). \n",th0, th1, th1 );
   //   printf("        charge = %f, relative_permittivity = %f. \n", charge, epsilon_in );
   //   printf("        MPI_rank = %d\n",this->comm().rank() );
   //   printf("        # of elements is %u\n",particle->mesh().n_elem() );
   //   printf("        # of nodes is %u\n",particle->mesh().n_nodes() );
   // }
  } // end for i-loop
  // Close the file and end the function
  infile.close();
  this->comm().barrier();
  std::cout << "Reading particle data from "<<filename<<" is completed!" << std::endl << std::endl; 
}   

// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::read_particles_data_restart(const std::string& filename,
                                                      const std::string& particle_mesh_type)
{
  std::cout << std::endl << "### Read rigid particle from = " << filename << std::endl
                         << "### Read surface nodes of each particle from 'restart_surface_mesh_particle_*.e"<< std::endl;
  // Check the existence of the particle input file
  std::ifstream infile;
  infile.open (filename, std::ios_base::in);
  if( !infile.good() )
  {
    std::cout << "***error in read_particles_data(): "<<filename<<" does NOT exist!" << std::endl;
    libmesh_error();
  }
  // init variables
  const std::size_t dim = _mesh.mesh_dimension(); // fluid mesh dimension
  unsigned int p_id;
  unsigned  p_type;
  Point p_center(0.);
  Point mag(1.);
  Point rot(1.);
  Real charge=0., epsilon_in=0.;  // parameters for Electrostatics
  Point sedimentation_body_force_density;
  std::size_t n1;
  Real x1;
  // Read file line by line
  std::string line_str, str_tmpt;
  std::getline(infile, line_str); // 0. Header line
  infile >> _n_rigid_particles >> line_str; // # number of particles
  infile >> _n_rigid_particle_types >> line_str >> str_tmpt; // # number of particle type ( == number of mesh types == number of mesh files)
  // init
  _mass.resize(_n_rigid_particle_types);
  _particles.resize(_n_rigid_particles);
  _rigid_particle_mesh_files.resize(_n_rigid_particles);
  // Read the particle mass
  while (std::getline(infile, line_str)){
    if(line_str=="Masses") break;
  }
  std::getline(infile, line_str); // skip the empty line
  for(std::size_t i=0; i<_n_rigid_particle_types; ++i)
  {
    infile >> n1 >> _mass[i];
  }
  // Read mesh files from mesh files saved for restart
  for (std::size_t i=0; i<_n_rigid_particles; ++i){
    std::ostringstream particle_mesh_oss;
    particle_mesh_oss <<"restart_surface_mesh_particle_"<<i+1<<".e";
    _rigid_particle_mesh_files[i] = particle_mesh_oss.str();
  }
  // Read the rigid particles
  while (std::getline(infile, line_str)){
    if(line_str=="Particles") break;
  }
  std::getline(infile, line_str); // skip the empty line
  for (std::size_t i=0; i<_n_rigid_particles; ++i)
  {   
    infile >> p_id
           >> p_type 
           >> x1 >> x1 >> x1 
           >> x1 >> x1 >> x1
           >> x1 >> x1 >> x1 
           >> charge 
           >> epsilon_in 
           >> sedimentation_body_force_density(0) >> sedimentation_body_force_density(1) >>sedimentation_body_force_density(2);
    RigidParticle* particle = new RigidParticle(p_id-1, p_type-1, p_center, mag, rot, charge, epsilon_in, sedimentation_body_force_density, this->comm());
    // When particles have different shapes, read mesh from different mesh files!
      //if( (!smesh_exist[p_type]) || i==0 )  particle->extract_surface_mesh( vmesh_file[p_type], smesh_file[p_type] );
    if(particle_mesh_type == "surface_mesh"){
      const std::string& particle_mesh_filename = _rigid_particle_mesh_files[i];
      bool particle_mesh_exist = PMToolBox::file_exist(particle_mesh_filename);
      if( !particle_mesh_exist )
      {
        std::cout << "***error in read_particles_data_restart(): particle mesh file '" << particle_mesh_filename << "' does NOT exist!" << std::endl;
        libmesh_error();
      }
      else{
        particle->read_mesh(particle_mesh_filename, particle_mesh_type); 
      } 
    } // end if
    else
    {
      std::cout << "***error in read_particles_data(): invalid mesh type!" << std::endl;
      libmesh_error();
    } // end else
    // Assignment
    _particles[i] = particle; 
   // --------------------- test ------------------------------------------
   // if(this->comm().rank()==0 )
   // {
   //   printf("ParticleMesh::read_particles_data: x = %f, y = %f, z = %f. \n",x, y, z );
   //   printf("        radius = %f, relative density = %f. \n",r, den );
   //   printf("        mgx = %f, mgy = %f, mgz = %f. \n",mgx, mgy, mgz);
   //   printf("        rotation angle = (%f, %f, %f). \n",th0, th1, th1 );
   //   printf("        charge = %f, relative_permittivity = %f. \n", charge, epsilon_in );
   //   printf("        MPI_rank = %d\n",this->comm().rank() );
   //   printf("        # of elements is %u\n",particle->mesh().n_elem() );
   //   printf("        # of nodes is %u\n",particle->mesh().n_nodes() );
   // }
  } // end for i-loop
  // Close the file and end the function
  infile.close();
  this->comm().barrier();
  std::cout << "Reading particle data for restart is completed!" << std::endl << std::endl; 
}   
  
// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::generate_random_particles(const std::size_t N,
                                                    const Real bbox_XA, const Real bbox_XB,
                                                    const Real bbox_YA, const Real bbox_YB,
                                                    const Real bbox_ZA, const Real bbox_ZB)
{
  // problem dimension and domain size
  const std::size_t dim = _mesh.mesh_dimension();
  const Real max_range_x = bbox_XB - bbox_XA;
  const Real max_range_y = bbox_YB - bbox_YA;
  const Real max_range_z = bbox_ZB - bbox_ZA;
  const Real r = 1.0, den = 1.0;
  
  // generate random particle coordinates inside the domain, and write out the file
  // We only use rank=0 proccessor to avoid generating different random numbers on each processor.
  if( this->comm().rank()==0 )
  {
    printf("---> test in generate_random_particles: Generating %lu random particles ...\n",N);
    
    // write the particle coordinates into a file
    std::string filename = "particle_random_file.txt";
    int o_width = 5, o_precision = 9;
    
    std::ofstream outfile;
    outfile.open(filename,std::ios_base::out);
    outfile << N << "\n";
    for (size_t i=0;i<N;i++)
    {
      // generate random coordinates
      Real x = max_range_x * (std::rand() % 1000) / 1000 + bbox_XA;
      Real y = max_range_y * (std::rand() % 1000) / 1000 + bbox_YA;
      Real z = max_range_z * (std::rand() % 1000) / 1000 + bbox_ZA;
      if (KDDim==2 || dim==2) z = 0.0;
    
      outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
      outfile.precision(o_precision);   outfile.width(o_width);
      outfile << x << "  " << y << "  " << z << "  " << r << "  "<< den << "  \n";
    } // end loop-i
    
    // close the file
    outfile.close();
    printf("---> test in generate_random_particles: random particle file is created!\n");
  }
  
  this->comm().barrier();
}
  

  
// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>:: generate_random_particles(const std::size_t N,
                                                     const Point& bbox_min,
                                                     const Point& bbox_max)
{
  this->generate_random_particles(N,
                                  bbox_min(0),bbox_max(0),
                                  bbox_min(1),bbox_max(1),
                                  bbox_min(2),bbox_max(2) );
}
  


  
// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::reinit()
{
  START_LOG ("reinit()", "ParticleMesh<KDDim>");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if the KD tree is already built, we need to re-construct the KD tree in reinit()
   This is necessary at each time step when particles move!
   FIXME: is this necessary when the total number of particles does not change?
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if ( _kd_tree.get() )
    this->clear_kd_tree();
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Then we re-construct the kd-tree after clearing it!
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  this->construct_kd_tree();
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Construct the particle-particle, elem-particle neighbor list and particle force.
   Loop over ALL particles to calculate their neighbor lists and forces.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for (std::size_t j=0; j<_particles.size(); ++j)
  {
   //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     //get the neighbor indices & distance values
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    const Point &tgt( _particles[j]->get_centroid() );
    std::vector<std::pair<std::size_t,Real> > IndicesDists0, IndicesDists;
    this->build_particle_neighbor_list(tgt, _is_sorted, IndicesDists);
    /* IndicesDists above returns <particle_id, distance>! */
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Exclude the current particle itself, so that the particle's neighbor list
     does NOT contain itself.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    const std::size_t pid = _particles[j]->id();
    for (std::size_t i=0; i<IndicesDists.size(); ++i)
      if ( IndicesDists[i].first!= pid)
        IndicesDists0.push_back( IndicesDists[i] );
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     set the neighbor list of the j-th particle
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    _particles[j]->set_neighbor_list (IndicesDists0);
    
  } // end for j-loop over particles
  //printf("--->debug: ParticleMesh::reinit() I am here 000 \n");
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Construct the element-particle neighbor list and particle-element id map
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //this->build_elem_neighbor_list();
  
  
  STOP_LOG ("reinit()", "ParticleMesh<KDDim>");
}
  
  
  
// ======================================================================
template <unsigned int KDDim>
Real ParticleMesh<KDDim>::search_radius(const std::string & p_e) const
{
  if (p_e=="p") return _search_radius_p;
  else if (p_e=="e") return _search_radius_e;
  else
  {
    printf("ParticleMesh::search_radius(): the input option must be either p or e!\n");
    libmesh_error();
  }
}
  
  
  
// ======================================================================
template <unsigned int KDDim>
std::size_t ParticleMesh<KDDim>::num_mesh_points() const
{
  START_LOG ("num_mesh_points()", "ParticleMesh<KDDim>");
  
  std::size_t n_points = 0;
  for (std::size_t i=0; i<_particles.size(); ++i)
    n_points += _particles[i]->num_mesh_nodes();
  
  STOP_LOG ("num_mesh_points()", "ParticleMesh<KDDim>");
  return n_points;
}

  
  
// ======================================================================
// =========================== print functions ==========================
// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::print_particle_info() const
{
  printf("======================= printing the particle information: =======================\n\n");
  printf("There are totally %lu particles (search radius = %f)\n",_particles.size(),_search_radius_p);
  for (std::size_t j=0; j<_particles.size(); ++j)
    _particles[j]->print_info();
  
  printf("========================= end of the particle information =========================\n\n");
} 
  
  

  
  
// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::update_mesh(const std::vector<Point>& nodal_pos,
                                      const std::vector<Point>& nodal_vel)
{
  START_LOG ("update_mesh()", "ParticleMesh<KDDim>");
  // loop over all rigid particles
  std::size_t node_start_id = 0;
  for (int particle_id = 0; particle_id < _n_rigid_particles; particle_id++){
    std::size_t num_mesh_nodes = _particles[particle_id]->num_mesh_nodes();
    std::vector<Point> nodal_points(num_mesh_nodes);
    std::vector<Point> nodal_velocities(num_mesh_nodes);
    // extract positions from nodal_pos to nodal_points vector
    for (int node_id = 0; node_id < num_mesh_nodes; node_id++){
      nodal_points[node_id] = nodal_pos[node_start_id+node_id];
      nodal_velocities[node_id] = nodal_vel[node_start_id+node_id];
    } // end loop over node id
    _particles[particle_id]->update_mesh(nodal_points, nodal_velocities);
    node_start_id+=num_mesh_nodes;
  }//end loop over particle id

  STOP_LOG ("update_mesh()", "ParticleMesh<KDDim>");
}
  
// // ======================================================================
// template <unsigned int KDDim>
// void ParticleMesh<KDDim>::update_point_mesh(PointMesh<KDDim>* point_mesh) const
// {
//   START_LOG ("update_point_mesh()", "ParticleMesh<KDDim>");
  
//   // Loop over each Particle
//   const std::size_t n_particles = this->num_particles();
//   std::size_t start_id = 0;
//   for(std::size_t i=0; i<n_particles; ++i)
//   {
//     // Extract the point(node) coordiantes from particle_mesh
//     std::vector<Point> node_xyz;
//     _particles[i]->extract_nodes(node_xyz);
    
//     // Assign the values to the point_mesh
//     const std::size_t n_points = node_xyz.size();
//     for(std::size_t j=0; j<n_points; ++j)
//     {
//       point_mesh->particles()[start_id+j]->point() = node_xyz[j];
//     }
    
//     // Update the start id
//     start_id += n_points;
//   }
  
  
//   STOP_LOG ("update_point_mesh()", "ParticleMesh<KDDim>");
// }
  
  
  
// // ======================================================================
// template <unsigned int KDDim>
// void ParticleMesh<KDDim>::update_particle_mesh(const PointMesh<KDDim>* point_mesh)
// {
//   START_LOG ("update_particle_mesh()", "ParticleMesh<KDDim>");
  
//   // Extract the point coordinates from the point_mesh object
//   const std::size_t  n_points = point_mesh->num_particles();
//   std::vector<Point> nodal_vec(n_points);
//   for(std::size_t i=0; i<n_points; ++i){
//     nodal_vec[i] = point_mesh->particles()[i]->point();
//   }
  
//   // Update the mesh nodes using the extracted values
//   this->update_mesh(nodal_vec);
  
//   STOP_LOG ("update_particle_mesh()", "ParticleMesh<KDDim>");
// }


// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::zero_node_force()
{
  START_LOG ("zero_particle_force_density", "ParticleMesh<KDDim>");
  for (int i=0; i<_n_rigid_particles;i++) _particles[i]->zero_node_force();
  
  STOP_LOG ("zero_particle_force_density", "ParticleMesh<KDDim>");
}
  
  
  
// ======================================================================
template <unsigned int KDDim>
std::vector<Real> ParticleMesh<KDDim>::mesh_size() const
{
  START_LOG ("mesh_size()", "ParticleMesh<KDDim>");
  
  const std::size_t np = this->num_particles();
  std::vector<Real> hmin_max(2);  // store hmin and hmax
  for(std::size_t i=0; i<np; ++i)
  {
    if(i==0){
      hmin_max = _particles[i]->mesh_size();
    }
    else{
      const std::vector<Real> vals = _particles[i]->mesh_size();
      if( vals[0]<hmin_max[0] ) hmin_max[0] = vals[0];
      if( vals[1]>hmin_max[1] ) hmin_max[1] = vals[1];
    }
  }
  
  STOP_LOG ("mesh_size()", "ParticleMesh<KDDim>");
  return hmin_max;
}
  
  
  
  
// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::volume_conservation()
{
  START_LOG ("volume_conservation()", "ParticleMesh<KDDim>");
  
  // Loop over all the particles, and update their nodal values
  for (std::size_t i=0; i<_particles.size(); ++i)
  {
    _particles[i]->volume_conservation();
  }
  
  STOP_LOG ("volume_conservation()", "ParticleMesh<KDDim>");
}
  

// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::initial_particle_center_of_mass(std::vector<Point>& center0) const
{
  START_LOG("initial_center_of_mass()", "ParticleMesh<KDDim>");
  center0.resize(_n_rigid_particles);
  for (std::size_t i=0; i<_n_rigid_particles; i++){
    center0[i] = _particles[i]->get_centroid0();
  }
  STOP_LOG("initial_center_of_mass()", "ParticleMesh<KDDim>");
}

// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::write_particle(const unsigned int& step_id,
                                          const unsigned int& o_step,
                                          const Real& real_time,
                                          const std::vector<std::string>& output_file,
                                          unsigned int comm_in_rank) const
{
  START_LOG("write_particle()", "ParticleMesh<KDDim>");
  this->write_time(step_id, o_step, real_time, comm_in_rank);
  this->write_particle_mesh_restart();
  for (int i = 0; i < output_file.size(); i++){
    if(output_file[i] == "equation_systems") {
      // this output has been written in Copss.C
    }
    else if(output_file[i] == "trajectory") this -> write_particle_trajectory(o_step, comm_in_rank);
    else if (output_file[i] == "particle_mesh") this->write_particle_mesh(o_step);
    else if(output_file[i] == "mean_square_displacement"){
        std::cout <<"Error: there is difficulty to calculate msd of rigid particles from ROUT, fix it before output msd" << std::endl;
        libmesh_error();
        // this -> write_particle_msd(step_id, o_step,center0, lvec, comm_in_rank);
      
    }
    else {
      std::cout << "unsupported output_file content: (" << output_file[i] <<")" << std::endl; 
      libmesh_error();
    }
  }  

  STOP_LOG("write_particle()", "ParticleMesh<KDDim>");
}

// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::write_time(const unsigned int& step_id,
                              const unsigned int& o_step,
                              const Real& real_time,
                              unsigned int comm_in_rank) const
{
  // write step and real time
  std::ofstream out_file;  
  if(comm_in_rank == 0){
    if(step_id == 0){
      out_file.open("out.time", std::ios_base::out);
      out_file << "step_id" <<"  " <<"o_step" <<"  " << "real_time" << "\n";      
    }
    else{
      out_file.open("out.time", std::ios_base::app);
    }
    out_file.precision(o_precision);
    out_file <<step_id <<"  "<< o_step <<"  " <<real_time <<"\n";
    out_file.close();  
  }
}

// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::write_particle_trajectory(const unsigned int& o_step,
                                                    unsigned int comm_in_rank) const
{
  START_LOG ("write_particle_trajectory()", "ParticleMesh<KDDim>");
  std::ostringstream oss;
  oss << "output_particle_" << o_step <<".csv";
  std::ofstream out_file;
  if (comm_in_rank == 0){
    out_file.open(oss.str(), std::ios_base::out);
    // write out the csv file  
    // POINT data
    out_file <<"scalar x_coord y_coord z_coord x_vel y_vel z_vel x_force y_force z_force\n";
    out_file.precision(o_precision);
    for(std::size_t i=0; i<_n_rigid_particles; ++i)
    {
      out_file << i << " ";
      // write position
      for(std::size_t j=0; j<_dim; ++j){
        out_file << _particles[i]->get_centroid()(j) << " ";
      }
      // write velocity
      for (std::size_t j=0; j<_dim; ++j){
        out_file <<_particles[i]->centroid_velocity()(j) << " ";
      }
      // write force
      for (std::size_t j=0; j<_dim; ++j){
        out_file <<_particles[i]->centroid_force()(j) <<" ";
      }
      out_file <<"\n";
    }
    out_file << "\n";
    out_file.close();
  } // end if comm_in_rank == 0

  STOP_LOG ("write_particle_trajectory()", "ParticleMesh<KDDim>");
}

// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::write_particle_mesh(const unsigned int& o_step) const
{
  START_LOG ("write_particle_mesh()", "ParticleMesh<KDDim>");
  this->comm().barrier();
  for (std::size_t particle_id = 0; particle_id < _n_rigid_particles; particle_id++){
    std::ostringstream surface_mesh_filename;
    surface_mesh_filename << "surface_mesh_particle_"<<particle_id+1<<"_step_"<<o_step<<".e";
    _particles[particle_id]->mesh().write(surface_mesh_filename.str());
  }
  STOP_LOG ("write_particle_mesh()", "ParticleMesh<KDDim>");
}

// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::write_particle_mesh_restart() const
{
  START_LOG ("write_particle_mesh_restart()", "ParticleMesh<KDDim>");
  this->comm().barrier();
  for (std::size_t particle_id = 0; particle_id < _n_rigid_particles; particle_id++){
    std::ostringstream surface_mesh_filename;
    surface_mesh_filename << "restart_surface_mesh_particle_"<<particle_id+1 << ".e";
    _particles[particle_id]->mesh().write(surface_mesh_filename.str());
  }
  STOP_LOG ("write_particle_mesh_restart()", "ParticleMesh<KDDim>");
}
  

// ======================================================================
template <unsigned int KDDim>
SerialMesh& ParticleMesh<KDDim>::stitched_mesh()
{
  START_LOG ("merge_particle_mesh()", "ParticleMesh<KDDim>");

  SerialMesh& _particles_mesh (_particles[0]->mesh());

  // assign subdomain_id of 1st particle to 0
  MeshBase::element_iterator     el = _particles_mesh.elements_begin();
  MeshBase::element_iterator end_el = _particles_mesh.elements_end();
  for ( ; el != end_el; ++el)
  {
     Elem* elem = *el;
     elem->subdomain_id() = 0;
  }

  const std::size_t n_particles = this->num_particles();
  if(n_particles>1)
  {
    for(std::size_t i=1; i<this->num_particles(); ++i)
    {
      // intermediate mesh to store each particle's mesh
      SerialMesh mesh(_particles[i]->mesh());

      // assign subdomain_id to particle_id
      el     = mesh.elements_begin();
      end_el = mesh.elements_end();

      for ( ; el != end_el; ++el)
      {
         Elem* elem = *el;
         elem->subdomain_id() = i;
      }

      _particles_mesh.stitch_meshes(mesh,0,0);
    }
  }

  STOP_LOG ("merge_particle_mesh()", "ParticleMesh<KDDim>");

  // return the stitched mesh
  return _particles_mesh;
}
  

/********************************************************************************
    Uncomment play around the section below to construct neighbor_list for 
    rigidParticle-rigidParticle and rigidParticle-elem
********************************************************************************/    
   
   
// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::construct_kd_tree ()
{
#ifdef LIBMESH_HAVE_NANOFLANN
  
  START_LOG ("construct_kd_tree()", "ParticleMesh<KDDim>");
  
  // Initialize underlying KD tree if this is not constructed.
  if (_kd_tree.get() == NULL)
    _kd_tree.reset(new kd_tree_t(KDDim,
                                 _point_list_adaptor,
                                 nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */)));
  
  libmesh_assert (_kd_tree.get() != NULL);
  
  _kd_tree->buildIndex();
  
  STOP_LOG ("construct_kd_tree()", "ParticleMesh<KDDim>");
#endif
}

  

// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::clear_kd_tree()
{
#ifdef LIBMESH_HAVE_NANOFLANN
  if (_kd_tree.get())   // If exist, delete the KD Tree and start fresh
    _kd_tree.reset (NULL);
#endif
}

  

// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::build_particle_neighbor_list(const Point &tgt,
                                                       const bool is_sorted,
                                                       std::vector<std::pair<std::size_t,Real> >& IndicesDists)
{
#ifdef LIBMESH_HAVE_NANOFLANN
  START_LOG ("build_particle_neighbor_list(point)", "ParticleMesh<KDDim>");
  
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
  for (std::size_t j=0; j<IndicesDists.size(); ++j)
    IndicesDists[j].second = std::sqrt( IndicesDists[j].second );
  
  
  /* ------------------------------------------------------------------------
   * if the periodic boundary condition is applied, we must find the neighbor
   * list around its image particles for computing the interaction forces.
   *
   * In the following implementation, we assume that the domain size MUST
   * be larger than 4X search_radius so that only one image in this direction
   * needs to be considered. This is typically reasonable in realistic simulations!
   * ------------------------------------------------------------------------*/
  // loop over each direction to find its images
  std::size_t  NImage = 0;
  if(KDDim==2) NImage = 3;
  if(KDDim==3) NImage = 7;
  for (std::size_t i=0; i<NImage; ++i)
  {
    Point im_pt;
    bool has_image = _periodic_boundary->get_image_point(tgt,_search_radius_p,i,im_pt);
    
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
//      if (this->comm().rank()==0 && IndicesDists_image.size()>0)
//      {
//        printf("--->test in build_particle_neighbor_list() i=%lu:\n",i);
//        printf("    Original point = (%f %f %f) \n", tgt(0), tgt(1), tgt(2) );
//        printf("    Image point    = (%f %f %f) has the neighbors:\n", im_pt(0),im_pt(1),im_pt(2));
//        for (std::size_t j=0; j<IndicesDists_image.size(); ++j)
//          printf("      particle id %lu, distance to the image point is %f\n",
//                 IndicesDists_image[j].first, IndicesDists_image[j].second);
//      }
      // -------------------------- test --------------------------
      
    } // end if

  } // end for i-loop
  /* -----------------------------------------------------------------------*/
  
  STOP_LOG ("build_particle_neighbor_list(point)", "ParticleMesh<KDDim>");
#endif
}


  
// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::build_particle_neighbor_list()
{
#ifdef LIBMESH_HAVE_NANOFLANN

  // do a radius search to build the neighbor list for ALL particles.
  for (std::size_t j=0; j<_particles.size(); ++j)
  {
    // get the neighbor indices & distance values
    const Point &tgt( _particles[j]->get_centroid() );
    std::vector<std::pair<std::size_t,Real> > IndicesDists0, IndicesDists;
    this->build_particle_neighbor_list(tgt, _is_sorted, IndicesDists);
    /* IndicesDists above returns <particle_id, distance>! */
    
    // Exclude the current particle itself!
    const std::size_t pid = _particles[j]->id();
    for (std::size_t i=0; i<IndicesDists.size(); ++i)
      if ( IndicesDists[i].first!= pid)
        IndicesDists0.push_back( IndicesDists[i] );
    
    // set the neighbor list of the j-th particle
    _particles[j]->set_neighbor_list (IndicesDists0);
  } // end for j-loop over particles
  
#endif
}

  
// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::build_particle_neighbor_list_naively()
{
  START_LOG ("build_particle_neighbor_list_naively()", "ParticleMesh<KDDim>");
  
  // do a radius search to build the neighbor list for ALL particles.
  for (std::size_t i=0; i<_particles.size(); ++i)
  {
    std::vector<std::pair<std::size_t,Real> > IndicesDists;
    
    const Point pt0( _particles[i]->get_centroid() );
    for (std::size_t j=0; j<_particles.size(); ++j)
    {
      const Point ptj( _particles[j]->get_centroid() );
      const Point pt_ij = ptj - pt0;
      const Real dist = pt_ij.norm(); // the real distance
      
      if( dist < _search_radius_p )
        IndicesDists.push_back ( std::make_pair(_particles[j]->id(), dist) );
    }
    
    // if needed, sort the particle neighbor list.
    if (_is_sorted)
    {
      // ***sort the particle neighbor list (NOT implemented here)
      printf("***warning: the particle neighbor list is not sorted in this function!");
    }
    
    // set the neighbor list of the j-th particle (the result is not sorted)
    _particles[i]->set_neighbor_list (IndicesDists);
  }
  
  STOP_LOG ("build_particle_neighbor_list_naively()", "ParticleMesh<KDDim>");
}
  
  
  
// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>:: build_elem_neighbor_list(const Elem* elem,
                                                    const bool is_sorted,
                                                    std::vector<std::size_t>& n_list)
{
  // This function must be called on every processor.
  parallel_object_only();
  
#ifdef LIBMESH_HAVE_NANOFLANN
  START_LOG ("build_elem_neighbor_list(elem)", "ParticleMesh<KDDim>");
  
  
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
  // loop over each direction to find its images
  std::size_t  NImage = 0;
  if(KDDim==2) NImage = 3;
  if(KDDim==3) NImage = 7;
  for (std::size_t i=0; i<NImage; ++i)
  {
    if ( _periodic_boundary )
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
        for (std::size_t j=0; j<IndicesDists_image.size(); ++j)
        IndicesDists.push_back( IndicesDists_image[j] );
      } // end if
    } // end if ( _periodic_boundary )
    
  } // end for i-loop
  /* -----------------------------------------------------------------------*/
  
  
  // The distance is L2 form, which is the distance square, so we should take sqrt().
  // However, we don't need to store the distance values in this list, so it is ignored
  const std::size_t np = IndicesDists.size();
  if(np > 0)
  {
    n_list.resize (np);
    for (std::size_t j=0; j<np; ++j)
      n_list[j] = IndicesDists[j].first;
  }
  
  STOP_LOG ("build_elem_neighbor_list(elem)", "ParticleMesh<KDDim>");
#endif
}
  
  
  
// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::build_elem_neighbor_list()
{
  START_LOG ("build_elem_neighbor_list()", "ParticleMesh<KDDim>");
  
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
  std::vector<std::size_t> particle_id_send_list_vec, element_id_send_list_vec;
  
  
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
//    printf("--->debug: ParticleMesh::build_elem_neighbor_list()11 elem_id = %lu \n",elem_id);
    std::vector<std::size_t> n_list;
    this->build_elem_neighbor_list(elem, _is_sorted, n_list);
    call_count++;
//    printf("--->debug: ParticleMesh::build_elem_neighbor_list()22 elem_id = %lu \n",elem_id);
    
    
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
    for (std::size_t j=0; j<n_list.size(); ++j)
    {
      const std::size_t pid = n_list[j];
      const Point pt        = _particles[pid]->get_centroid();
      const bool  inside    = elem->contains_point(pt);
//      printf("--->debug: ParticleMesh::build_elem_neighbor_list()1 p_id = %lu, pt = (%f,%f,%f)\n",
//             elem_id, pt(0), pt(1), pt(2) );
      
      if (inside )
      {
        particle_id_send_list_vec.push_back(pid);
        element_id_send_list_vec. push_back(elem_id);
        //printf("*****pid = %lu, eid = %lu, processor_id = %i\n",pid,elem_id,elem->processor_id() );
        //printf("*****point xyz = (%f %f %f)"\n\n",pt(0),pt(1),pt(2)); //elem->print_info();
      } // end if
      
    } // end for j-loop
    
    
  } // end for elem-loop
  
  this->comm().barrier();
  printf("--->test in build_elem_neighbor_list() call_count = %lu on rank = %d\n",
         call_count,this->comm().rank() );
  printf("--->debug: ParticleMesh::build_elem_neighbor_list()3 I am here 0000 \n");
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Allgather the particle_send_list and element_send_list
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  this->comm().allgather(particle_id_send_list_vec);
  this->comm().allgather(element_id_send_list_vec);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set the element id for each particle. *** No elem id for finite size particles!
   This operation is on all the processors
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//  for(std::size_t i=0; i<particle_id_send_list_vec.size(); ++i)
//  {
//    const std::size_t pid = particle_id_send_list_vec[i];
//    const std::size_t eid =  element_id_send_list_vec[i];
//    _particles[pid]->set_elem_id(eid);
//  }
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
  
  
  STOP_LOG ("build_elem_neighbor_list()", "ParticleMesh<KDDim>");
}
  
// ======================================================================
template <unsigned int KDDim>
void ParticleMesh<KDDim>::print_elem_neighbor_list(std::ostream &out) const
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
      printf("========== There are %lu neighbor particles around the element %lu, output rank = %u :===\n",
             n_list.size(), elem_id, this->comm().rank() );
      printf("element centroid = (%f, %f, %f), and hmax = %f \n",
             center_pt(0), center_pt(1), center_pt(2), hmax);
      
      if (!n_list.empty())
        for (std::size_t i=0; i<n_list.size(); ++i)  //printf("%lu    ", n_list[i]);
          _particles[ n_list[i] ]->print_info();
      else
        out << "There is no neighboring particle around this element! \n";
      printf("\n");
    } // end if
    
  } // end for p-loop
  printf("======================= end of the element neighbor list ======================\n\n");
}
   

// ------------------------------------------------------------
// Explicit Instantiations
template class ParticleMesh<1>;
template class ParticleMesh<2>;
template class ParticleMesh<3>;

} // end of namespace

