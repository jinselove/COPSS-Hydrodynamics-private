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
#include <iomanip> // providing parametric manipulators: setw
#include <fstream>
#include <vector>

#include "pm_toolbox.h"
#include "rigid_particle.h"
#include "point_mesh.h"
#include "particle_mesh.h"


namespace libMesh {
// ======================================================================
template<unsigned int KDDim>
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
template<unsigned int KDDim>
ParticleMesh<KDDim>::ParticleMesh(MeshBase  & mesh,
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
template<unsigned int KDDim>
ParticleMesh<KDDim>::ParticleMesh(MeshBase          & mesh,
                                  PMPeriodicBoundary& pmpb,
                                  const Real        & search_radius_p,
                                  const Real        & search_radius_e)
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
template<unsigned int KDDim>
ParticleMesh<KDDim>::~ParticleMesh()
{
  // delete the particle pointers
  for (std::size_t i = 0; i < _particles.size(); ++i) delete _particles[i];

  _particles.clear();
}

// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::read_particles_data(const std::string& filename,
                                              const std::string& particle_mesh_type)
{
  std::cout << std::endl << "###particle coordinate filename = " << filename <<
    std::endl;

  // Check the existence of the particle input file
  std::ifstream infile;

  infile.open(filename, std::ios_base::in);

  if (!infile.good())
  {
    std::cout <<
      "***error in read_particles_data(): particle coordinate file does NOT exist!"
              << std::endl;
    libmesh_error();
  }

  // init variables
  const std::size_t dim = _mesh.mesh_dimension(); // fluid mesh dimension
  unsigned int  p_id;
  unsigned  int p_type;
  unsigned int  mesh_id;
  Point p_center;
  Point mag;
  Point rot;
  Real  charge = 0., epsilon_in = 0.; // parameters for Electrostatics
  Point sedimentation_body_force_density;
  std::size_t n1;

  // Read file line by line
  std::string line_str, str_tmpt;
  std::getline(infile, line_str);                                 // 0. Header
                                                                  // line
  infile >> _n_rigid_particles >> line_str;                       // # number of
                                                                  // particles
  infile >> _n_rigid_particle_types >> line_str >> str_tmpt;      // # number of
                                                                  // particle
                                                                  // type
  infile >> _n_rigid_particle_mesh_types >> line_str >> str_tmpt; // # number of
                                                                  // particle
                                                                  // mesh type

  // init
  _mass.resize(_n_rigid_particle_types);
  _particles.resize(_n_rigid_particles);
  _rigid_particle_mesh_files.resize(_n_rigid_particle_types);

  // Read the particle mass
  while (std::getline(infile, line_str)) {
    if (line_str == "Masses") break;
  }
  std::getline(infile, line_str); // skip the empty line

  for (std::size_t i = 0; i < _n_rigid_particle_types; ++i)
  {
    infile >> n1 >> _mass[i];
  }

  // Read mesh files
  while (std::getline(infile, line_str)) {
    if (line_str == "Particle Meshes") break;
  }
  std::getline(infile, line_str); // skip the empty line

  for (std::size_t i = 0; i < _n_rigid_particle_mesh_types; ++i)
  {
    infile >> n1 >> _rigid_particle_mesh_files[i];
  }

  // Read the rigid particles
  while (std::getline(infile, line_str)) {
    if (line_str == "Particles") break;
  }
  std::getline(infile, line_str); // skip the empty line

  for (std::size_t i = 0; i < _n_rigid_particles; ++i)
  {
    infile >> p_id
    >> p_type
    >> mesh_id
    >> p_center(0) >> p_center(1) >> p_center(2)
    >> mag(0) >> mag(1) >> mag(2)
    >> rot(0) >> rot(1) >> rot(2)
    >> charge
    >> epsilon_in
    >> sedimentation_body_force_density(0) >>
    sedimentation_body_force_density(1) >> sedimentation_body_force_density(2);
    RigidParticle *particle = new RigidParticle(p_id - 1,
                                                p_type - 1,
                                                p_center,
                                                mag,
                                                rot,
                                                charge,
                                                epsilon_in,
                                                sedimentation_body_force_density,
                                                this->comm());

    // When particles have different shapes, read mesh from different mesh
    // files!
    // if( (!smesh_exist[p_type]) || i==0 )  particle->extract_surface_mesh(
    // vmesh_file[p_type], smesh_file[p_type] );
    if (particle_mesh_type == "surface_mesh" or particle_mesh_type ==
        "volume_mesh") {
      const std::string& particle_mesh_filename =
        _rigid_particle_mesh_files[mesh_id - 1];
      bool particle_mesh_exist = PMToolBox::file_exist(
        particle_mesh_filename);

      if (!particle_mesh_exist)
      {
        std::cout << "***error in read_particles_data(): particle mesh file '" <<
          particle_mesh_filename << "' does NOT exist!" << std::endl;
        libmesh_error();
      }
      else {
        particle->read_mesh(particle_mesh_filename, particle_mesh_type);
      }
    } // end if
    else
    {
      std::cout << "***error in read_particles_data(): invalid mesh type!" <<
        std::endl;
      libmesh_error();
    } // end else
    // Assignment
    _particles[i] = particle;

    // --------------------- test ------------------------------------------
    // if(this->comm().rank()==0 )
    // {
    //   printf("ParticleMesh::read_particles_data: x = %f, y = %f, z = %f.
    // \n",x, y, z );
    //   printf("        radius = %f, relative density = %f. \n",r, den );
    //   printf("        mgx = %f, mgy = %f, mgz = %f. \n",mgx, mgy, mgz);
    //   printf("        rotation angle = (%f, %f, %f). \n",th0, th1, th1 );
    //   printf("        charge = %f, relative_permittivity = %f. \n", charge,
    // epsilon_in );
    //   printf("        MPI_rank = %d\n",this->comm().rank() );
    //   printf("        # of elements is %u\n",particle->mesh().n_elem() );
    //   printf("        # of nodes is %u\n",particle->mesh().n_nodes() );
    // }
  } // end for i-loop
  // Close the file and end the function
  infile.close();
  this->comm().barrier();
  std::cout << "Reading particle data from " << filename << " is completed!" <<
    std::endl << std::endl;
}

// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::read_particles_data_restart(
  const std::string & filename,
  const unsigned int& o_step)
{
  std::cout << std::endl << "### Read rigid particle from = " << filename <<
    ", restart step = " << o_step << std::endl;

  // Check the existence of the particle input file
  std::ifstream infile;

  infile.open(filename, std::ios_base::in);

  if (!infile.good())
  {
    std::cout << "***error in read_particles_data_restart(): " << filename <<
      " does NOT exist!" << std::endl;
    libmesh_error();
  }
  std::size_t tmp_n;
  std::string line_str, str_tmpt;

  // std::getline(infile, line_str); // skip header line
  // ignore lines until "#output_step_id,o_step"
  std::ostringstream oss;
  oss << "#output_step_id," << o_step;

  while (std::getline(infile, line_str)) {
    if (line_str == oss.str()) break;
  }

  // read surface node
  for (std::size_t i = 0; i < _n_rigid_particles; ++i)
  {
    std::vector<Point> node_pos;
    std::size_t num_mesh_nodes = _particles[i]->num_mesh_nodes();
    node_pos.resize(num_mesh_nodes);
    Real node_center_dist;

    for (std::size_t j = 0; j < num_mesh_nodes; j++) {
      infile >> tmp_n
      >> tmp_n
      >> node_pos[j](0) >> node_pos[j](1) >> node_pos[j](2)
      >> node_center_dist;
    }
    _particles[i]->update_mesh(node_pos);
  }
  infile.close();
  this->comm().barrier();
  std::cout << "Update particle surface node position from " << filename <<
    ", restart step = " << o_step << " is completed\n";
}

// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::read_particles_data_restart(const std::string& filename,
                                                      const std::string& particle_mesh_type)
{
  std::cout << std::endl << "### Read rigid particle from = " << filename <<
    std::endl
            <<
    "### Read surface nodes of each particle from 'restart_surface_mesh_particle_*.e"
            << std::endl;

  // Check the existence of the particle input file
  std::ifstream infile;

  infile.open(filename, std::ios_base::in);

  if (!infile.good())
  {
    std::cout << "***error in read_particles_data_restart(): " << filename <<
      " does NOT exist!" << std::endl;
    libmesh_error();
  }

  // init variables
  const std::size_t dim = _mesh.mesh_dimension(); // fluid mesh dimension
  unsigned int p_id;
  unsigned     p_type;
  Point p_center(0.);
  Point mag(1.);
  Point rot(1.);
  Real  charge = 0., epsilon_in = 0.; // parameters for Electrostatics
  Point sedimentation_body_force_density;
  std::size_t n1;
  Real x1;

  // Read file line by line
  std::string line_str, str_tmpt;
  std::getline(infile, line_str);                            // 0. Header line
  infile >> _n_rigid_particles >> line_str;                  // # number of
                                                             // particles
  infile >> _n_rigid_particle_types >> line_str >> str_tmpt; // # number of
                                                             // particle type (
                                                             // == number of
                                                             // mesh types ==
                                                             // number of mesh
                                                             // files)
  // init
  _mass.resize(_n_rigid_particle_types);
  _particles.resize(_n_rigid_particles);
  _rigid_particle_mesh_files.resize(_n_rigid_particles);

  // Read the particle mass
  while (std::getline(infile, line_str)) {
    if (line_str == "Masses") break;
  }
  std::getline(infile, line_str); // skip the empty line

  for (std::size_t i = 0; i < _n_rigid_particle_types; ++i)
  {
    infile >> n1 >> _mass[i];
  }

  // Read mesh files from mesh files saved for restart
  for (std::size_t i = 0; i < _n_rigid_particles; ++i) {
    std::ostringstream particle_mesh_oss;
    particle_mesh_oss << "restart_surface_mesh_particle_" << i + 1 << ".e";
    _rigid_particle_mesh_files[i] = particle_mesh_oss.str();
  }

  // Read the rigid particles
  while (std::getline(infile, line_str)) {
    if (line_str == "Particles") break;
  }
  std::getline(infile, line_str); // skip the empty line

  for (std::size_t i = 0; i < _n_rigid_particles; ++i)
  {
    infile >> p_id
    >> p_type
    >> x1 >> x1 >> x1
    >> x1 >> x1 >> x1
    >> x1 >> x1 >> x1
    >> charge
    >> epsilon_in
    >> sedimentation_body_force_density(0) >>
    sedimentation_body_force_density(1) >> sedimentation_body_force_density(2);
    RigidParticle *particle = new RigidParticle(p_id - 1,
                                                p_type - 1,
                                                p_center,
                                                mag,
                                                rot,
                                                charge,
                                                epsilon_in,
                                                sedimentation_body_force_density,
                                                this->comm());

    // When particles have different shapes, read mesh from different mesh
    // files!
    // if( (!smesh_exist[p_type]) || i==0 )  particle->extract_surface_mesh(
    // vmesh_file[p_type], smesh_file[p_type] );
    if (particle_mesh_type == "surface_mesh") {
      const std::string& particle_mesh_filename = _rigid_particle_mesh_files[i];
      bool particle_mesh_exist                  = PMToolBox::file_exist(
        particle_mesh_filename);

      if (!particle_mesh_exist)
      {
        std::cout <<
          "***error in read_particles_data_restart(): particle mesh file '" <<
          particle_mesh_filename << "' does NOT exist!" << std::endl;
        libmesh_error();
      }
      else {
        particle->read_mesh(particle_mesh_filename, particle_mesh_type);
      }
    } // end if
    else
    {
      std::cout << "***error in read_particles_data(): invalid mesh type!" <<
        std::endl;
      libmesh_error();
    } // end else
    // Assignment
    _particles[i] = particle;

    // --------------------- test ------------------------------------------
    // if(this->comm().rank()==0 )
    // {
    //   printf("ParticleMesh::read_particles_data: x = %f, y = %f, z = %f.
    // \n",x, y, z );
    //   printf("        radius = %f, relative density = %f. \n",r, den );
    //   printf("        mgx = %f, mgy = %f, mgz = %f. \n",mgx, mgy, mgz);
    //   printf("        rotation angle = (%f, %f, %f). \n",th0, th1, th1 );
    //   printf("        charge = %f, relative_permittivity = %f. \n", charge,
    // epsilon_in );
    //   printf("        MPI_rank = %d\n",this->comm().rank() );
    //   printf("        # of elements is %u\n",particle->mesh().n_elem() );
    //   printf("        # of nodes is %u\n",particle->mesh().n_nodes() );
    // }
  } // end for i-loop
  // Close the file and end the function
  infile.close();
  this->comm().barrier();
  std::cout << "Reading particle data for restart is completed!" << std::endl <<
    std::endl;
}

// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::generate_random_particles(const std::size_t N,
                                                    const Real        bbox_XA,
                                                    const Real        bbox_XB,
                                                    const Real        bbox_YA,
                                                    const Real        bbox_YB,
                                                    const Real        bbox_ZA,
                                                    const Real        bbox_ZB)
{
  // problem dimension and domain size
  const std::size_t dim = _mesh.mesh_dimension();
  const Real max_range_x = bbox_XB - bbox_XA;
  const Real max_range_y = bbox_YB - bbox_YA;
  const Real max_range_z = bbox_ZB - bbox_ZA;
  const Real r = 1.0, den = 1.0;

  // generate random particle coordinates inside the domain, and write out the
  // file
  // We only use rank=0 proccessor to avoid generating different random numbers
  // on each processor.
  if (this->comm().rank() == 0)
  {
    printf(
      "---> test in generate_random_particles: Generating %lu random particles ...\n",
      N);

    // write the particle coordinates into a file
    std::string filename = "particle_random_file.txt";
    int o_width = 5, o_precision = 9;

    std::ofstream outfile;
    outfile.open(filename, std::ios_base::out);
    outfile << N << "\n";

    for (size_t i = 0; i < N; i++)
    {
      // generate random coordinates
      Real x = max_range_x * (std::rand() % 1000) / 1000 + bbox_XA;
      Real y = max_range_y * (std::rand() % 1000) / 1000 + bbox_YA;
      Real z = max_range_z * (std::rand() % 1000) / 1000 + bbox_ZA;

      if ((KDDim == 2) || (dim == 2)) z = 0.0;

      outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
      outfile.precision(o_precision);   outfile.width(o_width);
      outfile << x << "  " << y << "  " << z << "  " << r << "  " << den <<
        "  \n";
    } // end loop-i

    // close the file
    outfile.close();
    printf(
      "---> test in generate_random_particles: random particle file is created!\n");
  }

  this->comm().barrier();
}

// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::generate_random_particles(const std::size_t N,
                                                    const Point     & bbox_min,
                                                    const Point     & bbox_max)
{
  this->generate_random_particles(N,
                                  bbox_min(0), bbox_max(0),
                                  bbox_min(1), bbox_max(1),
                                  bbox_min(2), bbox_max(2));
}

// ======================================================================
template<unsigned int KDDim>
Real ParticleMesh<KDDim>::search_radius(const std::string& p_e) const
{
  if (p_e == "p") return _search_radius_p;
  else if (p_e == "e") return _search_radius_e;
  else
  {
    printf(
      "ParticleMesh::search_radius(): the input option must be either p or e!\n");
    libmesh_error();
  }
}

// ======================================================================
template<unsigned int KDDim>
std::size_t ParticleMesh<KDDim>::num_mesh_points() const
{
  START_LOG("num_mesh_points()", "ParticleMesh<KDDim>");

  std::size_t n_points = 0;

  for (std::size_t i = 0; i < _particles.size();
       ++i) n_points += _particles[i]->num_mesh_nodes();

  STOP_LOG("num_mesh_points()", "ParticleMesh<KDDim>");
  return n_points;
}

// ======================================================================
// =========================== print functions ==========================
// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::print_particle_info() const
{
  printf(
         "======================= printing the particle information: =======================\n\n");
  printf("There are totally %lu particles (search radius = %f)\n",
         _particles.size(),
         _search_radius_p);

  for (std::size_t j = 0; j < _particles.size(); ++j) _particles[j]->print_info();

  printf(
    "========================= end of the particle information =========================\n\n");
}

// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::update_mesh(const std::vector<Point>& nodal_pos,
                                      const std::vector<Point>& nodal_vel)
{
  START_LOG("update_mesh()", "ParticleMesh<KDDim>");

  // loop over all rigid particles
  std::size_t node_start_id = 0;

  for (int particle_id = 0; particle_id < _n_rigid_particles; particle_id++) {
    std::size_t num_mesh_nodes = _particles[particle_id]->num_mesh_nodes();
    std::vector<Point> nodal_points(num_mesh_nodes);
    std::vector<Point> nodal_velocities(num_mesh_nodes);

    // extract positions from nodal_pos to nodal_points vector
    for (int node_id = 0; node_id < num_mesh_nodes; node_id++) {
      nodal_points[node_id]     = nodal_pos[node_start_id + node_id];
      nodal_velocities[node_id] = nodal_vel[node_start_id + node_id];
    } // end loop over node id
    _particles[particle_id]->update_mesh(nodal_points, nodal_velocities);
    node_start_id += num_mesh_nodes;
  }   // end loop over particle id

  STOP_LOG("update_mesh()", "ParticleMesh<KDDim>");
}

// // ======================================================================
// template <unsigned int KDDim>
// void ParticleMesh<KDDim>::update_point_mesh(PointMesh<KDDim>* point_mesh)
// const
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
// void ParticleMesh<KDDim>::update_particle_mesh(const PointMesh<KDDim>*
// point_mesh)
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
template<unsigned int KDDim>
void ParticleMesh<KDDim>::zero_node_force()
{
  START_LOG("zero_particle_force_density", "ParticleMesh<KDDim>");

  for (int i = 0; i < _n_rigid_particles; i++) _particles[i]->zero_node_force();

  STOP_LOG("zero_particle_force_density", "ParticleMesh<KDDim>");
}

// ======================================================================
template<unsigned int KDDim>
std::vector<Real>ParticleMesh<KDDim>::mesh_size() const
{
  START_LOG("mesh_size()", "ParticleMesh<KDDim>");

  const std::size_t np = this->num_particles();
  std::vector<Real> hmin_max(2); // store hmin and hmax

  for (std::size_t i = 0; i < np; ++i)
  {
    if (i == 0) {
      hmin_max = _particles[i]->mesh_size();
    }
    else {
      const std::vector<Real> vals = _particles[i]->mesh_size();

      if (vals[0] < hmin_max[0]) hmin_max[0] = vals[0];

      if (vals[1] > hmin_max[1]) hmin_max[1] = vals[1];
    }
  }

  STOP_LOG("mesh_size()", "ParticleMesh<KDDim>");
  return hmin_max;
}

// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::volume_conservation()
{
  START_LOG("volume_conservation()", "ParticleMesh<KDDim>");

  // Loop over all the particles, and update their nodal values
  for (std::size_t i = 0; i < _particles.size(); ++i)
  {
    _particles[i]->volume_conservation();
  }

  STOP_LOG("volume_conservation()", "ParticleMesh<KDDim>");
}

// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::initial_particle_center_of_mass(
  std::vector<Point>& center0) const
{
  START_LOG("initial_center_of_mass()", "ParticleMesh<KDDim>");
  center0.resize(_n_rigid_particles);

  for (std::size_t i = 0; i < _n_rigid_particles; i++) {
    center0[i] = _particles[i]->get_centroid0();
  }
  STOP_LOG("initial_center_of_mass()", "ParticleMesh<KDDim>");
}

// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::write_particle(const unsigned int            & step_id,
                                         const unsigned int            & o_step,
                                         const Real                    & real_time,
                                         const std::vector<std::string>& output_file,
                                         unsigned int                    comm_in_rank)
const
{
  START_LOG("write_particle()", "ParticleMesh<KDDim>");
  this->write_time(step_id, o_step, real_time, comm_in_rank);

  for (int i = 0; i < output_file.size(); i++) {
    if (output_file[i] == "equation_systems") {
      // this output has been written in Copss.C
    }
    else if (output_file[i] == "trajectory") this->write_particle_trajectory(
        o_step,
        comm_in_rank);
    else if (output_file[i] == "particle_mesh") this->write_particle_mesh(o_step);

    // else if (output_file[i] == "restart_file/particle_mesh") this ->
    // write_particle_mesh_restart();
    else if (output_file[i] == "surface_node") this->write_surface_node(o_step,
                                                                        comm_in_rank);

    // else if (output_file[i] == "restart_file/surface_node") this->
    // write_surface_node_restart();
    else if (output_file[i] == "mean_square_displacement") {
      std::cout <<
        "Error: there is difficulty to calculate msd of rigid particles from ROUT, fix it before output msd"
                << std::endl;
      libmesh_error();

      // this -> write_particle_msd(step_id, o_step,center0, lvec,
      // comm_in_rank);
    }
    else {
      std::cout << "unsupported output_file content: (" << output_file[i] <<
        ")" << std::endl;
      libmesh_error();
    }
  }

  STOP_LOG("write_particle()", "ParticleMesh<KDDim>");
}

// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::write_time(const unsigned int& step_id,
                                     const unsigned int& o_step,
                                     const Real        & real_time,
                                     unsigned int        comm_in_rank) const
{
  // write step and real time
  std::ofstream out_file;

  if (comm_in_rank == 0) {
    if (step_id == 0) {
      out_file.open("out.time", std::ios_base::out);
      out_file << "step_id" << "  " << "o_step" << "  " << "real_time" << "\n";
    }
    else {
      out_file.open("out.time", std::ios_base::app);
    }
    out_file.precision(o_precision);
    out_file << step_id << "  " << o_step << "  " << real_time << "\n";
    out_file.close();
  }
}

// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::write_particle_trajectory(const unsigned int& o_step,
                                                    unsigned int        comm_in_rank)
const
{
  START_LOG("write_particle_trajectory()", "ParticleMesh<KDDim>");
  std::ostringstream oss;
  oss << "output_particle_trajectory.csv";
  std::ofstream out_file;
  out_file.precision(o_precision);
  out_file.setf(std::ios::fixed);

  if (comm_in_rank == 0) {
    if (o_step == 0) {
      out_file.open(oss.str(), std::ios_base::out);
      out_file <<
        "scalar x_coord y_coord z_coord x_vel y_vel z_vel x_force y_force z_force \n";
    }
    else {
      out_file.open(oss.str(), std::ios_base::app);
    }
    out_file << "#output_step_id," << o_step << " \n";

    for (std::size_t i = 0; i < _n_rigid_particles; ++i)
    {
      out_file << i << " ";

      // write position
      for (std::size_t j = 0; j < _dim; ++j) {
        out_file << _particles[i]->get_centroid()(j) << " ";
      }

      // write velocity
      for (std::size_t j = 0; j < _dim; ++j) {
        out_file << _particles[i]->centroid_velocity()(j) << " ";
      }

      // write force
      for (std::size_t j = 0; j < _dim; ++j) {
        out_file << _particles[i]->centroid_force()(j) << " ";
      }
      out_file << "\n";
    }
    out_file.close();
  } // end if comm_in_rank == 0
  STOP_LOG("write_particle_trajectory()", "ParticleMesh<KDDim>");
}

// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::write_surface_node(const unsigned int& o_step,
                                             unsigned int        comm_in_rank)
const
{
  std::ostringstream oss;

  oss << "output_surface_node.csv";
  std::ofstream out_file;
  out_file.precision(o_precision);
  out_file.setf(std::ios::fixed);

  if (comm_in_rank == 0) {
    if (o_step == 0) {
      out_file.open(oss.str(), std::ios_base::out);
      out_file <<
        "rigid_particle_id node_id x_coord y_coord z_coord node_center_distance \n";
    }
    else {
      out_file.open(oss.str(), std::ios_base::app);
    }
    out_file << "#output_step_id," << o_step << "\n";

    for (std::size_t i = 0; i < _n_rigid_particles; i++) {
      for (std::size_t j = 0; j < _particles[i]->num_mesh_nodes(); ++j) {
        const Point node_pos         = _particles[i]->mesh_point(j);
        const Real  node_center_dist = _particles[i]->node_center_dist(j);
        out_file << i << " " << j << " ";

        for (std::size_t k = 0; k < KDDim; k++) {
          out_file << node_pos(k) << " ";
        } // end loop over all dimensions
        out_file << node_center_dist << " \n";
      }   // end loop over all mesh points of one rigid particle
    }     // end loop over all rigid particles
    out_file.close();
  }       // end if comm_in_rank
}

// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::write_particle_mesh(const unsigned int& o_step) const
{
  START_LOG("write_particle_mesh()", "ParticleMesh<KDDim>");
  this->comm().barrier();

  for (std::size_t particle_id = 0; particle_id < _n_rigid_particles;
       particle_id++) {
    std::ostringstream surface_mesh_filename;
    surface_mesh_filename << "out_particle_" << particle_id + 1
                          << ".e-s."
                          << std::setw(7)
                          << std::setfill('0')
                          << std::right
                          << o_step;
    ExodusII_IO(_particles[particle_id]->mesh()).write(surface_mesh_filename.str());
  }
  STOP_LOG("write_particle_mesh()", "ParticleMesh<KDDim>");
}

// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::write_particle_mesh_restart() const
{
  START_LOG("write_particle_mesh_restart()", "ParticleMesh<KDDim>");
  this->comm().barrier();

  for (std::size_t particle_id = 0; particle_id < _n_rigid_particles;
       particle_id++) {
    std::ostringstream surface_mesh_filename;
    surface_mesh_filename << "restart_surface_mesh_particle_" << particle_id +
      1 << ".e";
    _particles[particle_id]->mesh().write(surface_mesh_filename.str());
  }
  STOP_LOG("write_particle_mesh_restart()", "ParticleMesh<KDDim>");
}

// ======================================================================
template<unsigned int KDDim>
SerialMesh& ParticleMesh<KDDim>::stitched_mesh()
{
  START_LOG("merge_particle_mesh()", "ParticleMesh<KDDim>");

  SerialMesh& _particles_mesh(_particles[0]->mesh());

  // assign subdomain_id of 1st particle to 0
  MeshBase::element_iterator el     = _particles_mesh.elements_begin();
  MeshBase::element_iterator end_el = _particles_mesh.elements_end();

  for (; el != end_el; ++el)
  {
    Elem *elem = *el;
    elem->subdomain_id() = 0;
  }

  const std::size_t n_particles = this->num_particles();

  if (n_particles > 1)
  {
    for (std::size_t i = 1; i < this->num_particles(); ++i)
    {
      // intermediate mesh to store each particle's mesh
      SerialMesh mesh(_particles[i]->mesh());

      // assign subdomain_id to particle_id
      el     = mesh.elements_begin();
      end_el = mesh.elements_end();

      for (; el != end_el; ++el)
      {
        Elem *elem = *el;
        elem->subdomain_id() = i;
      }

      _particles_mesh.stitch_meshes(mesh, 0, 0);
    }
  }

  STOP_LOG("merge_particle_mesh()", "ParticleMesh<KDDim>");

  // return the stitched mesh
  return _particles_mesh;
}

// ======================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::reinit()
{
  START_LOG("reinit()", "ParticleMesh<KDDim>");

  PMToolBox::output_message("===> Updating neighbor list for rigid particles "
                            "....", this->comm());
  // --> 1. build elem-elem neighbor list. Notice that this neighbor list
  // should and should only be built once since element locations don't
  // change during the simulation
  if (_elem_elem_neighbor_list.size()==0)
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

  STOP_LOG("reinit()", "ParticleMesh<KDDim>");
}

// =====================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::build_elem_elem_neighbor_list()
{
  // let's be sure we properly initialize on every processor at once: this
  // function is only executed in debug mode
  parallel_object_only();
  START_LOG("build_elem_elem_neighbor_list()", "ParticleMesh<KDDim>");

  std::ostringstream oss;
  oss << "\n>>>>>>>>>>> Building elem-elem neighbor list";
  PMToolBox::output_message(oss, this->comm());
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
  oss << "Done!!!";
  PMToolBox::output_message(oss, this->comm());
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

  STOP_LOG("build_elem_elem_neighbor_list()", "ParticleMesh<KDDim>");
}

// =====================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::build_elem_point_containing_list()
{
  START_LOG("build_elem_point_containing_list()", "ParticleMesh<KDDim>");

  std::ostringstream oss;
  oss << "\n>>>>>>>>>>> Building elem-RigidParticle containing list";
  PMToolBox::output_message(oss, this->comm());

  // clear the mapping
  _elem_point_containing_list.clear();

  for (std::size_t p_id=0; p_id<_particles.size(); p_id++)
  {
    bool use_locator = false;

    // if current element has not been assigned, we use locator
    if (_particles[p_id]->elem_id()==-1)
      use_locator = true;
      // if the current element doesn't contain this point anymore, we use locator
    else if (!_mesh.elem(_particles[p_id]->elem_id())->contains_point
      (_particles[p_id]->get_centroid()))
      use_locator = true;

    // use locator. If found, reset point elem id
    if (use_locator)
    {
      const Elem* elem = _mesh.point_locator().operator()
        (_particles[p_id]->get_centroid());
      if (elem== nullptr)
      {
        std::ostringstream oss;
        oss << "Error: cannot find the element that current particle is "
               "within : "
            <<"particle id = " << p_id <<", "
            <<"particle location = (" << _particles[p_id]->get_centroid()(0) <<
            ", "
            << _particles[p_id]->get_centroid()(1) << ", "
            << _particles[p_id]->get_centroid()(2) << "). Exiting ...";
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

  oss << "Done!!!";
  PMToolBox::output_message(oss, this->comm());
//  for (dof_id_type elem_id=0; elem_id<_mesh.n_elem(); elem_id++) {
////    if (_elem_point_containing_list.count(elem_id)){
//    oss << "Elem id " << elem_id << " contains RigidParticle: ";
//    const std::vector <dof_id_type> &contained_points =
//      this->get_elem_point_containing_list(elem_id);
//    for (dof_id_type i = 0; i < contained_points.size(); i++)
//      oss << contained_points[i] << ", ";
//    oss << "\n";
////    }
//  }
//  PMToolBox::output_message(oss, this->comm());

  STOP_LOG("build_elem_point_containing_list()", "ParticleMesh<KDDim>");
}

// =====================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::build_elem_point_neighbor_list()
{
  // let's be sure we properly initialize on every processor at once: this
  // function is only executed in debug mode
  parallel_object_only();
  START_LOG("build_elem_point_neighbor_list", "ParticleMesh<KDDim>");

  std::ostringstream oss;
  oss << "\n>>>>>>>>>>> Building elem-RigidParticle neighbor list (using "
         "information in elem-elem neighbor list and elem-RigidParticle "
         "containing list)";
  PMToolBox::output_message(oss, this->comm());
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
  oss << "Done!!!";
  PMToolBox::output_message(oss, this->comm());
//  for (dof_id_type elem_id=0; elem_id<n_elem; elem_id++)
//  {
//    oss << "elem id " << elem_id << ", neighbor RigidParticle list: ";
//    const std::vector<dof_id_type>& point_nb_list =
//      this->get_elem_point_neighbor_list(elem_id);
//    for (dof_id_type i=0; i<point_nb_list.size(); i++)
//      oss << point_nb_list[i] << ", ";
//    oss << "\n";
//  }
//  PMToolBox::output_message(oss, this->comm());


  STOP_LOG("build_elem_point_neighbor_list", "ParticleMesh<KDDim>");
}

// =====================================================================
template<unsigned int KDDim>
void ParticleMesh<KDDim>::build_point_point_neighbor_list()
{
  START_LOG("build_point_point_neighbor_list()", "ParticleMesh<KDDim>");

  std::ostringstream oss;
  oss << "\n>>>>>>>>>>> Building RigidParticle-RigidParticle neighbor list "
         "(first get the element id of this RigidParticle and then find all "
         "neighbor RigidParticle of the element)";
  PMToolBox::output_message(oss, this->comm());

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

  oss << "Done!!!";
  PMToolBox::output_message(oss, this->comm());
//  for (dof_id_type p_id=0; p_id<_particles.size(); p_id++) {
//    oss << "RigidParticle id " << p_id <<" (elem_id = " <<_particles[p_id]->elem_id()
//        <<"), neighbor RigidParticle: ";
//    const std::vector <dof_id_type> &neighbor_points =
//      _particles[p_id]->neighbor_list();
//    for (dof_id_type i = 0; i<neighbor_points.size(); i++)
//      oss << neighbor_points[i] << ", ";
//    oss << "\n";
//  }
//  PMToolBox::output_message(oss, this->comm());

  STOP_LOG("build_point_point_neighbor_list()", "ParticleMesh<KDDim>");
}

// ======================================================================
//template<unsigned int KDDim>
//void ParticleMesh<KDDim>::print_elem_neighbor_list(std::ostream& out) const
//{
//  printf(
//    "======================= printing the element neighbor list: ========================\n\n");
//  std::map<const dof_id_type, std::vector<dof_id_type> >::const_iterator p;
//
//  for (p = _elem_point_neighbor_list.begin(); p != _elem_point_neighbor_list
//  .end(); ++p)
//  {
//    const dof_id_type& elem_id             = p->first;
//    const std::vector<dof_id_type>& n_list = p->second;
//    const Elem *elem                      = _mesh.elem(elem_id);
//    const Point center_pt                 = elem->centroid();
//    const Real  hmax                      = elem->hmax();
//
//    // ----- Scheme 2: print out information using printf (from all the
//    // processors) -----
//    if (n_list.size() > 0)
//    {
//      printf(
//             "========== There are %lu neighbor particles around the element %lu, output rank = %u :===\n",
//        n_list.size(),
//        elem_id,
//        this->comm().rank());
//      printf("element centroid = (%f, %f, %f), and hmax = %f \n",
//             center_pt(0), center_pt(1), center_pt(2), hmax);
//
//      if (!n_list.empty())
//        for (std::size_t i = 0; i < n_list.size(); ++i) // printf("%lu    ",
//                                                        // n_list[i]);
//          _particles[n_list[i]]->print_info();
//      else out << "There is no neighboring particle around this element! \n";
//      printf("\n");
//    } // end if
//  }   // end for p-loop
//  printf(
//    "======================= end of the element neighbor list ======================\n\n");
//}

// ------------------------------------------------------------
// Explicit Instantiations
template class ParticleMesh<1>;
template class ParticleMesh<2>;
template class ParticleMesh<3>;
} // end of namespace
