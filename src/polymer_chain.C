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




#include <fstream>

#include "libmesh/libmesh_common.h"

#include "pm_periodic_boundary.h"
#include "polymer_chain.h"



namespace libMesh
{


// ======================================================================
PolymerChain::PolymerChain(const std::size_t chain_id)
  : _chain_id(chain_id),
    _periodic_boundary(NULL)
{
  // Do nothing
}


  
// ======================================================================
PolymerChain::PolymerChain(const std::size_t chain_id,
                           PMPeriodicBoundary& pm_pb)
  : _chain_id(chain_id),
  _periodic_boundary(&pm_pb)
{
  // Do nothing
}
  
// ======================================================================
PolymerChain::~PolymerChain()
{
  // Free the memory
  for (std::size_t i=0; i<this->n_beads(); ++i)
  {
    if(_beads[i]){
      delete _beads[i];
    }
  }
  
}


// ======================================================================
void PolymerChain::read_data(const std::string& filename)
{
  START_LOG ("read_data()", "PolymerChain");
  
  // Open the local file and check the existance
  std::cout <<"\n###Polymer chain filename = "<<filename <<std::endl;
  std::ifstream infile;
  infile.open (filename, std::ios_base::in);
  if( !infile.good() )
  {
    printf("***warning: read_data() can NOT read the polymer chain data!");
    libmesh_error();
  }
  
  // init variables:
  // point_type:  0 - polymer bead point; 1 - tracking point; or user-defined type
  const PointType point_type = POLYMER_BEAD;
  Real x=0., y=0., z=0.;            // initialize bead coords
  std::size_t n_beads, b_id, c_id, bead_type;   //
  std::vector<Real> rot_vec(4); // rotation vector (a,b,c) + theta.
  
  // read particle data
  infile >> n_beads;            // total number of beads
  _beads.resize(n_beads);
  for (std::size_t i=0; i<n_beads; ++i)
  {
    // We read, but don't use the bead_type when construct PointParticle.
    // because the file numbers the bead and its type from 1, while C++ prefers
    // to start the number from 0.
    infile >> b_id >> c_id >> bead_type >> x >> y >> z
    >> rot_vec[0] >> rot_vec[1] >>rot_vec[2] ;
    Point pt(x,y,z);
    PointParticle* particle = new PointParticle(pt, i, point_type, rot_vec);
    particle->set_parent_id(c_id);  // parent id = chain id
    
    // add to the beads list
    _beads[i] = particle;
    
  } // end for i-loop
  
  // change chain id
  _chain_id = c_id;
  
  // Finish and close the file
  infile.close();
  std::cout << "Reading polymer chain data from "<<filename<<" is completed!\n\n";
  
  STOP_LOG ("read_data()", "PolymerChain");
}


  
// ======================================================================
void PolymerChain::read_data_pizza(const std::string& filename)
{
  START_LOG ("read_data_pizza()", "PolymerChain");
  
  // Open the local file and check the existance
//  std::cout <<"-----------------> Read polymer chain in pizza format from "<<filename <<std::endl;
  std::ifstream infile;
  infile.open (filename, std::ios_base::in);
  if( !infile.good() )
  {
    printf("***warning: read_data_pizza() can NOT read the polymer chain data!\n");
    libmesh_error();
  }
  
  // init variables:
  // point_type:  0 - polymer bead point; 1 - tracking point; or user-defined type
  const PointType point_type = POLYMER_BEAD;
  Real x=0., y=0., z=0.;            // initialize bead coords
  int bead_id, chain_id, bead_type; //
  std::vector<Real> rot_vec(4); // rotation vector (a,b,c) + theta.
  std::size_t n1, n2, n3;
  
  
  // Read file line by line
  std::string line_str, str_tmpt;
  std::getline(infile, line_str); // 0. Header line
  infile >> _n_beads >> line_str; // 1. # of beads
  infile >> _n_bonds >> line_str; // 2. # of bonds
  infile >> _n_bead_types >> line_str >> str_tmpt; // 3. # of bead types
  infile >> _n_bond_types >> line_str >> str_tmpt; // 4. # of bond types
  // init
  _mass.resize(_n_bead_types);
  _beads.resize(_n_beads);
  _bonds.resize(_n_bonds);
  for(std::size_t i=0; i<_n_bonds; ++i) _bonds[i].resize(3);
  
  // Read the bead mass
  while (std::getline(infile, line_str)){
    if(line_str=="Masses") break;
  }
  std::getline(infile, line_str); // skip the empty line
  for(std::size_t i=0; i<_n_bead_types; ++i)
  {
    infile >> n1 >> _mass[i];
  }
  
  
  // Read the beads
  while (std::getline(infile, line_str)){
    if(line_str=="Atoms") break;
  }
  std::getline(infile, line_str); // skip the empty line
  for(std::size_t i=0; i<_n_beads; ++i)
  {
    // read data in this line
    infile >> bead_id >> chain_id >> bead_type >> x >> y >> z
           >> rot_vec[0] >> rot_vec[1] >>rot_vec[2];
    
    // create PointParticle
    Point pt(x,y,z);
    PointParticle* particle = new PointParticle(pt, bead_id-1, point_type, rot_vec);
    particle->set_parent_id(chain_id-1);  // parent id = chain id
    
    // add to the beads list
    _beads[i] = particle;
  }
  
  // Use id of the last chain to set number of chains, and initialize vector
  _n_chains = chain_id;
  _n_beads_per_chain.resize(_n_chains, 0);
 
  // Read the bonds
  while (std::getline(infile, line_str)){
    if(line_str=="Bonds") break;
  }
  std::getline(infile, line_str); // skip the empty line
  for(std::size_t i=0; i<_n_bonds; ++i)
  {
    // read data in this line
    infile >> bead_id >> n1 >> n2 >> n3;
    
    // create PointParticle
    _bonds[i][0] = n1-1;  // bond type
    _bonds[i][1] = n2-1;  // id of connected bead 1
    _bonds[i][2] = n3-1;  // id of connected bead 2
  }

  
  // Finish and close the file
  infile.close();
  
  // set number of beads on each chain 
  for(std::size_t i=0; i<_n_beads; ++i)
  {
    // chain id
    chain_id = _beads[i]->parent_id();
    _n_beads_per_chain[chain_id] += 1; 
  }

  STOP_LOG ("read_data_pizza()", "PolymerChain");
}

  

// ======================================================================
void PolymerChain::read_data_vtk(const std::string& filename)
{
  START_LOG ("read_data_vtk()", "PolymerChain");
  
  // Open the local file and check the existance
  std::cout <<"\n###Read polymer chain with vtk format"<<std::endl;
  std::cout <<"   filename = "<<filename <<std::endl;
  std::ifstream infile;
  infile.open (filename, std::ios_base::in);
  if( !infile.good() )
  {
    printf("***warning: read_data_vtk() can NOT read the polymer chain data!");
    libmesh_error();
  }
  
  // init variables:
  // point_type:  0 - polymer bead point; 1 - tracking point; or user-defined type
  const PointType point_type = POLYMER_BEAD;
  Real x=0., y=0., z=0.;            // initialize bead coords
  int bead_id, chain_id, bead_type; //
  std::vector<Real> rot_vec(4); // rotation vector (a,b,c) + theta.
  std::size_t n1, n2, n3;
  
  
  // Read file line by line
  std::string line_str, str_tmpt;
  std::getline(infile, line_str); // 0. Header line
  std::getline(infile, line_str); // 1. Header line
  std::getline(infile, line_str); // 2. Header line
  std::getline(infile, line_str); // 3. Header line
  infile >> line_str >> _n_beads >> str_tmpt; // 4. # of beads

  // Read beads data
  _beads.resize(_n_beads);
  for(std::size_t i=0; i<_n_beads; ++i)
  {
    // read data in this line: coordinates
    infile >> x >> y >> z;
    
    // create PointParticle
    Point pt(x,y,z);
    bead_type = 0;    // bead type = 0 for polymer chain.
    PointParticle* particle = new PointParticle(pt, i, point_type, rot_vec);
    
    // add to the beads list
    if( _beads[i] ) delete _beads[i];
    _beads[i] = particle;
  }
  
  // Read bonds data
  std::getline(infile, line_str);       // skip the empty line
  infile >> line_str >> _n_bonds >> n1; // # of bonds
  _bonds.resize(_n_bonds);
  for(std::size_t i=0; i<_n_bonds; ++i) _bonds[i].resize(3);
  for(std::size_t i=0; i<_n_bonds; ++i)
  {
    // read data in this line
    infile >> n1 >> n2 >> n3;
    
    // create PointParticle
    _bonds[i][0] = 0;   // bond type is set 0
    _bonds[i][1] = n2;  // connect bead 1
    _bonds[i][2] = n3;  // connect bond 2
  }
  
  // Read parent id
  std::getline(infile, line_str);   // skip(this line can't be removed)
  std::getline(infile, str_tmpt);   // skip
  std::getline(infile, str_tmpt);   // skip
  std::getline(infile, str_tmpt);   // skip
  for(std::size_t i=0; i<_n_beads; ++i)
  {
    infile >> chain_id;
    _beads[i]->set_parent_id(chain_id);
    //printf("---> chain id = %d\n",chain_id);
  }
  // Use id of the last chain to set number of chains, and initialize vector
  _n_chains = chain_id + 1; // chain ID start from 0 in .vtk
  _n_beads_per_chain.resize(_n_chains, 0);
  
  // Finish and close the file
  infile.close();

  // set number of beads on each chain
  for (std::size_t i=0; i<_n_beads; ++i){
    _n_beads_per_chain[_beads[i]->parent_id()] += 1;
  }
  std::cout << "Reading polymer chain data from "<<filename<<" is completed!\n\n";
  STOP_LOG ("read_data_vtk()", "PolymerChain");
}

//=======================================================================
void PolymerChain::read_data_csv(const std::string& filename)
{
  START_LOG ("read_data_csv()", "PolymerChain");
  
  // Open the local file and check the existance
  std::cout <<"\n###Read bead with csv format"<<std::endl;
  std::cout <<"   filename = "<<filename <<std::endl;
  std::ifstream infile;
  infile.open (filename, std::ios_base::in);
  if( !infile.good() )
  {
    printf("***warning: read_data_csv() can NOT read the bead data!");
    libmesh_error();
  }
  
  // init variables:
  // point_type:  0 - polymer bead point; 1 - tracking point; or user-defined type
  const PointType point_type = POLYMER_BEAD;
  Real x=0., y=0., z=0.;            // initialize bead coords
  Real vx=0., vy=0., vz=0.;         // initialize bead velocities
  Real fx=0., fy=0., fz=0.;         // initialize bead forces
  int bead_id, chain_id, bead_type; //
  std::vector<Real> rot_vec(4); // rotation vector (a,b,c) + theta.
  _n_chains = 1;
  // Read file line by line
  std::string line_str, str_tmpt;
  std::getline(infile, line_str); // 0. Header line
  // Read beads data
  std::cout << line_str << std::endl;
  int i = 0;
  int bead_old = 0;
  while (!infile.eof()) {
    infile >> bead_id >> x >> y >> z >> vx >> vy >> vz >> fx >> fy >> fz;
    // create PointParticle
    Point pt(x,y,z);
    bead_type = 0;    // bead type = 0 for polymer chain.
    if (bead_id == 0 or !(bead_id == bead_old)){
      PointParticle* particle = new PointParticle(pt, i, point_type, rot_vec);
      // add to the beads list
      particle->set_parent_id(0);
      _beads.push_back(particle);
    }
    bead_old = bead_id; // avoid repeating last line
    i+=1;
  }
  _n_beads = bead_id + 1; // number of beads = bead id of last bead read from csv
  _beads.resize(_n_beads);
  _n_beads_per_chain.resize(_n_chains,0);
  // Finish and close the file
  infile.close();
  std::cout <<"_n_beads = "<<_n_beads <<std::endl;
  // set number of beads on each chain
  for (std::size_t i=0; i<_n_beads; ++i){
   _n_beads_per_chain[_beads[i]->parent_id()] += 1;
  }
  std::cout << "Reading bead data from "<<filename<<" is completed!\n\n";
  STOP_LOG ("read_data_csv()", "PolymerChain");

}
  
  
// ======================================================================
void PolymerChain::read_pbc_count()
{
  START_LOG ("read_pbc_count()", "PolymerChain");
  
  // Open the local file and check the existance
  std::ifstream infile;
  infile.open ("pbc_count.dat", std::ios_base::in);
  if( !infile.good() )
  {
    printf("***warning: read_pbc_count() can NOT read the PBC count data!");
    libmesh_error();
  }

  // init variables:
  int nx=0, ny=0, nz=0;
 
  // Read pbc counter data
  for(std::size_t i=0; i<_n_beads; ++i)
  {
    // read data in this line: counters in x, y, z direction
    infile >> nx >> ny >> nz;
 
    // pass counter to PointParticle
    _beads[i]->set_pbc_counter(nx, ny, nz);
  }
 
  // Finish and close the file
  infile.close();
 
  STOP_LOG ("read_pbc_count()", "PolymerChain");
}



// ======================================================================
void PolymerChain::add_bead(const Point& pt,
                            const PointType point_type,
                            const int parent_id,
                            const std::vector<Real>& rot_vec)
{
  START_LOG ("add_bead()", "PolymerChain");
  
  // bead id = the last id + 1
  const std::size_t bead_id = _beads.size();
  PointParticle* particle = new PointParticle(pt,bead_id,point_type,rot_vec);
  particle->set_parent_id(parent_id);
  
  // add to the beads list
  _beads.push_back(particle);

  // build the bond
  std::vector<std::size_t> new_bond(3);
  new_bond[0] = 0;
  new_bond[1] = bead_id - 1;
  new_bond[2] = bead_id;
  _bonds.push_back(new_bond);
  
  STOP_LOG ("add_bead()", "PolymerChain");
}
  
  
// ======================================================================
void PolymerChain::generate_polymer_chains(const unsigned int n_beads,
                                          const unsigned int n_chains,
                                          const Real init_bondLength ,
                                          const Point init_bbox_min,
                                          const Point init_bbox_max,
                                          const std::string& filename,
					  unsigned int comm_in_rank)
{
  // problem dimension and domain size
  const std::size_t dim = 3;
  const Real r = 1.0, den = 1.0;
  const std::size_t n_beads_perChain = n_beads / n_chains;
  const std::size_t n_bonds_perChain = n_beads_perChain - 1;  // number of springs
  const std::size_t n_bonds = n_chains * n_bonds_perChain;
  // generate random particle coordinates inside the domain, and write out the file
  if( comm_in_rank==0 )
   { 
   
    std::cout << "\n-------------------------------------------------\n";
    std::cout << "generating polymer chains... ---> "<<n_chains <<" chains, "<<n_beads_perChain<<" beads per chain, maximum bond length = "<<init_bondLength <<"\n";
    std::cout << "within box ----> box min = ( "<<init_bbox_min(0) << ", " <<init_bbox_min(1) <<", "<<init_bbox_min(2)<<" );\n";
    std::cout << "                 box max = ( "<<init_bbox_max(0) << ", " <<init_bbox_max(1) <<", "<<init_bbox_max(2)<<" );\n";
    // write the particle coordinates into a file
    //std::string filename = "polymer_data.in";
    std::ofstream outfile;
    outfile.open(filename, std::ios_base::out);
    // header line
    outfile << "COPSS FENE chain data file\n";
    // # of beads
    outfile << n_beads << " atoms\n";
    // # of bonds
    outfile << n_bonds <<" bonds\n";
    // atoms types (default = 1)
    outfile << "1 atom types\n";
    // bond types (default = 1)
    outfile << "1 bond types\n";
    // box info
    outfile << init_bbox_min(0) << " " << init_bbox_max(0) << " xlo xhi\n";
    outfile << init_bbox_min(1) << " " << init_bbox_max(1) << " ylo yhi\n";
    outfile << init_bbox_min(2) << " " << init_bbox_max(2) << " zlo zhi\n";
    // masses (default only one kind)
    outfile <<"\nMasses\n";
    outfile <<"\n1 1.0\n";
    // beads info
    outfile<<"\nAtoms\n\n";
    
    Point pt1, nv, bv;    
    Point pt0 = init_bbox_max + init_bbox_min;
    Real dr, theta, phi;
    //generate random coordinates for each bead
    for(std::size_t i=0; i<n_beads; ++i)
    {
      
      bool in_domain = false;
      // We will generate a bead, whose coordinate is in the domain
      std::size_t count = 0;
      while(in_domain==false)
      {
        // generate r, theta, phi in spherical coordinates and then convert it to cartesian corrdinates
	// bond length is between 2 and init_bondLength
	dr = (std::rand() % 1000 / 1000.) * (1.- 2./ init_bondLength) + 2. / init_bondLength;
	theta = (std::rand() % 1000) / 1000.;
	phi = 2.* (std::rand() % 1000) / 1000.;
	nv(0) = dr * std::sin(theta * M_PI) * std::cos(phi * M_PI);
	nv(1) = dr * std::sin(theta * M_PI) * std::sin(phi * M_PI);
	nv(2) = dr * std::cos(theta * M_PI);
        // bond (spring) vector: lenght and direction
        for (std::size_t j=0; j<3; ++j)   bv(j) = init_bondLength * nv(j);
        // position of the new bead.
        if (i==0)
          pt1 = pt0;
        else
          pt1 += bv;
        // end if-else
        
        // check if the new bead is in the assigned domain.
        for(std::size_t j=0; j<3; ++j)
        {
          if( (pt1(j)<=init_bbox_min(j))  || (pt1(j)>=init_bbox_max(j)) )
          {
            in_domain = false;
            pt1 -= bv;
            break;    // if the new bead is outside the domain, break j-loop
          }
          else
            in_domain = true;
          // end if-else
        } // end for j-loop
        
       // count++;
       // printf("---> test in generate_polymer_chain: bead %lu, count = %lu, vn = (%f,%f,%f)\n",
         //      i, count, nv(0), nv(1), nv(2) );
        
      } // end while loop
      printf("---> test in generate_polymer_chain: new bead xyz = (%f,%f,%f)\n",pt1(0),pt1(1),pt1(2));
      
      // write out the information
      // bead_id chain_id atom_type x y z rot_x rot_y rot z
      // default: rot_x = 0, rot_y = 0, rot_z = 0
      // default: atom_type = 1
      outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
      outfile.precision(o_precision); 
      // compute chain_id_i, starting from 0
      std::size_t chain_id_i = i / n_beads_perChain;
      // bead_id starting from 1; chain_id starts from 1
      outfile << i+1 <<"  "<<chain_id_i + 1<<"  1"<<"  "<< pt1(0) << "  " << pt1(1) << "  " << pt1(2) << "  0  0  0\n";
    } // end loop-i

    //output bonds
    outfile<<"\nBonds\n\n";
    
    for(std::size_t i=0; i<n_bonds; ++i){
    // bond_id bond_type bondleft bond_right
    // bond_type = 1(default)

    // chain_id_i starts from 0
    std::size_t chain_id_i = i / n_bonds_perChain;
    outfile <<i+1 <<"  1  "<< i + 1 + chain_id_i<<"  "<< i + 2 + chain_id_i<<"\n";
    }
    outfile.close();
    std::cout << "\n\nfinished!!generating polymer chains...\n\n";
  } 
  return;
}

 
// ======================================================================
void PolymerChain::generate_polymer_chain(const Point pt0,
                                          const std::size_t Nb,
                                          const Real Ls,
                                          const Point& bbox_min,
                                          const Point& bbox_max,
                                          const std::string& filename)
{
  // problem dimension and domain size
  const std::size_t dim = 3;
  const Real r = 1.0, den = 1.0;
  const std::size_t Ns = Nb - 1;  // number of springs
  
  // generate random particle coordinates inside the domain, and write out the file
//  if( this->comm().rank()==0 )
  {
    printf("---> test in generate_polymer_chain: Generating %lu beads and %lu springs...\n",Nb,Ns);
    
    // write the particle coordinates into a file
    //std::string filename = "polymer_data.in";
    std::ofstream outfile;
    outfile.open(filename, std::ios_base::out);
    outfile << Nb << "\n";
    Point pt1, nv, bv;
    for(std::size_t i=0; i<Nb; ++i)
    {
      
      bool in_domain = false;
      
      // We will generate a bead, whose coordinate is in the domain
      std::size_t count = 0;
      while(in_domain==false)
      {
        // generate random vector vn,  -1 < vn(j) < +1
        for (std::size_t j=0; j<3; ++j)
          nv(j) = 2.0*(std::rand() % 1000 ) / 1000. - 1.0;
        
        if (dim==2) nv(2) = 0.0;
        
        // normalize the vector, and let |vn| = 1
        const Real sr = nv.size();
        for (std::size_t j=0; j<3; ++j)   nv(j) /= sr;
        
        // bond (spring) vector: lenght and direction
        for (std::size_t j=0; j<3; ++j)   bv(j) = Ls*nv(j);
        
        // position of the new bead.
        if (i==0)
          pt1 = pt0;
        else
          pt1 += bv;
        // end if-else
        
        // check if the new bead is in the assigned domain.
        for(std::size_t j=0; j<3; ++j)
        {
          if( (pt1(j)<=bbox_min(j))  || (pt1(j)>=bbox_max(j)) )
          {
            in_domain = false;
            pt1 -= bv;
            break;    // if the new bead is outside the domain, break j-loop
          }
          else
            in_domain = true;
          // end if-else
        } // end for j-loop
        
        count++;
        printf("---> test in generate_polymer_chain: bead %lu, count = %lu, vn = (%f,%f,%f)\n",
               i, count, nv(0), nv(1), nv(2) );
        
      } // end while loop
      printf("---> test in generate_polymer_chain: new bead xyz = (%f,%f,%f)\n",pt1(0),pt1(1),pt1(2));
      
      // write out the coordinates
      outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
      outfile.precision(o_precision);
      outfile << pt1(0) << "  " << pt1(1) << "  " << pt1(2) << "  " << r << "  "<< den << "  \n";
    } // end loop-i
    
    outfile.close();
    printf("---> test in generate_polymer_chain: polymer chain file is created!\n");
  }
  
//  this->comm().barrier();
  return;
}



// ======================================================================
void PolymerChain::print_info() const
{
  START_LOG ("print_info()", "PolymerChain");
  
  
  printf("============================================================================\n");
  printf("*** printing the information of polymer chain: \n");
  printf("============================================================================\n");
  
  printf("There are totally %lu beads and %lu bonds in the polymer chain\n",_n_beads, _n_bonds);
  printf("   %lu bead types\n",_n_bead_types);
  printf("   %lu bond types\n",_n_bond_types);
  
  
  // Masses
  printf("\nMasses:\n");
  for(std::size_t i=0; i<_n_bead_types; ++i) {
    printf("   %lu  %f\n",i, _mass[i]);
  }
  
  // beads
  printf("\nBeads:\n");
  for(std::size_t i=0; i<_n_beads; ++i) {
    printf("   %d %d %d  %f %f %f  %f %f %f\n",
           _beads[i]->id(),_beads[i]->parent_id(),_beads[i]->point_type(),
           _beads[i]->center()(0), _beads[i]->center()(1), _beads[i]->center()(2),
           _beads[i]->orientation()[0], _beads[i]->orientation()[1], _beads[i]->orientation()[2]);
  }
  
  // Bonds
  printf("\nBonds:\n");
  for(std::size_t i=0; i<_n_bonds; ++i) {
    printf("   %lu  %lu %lu %lu\n", i, _bonds[i][0], _bonds[i][1], _bonds[i][2]);
  }
  printf("\n");
  
  
  // Bead info (Note there is no neighbor list info before initialized)
  printf("There are totally %lu beads in the polymer chain:\n",_n_beads);
  for (std::size_t j=0; j<_n_beads; ++j) {
    _beads[j]->print_info();
  }
  printf("======================= end of the polymer information =======================\n\n");
  
  STOP_LOG ("print_info()", "PolymerChain");
}

  
  
// ======================================================================
void PolymerChain::write_polymer_chain(const unsigned int& step_id,
                                       const unsigned int& o_step,
                                       const Real& real_time,
                                       const std::vector<Point>& center0,
                                       const std::vector<Real>& lvec,
                                       std::vector<std::string>& output_file,
                                       unsigned int comm_in_rank) const
{
  START_LOG ("write_polymer_chain()", "PolymerChain");
  this -> write_time(step_id, o_step, real_time, comm_in_rank);
  for (int i = 0; i < output_file.size(); i++){
    if(output_file[i] == "equation_systems") {
      // this output has been written in COPSS.C
    }
    else if(output_file[i] == "trajectory") this -> write_polymer_trajectory(o_step, comm_in_rank);
    else if(output_file[i] == "center_of_mass") this -> write_polymer_com(step_id, o_step,lvec, comm_in_rank);
    else if(output_file[i] == "chain_stretch") this -> write_polymer_stretch(step_id,o_step, lvec, comm_in_rank);
    else if(output_file[i] == "radius_of_gyration") this -> write_polymer_rog(step_id, o_step, lvec, comm_in_rank);
    else if(output_file[i] == "mean_square_displacement") this -> write_polymer_msd(step_id, o_step, center0, lvec, comm_in_rank);
    else {
      std::cout << "unsupported output_file content: (" << output_file[i] <<")" << std::endl; 
      libmesh_error();
    }
  }
  STOP_LOG ("write_polymer_chain()", "PolymerChain");
}

// ========================================================================
void PolymerChain::write_polymer_trajectory(const unsigned int& o_step,
                                            unsigned int comm_in_rank) const
{
  // OFSTEAM
  std::ostringstream oss;
  oss << "output_polymer_" << o_step <<".vtk";
  std::ofstream outfile;
  if(comm_in_rank == 0){
    outfile.open(oss.str(), std::ios_base::out);
//    outfile.precision(o_precision);
    const std::size_t n_beads = _beads.size();  
    // write out the VTK file
    outfile << "# vtk DataFile Version 4.0\n";
    outfile << "Polymer chain\n";
    outfile << "ASCII\n";
    outfile << "DATASET POLYDATA\n";
    
    // POINT data
    outfile << "POINTS " << n_beads << " float\n";
    for(std::size_t i=0; i<n_beads; ++i)
    {
      for(std::size_t j=0; j<3; ++j){
        outfile << _beads[i]->center()(j) << " ";
      }
      outfile << "\n";
    }
    outfile << "\n";
    
    // LINE(spring) data, which can be visualized by connectivity in Paraview
    const std::size_t Ns = _bonds.size();
    outfile << "LINES " << Ns << " " << Ns*3 << "\n";
    for(std::size_t i=0; i<Ns; ++i)
    {
      outfile << 2 << " " << _bonds[i][1] << " " << _bonds[i][2] << "\n";
    }
    
    // Define the bead types
    outfile << "POINT_DATA " << n_beads << "\n";
    outfile << "SCALARS " << "bead_type "<< "int" << "\n";
    outfile << "LOOKUP_TABLE default \n";
    for(std::size_t i=0; i<n_beads; ++i)
    {
      const int parent_id = _beads[i]->parent_id();
      //const int b_type = (parent_id>0)? parent_id : 1;
      const int b_type = parent_id;
      outfile << b_type << " \n";
    }
    outfile.close();    
  } // end if comm_in_rank
}

// =============================================================================
void PolymerChain::write_polymer_com(const unsigned int& step_id,
                                     const unsigned int& o_step,
                                     const std::vector<Real>& lvec,
                                     unsigned int comm_in_rank) const
{
  // compute center of mass at step i
  std::vector<Point> centeri(_n_chains);
  compute_center_of_mass(lvec, centeri);
  // write center of mass of every chain to out.center_of_mass
  std::ofstream out_file;  
  if(comm_in_rank == 0){
    if(step_id == 0){
      out_file.open("out.center_of_mass", std::ios_base::out);
      out_file <<"o_step" <<"  " <<"center_of_mass\n"; 
    }
    else{
      out_file.open("out.center_of_mass", std::ios_base::app);
    }
    out_file.precision(o_precision);
    out_file << o_step <<"  "; 
    for (std::size_t j = 0; j < _n_chains; j++){
            for (std::size_t k = 0; k < _dim; k++){
                    out_file << centeri[j](k) << " ";
            } // end k loop
    } // end j loop
    out_file <<"\n";
    out_file.close();
  }// end if comm_in_rank == 0
}

//==============================================================================
void PolymerChain::write_polymer_stretch(const unsigned int& step_id,
                                         const unsigned int& o_step,
                                         const std::vector<Real>& lvec,
                                         unsigned int comm_in_rank) const
{
  // compute chain strech at step i
  std::vector<Point> stretch(_n_chains);
  this -> compute_chain_stretch(lvec,stretch);
  // write strech of each chain to out.chain_strech
  std::ofstream out_file;  
  if(comm_in_rank == 0){
    if(step_id == 0){
      out_file.open("out.chain_strech", std::ios_base::out);
      out_file <<"o_step" <<"  " <<"chain_strech\n";
    }
    else{
      out_file.open("out.chain_strech", std::ios_base::app);
    }
    out_file.precision(o_precision);    
    out_file << o_step <<"  ";     
    for (std::size_t j = 0; j < _n_chains; j++){
      for (std::size_t k = 0; k < _dim; k++){
              out_file << stretch[j](k) << " ";
      } // end k loop
    }// end j loop
    out_file <<"\n";
    out_file.close();
  }// end if comm_in_rank == 0
} // end function

//==========================================================================
void PolymerChain::write_polymer_rog(const unsigned int& step_id,
                                     const unsigned int& o_step,
                                     const std::vector<Real>& lvec,
                                     unsigned int comm_in_rank) const
{
  //compute center of mass of each chain at step i
  std::vector<Point> centeri(_n_chains);
  this -> compute_center_of_mass(lvec, centeri);
  //compute radius of gyration:sqrt(Rg^2)
  std::vector<Real> RgSqrt(_n_chains);
  this -> compute_radius_of_gyration(lvec, centeri, RgSqrt);
  // write radius of gyration of each chain to out.radius_of_gyration
  std::ofstream out_file;  
  if(comm_in_rank == 0){
    if(step_id == 0){
      out_file.open("out.radius_of_gyration", std::ios_base::out);
      out_file <<"o_step" <<"  " << "radius_of_gyration\n";    
    }
    else{
      out_file.open("out.radius_of_gyration", std::ios_base::app);
    }
    out_file.precision(o_precision);
    out_file <<o_step  <<"  ";     
    for (std::size_t j = 0; j < _n_chains; j++){
      out_file << RgSqrt[j] <<"  ";
    }// end j loop 
    out_file <<"\n";
    out_file.close(); 
  }// end if comm_in_rank == 0
}

//==========================================================================
void PolymerChain::write_polymer_msd(const unsigned int& step_id,
                                     const unsigned int& o_step,
                                     const std::vector<Point>& center0,
                                     const std::vector<Real>& lvec,
                                     unsigned int comm_in_rank) const
{
  //compute center of mass of each chain at step i
  std::vector<Point> centeri(_n_chains);
  this -> compute_center_of_mass(lvec, centeri);
  // compute mean square displacement
  Point msd;
  this->compute_mean_square_displacement(center0, centeri, msd);
  // write mean square displacment to out.mean_square_displacement
  std::ofstream out_file;  
  if(comm_in_rank == 1){
    if(step_id ==0){
      out_file.open("out.mean_square_displacement", std::ios_base::out);
      out_file << "o_step" <<"  " << "msd" << "\n";      
    }
    else{
      out_file.open("out.mean_square_displacement", std::ios_base::app);
    }
    out_file.precision(o_precision);    
    out_file << o_step <<"  "; 
    for (int i = 0; i <_dim; i++){
      out_file <<msd(i) <<"  ";
    }
    out_file <<"\n";
    out_file.close(); 
  }// end if comm_in_rank == 0  
}

  
// ======================================================================
void PolymerChain::write_bead(const unsigned int& step_id,
                              const unsigned int& o_step,
                              const Real& real_time,
                              const std::vector<Point>& center0,
                              const std::vector<Real>& lvec,
                              const std::vector<std::string>& output_file,
                              unsigned int comm_in_rank) const
{
  START_LOG("write_bead()", "PolymerChain");
  this->write_time(step_id, o_step, real_time, comm_in_rank);
  for (int i = 0; i < output_file.size(); i++){
    if(output_file[i] == "equation_systems") {
      // this output has been written in Copss.C
    }
    else if(output_file[i] == "trajectory") this -> write_bead_trajectory(o_step, comm_in_rank);
    else if(output_file[i] == "center_of_mass") this -> write_bead_com(step_id, o_step,lvec, comm_in_rank);
    else if(output_file[i] == "mean_square_displacement") this -> write_bead_msd(step_id, o_step,center0, lvec, comm_in_rank);
    else {
      std::cout << "unsupported output_file content: (" << output_file[i] <<")" << std::endl; 
      libmesh_error();
    }
  }
  STOP_LOG("write_bead()", "PolymerChain");
}


void PolymerChain::write_bead_trajectory(const unsigned int& o_step,
                                         unsigned int comm_in_rank) const
{
  START_LOG ("write_bead_trajectory()", "PolymerChain");
  std::ostringstream oss;
  oss << "output_bead_" << o_step <<".csv";
  std::ofstream out_file;
  if (comm_in_rank == 0){
    out_file.open(oss.str(), std::ios_base::out);
    // write out the csv file  
    // POINT data
    out_file <<"scalar x_coord y_coord z_coord x_vel y_vel z_vel x_force y_force z_force\n";
    out_file.precision(o_precision);
    for(std::size_t i=0; i<_n_beads; ++i)
    {
      out_file << i << " ";
      // write position
      for(std::size_t j=0; j<_dim; ++j){
        out_file << _beads[i]->center()(j) << " ";
      }
      // write velocity
      for (std::size_t j=0; j<_dim; ++j){
        out_file <<_beads[i]->particle_velocity()[j] << " ";
      }
      // write force
      for (std::size_t j=0; j<_dim; ++j){
        out_file <<_beads[i]->particle_force()[j] <<" ";
      }
      out_file <<"\n";
    }
    out_file << "\n";
    out_file.close();
  } // end if comm_in_rank == 0
  STOP_LOG ("write_bead_trajectory()", "PolymerChain");
}

void PolymerChain::write_bead_com(const unsigned int& step_id,
                                  const unsigned int& o_step,
                                  const std::vector<Real>& lvec,
                                  unsigned int comm_in_rank) const
{
  START_LOG ("write_bead_com()", "PolymerChain");
  // compute center of mass at step i
  std::vector<Point> centeri(_n_chains);
  this->compute_center_of_mass(lvec, centeri);  
  // write center of mass of every chain to out.center_of_mass
  std::ofstream out_file;  
  if(comm_in_rank == 0){
    if (step_id == 0){
      out_file.open("out.center_of_mass", std::ios_base::out);
      out_file << "o_step" <<"  " << "center_of_mass";       
    }
    else{
      out_file.open("out.center_of_mass", std::ios_base::app);
    }
    out_file.precision(o_precision);
    out_file << o_step <<"  "; 
    for (std::size_t j = 0; j < _n_chains; j++){
            for (std::size_t k = 0; k < _dim; k++){
                    out_file << centeri[j](k) << " ";
            } // end k loop
    } // end j loop      
    out_file <<"\n";
    out_file.close();
  }// end if comm_in_rank == 0
  STOP_LOG ("write_bead_com()", "PolymerChain");
}

void PolymerChain::write_bead_msd(const unsigned int& step_id,
                                  const unsigned int& o_step,
                                  const std::vector<Point>& center0,
                                  const std::vector<Real>& lvec,
                                  unsigned int comm_in_rank) const
{
  START_LOG("write_bead_msd()", "PolymerChain");
  //get centers of all beads
  std::vector<Point> centeri(_n_beads);
  for (unsigned int i = 0; i < _n_beads; i++){
    for (int j = 0; j <_dim; j++){
      centeri[i](j) = lvec[i*_dim+j];
    }
  }
  // compute mean square displacement
  Point msd;
  this->compute_mean_square_displacement(center0, centeri, msd);
  // write mean square displacment to out.mean_square_displacement
  std::ofstream out_file;  
  if(comm_in_rank == 0){
    if(step_id == 0){
      out_file.open("out.mean_square_displacement", std::ios_base::out);
      out_file << "o_step" <<"  " << "mean_square_displacement" << "\n";      
    }
    else{
      out_file.open("out.mean_square_displacement", std::ios_base::app);
    }
    out_file.precision(o_precision);
    out_file << o_step <<"  ";
    for (int i = 0; i < _dim; i++){
      out_file <<msd(i) <<"  ";
    }
    out_file<<"\n";
    out_file.close();       
  }// end if comm_in_rank == 0  
  STOP_LOG("write_bead_msd()", "PolymerChain");
}

// ======================================================================
void PolymerChain::write_time(const unsigned int& step_id,
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
void PolymerChain::write_unfolded_polymer_chain(const std::string& filename) const
{
  START_LOG ("write_unfolded_polymer_chain()", "PolymerChain");
  
  // OFSTEAM
  std::ofstream outfile, cntfile;
  outfile.open(filename, std::ios_base::out);
  cntfile.open("pbc_count.dat", std::ios_base::out);

  const std::size_t n_beads = _beads.size();
  Point box_len = _periodic_boundary->box_length();
 
  // POINT data
  outfile << n_beads << " 0\n\n";
  for(std::size_t i=0; i<n_beads; ++i)
  {
    outfile << "16 ";

    std::vector<int> count = _beads[i]->counter();

    for(std::size_t j=0; j<3; ++j){
      outfile << _beads[i]->center()(j) + count[j] * box_len(j) << " ";
      cntfile << count[j] << " ";
    }
    outfile << "\n";
    cntfile << "\n";
  }
  
  outfile.close();
  cntfile.close();
 
  STOP_LOG ("write_unfolded_polymer_chain()", "PolymerChain");
}



// ======================================================================
void PolymerChain::write_unfolded_com(const std::string& filename) const
{
  START_LOG ("write_unfolded_com()", "PolymerChain");
 
  // OFSTEAM
  std::ofstream outfile, cntfile;
  outfile.open(filename, std::ios_base::out);
  cntfile.open("pbc_count.dat", std::ios_base::out);

  Point box_len = _periodic_boundary->box_length();
 
  // Calculate unfolded center of mass
  std::vector<std::vector<Real>> com_chains(_n_chains, std::vector<Real>(3,0.));

  for(std::size_t i=0; i<_n_beads; ++i)
  {
    int chain_id = _beads[i]->parent_id();

    std::vector<int> count = _beads[i]->counter();

    for(std::size_t j=0; j<3; ++j){
      com_chains[chain_id][j] += _beads[i]->center()(j) + count[j] * box_len(j);
      cntfile << count[j] << " ";
    }
    cntfile << "\n";
  }
  cntfile.close();

  // Output center of mass for each chain
  outfile << _n_chains << " 0\n\n";
  for(std::size_t i=0; i<_n_chains; ++i) 
  {
    outfile << "16 ";

    for(std::size_t j=0; j<3; ++j){
      outfile << com_chains[i][j]/_n_beads_per_chain[i] << " ";
    }
    outfile << "\n";
  }
  outfile.close();
 
  STOP_LOG ("write_unfolded_com()", "PolymerChain");
}



// ======================================================================
Point PolymerChain::spring_vector(const unsigned int i) const
{
  START_LOG ("spring_vector()", "PolymerChain");
  
  // Check: i cannot exceed the largest number of spring
  libmesh_assert_less(i, this->n_beads()-1);
  
  // Get the two connected beads
  const Point& b0 = _beads[i]->center();
  const Point& b1 = _beads[i+1]->center();
  
  // Compute the spring vector b1 - b0
  const Point b01 = this->bead_vector(b0,b1);
  
  STOP_LOG ("spring_vector()", "PolymerChain");
  return b01;
}

  

// ======================================================================
Point PolymerChain::end_to_end_vector() const
{
  START_LOG ("end_to_end_vector()", "PolymerChain");
  
  // R_end - R_start
  const std::size_t end_id = _beads.size() - 1;
  
  // Check: the polymer at least contains two beads!
  libmesh_assert_less(0, end_id); // end_id > 0 => 1,2,...
  
  // Get the two connected beads
  const Point& b0 = _beads[0]->center();
  const Point& b1 = _beads[end_id]->center();
  
  // Compute the spring vector b1 - b0
  const Point b01 = this->bead_vector(b0,b1);
  
  STOP_LOG ("end_to_end_vector()", "PolymerChain");
  return b01;
}

  
  
// ======================================================================
Point PolymerChain::end_to_end_vector_square() const
{
  START_LOG ("end_to_end_vector_square()", "PolymerChain");
  
  Real val = 0.;
  
  // # of springs; and loop over each spring
  const std::size_t Ns = _beads.size() - 1;
  for(std::size_t m=0; m<Ns; ++m)
  {
    // bond vector r_m
    const Point& b0 = _beads[m]->center();
    const Point& b1 = _beads[m+1]->center();
    const Point bm = this->bead_vector(b0,b1);
    
    for(std::size_t n=0; n<Ns; ++n)
    {
      const Point& b2 = _beads[n]->center();
      const Point& b3 = _beads[n+1]->center();
      const Point bn = this->bead_vector(b2,b3);
      
      val += bm(0)*bn(0) + bm(1)*bn(1) +bm(2)*bn(2);
    } // end for n
  } // end for m
  
  
  STOP_LOG ("end_to_end_vector_square()", "PolymerChain");
  return val;
}
  
  
  
// ======================================================================
Real PolymerChain::compute_chain_length()
{
  START_LOG ("compute_chain_length()", "PolymerChain");
  
  Real len = 0.0;
  
  // # of springs; and loop over each spring
  const std::size_t Ns = _beads.size() - 1;
  for(std::size_t i=0; i<Ns; ++i)
  {
    //
    const Point& b0 = _beads[i]->center();
    const Point& b1 = _beads[i+1]->center();
    
    // Evaluate the bead distance (spring length)
    const Point dpt = this->bead_vector(b0,b1);
    
    len += dpt.size();
  }
  
  STOP_LOG ("compute_chain_length()", "PolymerChain");
  return len;
}

// ======================================================================
void PolymerChain::compute_center_of_mass(const std::vector<Real>& lvec,
                                          std::vector<Point>& center) const
{
  START_LOG("compute_center_of_mass()", "PolymerChain");
  //loop over each chain
  std::size_t start_id = 0;
  for(std::size_t i = 0; i < _n_chains; i++){
   // loop over each bead on this chain
    for(std::size_t j = 0; j < _n_beads_per_chain[i]; j++){
      for (std::size_t k = 0; k < _dim; k++){
        center[i](k) += lvec[(start_id + j) * _dim + k];
      }// end k loop
    }// end j loop over each beads of a chain
    center[i] /= Real(_n_beads_per_chain[i]);
    start_id += _n_beads_per_chain[i]; 
  }//end i loop over each chain
  STOP_LOG("compute_center_of_mass()", "PolymerChain");
}

//=======================================================================
void PolymerChain::initial_bead_center_of_mass(std::vector<Point>& center0) const
{
  START_LOG("initial_bead_center_of_mass()", "PolymerChain");
  center0.resize(_n_beads);
  for (std::size_t i = 0; i < _n_beads; i++){
    center0[i] = _beads[i] -> point();
  }
  STOP_LOG("initial_bead_center_of_mass()", "PolymerChain");
}

// =====================================================================
void PolymerChain::initial_chain_center_of_mass(std::vector<Point>& center0) const
{
  START_LOG("initial_chain_center_of_mass()", "PolymerChain");
  center0.resize(_n_chains);
  std::size_t start_id = 0;
  for (std::size_t i = 0; i < _n_chains; i ++){
    for (std::size_t j = 0; j < _n_beads_per_chain[i]; j++){
      center0[i] += _beads[start_id + j]->point();
    } // end j loop
    center0[i] /= Real(_n_beads_per_chain[i]);
    start_id += _n_beads_per_chain[i];
  }// end i loop
  STOP_LOG("initial_chain_center_of_mass()", "PolymerChain");
}

// ======================================================================
void PolymerChain::compute_chain_stretch(const std::vector<Real>& lvec,
                                         std::vector<Point>& stretch) const
{
  START_LOG("compute_chain_stretch()", "PolymerChain");
  // loop over each particle and find the max/min (X,Y,Z) 
  std::vector<Point> max_xyz(_n_chains);
  std::vector<Point> min_xyz(_n_chains);
  Point pti;
  //loop over each chain
  std::size_t start_id  = 0;
  for (std::size_t i = 0; i < _n_chains; i++){
    // loop over each bead in this chain
    for (std::size_t j=0; j < _n_beads_per_chain[i]; j++)
    {
     // loop over each direction
     for (std::size_t k = 0;  k < _dim; k++){
       pti(k) = lvec[(start_id + j)*_dim + k];
     }// end k loop
     // compute min_xyz and max_xyz of this chain
     if (j==0){
       min_xyz[i] = pti;
       max_xyz[i] = pti;
     } // end if
     else{
       for (std::size_t k = 0; k < _dim; k++){
         if (pti(k) < min_xyz[i](k)) {min_xyz[i](k) = pti(k);}
         if(pti(k) > max_xyz[i](k)) {max_xyz[i](k) = pti(k);}
       }// end loop k
     }// end else
    }// end loop j
    start_id += _n_beads_per_chain[i];
  }// end loop i 
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Compute the stretch in different directions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  for(std::size_t i =0; i < _n_chains; i++) {
    for(std::size_t k=0; k<_dim; k++ ){
       stretch[i](k) = max_xyz[i](k)-min_xyz[i](k);
    }// end dim
  }//end chain_id

  STOP_LOG("compute_chain_stretch()", "PolymerChain");
}

void PolymerChain::compute_radius_of_gyration(const std::vector<Real>& lvec,
                                              const std::vector<Point>& center,
                                              std::vector<Real>& RgSqrt) const
{
  START_LOG("compute_radius_of_gyration()", "PolymerChain");
  std::size_t start_id = 0;
  // loop over each chain
  for (std::size_t i = 0; i < _n_chains; i++)
  {
    Real Rg = 0.;
    // loop over each bead on this chain
    for (std::size_t j = 0; j < _n_beads_per_chain[i] ; j++){
     Point pti;
     for (std::size_t k = 0; k < _dim; k++){
       pti(k) = lvec [ (start_id + j) *_dim + k] - center[i](k);
     } // end k loop
     Rg += pti.norm_sq() / _n_beads_per_chain[i];
    } // end j loop
    RgSqrt[i] = std::sqrt(Rg);
  }
  STOP_LOG("compute_radius_of_gyration()", "PolymerChain");
}

void PolymerChain::compute_mean_square_displacement(const std::vector<Point>& R0,
                                                    const std::vector<Point>& Rt,
                                                    Point& msd) const
{
  START_LOG("compute_mean_square_displacement()", "PolymerChain");
  if(R0.size() != Rt.size()) {
    std::cout <<"Warning in PolymerChain::compute_mean_square_displacement(R0, Rt): size of R0 (="<<R0.size()
    <<") and Rt (=" << Rt.size() 
    <<") are not the same" << std::endl;
    libmesh_error();
  }
  unsigned int num_center = R0.size();
  for (unsigned int i = 0; i < num_center; i++){
    Point displacement = Rt[i] - R0[i];
    for (int j = 0; j < _dim; j++){
      msd(j) += displacement(j) * displacement(j);
    }
  }
  msd /= Real(num_center);
  STOP_LOG("compute_mean_square_displacement()", "PolymerChain");
}

// ======================================================================
bool PolymerChain::check_chain(const Real& Ls)
{
  START_LOG ("check_chain()", "PolymerChain");
  
  bool chain_broken = false;
  
  // loop over each bond
  for(std::size_t i=0; i<_n_bonds; ++i)
  {
    // IDs of two connected beads
    std::size_t id_bead_1 = _bonds[i][1];
    std::size_t id_bead_2 = _bonds[i][2];

    const Point& b0 = _beads[id_bead_1]->center();
          Point& b1 = _beads[id_bead_2]->center();
 
    // Evaluate the bond length
    const Point dpt = this->bead_vector(b0,b1);

    // Check if the length is larger than the maximum length
    Real len = dpt.norm();

    if( len >= Ls ){
      chain_broken = true;

      // Move bead (b1) to make the chain length less than maximum length
      Point new_b1 (b0);
      new_b1.add_scaled(dpt, Ls * 0.9 / len); // Scale the broken chain length to be 0.9*max_chain_length
      b1 = new_b1;

      // Move bead (b1) into the box when considering PBC
      if(_periodic_boundary!=NULL){
        _periodic_boundary->correct_position(b1); // Move new position into box when considering PBC
      }
    }
  }

  STOP_LOG ("check_chain()", "PolymerChain");
  return chain_broken;
}




// ======================================================================
Point PolymerChain::bead_vector(const Point& bead0,
                                const Point& bead1) const
{
  // Evaluate the bead distance (spring length)
  Point dpt;
  if(_periodic_boundary==NULL) {
    dpt = bead1 - bead0;
  }
  else {
    dpt = _periodic_boundary->point_vector(bead0,bead1);
  }
  
  // return
  return dpt;
}
 

} // end of the namespace
