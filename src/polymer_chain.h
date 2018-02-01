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
#include "math.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/parallel_object.h"


#include "point_particle.h"
#include "pm_toolbox.h"




namespace libMesh
{

/*
 * This class defines a polymer chain for bead-spring model,
 * Note that this is only used for a single chain.
 */

  class PMPeriodicBoundary;
  
  
class PolymerChain :  public ReferenceCountedObject<PolymerChain>
//                      public ParallelObject
{
public:

  // Constructor
  PolymerChain(const std::size_t chain_id);
  
  
  // Constructor with periodic boundary
  PolymerChain(const std::size_t chain_id,
               PMPeriodicBoundary& pm_pb);
  
  
  // Destructor
  ~PolymerChain();
  
  
  /* 
   * Read the data of chain from local file, whose data structure is:
   * Nb
   * global#  local#  bead_type x0 y0 z0 a0 b0 c0
   * global#  local#  bead_type x1 y1 z1 a1 b1 c1
   * ......
   * Generated by Pizza.py chain tool.
   * Note this is only used to read one chain data.
   */
  void read_data(const std::string& filename);      // particle xyz file);
  
  
  /*
   * Read the data of chain from local file 
   * Generated by Pizza.py chain tool.
   */
  void read_data_pizza(const std::string& filename);
  
  
  /*
   * Read the data of chain from local file with vtk format.
   * This will be useful to restart the simulation
   */
  void read_data_vtk(const std::string& filename);
  
  /*
   * Read the data of chain from local file with csv format.
   * This will be useful to restart the simulation
   */
  void read_data_csv(const std::string& filename); 


   /*
   * Read the data of PBC counter from local file.
   * This is useful to restart the simulation
   */
  void read_pbc_count();


  /*
   * Add a bead to the polymer chain, bonded with \bond_bead
   */
  void add_bead(const Point& pt,
                const PointType point_type,
                const int parent_id,
                const std::vector<Real>& rot_vec);
  
 
/* generate multiple polymer_chains for simulation
 
 ========================================================================================================= 
 * !n_chains  ---> number of polymer chains you want to generate
 * !n_beads_perChain ---> number of beads per chain (all chains have the same number of beads)
 * !LS ---> the maximum bond length
 * !init_bbox_min   ---> min of the box you want to generate your initial polymers in.
 * !init_bbox_max   ---> mat of the box you wnat to generate your initial polymers in.
 * !filename	    ---> where you want to output the polymer data you generated.
 */
 void generate_polymer_chains(const unsigned int n_beads,
                             const unsigned int n_chains,
                             const Real init_bondLength, 
                             const Point init_bbox_min,
                             const Point init_bbox_max,
                             const std::string& filename,
                             unsigned int comm_in_rank);


 
  
  /************ NOT ready for use right now!
   *
   * Generate a polymer chain within the bounding box with random walk.
   * pt0 defines the location of the first bead;
   * Number of beads Nb, and number of spring Ns = Nb - 1;
   * Equilibrium Length of spring is Ls.
   */
  void generate_polymer_chain(const Point pt0,
                              const std::size_t Nb, const Real Ls,
                              const Point& bbox_min,
                              const Point& bbox_max,
                              const std::string& filename);
  
  
  
  /*
   * Print out the polymer info
   */
  void print_info() const;
  
  
  /*
   * Write out the polymer chain to a local file
   */
  void write_polymer_chain(const unsigned int& step_id,
                           const unsigned int& o_step,
                           const Real& real_time,
                           const std::vector<Point>& center0,
                           const std::vector<Real>& lvec,
                           std::vector<std::string>& output_file,
                           unsigned int comm_in_rank) const;
  /*
   * Write out trajectories of all polymers (in vtk format)
   * "output_polymer_o_step.vtk"
   */
  void write_polymer_trajectory (const unsigned int& o_step,
                                 unsigned int comm_in_rank) const;

  /*
   * Write out center of mass of each chain (in vtk format)
   * "out.center_of_mass"
   */
  void write_polymer_com(const unsigned int& step_id,
                         const unsigned int& o_step,
                         const std::vector<Real>& lvec,
                         unsigned int comm_in_rank) const;
  /*
   * Write out stretch of each polymer chain
   * "out.stretch"
   */
  void write_polymer_stretch(const unsigned int& step_id,
                             const unsigned int& o_step,
                             const std::vector<Real>& lvec,
                             unsigned int comm_in_rank) const;
  /*
   * Write out radius of gyration of each polymer chain
   * "out.radius_of_gyration"
   */
  void write_polymer_rog(const unsigned int& step_id,
                         const unsigned int& o_step,
                         const std::vector<Real>& lvec,
                         unsigned int comm_in_rank) const;

  /*
   * Write out mean square displacement
   * Average over all chains
   * "out.mean_square_displacement"
   */
  void write_polymer_msd(const unsigned int& step_id,
                         const unsigned int& o_step,
                         const std::vector<Point>& center0,
                         const std::vector<Real>& lvec,
                         unsigned int comm_in_rank) const;

  /*
   * Write out the bead to a local file
   */
  void write_bead(const unsigned int& step_id,
                  const unsigned int& o_step,
                  const Real& real_time,
                  const std::vector<Point>& center0,
                  const std::vector<Real>& lvec,
                  const std::vector<std::string>& output_file,
                  unsigned int comm_in_rank) const;

  /*
   * Write out bead trajectories && forces && velocity to csv files
   * "output_bead_o_step.csv"
   */
  void write_bead_trajectory(const unsigned int& o_step,
                             unsigned int comm_in_rank) const;

  /*
   * Write out center of mass 
   * average over all beads, output is a point
   * this center of mass usually is of no use, but I implemented it just in case
   * "out.center_of_mass"
   */
  void write_bead_com(const unsigned int& step_id,
                      const unsigned int& o_step,
                      const std::vector<Real>& lvec,
                      unsigned int comm_in_rank) const; 

  /*
   * Write out mean square displacement of all beads
   * Average over all beads
   * "out.center_of_mass"
   */
  void write_bead_msd(const unsigned int& step_id,
                      const unsigned int& o_step,
                      const std::vector<Point>& center0,
                      const std::vector<Real>& lvec,
                      unsigned int comm_in_rank) const;
  /*
   * write out step and real time
   */
  void write_time(const unsigned int& step_id,
                  const unsigned int& o_step,
                  const Real& real_time,
                  unsigned int comm_in_rank) const;

  /*
   * Write out the unfolded polymer chain to a local file in .xyz format
   */
  void write_unfolded_polymer_chain(const std::string& filename) const;


  /*
   * Write out the center of mass of unfolded polymer chain
   * to a local file in .xyz format
   */ 
  void write_unfolded_com(const std::string& filename) const;


  /*
   * Return the chain id
   */
  std::size_t chain_id() const {  return _chain_id; }
 

  std::size_t n_chains() const { return _n_chains;} 
  
  /*
   * Return the total number of beads
   */
  std::size_t n_beads() const {  return _beads.size(); }


  /*
   * Return the total number of bonds
   */
  std::size_t n_bonds() const {  return _n_bonds; }


  /*
   * Return the total number of springs
   */
  std::size_t n_springs() const {  return _beads.size()-1; }
  
  
  /*
   * The beads vector for both const Ref and non-const Ref
   */
  const std::vector<PointParticle*>& beads() { return _beads;  };
 
 
  /*
   * The bonds vector for both const Ref and non-const Ref
   */
  const std::vector<std::vector<std::size_t> >& bonds() { return _bonds;  };


  /*
   * The i-th spring vector of coarse grained polymer chain
   * R(i+1) - R(i), where R is position vector of the bead
   */
  Point spring_vector(const unsigned int i) const;
  
  
  /*
   * Return the end-to-end vector of a polymer chain
   * sum[n=1:N]r(n), where r(n) is the n-th bond vector
   */
  Point end_to_end_vector() const;
  
  
  /*
   * Return the average square end-to-end vector of a polymer chain
   * sum[m=1:N]sum[n=1:N]r(m)r(n)
   * for an ideal chain, <R2> = N*b^2
   */
  Point end_to_end_vector_square() const;
  

  /*
   * Compute the current length of the chain
   * Note the chain length is NOT the contour length of the chain
   */
  Real compute_chain_length();
  
  
  /*
   * Contour length of the chain
   */
  //const Real contour_length() const;
  
  /*
   * Compute center_of_mass of each chain
   */ 
  void compute_center_of_mass(const std::vector<Real>& lvec,
                              std::vector<Point>& center) const;

  /*
   * Compute center of mass at step 0: center0, for beads
   */
  void initial_bead_center_of_mass(std::vector<Point>& center0) const;

  /*
   * Compute center of mass at step 0: center0, for chains
   */
  void initial_chain_center_of_mass(std::vector<Point>& center0) const;

  /*
   * Compute chain_stretch of each chain
   */
  void compute_chain_stretch(const std::vector<Real>& lvec,
                             std::vector<Point>& stretch) const;

  /*
   * Compute radius_of_gyration of each chain
   */
  void compute_radius_of_gyration(const std::vector<Real>& lvec,
                                  const std::vector<Point>& center,
                                  std::vector<Real>& RgSqrt) const;

  
  void compute_mean_square_displacement(const std::vector<Point>& R0,
                                        const std::vector<Point>& Rt,
                                        Point& msd) const;
  /*
   * Check if the chain is broken. If so, repair it.
   * It loops over every bond, so the created chain MUST
   * have bond information.
   */
  bool check_chain(const Real& Ls);      // Maxium length of the spring.
                   
  
  /*
   * Return the direction vector of two beads
   */
  Point bead_vector(const Point& bead0,
                    const Point& bead1) const;
  

private:
  
  // -------------------- New data structure -------------------
  
  // number of atoms(beads),
  std::size_t _n_beads;
  
  // number of bonds(springs)
  std::size_t _n_bonds;
 
  // number of chains
  int _n_chains;
  
  // number of atom types
  std::size_t _n_bead_types;
  
  // number of bond types
  std::size_t _n_bond_types;
 
  // number of beads on each chain
  std::vector<size_t> _n_beads_per_chain;
 
  // mass values of atoms(beads)
  std::vector<Real> _mass;
  
  // atom id, chain id and atom types are stored in PointParticle
  
  // atoms(beads) info
  std::vector<PointParticle*> _beads;
  
  // bonds info
  std::vector<std::vector<std::size_t> > _bonds;
  // -----------------------------------------------------------
  
  
  // the identity of the chain
  std::size_t _chain_id;

  // dimension of the system
  const unsigned int _dim = 3;
  
  // Bead list: two neighboring beads are connected by a spring.
  //std::vector<PointParticle*> _beads;
  
  // Pointer to the Periodic boundary condition
  PMPeriodicBoundary* _periodic_boundary;

  // precision of output
  const int o_precision = 6;
}; // end of class
  
}  // end of namespace

