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
#include <map>
#include <cmath>

#include "pm_linear_implicit_system.h"
#include "elasticity_system.h"
#include "force_field_base.h"
#include "pm_toolbox.h"

namespace libMesh
{

  /*
   * The class is designed for computing the force field
   * of the ParticleMesh system, including
   * -- spring force
   * -- excluded volume force
   * -- gravity
   * -- other user defined forces
   */
  
  
class ForceField : public ForceFieldBase
{
public:
  
  typedef std::string type_force_name;
  
  typedef std::vector<Real> type_force_parameter;
  
  typedef std::pair <type_force_name, type_force_parameter> type_force;

  // Constructor for a system with elastic structures (finite size particles)
  ForceField(PMLinearImplicitSystem& pm_sys,
             ElasticitySystem& el_sys);

  
  // Constructor for a system with point particles
  ForceField(PMLinearImplicitSystem& pm_sys);
  
  
  // Destructor
  ~ForceField();
  
  
  typedef ForceFieldBase Parent;
  
  /*
   * Check force field parameters when build force field classes.
   * This function is only evaluated once at the beginning of a simulation
   * check the following things:
   * 1) if this "force field" is supported;
   * 2) if this "force field" is applicable to this "particle type"
   * 3) if correct number of parameters are assigned to this "force field"
  */

  void check_force_field(); 


  /* 
   * User-defined function for reinitializing the force field
   * at each time step
   * The velocities of all beads are given (for calculating friction)
   * This function is the same to reinit_force_field() unless velocities are given, 
   * i.e., friction need to be calculated
   */
  // virtual void reinit_force_field(const std::vector<Real>& v_beads);
  
   /* 
   * User-defined function for reinitializing the force field
   * at each time step
   * The velocities of all beads are given
   */
  void reinit_force_field() override;


   /* 
   * Attach surface spring force
   * spring constant k0 is needed
   */  
  void attach_p_surface_constraint(const std::vector<Real>& params);

  /*
   * compute the constraint forces for the i-th rigid particle
   * This is achieved by applying large stiff spring forces.
   * k0 is a large number for the stiffness constant.
   */
  void compute_surface_constraint_force(const std::size_t& i, // the i-th particle
                              const Real& k0,
                              std::vector<Point>& nodal_force);

 

  void attach_nodal(const std::vector<Real>& params); 
   /*
   * Attach the particle-particle Excluded volume Force (Gaussian form)
   * This force occurs between the particle and all its neighbors,
   * which are non-bonded interactions.
   */
  void attach_pp_ev_gaussian(const std::vector<Real>& params);

   /*
   * Attach the particle-particle Excluded volume Force (Gaussian form)
   * This force occurs between the particle and all its neighbors,
   * which are non-bonded interactions.
   * This gaussian force is specially for beads on polymer_chains
   */
  void attach_pp_ev_gaussian_polymerChain(const std::vector<Real>& params);
  
  /* Compute the particle-particle Excluded volume Force on a given particle
   * (Gaussian Form)
   * This force occurs between the particle and all its neighbors,
   * which are non-bonded interactions.
   * potential: u_gaussian  = 1/2 * c1 * exp( -c2*|r_ij|^2 )
   * where r_ij is dimensionless : r_ij = R_ij / a and r_ij = rj - ri 
   * where c1 > 0, c2 > 0
   * force    : f_ij = -du_gaussian/dr_ij = -c1*c2* exp( -c2*|r_ij|^2 ) * r_ij
   * there is a negative sign before the force because r_ij = r_j - r_i, which has
   * an opposite direction as f_ij
   */
  void compute_pp_ev_gaussian(Real& c1,
                              Real& c2,
                              const std::size_t&  p_id,
                              std::vector<Real>& pforce );

  /*
   * Attach particle-particle Excluded volume Force (Lennard-Jones-cut form)
   * This force occurs between the particle and all its neighbors,
   * which are non-bonded interactions.   
   */

  void attach_pp_ev_lj_cut(const std::vector<Real>& params);

  /*
   * Attach particle-particle Excluded volume Force (pure repulsive Lennard-Jones)
   * This force occurs between the particle and all its neighbors,
   * which are non-bonded interactions.   
   */
  void attach_pp_ev_lj_repulsive(const std::vector<Real>& params);

  /* Compute the particle-particle Excluded volume Force on a given particle
   * (lj form, only calculate the force when r_ij < rcut)
   * This force occurs between the particle and all its neighbors,
   * which are non-bonded interactions.
   * potential: u_lj = 4 * epsilon *[(sigma / |r_ij|)^12 - (sigma / |r_ij|)^6]
   * f_ij = -24 * epsilon * (2*(sigma/|r_ij|)^12 - (sigma/|r_ij|)^6 ) * r_ij / |r_ij|^2
   * where r_ij (vector) is dimensionless : r_ij = R_ij / a and r_ij = rj - ri
   * there is a negative sign before the force because r_ij = r_j - r_i, which has
   * an opposite direction as f_ij
   */
  void compute_pp_ev_lj_cut(const Real& epsilon,
                            const Real& sigma,
                            const Real& rcut,
                            const std::size_t&  p_id,
                            std::vector<Real>& pforce) const;

  /*
   * Attach particle-particle Excluded volume Force (Harmonic form)
   * This force occurs between the particle and all its neighbors,
   * which are non-bonded interactions.   
   */
  void attach_pp_ev_harmonic_repulsive(const std::vector<Real>& params);

  /*
   * Compute the particle-particle Excluded volume force on a given particle
   * Harmonic form, only repulsive, i.e., harmonic force = 0 if distance > sigma
   * This force occurs between the particle and all its neighbors
   * which are non-bonded interactions
   */
  void compute_pp_ev_harmonic_repulsive(const Real& k,
                                        const Real& r0,
                                        const std::size_t& p_id,
                                        std::vector<Real>& pforce) const;
   /*
   * Attach spring forces to multiple polymer chains. Assume the chain 
   * are connected by Worm-like springs (WLS).
   * The spring force only occurs between two connected beads.
   * This function requires information on bonds, that is read
   * from data files prepared by Pizza toolkit.
   */
  void attach_pp_wormLike_spring();
  
  /*
   * Compute the Spring Force on a given particle. Assume the chain 
   * are connected by Worm-like springs (WLS).
   * The spring force only occurs between two connected beads
   */
  // void compute_pp_wormLike_spring(const std::size_t  p_id,
  //                               std::vector<Real>& pforce ) const;

    /*
   * Compute the Spring Force on a given particle. Assume the chain
   * are connected by FENE springs.
   * The spring force only occurs between two connected beads
   */
  // void compute_pp_fene_spring(const std::size_t  p_id,
  //                                std::vector<Real>& pforce ) const;

  /*
   * attach constant force read from force file to all particles
   * Each particle has to be assigned a force vector in the force file
   */
  void attach_p_constant(const std::vector<Real>& params);

  /*
   * attach gravity force read from force file to all particles
   * Each particle has to be assigned a force vector in the force file
   * the params are body force density
   */
  void attach_p_gravity(const std::vector<Real>& params);




  /*
   * User-defined function that to calculate force between different objects
   */ 
  // void attach_pp_friction(const std::vector<Real>& vel_beads);
 
   /*
   * Attach the particle-wall Excluded volume Force
   * Empirical form
   * This force occurs between the particle and the wall,
   * which are non-bonded interactions.
   * This empirical force is specially for beads on polymer_chains
   * @Jendrejack, R. M., Schwartz, D. C., Graham, M. D., & de Pablo, J. J. (2003). Effect of confinement on DNA dynamics in microfluidic devices. The Journal of Chemical Physics, 119(2), 1165–10. http://doi.org/10.1063/1.1575200
   */
  void attach_pw_ev_empirical_polymerChain();

  /*
   * Compute the particle-wall Exclude volume Force on a given particle
   * Empirical form
   * This force occurs between the particle and wall, but
   * does NOT in the periodic boundary direction without walls!
   * The current algorithm only consider the parallel walls with simple geometries.
   */
  void compute_pw_ev_empirical_polymerChain(const std::size_t&  p_id,
                                   std::vector<Real>& pforce ) const;

  /*
   * Attach the particle-wall Exclude volume Force on a given particle
   * Lennard-Jones cut form
   * This force occurs between the particle and wall, but
   * does NOT in the periodic boundary direction without walls!
   * 
   * The current algorithm only consider the parallel walls with simple geometries.
   */
  void attach_pw_ev_lj_cut(const std::vector<Real>& params);

  /*
   * Attach particle-wall Excluded volume Force
   * (pure repulsive Lennard-Jones)
   * This force occurs between the particle and all its neighbors,
   * which are non-bonded interactions.   
   */

  void attach_pw_ev_lj_repulsive(const std::vector<Real>& params);

  /* Compute the particle-wall Excluded volume Force on a given particle
   * (lj form, only calculate the force when r_ij < rcut)
   * This force occurs between the particle and wall,
   * which are non-bonded interactions.
   */
  void compute_pw_ev_lj_cut(const Real& epsilon,
                            const Real& sigma,
                            const Real& rcut,
                            const std::size_t&  p_id,
                            std::vector<Real>& pforce) const;


  /*
   * Attach the particle-wall Exclude volume Force on a given particle
   * harmonic form
   * This force occurs between the particle and wall, but
   * does NOT in the periodic boundary direction without walls!
   * The current algorithm only consider the parallel walls with simple geometries.
   */
  void attach_pw_ev_harmonic_repulsive(const std::vector<Real>& params);


  /*
   * Compute the particle-wall Excluded volume Force on a given particle
   * Harmonic form, only calculate the force when r_i_wall < r0
   * This force occurs between the particle and the wall, 
   * which are non-bonded interactions
   */
  void compute_pw_ev_harmonic_repulsive(const Real& k,
                                        const Real& r0,
                                        const std::size_t& p_id,
                                        std::vector<Real>& pforce) const;
   /*
   * Re-init the force vector of each point (set to zeros)
   * This is required at the beginning of each time step when
   * re-compute the forces on each point.
   */
  void reinit_point_force_zeros();
  
  
  
  /*
   * User-defined function for reinitializing the force field
   * of a DNA(polymer) chain at each time step
   */
//  void compute_force_field_dna();
 
  
  /*
   * User-defined function that incorporates:
   * (a) excluded volume force for particles
   * (b) the friction force between different objects.
   */
  // void modify_force_field(const std::vector<Real>& vel_beads);

 

  /*
   * Correct the bead position if it moves out of the wall or periodic boundary.
   */
  void check_wall(const std::size_t& p_id);


  /* 
   * Correct the bead position if it moves out of the wall or periodic boundary, and
   * count how many times a point particle has cross the boundary, which is used to
   * calculate its unfolded position.
   */
  void check_wall_pbc_count(const std::size_t& p_id);


  /* 
   * Correct the CHAIN position if it moves out of the wall or periodic boundary.
   */ 
  void check_walls();
  

  /* 
   * Correct the CHAIN position if it moves out of the wall or periodic boundary, and
   * count how many times each bead has cross the boundary, which is used to
   * calculate its unfolded position.
   */ 
  void check_walls_pbc_count();

  

  /*
   * Return a constant force vector.
   * This is a general case of the above function gravity_force()
   */
  // void constant_force(std::vector<Real>& pforce ) const;
  
  
private:
  
  // the particle-mesh (linear implicit) system
  PMLinearImplicitSystem* _pm_system;
  
  // The elastic system for solids
  // Use pointer instead of const Ref. to avoid explicit initialization!
  ElasticitySystem* _elastic_system;

  // partiticle mesh for finite size particle
  ParticleMesh<3>* _particle_mesh;

  // point mesh for point particle
  PointMesh<3>* _point_mesh; 

  // system dimension
  unsigned int _dim;

  // count total number of times beads escape the simulation domain.
  unsigned int _out_domain_counter = 0;

  // boundary conditions
  std::string _wall_type;

  std::vector<Real> _wall_params;

  Point _box_min;

  Point _box_max;

  Point _box_len;

  std::vector<bool> _periodic;

  std::vector<bool> _inlet;  
  
  // particle type
  std::string _particle_type;

  // point particle model
  std::string _point_particle_model;

  // number of finite size particles
  std::size_t _num_particles = 0;

  // number of point particles
  std::size_t _num_points = 0;

  // particle-particle force types
  std::vector<std::string> _pp_force_types;

  int _num_pp_force_types;

  // particle-particle force types
  std::vector<std::string> _pw_force_types;
  int _num_pw_force_types;



  // enum force_types
  enum pp_forces {pp_ev_gaussian = 0, 
                  pp_ev_gaussian_polymerChain,
                  pp_ev_lj_cut,
                  pp_ev_lj_repulsive,
                  pp_ev_harmonic_repulsive, 
                  pp_wormLike_spring, 
                  p_constant,
                  p_gravity,
                  p_surface_constraint};

  enum pw_forces{pw_ev_empirical_polymerChain = 0,
                 pw_ev_lj_cut,
                 pw_ev_lj_repulsive,
                 pw_ev_harmonic_repulsive};

  // string to int map for force
  std::map <std::string,int> forceTypeMap;

  // Coefficients used to compute the non-dimensional forces (values will be updated for particle_type == "point_particle")
  Real _bead_r = 0;
  Real _Ss2 = 0;
  Real _bk = 0;
  Real _Nks = 0;
  Real _kBT = 0;
  
};  // end of class
  
} // end of namespace
