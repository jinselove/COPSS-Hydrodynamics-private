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

// Local Includes -----------------------------------
#include "assemble_stokes.h"
#include "solver_stokes.h"
#include "pm_linear_implicit_system.h"

namespace libMesh
{

/*
 * The PMSystemStokes is designed to solve the Stokes
 * equation with regularized Gaussian point forces
 * due to particles using mixed FEM.
 */

//class Fix;

class PMSystemStokes : public PMLinearImplicitSystem
{
public:

  /**
   * Constructor.
   */
  PMSystemStokes (EquationSystems& es,
                  const std::string& name,
                  const unsigned int number); // number of systems


  /**
   * Destructor.
   */
  virtual ~PMSystemStokes ();


  /**
   * The type of system.
   */
  typedef PMSystemStokes sys_type;


  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }


  /*
   * Re-init particle mesh, including:
   * (1) reinit() reinit point particles
   *              build_particle_neighbor_list()
   *              build_elem_neighbor_list()
   * (2) update the mesh of each finite sized particle if there are;
   * (3) compute particle force (by force field)
   *             modify the force field according to the vel_last_step.
   */
  void reinit_hi_system(bool& neighbor_list_update_flag);


  /*
   * Re-init free-draining system, including:
   * (1) reinit() reinit point particles
   *              build_particle_neighbor_list()
   * (2) compute particle force (by force field)
   *             modify the force field according to the vel_last_step.
   */
  void reinit_fd_system(bool& neighbor_list_update_flag);


  /**
   * Clear all the data structures associated with the system.
   */
  void clear ();


  /**
   * Assemble the system matrix.
   * option == "undisturbed" : undisturbed velocity field (without particles)
   * option == "disturbed"   : disturbed velocity field   (with particles)
   */
  void assemble_matrix (const std::string& system_name,
                        const std::string& option);


  /**
   * Assemble the system rhs.
   */
  void assemble_rhs (const std::string& system_name,
                     const std::string& option);


  /*
   * Solve the Stokes system.
   * option = "undisturbed", compute the undisturbed field of flow without particles
   * option = "disturbed",   compute the disturbed field of flow with particles
   * re_init = true => re-assemble the matrix and reinit the KSP solver.
   */
  void solve (const std::string& option,
              const bool& re_init);


  /*
   * Add the local solution to the global solution
   */
  void add_local_solution();


  /*
   * Compute the L2-error in an unbounded domain
   * This function will change the system solution vector by add local solution.
   */
  void test_l2_norm(bool& neighbor_list_update_flag);


  /*
   * Write out equation systems of Stokes. This requires combining the
   * local and global solution, and update the solution vectors.
   * output_format=="EXODUS" ; "VTK"; "GMV"
   */
  void write_equation_systems(const std::size_t time_step,
                              const std::string& output_filename,
                              const std::string& output_format);


  /*
   * Return the SolverStokes
   */
  SolverStokes& stokes_solver() { return _stokes_solver;  }


  /**
   * Compute the purturbed velocity vector of moving points
   * according to physical fields. which is a dim*N vector
   * NOTE that the "self-term" should be excluded, which is done in GGEMSystem.
   * See ref.
   * [1] Y Zhang, J J de Pablo, and M D Graham, J Chem Phys (2012) 014901.
   * [2] P Pranay, S G Anekal, J P Hernandez-Ortiz, Phys Fluids (2010)
   *
   * option = "disturbed" or "undistrubed"
   * pv = particle_velocity
   */
  void compute_point_velocity(const std::string& option, std::vector<Real>& pv);


  /**
   * Compute the unperturbed velocity vector at the location of particles
   */
  std::vector<Real> compute_unperturbed_point_velocity();


  /*
   * Obtain the velocity vector on the i-th point
   */
  std::vector<Real> point_velocity(const std::vector<Real>& vel_beads,
                                   const std::size_t i) const;


  /**
   * Local velocity of a fluid point in an unbounded space,
   * which is computed from Green's function
   * force_type: "regularized" or "smooth"
   */
  std::vector<Real> local_velocity_fluid(const Point &p,
                                         const std::string& force_type) const;


  /**
   * Local velocity of a fluid point in an unbounded space,
   * which is computed from Green's function
   * force_type: "regularized" or "smooth"
   */
  std::vector<Real> local_velocity_fluid(const Elem* elem,
                                         const Point &p,
                                         const std::string& force_type) const;


  /**
   * Local velocity of a bead in an unbounded space,
   * which is computed from Green's function.
   * force_type: "regularized"  or "smooth"
   */
  Point local_velocity_bead(const std::size_t& bead_id,
                                        const std::string& force_type) const;


  /*
   * Self-exclusion term for the velocity at the i-th bead
   */
  Point global_self_exclusion(const std::size_t p_id) const;


  /*
   * Test function. output velocity profile along x-y-z direction
   */
  void test_velocity_profile(bool& neighbor_list_update_flag);


  /*
   * Return the exact solution of Stokes eqn for any given
   * point \pt0 in an unbounded domain.
   */
  //const std::vector<Real> exact_solution(const Point& pt0) const;


  /*
   * Write nodal coordinates (x, y, z) and velocities (vx, vy, vz)
   * from system's solution vector to raw data file for step i.
   * If you want to write the total velocites, this function MUST ONLY
   * be called after you *add_local_solution* and *unperturbed* velocities
   * to the solution vector.
   */
  void write_fluid_velocity_data(const std::string& filename);



private:

  // Stokes solver
  SolverStokes _stokes_solver;

  // Assemble Stokes system
  AssembleStokes* _assemble_stokes;

};

} // end namespace libMesh
