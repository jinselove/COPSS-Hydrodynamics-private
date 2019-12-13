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
#include "analytical_solution_stokes.h"
#include "ggem_stokes.h"

namespace libMesh {
/*
 * The PMSystemStokes is designed to solve the Stokes
 * equation with regularized Gaussian point forces
 * due to particles using mixed FEM.
 */

class PMSystemStokes : public PMLinearImplicitSystem {
public:

  /**
   * Constructor.
   */
  PMSystemStokes(EquationSystems  & es,
                 const std::string& name,
                 const unsigned int number); // number of systems


  /**
   * Destructor.
   */
  virtual ~PMSystemStokes();


  /**
   * The type of system.
   */
  typedef PMSystemStokes sys_type;


  /**
   * @returns a clever pointer to the system.
   */
  sys_type& system() {
    return *this;
  };

  /**
   * override get_dt function in parent class since this system needs more
   * conditions on dt to make the system stable
   */
  Real get_dt();

  /*
   * Re-init particle mesh, including:
   * (1) reinit() reinit point particles
   *              build_particle_neighbor_list()
   * (2) update the mesh of each finite sized particle if there are;
   * (3) compute particle force (by force field)
   *             modify the force field according to the vel_last_step.
   * only if option=="disturbed", we calculate forces
   */
  void reinit_system(bool& neighbor_list_update_flag,
                     const std::string& option);


  /**
   * Clear all the data structures associated with the system.
   */
  void clear();


  /**
   * Assemble the system matrix.
   * option == "undisturbed" : undisturbed velocity field (without particles)
   * option == "disturbed"   : disturbed velocity field   (with particles)
   */
  void assemble_matrix(const std::string& system_name,
                       const std::string& option) override;


  /**
   * Assemble the system rhs.
   */
  void assemble_rhs(const std::string& system_name,
                    const std::string& option) override;


  /*
   * Solve the Stokes system.
   * option = "undisturbed", compute the undisturbed field of flow without
   *particles
   * option = "disturbed",   compute the disturbed field of flow with particles
   */
  void solve(const std::string& option) override;


  /*
   * Return the SolverStokes
   */
  SolverStokes& solver_stokes() {
    return _solver_stokes;
  }

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
  void compute_point_velocity(const std::string& option,
                              std::vector<Real>& pv);


  /**
   * Compute the unperturbed velocity vector at the location of particles
   */
  std::vector<Real>compute_unperturbed_point_velocity();


  /*
   * Obtain the velocity vector on the i-th point
   */
  std::vector<Real>point_velocity(const std::vector<Real>& vel_beads,
                                  const std::size_t        i) const;


  /**
   * Local velocity of a fluid point in an unbounded space,
   * which is computed from Green's function
   * force_type: "regularized" or "smooth"
   */
  std::vector<Real>local_velocity_fluid(const Point      & p,
                                        const std::string& force_type,
                                        dof_id_type p_elem_id=-1) const;


  /**
   * Local velocity of a fluid point in an unbounded space,
   * which is computed from Green's function
   * force_type: "regularized" or "smooth"
   */
  std::vector<Real>local_velocity_fluid(const Elem        *elem,
                                        const Point      & p,
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
  void test_velocity_profile();

  /**
   * Test function. Test solutions on all nodes
   */
  void test_l2_norm();


    /*
     * Return the exact solution of Stokes eqn for any given
     * point \pt0 in an unbounded domain.
     */

  // const std::vector<Real> exact_solution(const Point& pt0) const;


  /*
   * Write nodal coordinates (x, y, z) and velocities (vx, vy, vz)
   * from system's solution vector to raw data file for step i.
   * If you want to write the total velocites, this function MUST ONLY
   * be called after you *add_local_solution* and *unperturbed* velocities
   * to the solution vector.
   */
  void write_fluid_velocity_data(const std::string& filename);

  /**
   * Couple Stokes System (HI or FD) with Poisson System by adding electrostatic
   * force to all points (Point Particles or surface nodes on Rigid Particles.) 
   * This electrostatic force includes contributions from both local and global
   * solution of the Poisson system.
   *
   * If module_np is true, this coupling will couple the Nernst-Planck system
   * to Poisson system first by adding the contribution to the charge density
   * from ion cloud, and then couple the total electrostatic force to Stokes
   * system.
   */  
  void couple_poisson();

  /**
   * Couple Stokes System with Nernst-Planck System to create a
   * convection-diffusion system without electrostatics, i.e, Stokes + Fick's
   * second law
   * This function will first initialize all NP system and relax them with
   * Poisson on but Fluid off for some time. Once all NP systems are relaxed,
   * this function will just solve all NP systems for one step to get c(t+dt)
   *
   * PseudoCode:
         void couple_np()
         {
            if first np sys is relax (i.e., all system are relaxed):
              for each np_sys:
              {
                 real_time += dt;
                 solve this np_sys with poisson&stokes --> gives c(t=t+dt)
              }

            else:
            {
                if first np sys is initialized (i.e., all system are initialized):
                {

                    real_time += dt

                    for each np_sys
                    {
                        np_sys.solve("diffusion"+(module_poisson) ? ("electrostatics")
                        : ("")) --> gives c(t+dt)

                        if real_time > relax_t_final:
                        {
                            set relax = true;
                            set real_time = 0.;
                        }
                    }

                    poisson.solve(with ion concentration); --> gives phi(t+dt)
                }

                else:
                    for each np_sys
                    {
                        initialize this np_sys, i.e, set concentration to be 0 and then
                        is_initialized = true
                    }

                this->couple_np();
            }
          }
   */
  void couple_np(unsigned int relax_step_id=0,
                 const unsigned int output_interval=10);
  
  
   /**
    * Write out Total solution of equation systems
    * Total solution = local + global + undisturbed (if exists)
    */
  void write_equation_systems(const unsigned int& o_step,
                              const Real& real_time);


  /**
   *
   */
  void update_solution_to_total() override;

  /**
   *
   */
  void resume_solution_to_global() override;

  
   /**
    * Save a pointer to undisturbed solution
    */
   UniquePtr<NumericVector<Number>> undisturbed_solution;

private:

  // Stokes solver
  SolverStokes _solver_stokes;

  // Assemble Stokes system
  AssembleStokes *assemble_stokes = nullptr;

  // Get a pointer to AnalyticalSolutionStokes
  AnalyticalSolutionStokes *analytical_solution = nullptr;

  // Get a pointer to GGEMStokes
  GGEMStokes *ggem_stokes = nullptr;

  // output precision (defined in input file, default is 6)
  int o_precision;
  
}; // end class
} // end namespace libMesh
