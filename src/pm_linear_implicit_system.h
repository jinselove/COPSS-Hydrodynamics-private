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
#include "libmesh/libmesh_config.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/linear_implicit_system.h"

// Local Includes   -----------------------------------
#include "point_mesh.h"
#include "particle_mesh.h"
#include "elasticity_system.h"
#include "assemble_stokes.h"
#include "stokes_solver.h"

// C++ Includes   -------------------------------------
#include <stdio.h>
#include <cstddef>


namespace libMesh
{


/*
 * The PMLinearImplicitSystem is designed to
 * solve the stokes equation with regularized Gaussian
 * point forces due to particles using mixed FEM.
 */

class Fix;


class PMLinearImplicitSystem : public LinearImplicitSystem
{
public:

  /**
   * Constructor.
   */
  PMLinearImplicitSystem (EquationSystems& es,
                          const std::string& name,
                          const unsigned int number); // number of systems


  /**
   * Destructor.
   */
  virtual ~PMLinearImplicitSystem ();



  /**
   * The type of system.
   */
  typedef PMLinearImplicitSystem sys_type;


  /**
   * The type of the parent.
   */
  typedef LinearImplicitSystem Parent;


  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }



  /**
   * Clear all the data structures associated with the system.
   */
  virtual void clear ();



  /**
   * Assemble the system matrix.
   * option == "undisturbed" : undisturbed velocity field (without particles)
   * option == "disturbed"   : disturbed velocity field   (with particles)
   *
   * However, the global matrix here should be the same for both cases.
   * option only influences the RHS vector???
   */
  void assemble_matrix (const std::string& system_name,
                        const std::string& option);


  /**
   * Assemble the system rhs.
   */
  void assemble_rhs (const std::string& system_name,
                     const std::string& option);



  /*   - - - - - Recommended - - - - - -
   * Solve the Stokes system.
   * option = "undisturbed", compute the undisturbed field of flow without particles
   * option = "disturbed",   compute the disturbed field of flow with particles
   * re_init = true => re-assemble the matrix and reinit the KSP solver.
   */
  void solve_stokes (const std::string& option,
                     const bool& re_init);


  /*
   * Return the StokesSolver
   */
  StokesSolver& stokes_solver() { return _stokes_solver;  }


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


  /*
   * Attach PointMesh
   */
  void attach_point_mesh(PointMesh<3>* pm){ _point_mesh = pm;  };


  /*
   * return editable PointMesh ptr and const PointMesh ptr
   */
  PointMesh<3>* point_mesh(){  return _point_mesh;  };
  PointMesh<3>* point_mesh() const {  return _point_mesh;  };


  /*
   * Attach ParticleMesh
   */
  void attach_particle_mesh(ParticleMesh<3>* pm){ _particle_mesh = pm;  };


  /*
   * return editable ParticleMesh ptr and const ParticleMesh ptr
   */
  ParticleMesh<3>* particle_mesh(){  return _particle_mesh;  };
  ParticleMesh<3>* particle_mesh() const {  return _particle_mesh;  };


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



  /*
   * Attach the ForceField
   */
  void attach_fixes(std::vector<Fix*> fixes){ _fixes = fixes;  };
  std::vector<Fix*> fixes() { return _fixes; };
  std::vector<Fix*> fixes() const {  return _fixes; };


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
   * Compute the L2-error in an unbounded domain
   * This function will change the system solution vector by add local solution.
   */
  void test_l2_norm(bool& neighbor_list_update_flag);


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
   * write out single particle(point particle) info, including:
   * timestep, time, particle_xyz,  particle_velocity
   */
  void write_out_single_particle(const Point& coords,
                                 const std::vector<Real>& vel,
                                 const int i_step,
                                 const Real time) const;


  /*
   * Write out particle(point) coordinates represented by a PETSc vector
   */
  void write_out_point_coordinate(Vec* ROUT,
                                  const std::size_t istep,
                                  const Real& time,
                                  const std::string& f_name,
                                  const std::string& openmode) const;


  /*
   * Write out particle(point) coordinates
   * The output format is Comma Separated Variable(CSV) for ParaView.
   */
  void write_point_csv(const std::string& filename,
                       const std::vector<Real>& pv,
                       const bool write_velocity) const;

  PetscErrorCode write_point_csv(const std::string& filename,
                                 Vec * petsc_vector,
                                 const bool write_velocity) const;

  /*
   * Write out equation systems of Stokes. This requires combining the
   * local and global solution, and update the solution vectors.
   * output_format=="EXODUS" ; "VTK"; "GMV"
   */
  void write_equation_systems(const std::size_t time_step,
                              const std::string& output_filename,
                              const std::string& output_format);


  /*
   * Add the local solution to the global solution
   */
  void add_local_solution();

 
  /*
   * Write nodal coordinates (x, y, z) and velocities (vx, vy, vz)
   * from system's solution vector to raw data file for step i.
   * If you want to write the total velocites, this function MUST ONLY
   * be called after you *add_local_solution* and *unperturbed* velocities
   * to the solution vector.
   */
  void write_fluid_velocity_data(const std::string& filename);

private:

  // particle mesh pointer
  PointMesh<3>* _point_mesh;

  // particle mesh pointer
  ParticleMesh<3>* _particle_mesh;

  // force field for particle/points
  std::vector<Fix*> _fixes;

  // Stokes solver
  StokesSolver _stokes_solver;

  // Assemble Stokes system
  AssembleStokes* _assemble_stokes;

  // the type of particles in this PMSystem: "point_particle" or "rigid_particle"
//  std::string _particle_type;
};  // end of class PMLinearImplicitSystem



} // end namespace libMesh
