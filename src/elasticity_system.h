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
// #include "point_mesh.h"
#include "particle_mesh.h"
#include "mesh_spring_network.h"


// include PETSc KSP solver
EXTERN_C_FOR_PETSC_BEGIN
#include <petscksp.h>
EXTERN_C_FOR_PETSC_END


// C++ Includes   -------------------------------------
#include <stdio.h>
#include <cstddef>


namespace libMesh {
/*
 * The ElasticitySystem is designed to solve the elasticity
 * equation associated with immersed solid bodies using FEM.
 */


class ElasticitySystem : public LinearImplicitSystem {
public:

  /**
   * Constructor.
   */
  ElasticitySystem(EquationSystems  & es,
                   const std::string& name,
                   const unsigned int number); // number of systems

  /**
   * Destructor.
   */
  virtual ~ElasticitySystem();


  /**
   * The type of system.
   */
  typedef ElasticitySystem sys_type;


  /**
   * The type of the parent.
   */
  typedef LinearImplicitSystem Parent;


  /**
   * @returns a clever pointer to the system.
   */
  sys_type& system() {
    return *this;
  }

  /**
   * Clear all the data structures associated with the system.
   */
  virtual void clear();


  /**
   * Create KSP solver
   */
  void init_ksp_solver();


  /**
   * Destroy KSP solver
   */
  void destroy_ksp_solver();


  /**
   * return KSP solver
   */
  KSP& ksp_solver() {
    return _ksp;
  }

  /**
   * construct the nodal force vector due to gravity
   * f = (fx,fy,fz) is the force density.
   *
   * It is an area density for surface mesh, and volume density of volume mesh!
   */
  void build_nodal_force_gravity(const std::vector<Real>& f);


  /**
   * Return the nodal force vector
   */
  std::vector<Real>nodal_force() const {
    return _nodal_force;
  }

  /*
   * Attach ParticleMesh
   */
  void attach_particle_mesh(ParticleMesh<3> *pm) {
    _particle_mesh = pm;
  }

  /*
   * return editable ParticleMesh ptr and const ParticleMesh ptr
   */
  ParticleMesh<3>* particle_mesh() {
    return _particle_mesh;
  }

  ParticleMesh<3>* particle_mesh() const {
    return _particle_mesh;
  }

  /*
   * Attach MeshSpringNetwork
   */
  void attach_mesh_spring_network(MeshSpringNetwork *msn) {
    _mesh_spring_network = msn;
  }

  /*
   * return editable MeshSpringNetwork ptr and const MeshSpringNetwork ptr
   */
  MeshSpringNetwork* mesh_spring_network() {
    return _mesh_spring_network;
  }

  MeshSpringNetwork* mesh_spring_network() const {
    return _mesh_spring_network;
  }

  /*
   * Compute the min and max mesh size:
   */
  std::vector<Real>mesh_size() const;

private:

  // particle mesh pointer
  ParticleMesh<3> *_particle_mesh;


  // (surface or volume) mesh spring network pointer
  MeshSpringNetwork *_mesh_spring_network;


  // label whether the system matrix is assembled
  // if is assembled, there is no need to assemble it every time step
  bool _matrix_assembled;


  // Krylov Subspace
  KSP _ksp;


  // nodal force vector, which has a copy on each processor
  std::vector<Real>_nodal_force;
}; // end of class
} // end of namespace
