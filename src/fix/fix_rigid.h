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

#include "fix.h"
#include "../rigid_particle.h"
#include "../elasticity_system.h"

namespace libMesh
{

  /*
   * The class is designed for computing the force field
   * of the system
   */
  
  
class FixRigid: public Fix
{
public:

  //! Constructor for a system with elastic structures (finite size particles)
  FixRigid(PMLinearImplicitSystem& pm_sys,
             ElasticitySystem& el_sys);

  //! Constructor for a system with point particles
  FixRigid(PMLinearImplicitSystem& pm_sys);

  virtual ~FixRigid(){}

  /*! Prepare for running
   *
   * Check if particle_type == "rigid_particle"
  */
  void initParticleType();

  /*
   * 1) check if a rigid particle is on a periodic boundary, if so, rebuild the particlemesh (fix me, is this necessary)
   * 2) check if a particle is on or crosses a non-slip wall, if so, move the nodes that crosses the wall back 
   */
  void check_walls();

  /*
   * Not implemented yet
   */
  void check_walls_pbcCount();

  /*
   * check if particle is on periodic boundary, if so, rebuild particle mesh
   * this function is called once before all fix::compute(), no matter how many fixes we have 
   */
  void check_pbc_pre_fix();

  /*
   * restore particle mesh after fix::compute()
   */
  void check_pbc_post_fix();

  /*
   * attach nodal force on particle nodes once given nodal_force_density, which will be calculate in each sub-fix
   */
  void attach_nodal();

protected:
  // The elastic system for solids
  // Use pointer instead of const Ref. to avoid explicit initialization!
  ElasticitySystem* elastic_system;

  // partiticle mesh for finite size particle
  ParticleMesh<3>* particle_mesh;  

  // rigid particles
  std::vector<RigidParticle*> rigid_particles;

  // number of rigid particles
  std::size_t num_rigid_particles;
  
};  // end of class
  
} // end of namespace
