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

  void check_walls();

  void check_walls_pbcCount();

protected:
  // The elastic system for solids
  // Use pointer instead of const Ref. to avoid explicit initialization!
  ElasticitySystem* elastic_system;

  // partiticle mesh for finite size particle
  ParticleMesh<3>* particle_mesh;  

  // rigid particles
  std::vector<RigidParticle*> rigid_particles;

  // number of rigid particles
  std::size_t num_particles;

};  // end of class
  
} // end of namespace
