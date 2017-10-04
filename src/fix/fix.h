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

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/parallel_object.h"
#include "libmesh/point.h"
#include "libmesh/id_types.h"
#include "libmesh/mesh.h"
#include "libmesh/serial_mesh.h"

#include "../pm_linear_implicit_system.h"
#include "../elasticity_system.h"
#include "../pm_toolbox.h"
#include "../point_mesh.h"
#include "../pm_periodic_boundary.h"
#include "fix_base.h"


namespace libMesh
{

  /*
   * The class is designed for computing the force field
   * of the system
   */
class FixBase;
  
class Fix
{ 
public:  

  //! Constructor for a system with elastic structures (finite size particles)
  Fix(PMLinearImplicitSystem& pm_sys,
             ElasticitySystem& el_sys);

  
  //! Constructor for a system with point particles
  Fix(PMLinearImplicitSystem& pm_sys);

  void initialization();

  //! Destructor
  virtual ~Fix(){};

  typedef std::string type_force_name;
  
  typedef std::vector<Real> type_force_parameter;
  
  typedef std::pair <type_force_name, type_force_parameter> type_force;

  
  /*! Prepare ForceField for run
  /*!
   * Check force field parameters when build force field classes.
   * This function is only evaluated once at the beginning of a simulation
   * check the following things:
   * 1) if this "force field" is supported;
   * 2) if this "force field" is applicable to this "particle type"
   * 3) if correct number of parameters are assigned to this "force field"
  */

  virtual void initParticleType() {}


  //! Check if num of parameters are correct
  virtual void initParams() {}

  /*! Compute force field attach it to particles
  /*!
   * User-defined function for reinitializing the force field
   * at each time step
   * The velocities of all beads are given (for calculating friction)
   * This function is the same to reinit_force_field() unless velocities are given, 
   * i.e., friction need to be calculated
   */
  virtual void compute() {}
 
  /* 
   * Correct the CHAIN position if it moves out of the wall or periodic boundary.
   */ 
  virtual void check_walls(){}

  /* 
   * Correct the CHAIN position if it moves out of the wall or periodic boundary, and
   * count how many times each bead has cross the boundary, which is used to
   * calculate its unfolded position.
   */ 
  virtual void check_walls_pbcCount() {}


  // debug
  virtual void print_fix(){}

  /*! use this for rigid particle fixes
   * check if particle is on periodic boundary, if so, rebuild particle mesh
   * this function is called once before all fix::compute(), no matter how many fixes we have 
   */
  virtual void check_pbc_pre_fix(){}

  /* use this for rigid particle fixes
   * restore particle mesh after fix::compute()
   */
  virtual void check_pbc_post_fix() {}

  /* use this for rigid particle fixes
   * attach nodal force on particle nodes once given nodal_force_density, which will be calculate in each sub-fix
   */
  virtual void attach_nodal() {}


  // the particle-mesh (linear implicit) system
  PMLinearImplicitSystem* pm_system;
  
  // pointMesh object
  PointMesh<3>* point_mesh;

  // system dimension
  unsigned int dim;

  // count total number of times beads escape the simulation domain.
  unsigned int out_domain_counter = 0;

  // force field
  std::string force_type;
  std::vector<Real> force_params;
  FixBase fix_base;

  // particle type
  std::string particle_type;

  // bead 
  Real bead_r;
  std::size_t num_points;
  std::vector<PointParticle*> point_particles;

  // kBT
  Real kBT;

  // const 
  const Real PI  = libMesh::pi;  //3.1415926


  // boundaries
  PMPeriodicBoundary* pbc;  //pmPeriodicBoundary object
  std::string wall_type;
  std::vector<Real> wall_params;
  Point box_min;
  Point box_max;
  Point box_len;
  std::vector<bool> periodic;
  std::vector<bool> inlet;  

};  // end of class
  
} // end of namespace
