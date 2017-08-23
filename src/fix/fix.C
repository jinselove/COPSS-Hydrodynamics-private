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




#include "libmesh/libmesh_logging.h"
#include "libmesh/equation_systems.h"
#include "fix.h"


using namespace libMesh;


// ======================================================================
Fix::Fix(PMLinearImplicitSystem& pm_sys_,
                       ElasticitySystem& el_sys_)
: pm_system(&pm_sys_), elastic_system(&el_sys_)
{
  this -> initialization();
}



// ======================================================================
Fix::Fix(PMLinearImplicitSystem& pm_sys_)
: pm_system(&pm_sys_)
{
  this -> initialization();
}

void Fix::attach_system(PMLinearImplicitSystem& pm_system_)
{
  pm_system = &pm_system_;
  this -> initialization();
}

void Fix::initialization()
{
  point_mesh = pm_system->point_mesh();

  dim      = pm_system->get_mesh().mesh_dimension();
  
  particle_type = pm_system->get_equation_systems().parameters.get<std::string>("particle_type");  

  kBT      = pm_system->get_equation_systems().parameters.get<Real>("kBT");

  wall_type = pm_system->get_equation_systems().parameters.get<std::string>("wall_type");

  wall_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>>(wall_type);

  pbc = point_mesh->pm_periodic_boundary();

  box_min = pbc->box_min();
  
  box_max = pbc->box_max();

  box_len = pbc->box_length();
  
  periodic = pbc->periodic_direction();
  
  inlet = pbc->inlet_direction();

  num_points = point_mesh -> num_particles();
}