#include "fix_rigid.h"

FixRigid::FixRigid(PMLinearImplicitSystem& pm_sys_)
:Fix(pm_sys_)
{
  particle_type = pm_system->get_equation_systems().parameters.get<std::string>("particle_type");
}

//================================================================
void FixRigid::preSimulation()
{
  //
}


void FixRigid::check_walls()
{
  std::cout <<"Warning: check_walls() is not implemented yet for rigid particles" << std::endl;
}


void FixRigid::check_walls_pbcCount()
{
  // Fix this
}