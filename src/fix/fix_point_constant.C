#include "fix_point_constant.h"

FixPointConstant::FixPointConstant(PMLinearImplicitSystem& pm_sys_)
:FixPoint(pm_sys_)
{
  force_type = "p_constant";
  this -> initParams();
}

void FixPointConstant::initParams()
{
  START_LOG("FixPointConstant::initParams()", "FixPointConstant");
  force_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("p_constant");
  if(force_params.size()!=3){
    std::cout << std::endl << "********************Error message********************" << std::endl
              << "---------------> The force type 'p_constant' requires 3 parameter (fx, fy, fz) (dimensionless form)" << std::endl
              << "****************************************" << std::endl;
    libmesh_error();    
  }
  for (int i=0; i<3; i++){
      pforce(i) = force_params[i];
  }

  STOP_LOG("FixPointConstant::initParams()", "FixPointConstant");
}

void FixPointConstant::print_fix()
{
  START_LOG("FixPointConstant::print_fix()", "FixPointConstant");
  std::cout <<"this is FixPointConstant" << std::endl;
  STOP_LOG("FixPointConstant::print_fix()", "FixPointConstant");
}

void FixPointConstant::compute()
{
  START_LOG("FixPointConstant::compute()", "FixPointConstant");
  for(std::size_t p_id=0; p_id<num_points; ++p_id)
  {  
    point_particles[p_id]->add_particle_force(pforce);
  } 
  STOP_LOG("FixPointConstant::compute()", "FixPointConstant");
}
