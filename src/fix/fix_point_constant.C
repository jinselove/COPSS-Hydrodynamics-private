#include "fix_point_constant.h"

FixPointConstant::FixPointConstant(PMLinearImplicitSystem& pm_sys_)
  : FixPoint(pm_sys_)
{
  force_type = "p_constant";
  this->initParams();
}

void FixPointConstant::print_fix()
{
  START_LOG("FixPointConstant::print_fix()", "FixPointConstant");
  std::cout << "this is FixPointConstant" << std::endl;
  STOP_LOG("FixPointConstant::print_fix()", "FixPointConstant");
}

void FixPointConstant::compute()
{
  START_LOG("FixPointConstant::compute()", "FixPointConstant");

  for (std::size_t p_id = 0; p_id < num_points; ++p_id)
  {
    const Point pforce = point_particles[p_id]->constant_force();
    point_particles[p_id]->add_particle_force(pforce);
  }
  STOP_LOG("FixPointConstant::compute()", "FixPointConstant");
}
