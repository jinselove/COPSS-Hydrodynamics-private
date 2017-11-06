#include "fix_rigid_sedimentation.h"

FixRigidSedimentation::FixRigidSedimentation(PMLinearImplicitSystem& pm_sys_)
:FixRigid(pm_sys_)
{
  force_type = "sedimentation";
  this -> initParams();
}

void FixRigidSedimentation::initParams()
{
  force_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("sedimentation");
  if(force_params.size()!=0){
    std::cout << std::endl << "********************Error message********************" << std::endl
              << "---------------> The force type 'sedimentation' does not require parameters from control file, instead, it needs sedimentation body force density from data file" << std::endl
              << "****************************************" << std::endl;
    libmesh_error();    
  }
}

void FixRigidSedimentation::print_fix()
{
  START_LOG("FixRigidSedimentation::print_fix()", "FixRigidSedimentation");
  std::cout <<"this is FixRigidSedimentation" << std::endl;
  STOP_LOG("FixRigidSedimentation::print_fix()", "FixRigidSedimentation");
}

void FixRigidSedimentation::compute()
{
 START_LOG("FixRigidSedimentation::compute()", "FixRigidSedimentation");
 for (int i=0; i < num_rigid_particles; i++){
  rigid_particles[i]->add_sedimentation_body_force_density();
 }
 STOP_LOG("FixRigidSedimentation::compute()", "FixRigidSedimentation");  
}