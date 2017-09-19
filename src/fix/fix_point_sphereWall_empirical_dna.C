#include "fix_point_sphereWall_empirical_dna.h"

FixPointSphereWallEmpiricalDNA::FixPointSphereWallEmpiricalDNA(PMLinearImplicitSystem& pm_sys_)
:
 FixPoint(pm_sys_)
{
  this -> initParams();
}

void FixPointSphereWallEmpiricalDNA::initParams()
{
  wall_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>>("slit");
  for (int i = 0; i < dim; i++){
    if(box_min(i) != wall_params[2*i] or box_max(i) != wall_params[2*i+1]){
      std::cout << std::endl << "********************Error message********************" << std::endl
                << "for wall/empirical_dna (slit), wall_params has to be the same as box_min and box_max" << std::endl
                << "****************************************" << std::endl;
      libmesh_error();
    }
  }
  //Real c1 = ;
  //Real c2 = ;
  //c0 = ;
  //d0 = ;
}

void FixPointSphereWallEmpiricalDNA::print_fix()
{
  std::cout <<"this is FixPointSphereWallEmpiricalDNA" << std::endl;
}

void FixPointSphereWallEmpiricalDNA::compute()
{	
}

