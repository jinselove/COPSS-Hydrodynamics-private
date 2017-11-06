#include "fix_point_sphereWall_lj_cut.h"

FixPointSphereWallLJCut::FixPointSphereWallLJCut(PMLinearImplicitSystem& pm_sys_)
:
 FixPoint(pm_sys_)
{
  force_type = "wall/lj_cut";
  this -> initParams();
}

void FixPointSphereWallLJCut::initParams()
{
START_LOG("FixPointSphereWallLJCut::initParams()", "FixPointSphereWallLJCut");
  force_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("wall/lj_cut");
  wall_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>>("sphere");
  if(force_params.size()!=3){
    std::cout << std::endl << "********************Error message********************" << std::endl
              << "---------------> The force type 'wall/lj_cut' requires 3 parameter (epsilon, sigma, rcut) (dimensionless form)" << std::endl
              << "****************************************" << std::endl;
    libmesh_error();    
  }
  if(wall_params.size()!=1){
    std::cout << std::endl << "********************Error message********************" << std::endl
              << "---------------> The force type 'wall/lj_cut' (sphere) requires 1 wall parameter, i.e., radius (dimensionless form)" << std::endl
              << "****************************************" << std::endl;
    libmesh_error();    
  }
  epsilon = force_params[0];
  sigma = force_params[1];
  rCut = force_params[2];
STOP_LOG("FixPointSphereWallLJCut::initParams()", "FixPointSphereWallLJCut");
}

void FixPointSphereWallLJCut::print_fix()
{
START_LOG("FixPointSphereWallLJCut::print_fix()", "FixPointSphereWallLJCut");
  std::cout <<"this is FixPointSphereWallLJCut" << std::endl;
STOP_LOG("FixPointSphereWallLJCut::print_fix()", "FixPointSphereWallLJCut");
}


void FixPointSphereWallLJCut::compute()
{
START_LOG("FixPointSphereWallLJCut::compute()", "FixPointSphereWallLJCut");
  for(std::size_t i=0; i<num_points; ++i)
  {
    std::vector<Real> pforce(dim);
    const Point pti = point_particles[i] -> point();
    // Retrieve wall_info
    Real cavity_radius = wall_params[0];
    const Real  pti_norm = pti.norm();
    const Point pti_unit = pti.unit();
    // vector point from particle i to the sphere wall
    const Point r_i_sphereWall = pti_unit * (cavity_radius-pti_norm); 
    // bead-wall interaction force
    if(r_i_sphereWall.norm() <= rCut) pforce = fix_base.lj_force(r_i_sphereWall, epsilon, sigma);
    // attach this force to particle i
    point_particles[i]->add_particle_force(pforce);    
  } // end for i-loop 
STOP_LOG("FixPointSphereWallLJCut::compute()", "FixPointSphereWallLJCut");
}

