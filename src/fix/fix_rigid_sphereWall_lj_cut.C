#include "fix_rigid_sphereWall_lj_cut.h"

FixRigidSphereWallLJCut::FixRigidSphereWallLJCut(PMLinearImplicitSystem& pm_sys_)
:FixRigid(pm_sys_)
{
  force_type = "wall/lj_cut";
  this -> initParams();
}

void FixRigidSphereWallLJCut::initParams()
{
START_LOG("FixRigidSphereWallLJCut::initParams()", "FixRigidSphereWallLJCut");
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
  cavity_radius = wall_params[0];
  epsilon = force_params[0];
  sigma = force_params[1];
  rCut = force_params[2];
STOP_LOG("FixRigidSphereWallLJCut::initParams()", "FixRigidSphereWallLJCut");
}

void FixRigidSphereWallLJCut::print_fix()
{
START_LOG("FixRigidSphereWallLJCut::print_fix()", "FixRigidSphereWallLJCut");
  std::cout <<"this is FixRigidSphereWallLJCut" << std::endl;
STOP_LOG("FixRigidSphereWallLJCut::print_fix()", "FixRigidSphereWallLJCut");
}


void FixRigidSphereWallLJCut::compute()
{
START_LOG("FixRigidSphereWallLJCut::compute()", "FixRigidSphereWallLJCut");

 /*------------------------------------------------------
  Loop over each rigid particle
  -------------------------------------------------------*/
  std::size_t node_start_id = 0;
  std::size_t point_id;
  for(std::size_t i=0; i<num_rigid_particles; ++i){
    const std::size_t n_nodes = rigid_particles[i]->num_mesh_nodes();
    std::vector<Point> pforce(n_nodes);
    MeshBase& p_mesh = rigid_particles[i]->mesh();
    MeshBase::node_iterator nd = p_mesh.active_nodes_begin();
    const MeshBase::node_iterator end_nd = p_mesh.active_nodes_end();
    // loop over each node on this particle
    for (; nd != end_nd; ++nd){
      Node *node = *nd;
      const dof_id_type node_id = node->id();
      // get point id 
      point_id = node_start_id + node_id;
      // get point particle
      const Point pti = point_particles[point_id]->point();
      const Real  pti_norm = pti.norm();
      const Point pti_unit = pti.unit();
      const Point r_i_sphereWall = pti_unit * (cavity_radius - pti_norm);
      if(r_i_sphereWall.norm() <= rCut) pforce[node_id] = fix_base.lj_force(r_i_sphereWall, epsilon, sigma);
    } // loop over all nodes
    // add node forces to this particle
    rigid_particles[i]->add_node_force(pforce);
    node_start_id += n_nodes;
  } // loop over all rigid particles  
  
STOP_LOG("FixRigidSphereWallLJCut::compute()", "FixRigidSphereWallLJCut");
}

