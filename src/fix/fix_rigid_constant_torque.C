#include "fix_rigid_constant_torque.h"

FixRigidConstantTorque::FixRigidConstantTorque(PMLinearImplicitSystem& pm_sys_)
  : FixRigid(pm_sys_)
{
  force_type = "constant_torque";
  this->initParams();
}

void FixRigidConstantTorque::initParams()
{
  force_params =
    pm_system->get_equation_systems().parameters.get<std::vector<Real> >(
      "constant_torque");

  if (force_params.size() != 1) {
    std::cout << std::endl <<
    "********************Error message********************" << std::endl
    << "---------------> The force type 'constant_torque' needs one"
       "parameters 'torque_z' from the control file"
    << std::endl
    << "****************************************" << std::endl;
    libmesh_error();
  }
}

void FixRigidConstantTorque::print_fix()
{
  START_LOG("FixRigidConstantTorque::print_fix()", "FixRigidConstantTorque");
  std::cout << "this is FixRigidConstantTorque" << std::endl;
  STOP_LOG("FixRigidConstantTorque::print_fix()", "FixRigidConstantTorque");
}

void FixRigidConstantTorque::compute()
{
  START_LOG("FixRigidConstantTorque::compute()", "FixRigidConstantTorque");

  /*------------------------------------------------------
     Loop over each rigid particle
     -------------------------------------------------------*/
  std::size_t node_start_id = 0;
  std::size_t point_id;
  std::vector<std::size_t> node_list;
  node_list.push_back(0); node_list.push_back(7);

  for (std::size_t i = 0; i < num_rigid_particles; ++i) {
    const std::size_t  n_nodes = rigid_particles[i]->num_mesh_nodes();
    std::vector<Point> pforce(n_nodes);
    MeshBase& p_mesh = rigid_particles[i]->mesh();
    // only apply force on the nodes in the node list
    for (std::size_t j=0; j<node_list.size(); j++){
      Point r = p_mesh.point(node_list[j]) - rigid_particles[i]->get_centroid();
      Real r_norm_sq = r.norm_sq();
      if (std::abs(r(2))>1.e-2) {
        std::cout << "Error: the z-coordinate of node 0 is far from zero.";
        libmesh_error();
      }
      pforce[node_list[j]] = Point(-r(1), r(0), 0.) * (force_params[0] / r_norm_sq
        / node_list.size());
//      std::cout<<"pforce["<<node_list[j]<<")="<<pforce[node_list[j]]<<std::endl;
    }
    // add node forces to this particle
    rigid_particles[i]->add_node_force(pforce);
  } // loop over all rigid particles

  STOP_LOG("FixRigidConstantTorque::compute()", "FixRigidConstantTorque");
}
