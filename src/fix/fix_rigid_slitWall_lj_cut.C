#include "fix_rigid_slitWall_lj_cut.h"

FixRigidSlitWallLJCut::FixRigidSlitWallLJCut(PMLinearImplicitSystem& pm_sys_)
  :
  FixRigid(pm_sys_)
{
  force_type = "wall/lj_cut";
  this->initParams();
}

void FixRigidSlitWallLJCut::initParams()
{
  START_LOG("FixRigidSlitWallLJCut::initParams()", "FixRigidSlitWallLJCut");
  force_params =
    pm_system->get_equation_systems().parameters.get<std::vector<Real> >(
      "wall/lj_cut");
  wall_params =
    pm_system->get_equation_systems().parameters.get<std::vector<Real> >("slit");

  if (force_params.size() != 3) {
    std::cout << std::endl <<
      "********************Error message********************" << std::endl
              <<
      "---------------> The force type 'wall/lj_cut' requires 3 parameter (epsilon, sigma, rcut) (dimensionless form)"
              << std::endl
              << "****************************************" << std::endl;
    libmesh_error();
  }

  for (int i = 0; i < dim; i++) {
    if (box_min(i) != wall_params[2 * i] or box_max(i) !=
        wall_params[2 * i + 1]) {
      std::cout << std::endl <<
        "********************Error message********************" << std::endl
                <<
        "for wall/lj_cut (slit), wall_params has to be the same as box_min and box_max"
                << std::endl
                << "****************************************" << std::endl;
      libmesh_error();
    }
  }
  epsilon = force_params[0];
  sigma   = force_params[1];
  rCut    = force_params[2];
  STOP_LOG("FixRigidSlitWallLJCut::initParams()", "FixRigidSlitWallLJCut");
}

void FixRigidSlitWallLJCut::print_fix()
{
  START_LOG("FixRigidSlitWallLJCut::print_fix()", "FixRigidSlitWallLJCut");
  std::cout << "this is FixRigidSlitWallLJCut" << std::endl;
  STOP_LOG("FixRigidSlitWallLJCut::print_fix()", "FixRigidSlitWallLJCut");
}

void FixRigidSlitWallLJCut::compute()
{
  START_LOG("FixRigidSlitWallLJCut::compute()", "FixRigidSlitWallLJCut");

  /*------------------------------------------------------
     Loop over each rigid particle
     -------------------------------------------------------*/
  std::size_t node_start_id = 0;
  std::size_t point_id;

  for (std::size_t i = 0; i < num_rigid_particles; ++i) {
    const std::size_t  n_nodes = rigid_particles[i]->num_mesh_nodes();
    std::vector<Point> pforce(n_nodes);
    MeshBase& p_mesh                     = rigid_particles[i]->mesh();
    MeshBase::node_iterator nd           = p_mesh.active_nodes_begin();
    const MeshBase::node_iterator end_nd = p_mesh.active_nodes_end();

    // loop over each node on this particle
    for (; nd != end_nd; ++nd) {
      Node *node                = *nd;
      const dof_id_type node_id = node->id();

      // get point id
      point_id = node_start_id + node_id;

      // get point particle
      const Point pti = point_particles[point_id]->point();

      // loop over all dimensions
      for (std::size_t j = 0; j < dim; ++j) {
        Point r_i_lo, r_i_hi;

        if (periodic[j] == false and inlet[j] == false) {
          r_i_lo(j) = box_min(j) - pti(j);
          r_i_hi(j) = box_max(j) - pti(j);

          if (r_i_lo.norm() <= rCut) pforce[node_id](j) += fix_base.lj_force(
              r_i_lo,
              epsilon,
              sigma)(j);

          if (r_i_hi.norm() <= rCut) pforce[node_id](j) += fix_base.lj_force(
              r_i_hi,
              epsilon,
              sigma)(j);
        }
      }
    } // loop over all nodes
    // add node forces to this particle
    rigid_particles[i]->add_node_force(pforce);
    node_start_id += n_nodes;
  } // loop over all rigid particles

  STOP_LOG("FixRigidSlitWallLJCut::compute()", "FixRigidSlitWallLJCut");
}
