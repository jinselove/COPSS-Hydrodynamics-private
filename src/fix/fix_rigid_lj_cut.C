#include "fix_rigid_lj_cut.h"

FixRigidLJCut::FixRigidLJCut(PMLinearImplicitSystem& pm_sys_)
  : FixRigid(pm_sys_)
{
  force_type = "lj_cut";
  this->initParams();
}

void FixRigidLJCut::initParams()
{
  START_LOG("FixRigidLJCut::initParams()", "FixRigidLJCut");
  force_params =
    pm_system->get_equation_systems().parameters.get<std::vector<Real> >("lj_cut");

  if (force_params.size() != 3) {
    std::cout << std::endl <<
      "********************Error message********************" << std::endl
              <<
      "---------------> The force type 'lj_cut' requires 3 parameter (epsilon, sigma, rcut) (dimensionless form)"
              << std::endl
              << "****************************************" << std::endl;
    libmesh_error();
  }
  epsilon = force_params[0];
  sigma   = force_params[1];
  rCut    = force_params[2];
  STOP_LOG("FixRigidLJCut::initParams()", "FixRigidLJCut");
}

void FixRigidLJCut::print_fix()
{
  START_LOG("FixRigidLJCut::print_fix()", "FixRigidLJCut");
  std::cout << "this is FixRigidLJCut" << std::endl;
  STOP_LOG("FixRigidLJCut::print_fix()", "FixRigidLJCut");
}

void FixRigidLJCut::compute()
{
  START_LOG("FixRigidLJCut::compute()", "FixRigidLJCut");

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

      // get neighbor list of this point
      const std::vector<Point>& neighbor_vector =
        point_particles[point_id]->neighbor_vector();
      const std::vector<dof_id_type>& n_list =
        point_particles[point_id]->neighbor_list();

      // Loop over each neigbhor
      for (std::size_t j = 0; j < neighbor_vector.size(); ++j)
      {
        const std::size_t n_id = n_list[j];

        // if n_id is not on the same rigid particle and n_id != point_id
        //   printf("point_id = %d, n_id = %d, node_start_id = %d,
        // node_start_id+n_nodes = %d\n",point_id, n_id, node_start_id,
        // node_start_id+n_nodes );
        if (n_id<node_start_id or n_id>node_start_id + n_nodes) // make sure
                                                                // this bead and
                                                                // the
                                                                // neighboring
                                                                // bead are not
                                                                // the same
                                                                // bead.
        {
          // printf("something is wrong, point_id = %d, n_id = %d, node_start_id
          // = %d, node_start_id+n_nodes = %d\n",point_id, n_id, node_start_id,
          // node_start_id+n_nodes );
          if (neighbor_vector[j].norm() <= rCut) {
            pforce[node_id] += fix_base.lj_force(neighbor_vector[j],
                                                 epsilon,
                                                 sigma);
          }
        } // end if
      }   // end loop over neighbors
    }     // loop over all nodes
    // add node forces to this particle
    rigid_particles[i]->add_node_force(pforce);
    node_start_id += n_nodes;
  } // loop over all rigid particles
  STOP_LOG("FixRigidLJCut::compute()", "FixRigidLJCut");
}
