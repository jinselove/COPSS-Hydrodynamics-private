#include "fix_rigid_surface_constraint.h"

FixRigidSurfaceConstraint::FixRigidSurfaceConstraint(PMLinearImplicitSystem& pm_sys_)
:FixRigid(pm_sys_)
{
  force_type = "surface_constraint";
  this -> initParams();
}

void FixRigidSurfaceConstraint::initParams()
{
  START_LOG("FixRigidSurfaceConstraint::initParams()", "FixRigidSurfaceConstraint");
  force_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("surface_constraint");
  if (force_params.size() != 2) {
   std::cout << std::endl << "*******************Error message*********************" << std::endl
              << "'surface_constraint' needs two parameters: surface_spring_constant, center_spring_constant" << std::endl
              << "****************************************" << std::endl;
   libmesh_error();
  }
  else{
    surface_spring_constant = force_params[0];
    center_spring_constant = force_params[1];
  }
  STOP_LOG("FixRigidSurfaceConstraint::initParams()", "FixRigidSurfaceConstraint");
}

void FixRigidSurfaceConstraint::print_fix()
{
  START_LOG("FixRigidSurfaceConstraint::print_fix()", "FixRigidSurfaceConstraint");
  std::cout <<"this is FixRigidSurfaceConstraint" << std::endl;
  STOP_LOG("FixRigidSurfaceConstraint::print_fix()", "FixRigidSurfaceConstraint");
}

void FixRigidSurfaceConstraint::compute()
{
 START_LOG("FixRigidSurfaceConstraint::compute()", "FixRigidSurfaceConstraint");	
  std::size_t  point_start_id     = 0;
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Loop over each particle, and add the body force to each node.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  for(std::size_t i=0; i < num_rigid_particles; ++i)
  {
    // compute the rigid constraint force on each node
    std::vector<Point> constraint_force;
    this->compute_constraint_force(i,constraint_force);
    // add this constraint force to nodes
    rigid_particles[i]->add_node_force(constraint_force);
  } // end for i-loop
 STOP_LOG("FixRigidSurfaceConstraint::compute()", "FixRigidSurfaceConstraint");  
}

void FixRigidSurfaceConstraint::compute_constraint_force(const std::size_t& i,
																 std::vector<Point>& nodal_force)
{
	START_LOG("FixRigidSurfaceConstraint::compute_constraint_force()", "FixRigidSurfaceConstraint");
 /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   find the partilce center and number of tracking points on the surface.
   We need this center position to construct springs in the radial direction.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  MeshSpringNetwork* mesh_spring  = rigid_particles[i]->mesh_spring_network();
  const std::size_t n_nodes = rigid_particles[i]->num_mesh_nodes();
  const Point      p_centroid = rigid_particles[i]->get_centroid();
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   loop over each tracking point(node), and compute spring forces
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  nodal_force.resize(n_nodes);
  for(std::size_t j=0; j<n_nodes; ++j)
  {
    // 1. Get the coords and neighboring nodes of the j-th node
    std::vector< std::pair<std::size_t, Real> > node_neighbors;
    node_neighbors    = mesh_spring->nodes_neighbors(j);
    const Point& pt0  = rigid_particles[i]->mesh_point(j);
    // 2. Loop over each neighboring nodes and compute the SPRING forces
    Point spring_f(0.);
    for(std::size_t k=0; k<node_neighbors.size(); ++k)
    {
      const std::size_t& neigh_id = node_neighbors[k].first;
      const Real&        neigh_l0 = node_neighbors[k].second;
      const Point& ptk = rigid_particles[i]->mesh_point(neigh_id);
      const Point R_ij = particle_mesh->pm_periodic_boundary()->point_vector(pt0,ptk);
      Point sf  = fix_base.spring_force_lhs(R_ij,neigh_l0,surface_spring_constant);
      spring_f += sf;
    } // end for k-loop
    // 3. compute the Spring force between the node the the particle center
    const Real  lc0    = mesh_spring->node_center_equilibrium_dist(j);
    const Point Rc_ij  = particle_mesh->pm_periodic_boundary()->point_vector(pt0,p_centroid);
    Point sfc  = fix_base.spring_force_lhs(Rc_ij,lc0,center_spring_constant);
    
    //std::cout << "|Rc_ij| = " <<Rc_ij.size() <<"; lc0 = " <<lc0 << std::endl;
    spring_f += sfc;
    // 4. convert the point local id to the global id, and apply forces.
    nodal_force[j] = spring_f;
    /* - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - -
     * test: print out the spring force on each node.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   // if(pm_system->comm().rank()==0)
   // {
   //   printf("--->test in ForceField::rigid_constraint_force(), Particle %lu\n",i);
   //   printf("--->Spring force (surface+center=total) on node %lu is fs = (%f,%f,%f)\n",
   //          j,spring_f(0),spring_f(1),spring_f(2);
   //   printf("--->Spring force (center) on node %lu is               fs = (%f,%f,%f)\n",
   //          j,sfc(0),sfc(1),sfc(2);
   // }
    
  } // end for j-loop
  STOP_LOG("FixRigidSurfaceConstraint::compute_constraint_force()", "FixRigidSurfaceConstraint");
}
