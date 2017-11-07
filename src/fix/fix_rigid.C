#include "fix_rigid.h"

FixRigid::FixRigid(PMLinearImplicitSystem& pm_sys_,
          		   ElasticitySystem& el_sys_)
:Fix(pm_sys_),
 elastic_system(&el_sys_)
{
	this -> initParticleType();
	particle_mesh = pm_system->particle_mesh();
	rigid_particles = particle_mesh->particles();
}

FixRigid::FixRigid(PMLinearImplicitSystem& pm_sys_)
:Fix(pm_sys_)
{
	this -> initParticleType();
	particle_mesh = pm_system -> particle_mesh();
	num_rigid_particles = particle_mesh -> num_particles();
  rigid_particles = particle_mesh->particles();
}

//================================================================
void FixRigid::initParticleType()
{
  particle_type = pm_system->get_equation_systems().parameters.get<std::string>("particle_type");
  if(particle_type != "rigid_particle") {
	std::cout << std::endl << "*******************Error message*********************" << std::endl
              << "The force field : "<< force_type << " is only for rigid_particle, i.e.,"
              << "not applicable to particle type: " << particle_type << std::endl
              << "****************************************" << std::endl;
	libmesh_error();
  }   
}

//===================================================================
void FixRigid::check_walls(){
  START_LOG("FixRigid::check_walls()", "FixRigid");
  wall_type = pm_system->get_equation_systems().parameters.get<std::string>("wall_type");
   // for slit walls
  if(wall_type == "slit"){
    for(std::size_t i=0; i < num_rigid_particles; ++i)
    {
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        1) Check if this particle is on the periodic boundary. If so, move the nodes
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
      const bool on_pb = rigid_particles[i]->on_the_periodic_boundary();
      if( on_pb )
      {
        rigid_particles[i]->rebuild_periodic_mesh();
        std::cout<<" >>>>> Notice: particle " <<i << "is on the period boundary, particle mesh is rebuilt" <<std::endl;
      }
      // 2) check non-slip walls 
      MeshBase& p_mesh = rigid_particles[i]->mesh();
      MeshBase::node_iterator       nd     = p_mesh.active_nodes_begin();
      const MeshBase::node_iterator end_nd = p_mesh.active_nodes_end();
      // loops over each node
      for (; nd!=end_nd; ++nd){
        Node* node = *nd;
        for (std::size_t i=0; i<dim; i++)
        {
          if (!periodic[i]){
            if( (*node)(i) < box_min(i)) (*node)(i) = 2.*box_min(i) - (*node)(i);  
            if( (*node)(i) >= box_max(i)) (*node)(i) = 2.*box_max(i) - (*node)(i);  
          }
        }
      } // end for loop over nodes
    } // end loop over particles
  }//end if wall_type
  else{
    Real cavity_radius = wall_params[0];
    for(std::size_t i=0; i < num_rigid_particles; ++i)
    {
      // check non-slip walls 
      MeshBase& p_mesh = rigid_particles[i]->mesh();
      MeshBase::node_iterator       nd     = p_mesh.active_nodes_begin();
      const MeshBase::node_iterator end_nd = p_mesh.active_nodes_end();
      // loops over each node
      for (; nd!=end_nd; ++nd){
        Node* node = *nd;
      // Move particle back into spherical cavity
        if( node->size() > cavity_radius ) (*node) = node->unit()*(2.*cavity_radius) - (*node);
      } // end for loop over nodes
    } // end loop over particles
  } // end else
  STOP_LOG("FixRigid::check_walls()", "FixRigid");
} // end check walls



// fix this: when rigid particle crosses non-slip walls, we need to find the farest node
// outside the wall and then move the particle back.
//    

// void FixRigid::check_walls()
// {
//   wall_type = pm_system->get_equation_systems().parameters.get<std::string>("wall_type");
//   if(wall_type = "slit"){
//     for(std::size_t i=0; i < num_rigid_particles; ++i)
//     {
//       /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         1) Check if this particle is on the periodic boundary. If so, move the nodes
//        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
//       const bool on_pb = rigid_particles[i]->on_the_periodic_boundary();
//       if( on_pb )
//       {
//         rigid_particles[i]->rebuild_periodic_mesh();
//         std::cout<<" >>>>> Notice: particle " <<i << "is on the period boundary, particle mesh is rebuilt" <<std::endl;
//       }
//       /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         2) Check if this particle penetrates non-slip wall, If so, move center of mass and all nodes
//        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
//       SerialMesh& p_mesh = rigid_particles[i]->mesh();
//       MeshBase::node_iterator       nd     = p_mesh.active_nodes_begin();
//       const MeshBase::node_iterator end_nd = p_mesh.active_nodes_end();
//       for (; nd != end_nd; ++nd){
//         Node* node = *nd;
//         for (std::size_t i = )
//       }
//     }
//   }
//   else{

//   }
// }

//===================================================================
void FixRigid::check_walls_pbcCount()
{
  // Fix this
}

//===================================================================
void FixRigid::check_pbc_pre_fix()
{
  START_LOG("FixRigid::check_pbc_pre_fix()", "FixRigid");
  for(std::size_t i=0; i < num_rigid_particles; ++i)
  {
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Check if this particle is on the periodic boundary. If so, move the nodes
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    const bool on_pb = rigid_particles[i]->on_the_periodic_boundary();
    if( on_pb )
    {
      if(pm_system->comm().rank()==0){
        printf("--->check_pbc_pre_fix() before all fix::compute() \n");
        printf("         The particle %lu is on the periodic boundary!\n",i);
      }
      rigid_particles[i]->rebuild_periodic_mesh();
    }  
  }
  STOP_LOG("FixRigid::check_pbc_pre_fix()", "FixRigid");
}

//===================================================================
void FixRigid::check_pbc_post_fix()
{
  START_LOG("FixRigid::check_pbc_post_fix()", "FixRigid");
  for(std::size_t i=0; i < num_rigid_particles; ++i)
  {
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Check if this particle is on the periodic boundary. If so, move the nodes
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    const bool& on_pb = rigid_particles[i]->if_on_the_periodic_boundary();
    if( on_pb )
    {
      if(pm_system->comm().rank()==0){
        printf("--->restore_periodic_mesh() after all fix::compute()\n");
      }
      rigid_particles[i]->restore_periodic_mesh();
    }  
  }  
  STOP_LOG("FixRigid::check_pbc_post_fix()", "FixRigid");
}

//===================================================================
void FixRigid::attach_nodal()
{
  START_LOG("FixRigid::attach_nodal()", "FixRigid");
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over each particle, and add the body force to each node.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  std::size_t point_start_id = 0;
  for(std::size_t i=0; i < num_rigid_particles; ++i)
  {
    std::vector<Point> nodal_force;
    rigid_particles[i]->get_node_force(nodal_force); 
    /*
     * Loop over each node through node iterator
     * Compute the gravitational force vector on each node
     */
    const std::size_t n_nodes = particle_mesh->particles()[i]->num_mesh_nodes();
    MeshBase& p_mesh = particle_mesh->particles()[i]->mesh();
    MeshBase::node_iterator       nd     = p_mesh.active_nodes_begin();
    const MeshBase::node_iterator end_nd = p_mesh.active_nodes_end();
    for ( ; nd != end_nd; ++nd)
    {
      // Store a pointer to the element we are currently working on.
      Node* node = *nd;
      const dof_id_type node_id = node->id();
      // get the dof numbers at this node (only for force vector)
      std::vector<Real> gforce(dim,0.);
      for(std::size_t k=0; k<dim; ++k){
        gforce[k] = nodal_force[node_id](k);
      } // end k-loop
      point_mesh->particles()[point_start_id+node_id]->add_particle_force(gforce);
    //  ------------------ TEST: print out the nodal force ----------------------
      // if(_pm_system->comm().rank()==0){
      //   printf("--->TEST:reinit_force_field() gforce = (%f,%f,%f)\n",gforce[0],gforce[1],gforce[2]);
      //   printf("         nodal force = (%f,%f,%f)\n",
      //          nodal_force[node_id](0),nodal_force[node_id](1),nodal_force[node_id](2));
      //   printf("         rigid_nodal_force = (%f,%f,%f)\n",
      //          rigid_nodal_force[node_id](0),rigid_nodal_force[node_id](1),rigid_nodal_force[node_id](2));
      // }
      // -------------------------------------------------------------------------    
    } // end for nd-loop
    // Point test_centroid_force;
    // for (int test_i=0; test_i<point_mesh->num_particles(); test_i++){
    //   for(int test_dim = 0; test_dim<3; test_dim++) test_centroid_force(test_dim) += point_mesh->particles()[test_i]->particle_force()[test_dim];
    // }
    // printf("Debug: for particle %d, test_centroid_force = (%f, %f, %f)\n",i, test_centroid_force(0), test_centroid_force(1), test_centroid_force(2) );
     /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     update the point start id
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    point_start_id += n_nodes;
  } // end for i-loop
  STOP_LOG("FixRigid::attach_nodal()", "FixRigid");
}

