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
	num_particles = particle_mesh -> num_particles();
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


void FixRigid::check_walls()
{
  std::cout <<"Warning: check_walls() is not implemented yet for rigid particles" << std::endl;
}


void FixRigid::check_walls_pbcCount()
{
  // Fix this
}
