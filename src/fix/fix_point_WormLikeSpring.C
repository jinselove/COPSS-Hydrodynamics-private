#include "fix_point_WormLikeSpring.h"

FixPointWormLikeSpring::FixPointWormLikeSpring(PMLinearImplicitSystem& pm_sys_)
:FixPoint(pm_sys_)
{
  this -> initPointParticleType();
  this -> initParams();
}

void FixPointWormLikeSpring::initPointParticleType()
{
  point_particle_model = pm_system->get_equation_systems().parameters.get<std::string>("point_particle_model");
  if(point_particle_model != "polymer_chain") {
  std::cout << std::endl << "*******************Error message*********************" << std::endl
              << "The force field : "<< force_type << " is only for polymer_chain, i.e.,"
              << "not applicable to point_particle_model: " << point_particle_model << std::endl
              << "****************************************" << std::endl;
  libmesh_error();
  }   
}

void FixPointWormLikeSpring::initParams()
{
  n_bonds = point_mesh -> num_bonds();
  polymer_chain = point_mesh->polymer_chain();
  // Worm-like-spring force parameters
  Real Nks      = pm_system->get_equation_systems().parameters.get<Real>("Nks");
  Real bk       = pm_system->get_equation_systems().parameters.get<Real>("bk");
  c1  = bead_r/(2.0*bk);
  Ls  = Nks*bk/bead_r;

  force_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("point_WormLikeSpring");
  if(force_params.size()!=0){
    std::cout << std::endl << "********************Error message********************" << std::endl
              << "---------------> The force type 'point_WormLikeSpring' require 0 parameter" << std::endl
              << "****************************************" << std::endl;
    libmesh_error();    
  }	
}

void FixPointWormLikeSpring::print_fix()
{
  std::cout <<"this is FixPointWormLikeSpring" << std::endl;
}

void FixPointWormLikeSpring::compute()
{ 
  // get bonds from polymer_chain
  const std::vector<std::vector<std::size_t> >& bonds = polymer_chain->bonds();
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Loop each bond and apply forces
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for(std::size_t i=0; i<n_bonds; ++i)
  {
    // Connecting bead ids
    std::size_t bead_id_1 = bonds[i][1];  // connect bead 1
    std::size_t bead_id_2 = bonds[i][2];  // connect bead 2

    // Force on bead 1
    const Point pti   = point_mesh->particles()[bead_id_1]->point();
    const Point ptj   = point_mesh->particles()[bead_id_2]->point();
    const Point R_ij  = point_mesh->pm_periodic_boundary()->point_vector(pti,ptj);
    const std::vector<Real> F_ij = fix_base.spring_force_wls(R_ij,c1,Ls);

    // Force on bead 2
    std::vector<Real> F_ji (F_ij);
    for (std::size_t j=0; j<dim; ++j) F_ji[j] = -F_ij[j];

    // Add forces to beads
    point_mesh->particles()[bead_id_1]->add_particle_force(F_ij);
    point_mesh->particles()[bead_id_2]->add_particle_force(F_ji);
  }
}

