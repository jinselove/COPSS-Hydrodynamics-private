#include "fix_point_slitWall_empirical_dna.h"

FixPointSlitWallEmpiricalDNA::FixPointSlitWallEmpiricalDNA(PMLinearImplicitSystem& pm_sys_)
:
 FixPoint(pm_sys_)
{
  force_type = "wall/empirical_dna";
  this -> initPointParticleType();
  this -> initParams();
}

void FixPointSlitWallEmpiricalDNA::initPointParticleType()
{
START_LOG("FixPointSlitWallEmpiricalDNA::initPointParticleType()", "FixPointSlitWallEmpiricalDNA");
  point_particle_model = pm_system->get_equation_systems().parameters.get<std::string>("point_particle_model");
  if(point_particle_model != "polymer_chain") {
  std::cout << std::endl << "*******************Error message*********************" << std::endl
              << "The force field 'wall/empirical_dna' is only for polymer_chain, i.e.,"
              << "not applicable to point_particle_model: " << point_particle_model << std::endl
              << "****************************************" << std::endl;
  libmesh_error();
  }   
STOP_LOG("FixPointSlitWallEmpiricalDNA::initPointParticleType()", "FixPointSlitWallEmpiricalDNA");
}

void FixPointSlitWallEmpiricalDNA::initParams()
{
START_LOG("FixPointSlitWallEmpiricalDNA::initParams()", "FixPointSlitWallEmpiricalDNA");
  wall_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>>("slit");
  for (int i = 0; i < dim; i++){
    if(box_min(i) != wall_params[2*i] or box_max(i) != wall_params[2*i+1]){
      std::cout << std::endl << "********************Error message********************" << std::endl
                << "for force field 'wall/empirical_dna' (slit), wall_params has to be the same as box_min and box_max" << std::endl
                << "****************************************" << std::endl;
      libmesh_error();
    }
  }
  const Real bk       = pm_system->get_equation_systems().parameters.get<Real>("bk");
  const Real Nks      = pm_system->get_equation_systems().parameters.get<Real>("Nks");
  Real c1 = bead_r / bk;
  Real c2 = c1 / std::sqrt(Nks);
  d0 = 0.5/c2;
  c0 = 25.0*c1;
STOP_LOG("FixPointSlitWallEmpiricalDNA::initParams()", "FixPointSlitWallEmpiricalDNA");
}

void FixPointSlitWallEmpiricalDNA::print_fix()
{
START_LOG("FixPointSlitWallEmpiricalDNA::print_fix()", "FixPointSlitWallEmpiricalDNA");
  std::cout <<"this is FixPointSlitWallEmpiricalDNA" << std::endl;
STOP_LOG("FixPointSlitWallEmpiricalDNA::print_fix()", "FixPointSlitWallEmpiricalDNA");
}

void FixPointSlitWallEmpiricalDNA::compute()
{
START_LOG("FixPointSlitWallEmpiricalDNA::compute()", "FixPointSlitWallEmpiricalDNA");
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop each point and apply forces
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for(std::size_t i=0; i<num_points; ++i)
  {
    std::vector<Real> pforce(dim);
    //  retrieve bead position and the box boundary
    const Point pti = point_particles[i]->point();
    // compute the particle-wall interaction force
    for (std::size_t j=0; j<dim; ++j){  // loop over each direction
    // compute the distance to the wall if NO periodic boundary and no inlet 
    Point r_i_lo, r_i_hi;
    if ( periodic[j]==false and inlet[j]==false) {
      // distance to the parallel walls
      r_i_lo(j) = box_min(j) - pti(j) ; 
      r_i_hi(j) = box_max(j) - pti(j);       
      // bead-wall interaction force
      pforce[j] += fix_base.polymer_wall_empirical_force(r_i_lo, c0, d0)[j];
      pforce[j] += fix_base.polymer_wall_empirical_force(r_i_hi, c0, d0)[j];
    } // end if
  } // end for j-loop    
    point_particles[i]->add_particle_force(pforce);    
  } // end for i-loop	
STOP_LOG("FixPointSlitWallEmpiricalDNA::compute()", "FixPointSlitWallEmpiricalDNA");
}

