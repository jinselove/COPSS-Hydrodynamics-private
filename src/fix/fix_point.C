#include "fix_point.h"

FixPoint::FixPoint(PMLinearImplicitSystem& pm_sys_)
:Fix(pm_sys_)
{
 particle_type = pm_system->get_equation_systems().parameters.get<std::string>("particle_type");
 point_particle_model = pm_system->get_equation_systems().parameters.get<std::string>("point_particle_model");
}

//================================================================
void FixPoint::initParticleType()
{
 if(particle_type != "point_particle") {
	std::cout << std::endl << "*******************Error message*********************" << std::endl
              << "The force field = "<< force_type << " is only for point_particle, i.e., cannot be applied to particle type = " << particle_type << std::endl
              << "****************************************" << std::endl;
	libmesh_error();
  } 
}


void FixPoint::check_walls()
{
  wall_type = pm_system->get_equation_systems().parameters.get<std::string>("wall_type");
  // for slit walls
  if(wall_type == "slit"){
    for(std::size_t id=0; id<num_points; ++id){
      const std::size_t p_id = point_particles[id]->id();
      std::vector<int> & count = point_particles[p_id]->counter();
      Point& pt0     = point_particles[p_id]->point();
      // check all three directions of the box
      for(std::size_t i=0; i<dim; ++i)
      {
        if(periodic[i]) // periodic boundary
        {
          if( pt0(i) <  box_min(i) ) {pt0(i) += box_len(i); count[i] -= 1;}
          if( pt0(i) >= box_max(i) ) {pt0(i) -= box_len(i); count[i] += 1;}
        }
        else   // non-periodic boundary (inpenetrable wall)
        {
          if( pt0(i) <  box_min(i) ){
                  pt0(i) = 2.*box_min(i) - pt0(i);  
                  out_domain_counter += 1;
                  std::cout << "*** Warning: " << p_id << "-th bead is out of domain. Out of domain occurs "
                            << out_domain_counter << " times in total." << std::endl;
                }
          if( pt0(i) >=  box_max(i) ){
                  pt0(i) = 2.*box_max(i) - pt0(i);  
                  out_domain_counter += 1;
                  std::cout << "*** Warning: " << p_id << "-th bead is out of domain. Out of domain occurs "
                            << out_domain_counter << " times in total." << std::endl;
                }
        }
      } // end for i-loop
    } //end for id loop
  } // end wall_type == "slit"
  else{
    Real cavity_radius = wall_params[0];
    for(std::size_t id=0; id<num_points; ++id){
      const std::size_t p_id = point_particles[id]->id();
      std::vector<int> & count = point_particles[p_id]->counter();
      Point& pt0     = point_particles[p_id]->point();
      Real  pt0_norm = pt0.norm();
      Point pt0_unit = pt0.unit();
    // Move particle back into spherical cavity
      if( pt0_norm > cavity_radius ){ 
      pt0 = pt0_unit * (2.*cavity_radius) - pt0;
      out_domain_counter += 1;
      std::cout << "*** Warning: " << p_id << "-th bead ("<< pt0(0) << "," << pt0(1) << "," <<pt0(2) <<") is out of domain. Out of domain occurs "
                << out_domain_counter << " times in total." << std::endl;
      }
    } // end for loop 
  }  // end if wall_type =="sphere"
} // end function


void FixPoint::check_walls_pbcCount()
{
  // Fix this
}
