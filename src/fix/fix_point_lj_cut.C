#include "fix_point_lj_cut.h"

FixPointLJCut::FixPointLJCut(PMLinearImplicitSystem& pm_sys_)
:FixPoint(pm_sys_)
{
  force_type = "lj_cut";
  this -> initParams();
}

void FixPointLJCut::initParams()
{
  START_LOG("FixPointLJCut::initParams()", "FixPointLJCut");
  force_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("lj_cut");
  if(force_params.size()!=3){
    std::cout << std::endl << "********************Error message********************" << std::endl
              << "---------------> The force type 'lj_cut' requires 3 parameter (epsilon, sigma, rcut) (dimensionless form)" << std::endl
              << "****************************************" << std::endl;
    libmesh_error();    
  }
  epsilon = force_params[0];
  sigma = force_params[1];
  rCut = force_params[2];
  STOP_LOG("FixPointLJCut::initParams()", "FixPointLJCut");
}

void FixPointLJCut::print_fix()
{
  START_LOG("FixPointLJCut::print_fix()", "FixPointLJCut");
  std::cout <<"this is FixPointLJCut" << std::endl;
  STOP_LOG("FixPointLJCut::print_fix()", "FixPointLJCut");
}

void FixPointLJCut::compute()
{
  START_LOG("FixPointLJCut::compute()", "FixPointLJCut");
  for(std::size_t p_id=0; p_id<num_points; ++p_id)
  {  
    // apply lj force on particle i
    Point pforce;
    const std::vector<Point>& neighbor_vector = point_particles[p_id]->neighbor_vector();    
    // Loop over each neighbor
    for (std::size_t j=0; j<neighbor_vector.size(); ++j)
    {
      if(neighbor_vector[j].norm() <= rCut){
        pforce += fix_base.lj_force(neighbor_vector[j], epsilon, sigma);
      }
    } // end loop over neighbors 
    point_particles[p_id]->add_particle_force(pforce);
  } 
  STOP_LOG("FixPointLJCut::compute()", "FixPointLJCut");
}

