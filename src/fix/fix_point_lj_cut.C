#include "fix_point_lj_cut.h"

FixPointLJCut::FixPointLJCut(PMLinearImplicitSystem& pm_sys_)
:FixPoint(pm_sys_)
{
  this -> initParams();
}

void FixPointLJCut::initParams()
{
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
}

void FixPointLJCut::print_fix()
{
  std::cout <<"this is FixPointLJCut" << std::endl;
}

void FixPointLJCut::compute()
{	
  for(std::size_t p_id=0; p_id<num_points; ++p_id)
  {  
    // apply lj force on particle i
    std::vector<Real> pforce(dim);
    const Point pti   = point_mesh->particles()[p_id]->point();
    std::vector<std::pair<std::size_t,Real> > n_list = point_mesh->particles()[p_id]->neighbor_list();
    // Loop over each neigbhor
    for (std::size_t j=0; j<n_list.size(); ++j)
    {
     const std::size_t n_id  = n_list[j].first;
     if(p_id != n_id)  // make sure this bead and the neighboring bead are not the same bead.
     {
      const Point ptj   = point_mesh->particles()[n_id]->point();
      const Point r_ij  = point_mesh->pm_periodic_boundary()->point_vector(pti,ptj);
      if(r_ij.norm() <= rCut){
        std::vector<Real> f_ij = fix_base.lj_force(r_ij, epsilon, sigma);
        for (std::size_t _dim=0; _dim<dim; ++_dim) pforce[_dim] += f_ij[_dim];
      }
     } // end if
    } // end loop over neighbors 
    point_mesh->particles()[p_id]->add_particle_force(pforce);
  }  
}

