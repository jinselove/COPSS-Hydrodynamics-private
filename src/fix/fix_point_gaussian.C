#include "fix_point_gaussian.h"

FixPointGaussian::FixPointGaussian(PMLinearImplicitSystem& pm_sys_)
:
 FixPoint(pm_sys_)
{
  //force_type = "pp_ev_lj_cut";
	force_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("gaussian");
  this -> initParams();
}

void FixPointGaussian::initParams()
{
  if(force_params.size()!=2){
    std::cout << std::endl << "********************Error message********************" << std::endl
              << "---------------> The force type 'gaussian' requires 2 parameter (c1, c2) (dimensionless form)" << std::endl
              << "****************************************" << std::endl;
    libmesh_error();    
  }	
  c1 = force_params[0];
  c2 = force_params[1];
}

void FixPointGaussian::print_fix()
{
  std::cout <<"this is FixPointGaussian" << std::endl;
}

void FixPointGaussian::compute()
{	
  for(std::size_t p_id=0; p_id < num_points; ++p_id)
  {  
    // apply the excluded volume force to each particle i
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
        std::vector<Real> f_ij = fix_base.gaussian_force(r_ij, c1, c2);
        for (std::size_t _dim=0; _dim < dim; ++_dim) pforce[_dim] += f_ij[_dim];
      } // end if
    } // end for i-loop    
    point_mesh->particles()[p_id]->add_particle_force(pforce);
  }
}

