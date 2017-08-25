#include "fix_point_slitWall_lj_cut.h"

FixPointSlitWallLJCut::FixPointSlitWallLJCut(PMLinearImplicitSystem& pm_sys_)
:
 FixPoint(pm_sys_)
{
  this -> initParams();
}

void FixPointSlitWallLJCut::initParams()
{
  force_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("wall/lj_cut");
  wall_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>>("slit");
  if(force_params.size()!=3){
    std::cout << std::endl << "********************Error message********************" << std::endl
              << "---------------> The force type 'wall/lj_cut' requires 3 parameter (epsilon, sigma, rcut) (dimensionless form)" << std::endl
              << "****************************************" << std::endl;
    libmesh_error();    
  }
  for (int i = 0; i < dim; i++){
    if(box_min(i) != wall_params[2*i] or box_max(i) != wall_params[2*i+1]){
      std::cout << std::endl << "********************Error message********************" << std::endl
                << "for wall/lj_cut (slit), wall_params has to be the same as box_min and box_max" << std::endl
                << "****************************************" << std::endl;
      libmesh_error();
    }
  }
  epsilon = force_params[0];
  sigma = force_params[1];
  rCut = force_params[2];
}

void FixPointSlitWallLJCut::print_fix()
{
  std::cout <<"this is FixPointSlitWallLJCut" << std::endl;
}

void FixPointSlitWallLJCut::compute()
{	
  for(std::size_t i=0; i<num_points; ++i)
  {
    std::vector<Real> pforce(dim);
    const Point pti = point_mesh->particles()[i] -> point();
    for(std::size_t j = 0; j<dim;++j){
      Point r_i_lo, r_i_hi;
      if(periodic[j]==false and inlet[j]==false){
        r_i_lo(j) = box_min(j) - pti(j);// opposite sign compared to rij = rj - r_i 
        r_i_hi(j) = box_max(j)- pti(j);// opposite sign compared to rij = rj - r_i 
        if(r_i_lo.norm() <= rCut) pforce[j] += fix_base.lj_force(r_i_lo, epsilon,sigma)[j];
        if(r_i_hi.norm() <= rCut) pforce[j] += fix_base.lj_force(r_i_hi, epsilon,sigma)[j];
      }// end of
    } // end for (i<_dim)
    point_mesh->particles()[i]->add_particle_force(pforce);    
  } // end for i-loop 
}

