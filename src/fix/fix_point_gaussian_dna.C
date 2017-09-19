#include "fix_point_gaussian_dna.h"

FixPointGaussianDNA::FixPointGaussianDNA(PMLinearImplicitSystem& pm_sys_)
:FixPoint(pm_sys_)
{
  this -> initPointParticleType();
  this -> initParams();
}

void FixPointGaussianDNA::initPointParticleType()
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

void FixPointGaussianDNA::initParams()
{
  force_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("gaussian_dna");
  if(force_params.size()!=1){
    std::cout << std::endl << "********************Error message********************" << std::endl
              << "---------------> The force type 'gaussian_dna' require 1 parameter (ev) (dimensionless form)" << std::endl
              << "****************************************" << std::endl;
    libmesh_error();    
  }	
  Real ev = force_params[0] * bead_r * bead_r * bead_r;
  Real Ss2      = pm_system->get_equation_systems().parameters.get<Real>("Ss2");
  Real bk       = pm_system->get_equation_systems().parameters.get<Real>("bk");
  Real Nks      = pm_system->get_equation_systems().parameters.get<Real>("Nks");

  // calculate c1, c2 for polymer_chain
  Real c1_tmp = 3./(4.*PI*Ss2);
  Real c1_tmp3 = c1_tmp * c1_tmp * c1_tmp;
  c1 = ev*Nks*Nks*std::sqrt(c1_tmp3);
  c2   = 3.*bead_r*bead_r/(4.*Ss2);    
}

void FixPointGaussianDNA::print_fix()
{
  std::cout <<"this is FixPointGaussianDNA" << std::endl;
}

void FixPointGaussianDNA::compute()
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

