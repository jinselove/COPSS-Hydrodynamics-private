#include "fix_point_gaussian_dna.h"

FixPointGaussianDNA::FixPointGaussianDNA(PMLinearImplicitSystem& pm_sys_)
  : FixPoint(pm_sys_)
{
  force_type = "gaussian_dna";
  this->initPointParticleType();
  this->initParams();
}

void FixPointGaussianDNA::initPointParticleType()
{
  START_LOG("FixPointGaussianDNA::initPointParticleType()",
            "FixPointGaussianDNA");
  point_particle_model =
    pm_system->get_equation_systems().parameters.get<std::string>(
      "point_particle_model");

  if (point_particle_model != "polymer_chain") {
    std::cout << std::endl <<
    "*******************Error message*********************" << std::endl
              << "The force field 'gaussian_dna' is only for polymer_chain, i.e.,"
              << "not applicable to point_particle_model: " <<
    point_particle_model << std::endl
              << "****************************************" << std::endl;
    libmesh_error();
  }
  STOP_LOG("FixPointGaussianDNA::initPointParticleType()", "FixPointGaussianDNA");
}

void FixPointGaussianDNA::initParams()
{
  START_LOG("FixPointGaussianDNA::initParams()", "FixPointGaussianDNA");
  force_params =
    pm_system->get_equation_systems().parameters.get<std::vector<Real> >(
      "gaussian_dna");

  if (force_params.size() != 1) {
    std::cout << std::endl <<
    "********************Error message********************" << std::endl
              <<
    "---------------> The force type 'gaussian_dna' require 1 parameter (ev) (dimensionless form)"
              << std::endl
              << "****************************************" << std::endl;
    libmesh_error();
  }
  Real ev  = force_params[0] * bead_r * bead_r * bead_r;
  Real Ss2 = pm_system->get_equation_systems().parameters.get<Real>("Ss2");
  Real bk  = pm_system->get_equation_systems().parameters.get<Real>("bk");
  Real Nks = pm_system->get_equation_systems().parameters.get<Real>("Nks");

  // calculate c1, c2 for polymer_chain
  Real c1_tmp  = 3. / (4. * PI * Ss2);
  Real c1_tmp3 = c1_tmp * c1_tmp * c1_tmp;
  c1 = ev * Nks * Nks * std::sqrt(c1_tmp3);
  c2 = 3. * bead_r * bead_r / (4. * Ss2);
  STOP_LOG("FixPointGaussianDNA::initParams()", "FixPointGaussianDNA");
}

void FixPointGaussianDNA::print_fix()
{
  std::cout << "this is FixPointGaussianDNA" << std::endl;
}

void FixPointGaussianDNA::compute()
{
  START_LOG("FixPointGaussianDNA::compute()", "FixPointGaussianDNA");

  for (std::size_t p_id = 0; p_id < num_points; ++p_id)
  {
    // apply the excluded volume force to each particle i
    Point pforce(0.);
    const std::vector<Point>& neighbor_vector =
      point_particles[p_id]->neighbor_vector();

    // Loop over each neighbor
    for (std::size_t j = 0; j < neighbor_vector.size(); ++j)
    {
      pforce += fix_base.gaussian_force(neighbor_vector[j], c1, c2);
    } // end for i-loop
    point_particles[p_id]->add_particle_force(pforce);
  }
  START_LOG("FixPointGaussianDNA::compute()", "FixPointGaussianDNA");
}
