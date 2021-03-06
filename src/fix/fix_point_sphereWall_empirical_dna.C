#include "fix_point_sphereWall_empirical_dna.h"

FixPointSphereWallEmpiricalDNA::FixPointSphereWallEmpiricalDNA(
  PMLinearImplicitSystem& pm_sys_)
  :
  FixPoint(pm_sys_)
{
  force_type = "wall/empirical_dna";
  this->initPointParticleType();
  this->initParams();
}

void FixPointSphereWallEmpiricalDNA::initPointParticleType()
{
  START_LOG("FixPointSphereWallEmpiricalDNA::initPointParticleType()",
            "FixPointSphereWallEmpiricalDNA");
  point_particle_model =
    pm_system->get_equation_systems().parameters.get<std::string>(
      "point_particle_model");

  if (point_particle_model != "polymer_chain") {
    std::cout << std::endl <<
    "*******************Error message*********************" << std::endl
              <<
    "The force field 'wall/empirical_dna' is only for polymer_chain, i.e.,"
              << "not applicable to point_particle_model: " <<
    point_particle_model << std::endl
              << "****************************************" << std::endl;
    libmesh_error();
  }
  STOP_LOG("FixPointSphereWallEmpiricalDNA::initPointParticleType()",
           "FixPointSphereWallEmpiricalDNA");
}

void FixPointSphereWallEmpiricalDNA::initParams()
{
  START_LOG("FixPointSphereWallEmpiricalDNA::initParams()",
            "FixPointSphereWallEmpiricalDNA");
  wall_params =
    pm_system->get_equation_systems().parameters.get<std::vector<Real> >("sphere");

  if (wall_params.size() != 1) {
    std::cout << std::endl <<
    "********************Error message********************" << std::endl
              <<
    "---------------> The force type 'wall/empirical_dna' (sphere) requires 1 wall parameter, i.e., radius (dimensionless form)"
              << std::endl
              << "****************************************" << std::endl;
    libmesh_error();
  }
  sphereWallRadius = wall_params[0];
  const Real bk  = pm_system->get_equation_systems().parameters.get<Real>("bk");
  const Real Nks = pm_system->get_equation_systems().parameters.get<Real>("Nks");
  Real c1        = bead_r / bk;
  Real c2        = c1 / std::sqrt(Nks);
  d0 = 0.5 / c2;
  c0 = 25.0 * c1;
  STOP_LOG("FixPointSphereWallEmpiricalDNA::initParams()",
           "FixPointSphereWallEmpiricalDNA");
}

void FixPointSphereWallEmpiricalDNA::print_fix()
{
  START_LOG("FixPointSphereWallEmpiricalDNA::print_fix()",
            "FixPointSphereWallEmpiricalDNA");
  std::cout << "this is FixPointSphereWallEmpiricalDNA" << std::endl;
  STOP_LOG("FixPointSphereWallEmpiricalDNA::print_fix()",
           "FixPointSphereWallEmpiricalDNA");
}

void FixPointSphereWallEmpiricalDNA::compute()
{
  START_LOG("FixPointSphereWallEmpiricalDNA::compute()",
            "FixPointSphereWallEmpiricalDNA");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Loop each point and apply forces
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for (std::size_t i = 0; i < num_points; ++i)
  {
    //  retrieve bead position and the box boundary
    const Point pti            = point_particles[i]->point();
    const Point r_i_sphereWall = pti.unit() * (sphereWallRadius - pti.norm());

    // compute force
    Point pforce = fix_base.polymer_wall_empirical_force(r_i_sphereWall, c0, d0);

    // attach this force to particle i
    point_particles[i]->add_particle_force(pforce);
  } // end for i-loop
  STOP_LOG("FixPointSphereWallEmpiricalDNA::compute()",
           "FixPointSphereWallEmpiricalDNA");
}
