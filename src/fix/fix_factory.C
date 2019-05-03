#include "fix_factory.h"

// add else branch if you need to implment new force field
Fix * FixFactory::buildFix(std::string           & fix_name,
                           PMLinearImplicitSystem& pm_system)
{
  particle_type = pm_system.get_equation_systems().parameters.get<std::string>(
    "particle_type");
  wall_type = pm_system.get_equation_systems().parameters.get<std::string>(
    "wall_type");

  if (particle_type == "point_particle") {
    // particle-particle force
    if (fix_name == "lj_cut") return new FixPointLJCut(pm_system);
    else if (fix_name == "gaussian") return new FixPointGaussian(pm_system);
    else if (fix_name ==
             "gaussian_dna") return new FixPointGaussianDNA(pm_system);
    else if (fix_name == "worm_like_spring") return new FixPointWLS(pm_system);

    // particle wall force
    else if (fix_name == "wall/lj_cut")
    {
      if (wall_type == "slit") return new FixPointSlitWallLJCut(pm_system);
      else if (wall_type ==
               "sphere") return new FixPointSphereWallLJCut(pm_system);
      else {
        std::cout <<
          "Error: FixFactory::buildFix: undefined wall_type for point_particles: "
                  << fix_name << std::endl;
        libmesh_error();
      }
    }
    else if (fix_name == "discretized_wall/lj_cut") {
      return new FixPointDiscretizedWallLJCut(pm_system);
    }
    else if (fix_name == "wall/empirical_dna")
    {
      if (wall_type == "slit") return new FixPointSlitWallEmpiricalDNA(pm_system);
      else if (wall_type == "sphere") return new FixPointSphereWallEmpiricalDNA(
          pm_system);
      else {
        std::cout <<
          "Error: FixFactory::buildFix: undefined wall_type for point_particles: "
                  << fix_name << std::endl;
        libmesh_error();
      }
    }
    else if (fix_name == "p_constant") {
      return new FixPointConstant(pm_system);
    }
    else {
      std::cout << "Error: undefined force field for point_particles: " <<
        fix_name << std::endl;
      libmesh_error();
    }
  }
  else {
    // rigid particle surface constraint has to be added by default
    if (fix_name == "surface_constraint") return new FixRigidSurfaceConstraint(
        pm_system);
    else if (fix_name == "sedimentation") return new FixRigidSedimentation(
        pm_system);
    else if (fix_name == "lj_cut") return new FixRigidLJCut(pm_system);

    // particle wall force
    else if (fix_name == "wall/lj_cut")
    {
      if (wall_type == "slit") return new FixRigidSlitWallLJCut(pm_system);
      else if (wall_type ==
               "sphere") return new FixRigidSphereWallLJCut(pm_system);
      else {
        std::cout <<
        "Error: FixFactory::buildFix: undefined wall_type for rigid_particles: "
                  << fix_name << std::endl;
        libmesh_error();
      }
    }
    else {
      std::cout << "Error: undefined force field for rigid_particles: " <<
      fix_name << std::endl;
      libmesh_error();
    }
  }
}
