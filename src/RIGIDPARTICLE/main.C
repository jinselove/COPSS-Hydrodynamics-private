#include "copss_rigid_particle_system.h"

int main(int argc, char **argv) {
  // Init COPSS and libmesh
  CopssInit init(argc, argv);
  CopssRigidParticleSystem system(init);

  // Print out the current date and time
  time_t rawtime;
  struct tm *timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  system.start_time(timeinfo);

  // Init equation_systems
  EquationSystems equation_systems = system.init_system(
    "rigid_particle_control.in");

  // Run simulation
  system.run(equation_systems);

  // Print end time
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  system.end_time(timeinfo);

  // Return 0
  return 0;
}
