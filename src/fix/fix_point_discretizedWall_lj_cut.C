#include <fstream>
#include "fix_point_discretizedWall_lj_cut.h"

FixPointDiscretizedWallLJCut::FixPointDiscretizedWallLJCut(PMLinearImplicitSystem& pm_sys_)
:
 FixPoint(pm_sys_)
{
  force_type = "discretized_wall/lj_cut";
  this -> initParams();
}

void FixPointDiscretizedWallLJCut::initParams()
{
START_LOG("FixPointDiscretizedWallLJCut::initParams()", "FixPointDiscretizedWallLJCut");
  force_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("discretized_wall/lj_cut");
  if(force_params.size()!=4){
    std::cout << std::endl << "********************Error message********************" << std::endl
              << "---------------> The force type 'discretized_wall/lj_cut' requires 4 parameter (epsilon, sigma, rcut, num_wall_particles) (dimensionless form)" << std::endl
              << "****************************************" << std::endl;
    libmesh_error();    
  }
  epsilon = force_params[0];
  sigma = force_params[1];
  rCut = force_params[2];
  num_wall_particles = int(force_params[3]);
  this -> read_wall_particle(); 
STOP_LOG("FixPointDiscretizedWallLJCut::initParams()", "FixPointDiscretizedWallLJCut");
}

void FixPointDiscretizedWallLJCut::print_fix()
{
START_LOG("FixPointDiscretizedWallLJCut::print_fix()", "FixPointDiscretizedWallLJCut");
  std::cout <<"this is FixPointDiscretizedWallLJCut" << std::endl;
STOP_LOG("FixPointDiscretizedWallLJCut::print_fix()", "FixPointDiscretizedWallLJCut");
}


void FixPointDiscretizedWallLJCut::compute()
{
START_LOG("FixPointDiscretizedWallLJCut::compute()", "FixPointDiscretizedWallLJCut");
  for(std::size_t p_id=0; p_id<num_points; ++p_id)
  {  
    // apply lj force on particle i
    Point pforce(0.);
    const Point pti   = point_particles[p_id]->point();
    // Loop over each wall_particle
    for (std::size_t j=0; j<wall_particle_pos.size(); ++j)
    {
      const Point r_ij  = point_mesh->pm_periodic_boundary()->point_vector(pti,wall_particle_pos[j]);
      if(r_ij.norm() <= rCut){
        Point f_ij = fix_base.lj_force(r_ij, epsilon, sigma);
        pforce += f_ij;
      }
    } // end loop over neighbors 
    point_particles[p_id]->add_particle_force(pforce);
  } 
}

void FixPointDiscretizedWallLJCut::read_wall_particle()
{
START_LOG("FixPointDiscretizedWallLJCut::read_wall_particle()", "FixPointDiscretizedWallLJCut");
  wall_particle_pos.resize(num_wall_particles); 
  // Open the local file and check the existance
  std::cout <<"\n###Read wall particle positions"<<std::endl;
  std::cout <<"   filename = "<<wall_particle_filename <<std::endl;
  std::ifstream infile;
  infile.open (wall_particle_filename, std::ios_base::in);
  if( !infile.good() )
  {
    printf("***warning: FixPointDiscretizedWallLJCut::read_wall_particle() can NOT read the wall_particle data!");
    libmesh_error();
  }
  // Read beads data
  for (unsigned int i=0; i<num_wall_particles; i++){ 
   infile >> wall_particle_pos[i](0) >> wall_particle_pos[i](1) >> wall_particle_pos[i](2);
  }
 // Finish and close the file
  infile.close();
  std::cout << "Reading wall particle positions from "<<wall_particle_filename<<" is completed!\n\n";


STOP_LOG("FixPointDiscretizedWallLJCut::read_wall_particle()", "FixPointDiscretizedWallLJCut");
}
