#include "copss_rigid_particle_system.h"

using std::cout;
using std::endl;
using std::string;


namespace libMesh{
//==========================================================================
CopssRigidParticleSystem::CopssRigidParticleSystem(CopssInit& init)
:Copss(init)
{

}

CopssRigidParticleSystem::~CopssRigidParticleSystem(){
  delete _particle_mesh;
  _particle_mesh = NULL;
	delete point_mesh;
	point_mesh = NULL;
}


//==========================================================================
void CopssRigidParticleSystem::read_particle_info(){
	if (particle_type != "rigid_particle"){
		error_msg = "invalid particle type ("+particle_type+") defined\n";
		PMToolBox::output_message(error_msg,comm_in);
		libmesh_error();
	}

  // Particles' meshes
  _particle_mesh_type = input_file("particle_mesh_type", "surface_mesh");  
  _particle_mesh_file.resize( input_file.vector_variable_size(_particle_mesh_type) ); // set the length of input vector
  for (unsigned int i=0; i < _particle_mesh_file.size(); i++){
    _particle_mesh_file[i] = input_file(_particle_mesh_type, "particle_mesh.e", i);
  }

}// end read_particle_info()

//==========================================================================
void CopssRigidParticleSystem::create_object(){

  // initialize _particle_mesh
  _particle_mesh = new ParticleMesh<3> (*mesh,search_radius_p,search_radius_e);
  // add periodic boundary
  _particle_mesh -> add_periodic_boundary(*pm_periodic_boundary);
  //Read the particle data
  std::ostringstream pfilename;
    // if restart, read particle data from ...
  if(restart){

    // implement this for restart mode
    cout <<"Warning: restart mode for rigid particles has not been implemented yet" << endl;
    libmesh_error();
  }
    // if not restart, read particle data from "rigid_particle_data.in"
  else{
    pfilename <<"rigid_particle_data.in";
    _particle_mesh -> read_particles_data(pfilename.str(), _particle_mesh_type, _particle_mesh_file);    
  }

  // reinit _particle_mesh
  _particle_mesh -> reinit();
  _hsize_solid = _particle_mesh->mesh_size();// mesh size of solid
  _hmins = _hsize_solid[0];
  _hmaxs = _hsize_solid[1];
  _num_particles = _particle_mesh -> num_particles();

  // attach mesh spring network
  this -> attach_mesh_spring_network();

  // print out information
  cout<<"##########################################################\n"
           <<"#                  Particle Parameters                    \n"
           <<"##########################################################\n\n"
           <<"   particle type             : " << particle_type.c_str() << endl
           <<"   particle mesh type        : " << _particle_mesh_type.c_str() << endl
           <<"   minimum mesh size of particle surface: hmins  : " << _hmins << endl
           <<"   maximum mesh size of particle surface: hmaxs = " << _hmins << endl;
  cout <<"------------> The non-dimensional variables:\n"
       <<"   non-dimensional bead radius      a0     : " << 1.0 << "\n"
       <<"   non-dimensional ksi = sqrt(PI)/(3a0)    : " << std::sqrt(PI) / 3. <<"\n";

  pfilename.str("");
  pfilename.clear();
  comm_in.barrier();
}//end function

//==========================================================================
void CopssRigidParticleSystem::attach_mesh_spring_network()
{
  _mesh_spring_network.resize(_num_particles);
  //std::vector<MeshSpringNetwork*> (particle_mesh.num_particles());
  for(std::size_t i=0; i<_num_particles; ++i)
  {
    MeshBase& p_mesh = _particle_mesh->particles()[i]->mesh();
    const Point& center0 = _particle_mesh->particles()[i]->center();
    _mesh_spring_network[i] = new MeshSpringNetwork(p_mesh,*pm_periodic_boundary);
    _mesh_spring_network[i]->build_spring_network(center0);
    _particle_mesh->particles()[i]->attach_mesh_spring_network(_mesh_spring_network[i]);
  }

}

//=====================================================================
void CopssRigidParticleSystem::create_object_mesh(){
  // prepare domain and objects
  cout << "\n==>(1/4) Generate/Create domain Mesh\n";
  this -> create_domain_mesh();
  cout << "\n==>(2/4) Create periodic box \n ";
  this -> create_periodic_boundary();
  cout << "\n==>(3/4) Create particle mesh object\n";
  this -> create_object();

  cout << "\n==>(4/4) Create point_mesh object \n";
  // Create object mesh
  point_mesh = new PointMesh<3> (*_particle_mesh, search_radius_p, search_radius_e);

  // No need to add periodic boundary, which is already included in particle_mesh

  // Reinit point_mesh
  point_mesh->reinit();

  // finish point_mesh, print information

  cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
  cout << "### The particle-mesh and point-mesh info:\n";
  cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
  cout << "Total number of particles: " << _particle_mesh->num_particles() << "\n";
  cout << "Total number of points: " << _particle_mesh->num_mesh_points() << "\n";
  // iterate over all particles for debug purpose
  if (debug_info = true){
    for (std::size_t i=0; i<_particle_mesh->num_particles(); ++i) {
      Real vol  = _particle_mesh->particles()[i]->compute_volume();
      Real area = _particle_mesh->particles()[i]->compute_area();
      Point pc_i = _particle_mesh->particles()[i]->compute_centroid();
      const std::vector<Real> hs = _particle_mesh->particles()[i]->mesh_size();
      const std::size_t n_surf_elem = _particle_mesh->particles()[i]->num_mesh_elem();
      const std::size_t n_surf_node = _particle_mesh->particles()[i]->num_mesh_nodes();

      cout << "------> Volume of the "<< i <<"-th particle is "<< vol <<", and area is "<< area << endl;
      if(dim == 2){
        cout << "        The centroid = ("<< pc_i(0) <<","<< pc_i(1) <<")" <<endl;
      }
      else if(dim == 3){
        cout << "        The centroid = ("<< pc_i(0) <<","<< pc_i(1) <<","<< pc_i(2) <<")" <<endl;
      }
      cout << "        hmin = "<< hs[0] <<", hmax = "<< hs[1] << endl;
      cout << "  Particle "<<i<<": center " << _particle_mesh->particles()[i]->center() << endl;
      cout << "         " << n_surf_elem << " elements for the particle's mesh" << endl;
      cout << "         " << n_surf_node << " nodes  for the particle's mesh" << endl;
    }
  }
  cout << "search_radius_p = " << search_radius_p << ", search_radius_e = " << search_radius_e << endl<<endl<<endl;
} // end function


//==================================================================================
void CopssRigidParticleSystem::attach_object_mesh(PMLinearImplicitSystem& system)
{
  system.attach_particle_mesh(_particle_mesh);
	system.attach_point_mesh(point_mesh);
}

//======================================================================================
void CopssRigidParticleSystem::set_parameters(EquationSystems& equation_systems){
  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = max_linear_iterations;
  equation_systems.parameters.set<Real> ("linear solver rtol") = linear_solver_rtol;
  equation_systems.parameters.set<Real> ("linear solver atol") = linear_solver_atol;
  equation_systems.parameters.set<bool>    ("user_defined_pc") = user_defined_pc;
  equation_systems.parameters.set<bool>     ("schur_user_ksp") = schur_user_ksp;
  equation_systems.parameters.set<Real>("schur_user_ksp_rtol") = schur_user_ksp_rtol;
  equation_systems.parameters.set<Real>("schur_user_ksp_atol") = schur_user_ksp_atol;
  equation_systems.parameters.set<string>    ("schur_pc_type") = schur_pc_type;
  equation_systems.parameters.set<StokesSolverType> ("solver_type") = solver_type;
  equation_systems.parameters.set<Real>              ("alpha") = alpha;
  equation_systems.parameters.set<Real>         ("kBT")        = kBT;
  equation_systems.parameters.set<Real>   ("fluid mesh size")  = min_mesh_size;
  equation_systems.parameters.set<Real>    ("solid mesh size") = _hmins;
  equation_systems.parameters.set<Real>       ("viscosity_0")  = muc;
  equation_systems.parameters.set<Real>               ("br0")  = 1.0;
  equation_systems.parameters.set<Real>       ("bead radius")  = Rb;
  equation_systems.parameters.set<Real>              ("drag")  = drag_c;
  equation_systems.parameters.set<Real>                ("tc")  = tc;
  equation_systems.parameters.set<string> ("particle_type")  = particle_type;
  equation_systems.parameters.set<string> ("particle_mesh_type") = _particle_mesh_type;
  // Attach force fields
  equation_systems.parameters.set<std::vector<string>> ("pp_force_types") = pp_force_type;
  for (int i=0; i<num_pp_force; i++) equation_systems.parameters.set<std::vector<Real>> (pp_force[i].first) = pp_force[i].second;
  equation_systems.parameters.set<std::vector<string>> ("pw_force_types") = pw_force_type;
  for (int i=0; i<num_pw_force; i++) equation_systems.parameters.set<std::vector<Real>> (pw_force[i].first) = pw_force[i].second;
  equation_systems.parameters.set<string> ("test_name") = test_name;
  equation_systems.parameters.set<string> ("wall_type") = wall_type;
  equation_systems.parameters.set<std::vector<Real>> (wall_type) = wall_params;
}

void CopssRigidParticleSystem::update_object(std::string stage)
{
  point_mesh->update_particle_mesh(_particle_mesh);
  if(debug_info == true){
    std::cout << "particle meshes are updated " << stage << endl;
  }
}


void CopssRigidParticleSystem::run(EquationSystems& equation_systems){
  PerfLog perf_log("Copss-Hydrodynamics-RigidParticleSystem");
  cout<<endl<<"============================4. Start moving particles ============================"<<endl<<endl;
  // get stokes system from equation systems
  PMLinearImplicitSystem& system = equation_systems.get_system<PMLinearImplicitSystem> ("Stokes");
   
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Parameters for dynamic process
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  NP     = point_mesh->num_particles();
  n_vec  = dim*NP;
  hmin = std::min(_hmins, hminf);
  hmax = std::max(_hmaxs, hmaxf);
  _write_particle_mesh = input_file("write_particle_mesh", true);

  // Get a better conformation of polymer chains before simulation.
  this -> update_object("in initial data input");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Compute undisturbed velocity field without particles.
  NOTE: We MUST re-init particle-mesh before solving Stokes
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  cout<<"==>(1/3) Compute the undisturbed velocity field"<<endl;
  perf_log.push ("solve undisturbed_system");
  this -> solve_undisturbed_system(equation_systems); 
  perf_log.pop ("solve undisturbed_system");
  /* output particle data at the 0-th step in the VTK format */
  if(_write_particle_mesh==true and restart==false)
  {
    _particle_mesh -> write_particle_mesh(smesh_file_name);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create vectors and Shell Mat for use:
   U0:          particle velocity vector;
   R0/R_mid:    particle position vector;
   dw/dw_mid:   random vector;
   RIN/ROUT:    the initial and intermediate particle postion vector for msd output
   RIN will not change, and ROUT excludes pbc
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  cout<<"==>(2/3) Prepare RIN & ROUT and Brownian_system in binary format at step 0"<<endl;
  this -> create_brownian_system(equation_systems);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Advancing in time. Fixman Mid-Point algorithm
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  cout<<"==>(3/3) Start calculating dynamics and advancing time steps"<<endl;
  const unsigned int istart = restart_step + 1;
  const unsigned int iend   = restart_step + nstep;
  o_step = restart_step;  // output step
  vel0.resize(n_vec);
  vel1.resize(n_vec);
  real_time = restart_time;
  //start integration

  perf_log.push ("integration");
  for(unsigned int i=istart; i<=iend; ++i)
  {
    cout << "\nStarting Fixman Mid-Point algorithm at step "<< i << endl;
    // integrate particle movement using fixman's mid point scheme
    this -> fixman_integrate(equation_systems, i);  
    // Update the time
    if(i%write_interval==0){
      smesh_file_name_i <<smesh_file_name 
                        <<"-s."
                        << std::setw(8) 
                        << std::setfill('0') 
                        << std::right 
                        << o_step;

  //    brownian_sys->output_statistics_stepi(out_msd_flag, out_stretch_flag, out_gyration_flag, out_com_flag,
  //                                           i, real_time, center0, ROUT);
      /*----------------------------------------------------------------------------------------------------
       * Write out ROUT for restart mode at step i
      ----------------------------------------------------------------------------------------------------*/
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vector_ROUT.dat",FILE_MODE_WRITE,&viewer);
      VecView(ROUT,viewer);
    }

  } // end step integration
  
  perf_log.pop ("integration");
}

} // end of namespace
