#include "copss_rigid_particle_system.h"

using std::string;


namespace libMesh {
// ==========================================================================
CopssRigidParticleSystem::CopssRigidParticleSystem(CopssInit& init)
  : Copss(init)
{
  // nothing
}

CopssRigidParticleSystem::~CopssRigidParticleSystem() {
  delete particle_mesh; particle_mesh = nullptr;
  for (int i = 0; i < mesh_spring_network.size(); i++) 
  {
    delete mesh_spring_network[i]; mesh_spring_network[i] = nullptr;
  }
  mesh_spring_network.clear();
}

// ==========================================================================
void CopssRigidParticleSystem::read_ggem_info() {
  alpha    = input_file("alpha", 0.1);
  ibm_beta = input_file("ibm_beta", 0.75);
  ss << "##########################################################\n"
     << "#                 GGEM information                       \n"
     << "##########################################################\n"
     << "-----------> the smoothing parameter in GGEM alpha = " << alpha << "\n"
     << "-----------> recommend meshsize <= " << 1. / (std::sqrt(2) * alpha) << "\n"
     << "-----------> the ibm beta = " << ibm_beta <<"  (IBM-GGEM ksi = 1./(beta*hmins), the optimized value of this variable depends on number of surface nodes";
  PMToolBox::output_message(ss, *comm_in);
  ss.clear();
}

// ==========================================================================
void CopssRigidParticleSystem::read_particle_info() {
  particle_type = input_file("particle_type", "other");

  if (particle_type != "rigid_particle") {
    error_msg = "invalid particle type (" + particle_type + ") defined\n";
    PMToolBox::output_message(error_msg, *comm_in);
    libmesh_error();
  }

  // Particles' meshes
  particle_mesh_type = input_file("particle_mesh_type", "surface_mesh");
} // end read_particle_info()

// ==========================================================================
void CopssRigidParticleSystem::create_object() {
  // initialize _particle_mesh
  particle_mesh = new ParticleMesh<3>(*mesh, search_radius_p, search_radius_e);

  // add periodic boundary
  particle_mesh->add_periodic_boundary(*pm_periodic_boundary);

  // Read the particle data
  std::ostringstream pfilename;

  // if restart, read particle data from ...
  pfilename << "rigid_particle_data.in";

  // read particle data from "rigid_particle_data.in"
  particle_mesh->read_particles_data(pfilename.str(), particle_mesh_type);
  hsize_solid         = particle_mesh->mesh_size(); // mesh size of solid
  hmins               = hsize_solid[0];
  hmaxs               = hsize_solid[1];
  num_rigid_particles = particle_mesh->num_particles();

  // attach mesh spring network
  this->attach_mesh_spring_network();

  // initialize center0 for calculation of MSD
  particle_mesh->initial_particle_center_of_mass(center0);

  // if restart, update particle positions from previous trajectories
  if (restart) {
    pfilename.str("");
    pfilename.clear();
    pfilename << "output_surface_node.csv";
    particle_mesh->read_particles_data_restart(pfilename.str(), o_step);
    ss << "##################### Restart mode ##########################\n"
       << "----> read rigid_particle surface node from " << pfilename.str() << "\n"
       << "#############################################################";
    PMToolBox::output_message(ss, *comm_in);
  }
  // reinit _particle_mesh
  particle_mesh->reinit();
  // print out information
  ss << "##########################################################\n"
     << "#                  Particle Parameters                    \n"
     << "##########################################################\n\n"
     << "   particle type             : " << particle_type.c_str() << "\n"
     << "   particle mesh type        : " << particle_mesh_type.c_str() << "\n"
     << "   minimum mesh size of particle surface: hmins  : " << hmins << "\n"
     << "   maximum mesh size of particle surface: hmaxs = " << hmaxs << "\n"
     << "------------> The non-dimensional variables:\n"
     << "   non-dimensional bead radius      a0     : " << 1.0 << "\n"
     << "   non-dimensional ksi = sqrt(PI)/(3a0)    : " << std::sqrt(PI) / 3.;
  PMToolBox::output_message(ss, *comm_in);
  pfilename.str("");
  pfilename.clear();
  comm_in->barrier();
} // end function

// ==========================================================================
void CopssRigidParticleSystem::attach_mesh_spring_network()
{
  mesh_spring_network.resize(num_rigid_particles);

  // std::vector<MeshSpringNetwork*> (particle_mesh.num_particles());
  for (std::size_t i = 0; i < num_rigid_particles; ++i)
  {
    MeshBase& p_mesh       = particle_mesh->particles()[i]->mesh();
    const Point& centroid0 = particle_mesh->particles()[i]->get_centroid0();
    mesh_spring_network[i] = new MeshSpringNetwork(p_mesh, *pm_periodic_boundary);
    mesh_spring_network[i]->build_spring_network(centroid0);
    particle_mesh->particles()[i]->attach_mesh_spring_network(mesh_spring_network[
                                                                i]);
  }
}

// =====================================================================
void CopssRigidParticleSystem::create_object_mesh() {
  // prepare domain and objects
  PMToolBox::output_message("==>(1/4) Generate/Create domain Mesh", *comm_in);
  this->create_domain_mesh();
  PMToolBox::output_message("==>(2/4) Create periodic box", *comm_in);
  this->create_periodic_boundary();
  PMToolBox::output_message("==>(3/4) Create particle mesh object", *comm_in);
  this->create_object();
  PMToolBox::output_message("==>(4/4) Create point_mesh object", *comm_in);

  // Create object mesh
  point_mesh = new PointMesh<3>(*particle_mesh, search_radius_p,
    search_radius_e);

  // No need to add periodic boundary, which is already included in
  // particle_mesh
  // Reinit point_mesh
  point_mesh->reinit(neighbor_list_update_flag);

  // finish point_mesh, print information
  ss << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
     << "### The particle-mesh and point-mesh info:\n"
     << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
     << "Total number of particles: " << particle_mesh->num_particles() << "\n"
     << "Total number of points: " << particle_mesh->num_mesh_points() << "\n"
     << "search_radius_p = " << search_radius_p << "\n"
     << "search_radius_e = " << search_radius_e;
  PMToolBox::output_message(ss, *comm_in);
  ss.clear();
  // iterate over all particles for debug purpose
  if (debug_info) {
    const bool print_neighbor_list = false;

    for (std::size_t i = 0; i < particle_mesh->num_particles(); ++i) {
      particle_mesh->particles()[i]->print_info();
    }
  }
} // end function

// ==================================================================================
void CopssRigidParticleSystem::attach_object_mesh(PMLinearImplicitSystem& system)
{
  system.attach_particle_mesh(particle_mesh);
  system.attach_point_mesh(point_mesh);
}

// ======================================================================================
void CopssRigidParticleSystem::set_parameters(EquationSystems& equation_systems)
{
  equation_systems.parameters.set<Real>("real_time") = real_time;
  equation_systems.parameters.set<bool>("adaptive_dt") = adaptive_dt;
  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations")
    = max_linear_iterations;
  equation_systems.parameters.set<Real>("linear solver rtol")
    = linear_solver_rtol;
  equation_systems.parameters.set<Real>("linear solver atol")
    = linear_solver_atol;
  equation_systems.parameters.set<bool>("user_defined_pc")
    = user_defined_pc;
  equation_systems.parameters.set<bool>("schur_user_ksp")
    = schur_user_ksp;
  equation_systems.parameters.set<Real>(
    "schur_user_ksp_rtol") = schur_user_ksp_rtol;
  equation_systems.parameters.set<Real>(
    "schur_user_ksp_atol") = schur_user_ksp_atol;
  equation_systems.parameters.set<string>(          "schur_pc_type")
    = schur_pc_type;
  equation_systems.parameters.set<SystemSolverType>("solver_type_stokes")
    = solver_type_stokes;
  equation_systems.parameters.set<Real>("alpha")
    = alpha;
  equation_systems.parameters.set<Real>("kBT")
    = kBT;
  equation_systems.parameters.set<Real>("ibm_beta")
    = ibm_beta;
  equation_systems.parameters.set<Real>(
    "minimum fluid mesh size") = hminf;
  equation_systems.parameters.set<Real>(
    "minimum solid mesh size") = hmins;
  equation_systems.parameters.set<Real>("minimum mesh size")
    = hmin;
  equation_systems.parameters.set<Real>("viscosity_0")
    = muc;
  equation_systems.parameters.set<Real>("br0")
    = 1.0;
  equation_systems.parameters.set<Real>("bead radius")
    = Rb;
  equation_systems.parameters.set<Real>("drag")
    = drag_c;
  equation_systems.parameters.set<Real>("tc")
    = tc;
  equation_systems.parameters.set<string>("particle_type")
    = particle_type;
  equation_systems.parameters.set<string>("particle_mesh_type")
    = particle_mesh_type;
  equation_systems.parameters.set<std::vector<string> >("force_types")
    = forceTypes;

  for (int i = 0; i < numForceTypes;
       i++) equation_systems.parameters.set<std::vector<Real> >(forces[i].first) =
      forces[i].second;
  equation_systems.parameters.set<string>("simulation_name")
    = simulation_name;
  equation_systems.parameters.set<string>("wall_type")
    = wall_type;
  equation_systems.parameters.set<std::vector<Real> >(wall_type)
    = wall_params;
  equation_systems.parameters.set<std::vector<bool> >("shear")
    = shear;
  equation_systems.parameters.set<std::vector<Real> >("shear_rate")
    = shear_rate;
  equation_systems.parameters.set<std::vector<unsigned int> >("shear_direction")
    = shear_direction;
  equation_systems.parameters.set<bool> ("with_hi") = with_hi;
  // parameters for modules
  equation_systems.parameters.set<bool>("module_poisson") = module_poisson;
  equation_systems.parameters.set<bool>("module_np") = module_np;
  // parameters of Poisson system
  if (module_poisson)
  {
    equation_systems.parameters.set<SystemSolverType>("solver_type_poisson")
            = solver_type_poisson;
    equation_systems.parameters.set<Real>("phi0") = phi0;
    equation_systems.parameters.set<Real>("epsilon") = epsilon;
    equation_systems.parameters.set<std::vector<unsigned int> >(
            "boundary_id_dirichlet_poisson") = boundary_id_dirichlet_poisson;
    equation_systems.parameters.set<std::vector<unsigned int> >(
            "boundary_id_neumann_poisson") = boundary_id_neumann_poisson;
    equation_systems.parameters.set<std::vector<Real> >(
            "boundary_value_dirichlet_poisson") = boundary_value_dirichlet_poisson;
    equation_systems.parameters.set<std::vector<Real> >(
            "boundary_value_neumann_poisson") = boundary_value_neumann_poisson;
  }
  // parameters of NP system
  if (module_np)
  {
    equation_systems.parameters.set<SystemSolverType>("solver_type_np")
            = solver_type_np;
    equation_systems.parameters.set<Real>("c0")
            = c0;
    equation_systems.parameters.set<std::vector<std::string>>("ion_name")
            = ion_name;
    equation_systems.parameters.set<std::vector<Real>>("ion_diffusivity")
            = ion_diffusivity;
    equation_systems.parameters.set<std::vector<int>>("ion_valence")
            = ion_valence;
    equation_systems.parameters.set<Real>("bjerrum_length") = lambda_B;
    equation_systems.parameters.set<Real>("np_system_relaxation_time") =
      np_system_relaxation_time;
    equation_systems.parameters.set<unsigned int>
      ("np_system_relaxation_write_interval") =np_system_relaxation_write_interval;
    equation_systems.parameters.set<std::vector<unsigned int> >(
            "boundary_id_dirichlet_np") = boundary_id_dirichlet_poisson;
    equation_systems.parameters.set<std::vector<std::vector<Real>> >(
      "boundary_value_dirichlet_np") = boundary_value_dirichlet_np;
  }
  equation_systems.parameters.set<int> ("o_precision") = o_precision;
}

void CopssRigidParticleSystem::write_object(unsigned int step_id)
{
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Allgather the distributed vector ROUT to local vector lvec on all
       processors
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       */

  // std::vector<Real> lvec;
  // brownian_sys->vector_transform(lvec,&ROUT, "backward"); // ROUT -> lvec
  point_mesh->update_particle_mesh(particle_mesh);
  particle_mesh->write_particle(step_id,
                                o_step,
                                real_time,
                                output_file,
                                comm_in->rank());
}

void CopssRigidParticleSystem::run(EquationSystems& equation_systems) {
  PerfLog perf_log("Copss-Hydrodynamics-RigidParticleSystem");

  PMToolBox::output_message(
    "=====================4. Start moving particles ==========================\n",
    *comm_in);
    
  // get stokes system from equation systems
  PMSystemStokes& system = equation_systems.get_system<PMSystemStokes>("Stokes");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     Parameters for dynamic process
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  NP    = point_mesh->num_particles();
  n_vec = dim * NP;
  hmin  = std::min(hmins, hminf);
  hmax  = std::max(hmaxs, hmaxf);

  // only when with_hi==true and with_brownian ==false, we can use a larger time
  // step
  if (!with_brownian and with_hi) {
    for (int i = 1; i < max_dr_coeff.size(); i++) max_dr_coeff[i] *= hmin;  // only
                                                                            // magnify
                                                                            // max_dr_coeff[1]
                                                                            // and
                                                                            // max_dr_coeff[2]
  }
  // attach max_dr_coeff to equation systems for future use
  equation_systems.parameters.set<std::vector<Real> >("max_dr_coeff") = max_dr_coeff;

  if (update_neighbor_list_everyStep) {
    ss << "====> neighbor_list is updated at every time step (including half step of fixman if available)";
  }
  else {
    neighbor_list_update_interval =
      (real_time <
       max_dr_coeff[0]) ? int(search_radius_p / 2. /
                              max_dr_coeff[1]) : int(search_radius_p / 2. /
                                                     max_dr_coeff[2]);
    ss << "====> neighbor_list is updated every " << neighbor_list_update_interval << " steps\n"
       << "Warning: be careful of using this option. Although the difference"
       << "between results from updating neighborList every some steps and from"
       << "updating neighborList at each step seems tiny, but we have not fully validated it.\n";
  }
  PMToolBox::output_message(ss, *comm_in);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     Compute undisturbed velocity field without particles.
     NOTE: We MUST re-init particle-mesh before solving Stokes
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  PMToolBox::output_message("==>(1/3) Solve the undisturbed system", *comm_in);
  this->solve_undisturbed_system(equation_systems);

  PMToolBox::output_message(
    "==>(2/3) Prepare RIN & ROUT and Brownian_system in binary format at step 0",
    *comm_in);
  this->create_brownian_system(equation_systems);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     Advancing in time. Fixman Mid-Point algorithm
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  PMToolBox::output_message(
    "==>(3/3) Start calculating dynamics and advancing time steps",
    *comm_in);
  vel0.resize(n_vec);
  vel1.resize(n_vec);

  // start integration
  perf_log.push("integration");

  for (unsigned int i = istart; i <= istart + nstep; ++i)
  {
    if (with_hi) {
      // integrate particle movement using fixman's mid point scheme
      this->fixman_integrate(equation_systems, i);
    }
    else {
      // integrate particle movement using fixman's mid point scheme
      this->langevin_integrate(equation_systems, i);
    }
  }
  perf_log.pop("integration");

  // destroy objects after integration
  this->destroy();
}
} // end of namespace
