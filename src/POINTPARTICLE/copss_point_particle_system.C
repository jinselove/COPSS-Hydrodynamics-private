#include "copss_point_particle_system.h"

using std::cout;
using std::endl;
using std::string;


namespace libMesh {
// ==========================================================================
CopssPointParticleSystem::CopssPointParticleSystem(CopssInit& init)
  : Copss(init)
{
  // nothing
}

CopssPointParticleSystem::~CopssPointParticleSystem() 
{
  delete polymer_chain; polymer_chain = nullptr;
}

// ==========================================================================
void CopssPointParticleSystem::read_ggem_info() {
  alpha = input_file("alpha", 0.1);
  ss << "##########################################################\n" 
     << "#                 GGEM information                        \n" 
     << "##########################################################\n"
     << "-----------> the smoothing parameter in GGEM alpha = " << alpha << "\n" 
     << "-----------> recommend meshsize <= " << 1. / (std::sqrt(2)* alpha);
  PMToolBox::output_message(ss, *comm_in); 
}

// ==========================================================================
void CopssPointParticleSystem::read_particle_info() {
  particle_type = input_file("particle_type", "other");

  if (particle_type != "point_particle") {
    ss << "invalid particle type (" + particle_type + ") defined\n";
    PMToolBox::output_message(ss, *comm_in);
    libmesh_error();
  }
  point_particle_model = input_file("point_particle_model", "other");

  if (point_particle_model == "bead") {}
  else if (point_particle_model == "polymer_chain") {
    bk             = input_file("bk", 1E-6);  // Kuhn length(um)
    Nks            = input_file("Nks", 1E-6); // # of Kuhn length per spring
    Ss2            = Nks * bk * bk / 6.;      // (um^2)
    q0             = Nks * bk;                // maximum spring length (um)
    max_spring_len = q0 / Rb;                 // non-dimensional max spring
                                              // length
  }
  else {
    ss << "	Invalid point_particle_model !!!";
    PMToolBox::output_message(ss, *comm_in);
    libmesh_error();
  }
} // end read_particle_parameter()

// ==========================================================================
void CopssPointParticleSystem::create_object() {
  const unsigned int chain_id = 0;

  polymer_chain = new PolymerChain(chain_id, *pm_periodic_boundary);
  std::ostringstream pfilename;
  pfilename << "point_particle_data.in";
  ss << "--------------> skip generating datafile, will read in existed data file: "
     << pfilename.str();
  PMToolBox::output_message(ss, *comm_in);
  polymer_chain->read_particles_data(pfilename.str());
  if (restart)
  {
    pfilename.str("");
    pfilename.clear();
    if (point_particle_model == "polymer_chain") {
      pfilename << "output_polymer_" << o_step << ".vtk";
      polymer_chain->read_particles_data_restart_vtk(pfilename.str());
    }
    else if (point_particle_model == "bead") {
      pfilename << "output_bead_" << o_step << ".csv";
      polymer_chain->read_particles_data_restart_csv(pfilename.str());
    }
    ss << "##################### Restart mode ##########################\n"
       << "---> read point particle data from " << pfilename.str() << "\n"
       << "#############################################################";
    PMToolBox::output_message(ss, *comm_in);
  }

  // output data read from file
  if (point_particle_model == "bead") {
    Nb = polymer_chain->n_beads();
    Ns = Nb - 1;
    polymer_chain->initial_bead_center_of_mass(center0);
  }
  else if (point_particle_model == "polymer_chain") {
    Nb           = polymer_chain->n_beads();
    nChains      = polymer_chain->n_chains();
    nBonds       = polymer_chain->n_bonds();
    Ns           = nBonds / nChains;
    chain_length = Ns * q0;       // contour length of the spring (um)
    Dc           = Db / Real(Nb); // Diffusivity of the chain (um^2/s)
    polymer_chain->initial_chain_center_of_mass(center0);
  }
  // for particular models
  ss << "##########################################################\n"
     << "#                  Particle Parameters                    \n"
     << "##########################################################\n"
     << "   particle type             : " << particle_type.c_str() << "\n"
     << "   point Particle model      : " << point_particle_model.c_str() << "\n"
     << "   number of point particles Nb = " << Nb << "\n";
  if (point_particle_model == "polymer_chain") {
    ss << "   number of springs per Chain       Ns  = " << Ns << endl
       << "   number of Chains              nChains = " << nChains << endl
       << "   Kuhn length                      bk  = " << bk << " (um)\n"
       << "   # of Kuhn segment per spring      Nks = " << Nks << "\n"
       << "   second moment of polymer chain    Ss2 = " << Ss2 << " (um^2)\n"
       << "   maximum spring length             q0  = " << q0  << " (um)\n"
       << "   chain length of polymer           Lc  = " << chain_length <<" (um)\n"
       << "   chain diffusivity                 Dc  = " << Dc << " (um^2/s)\n";
  }
  ss << "------------> The non-dimensional variables:\n"
     << "   non-dimensional bead radius      a0     = " << 1.0 << "\n"
     << "   non-dimensional ksi = sqrt(PI)/(3a0)    = " << std::sqrt(PI) / 3. << "\n";
  if (point_particle_model == "polymer_chain") {
    ss << "   non-dimensional Kuhn length    bk/a     = " << bk / Rb << "\n"
       << "   non-dimensional spring length  q0/a     = " << q0 / Rb << "\n"
       << "   non-dimensional contour length Lc/a     = " << chain_length / Rb << "\n"
       << "   non-dimensional Ss/a = sqrt(Ss2/a^2)    = " << std::sqrt(Ss2 / Rb / Rb) << "\n"
       << "   non-dimensional ksi = sqrt(PI)/(3a0)    = " << std::sqrt(PI) / (3.) << "\n";
  }
  pfilename.str(""); pfilename.clear();
  PMToolBox::output_message(ss, *comm_in);
  comm_in->barrier();
} // end function

// =====================================================================
void CopssPointParticleSystem::create_object_mesh() {
  // prepare domain and objects
  ss << "==>(1/4) Generate/Create domain Mesh";
  PMToolBox::output_message(ss, *comm_in);
  this->create_domain_mesh();
  ss << "==>(2/4) Create periodic box";
  PMToolBox::output_message(ss, *comm_in);
  this->create_periodic_boundary();
  ss << "==>(3/4) Create polymer chain object (for beads or polymer_chain)";
  PMToolBox::output_message(ss, *comm_in);
  this->create_object();
  ss << "==>(4/4) Create point_mesh object";
  PMToolBox::output_message(ss, *comm_in);
  // PMToolBox::output_message("initialize pointmesh class", *comm_in);
  point_mesh = new PointMesh<3>(*mesh,
                                *polymer_chain,
                                search_radius_p,
                                search_radius_e);
  // PMToolBox::output_message("add pbc", *comm_in);
  point_mesh->add_periodic_boundary(*pm_periodic_boundary);

  // reinit point mesh (including particles and neighbor list)
  // PMToolBox::output_message("before reinit point_mesh", *comm_in);
  point_mesh->reinit(neighbor_list_update_flag);
  // PMToolBox::output_message("after reinit point_mesh", *comm_in);
  ss << "-------------> Reinit point mesh object, finished! \n"
     << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n"
     << "### The point-mesh info:\n"
     << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n"
     << "Total number of point particles:" << point_mesh->num_particles() <<"\n"
     << "search_radius_p = " << search_radius_p << "\n"
     << "search_radius_e = " << search_radius_e << "\n";
  PMToolBox::output_message(ss, *comm_in);
} // end function

// ==================================================================================
void CopssPointParticleSystem::attach_object_mesh(PMLinearImplicitSystem& system)
{
  system.attach_point_mesh(point_mesh);
}

// ======================================================================================
void CopssPointParticleSystem::set_parameters(EquationSystems& equation_systems)
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
  equation_systems.parameters.set<Real>(
    "minimum fluid mesh size") = hminf;
  equation_systems.parameters.set<Real>(  "viscosity_0")
    = muc;
  equation_systems.parameters.set<Real>(  "br0")
    = 1.0;
  equation_systems.parameters.set<Real>(  "bk")
    = bk;
  equation_systems.parameters.set<Real>(  "q0")
    = q0;
  equation_systems.parameters.set<Real>(  "bead radius")
    = Rb;
  equation_systems.parameters.set<Real>(  "drag")
    = drag_c;
  equation_systems.parameters.set<Real>(  "tc")
    = tc;
  equation_systems.parameters.set<Real>(  "Nks")
    = Nks;
  equation_systems.parameters.set<Real>(  "Ss2")
    = Ss2;

  equation_systems.parameters.set<string>("particle_type")
    = particle_type;
  equation_systems.parameters.set<string>(
    "point_particle_model") = point_particle_model;
  equation_systems.parameters.set<std::vector<string> >("force_types")
    = forceTypes;

  for (int i = 0; i < numForceTypes;
       i++) equation_systems.parameters.set<std::vector<Real> >(forces[i].first) =
      forces[i].second;
  equation_systems.parameters.set<string>(                    "simulation_name")
    = simulation_name;
  equation_systems.parameters.set<string>(                    "wall_type")
    = wall_type;
  equation_systems.parameters.set<std::vector<Real> >(        wall_type)
    = wall_params;
  equation_systems.parameters.set<std::vector<bool> >(        "shear")
    = shear;
  equation_systems.parameters.set<std::vector<Real> >(        "shear_rate")
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
    // calculate minimum dt required by all NP systems
    // this dt will be updated during integration process if dt required by
    // other systems are even smaller
    unsigned int n_sys = equation_systems.n_systems();
    std::vector<std::string> np_sys_names;
    for (unsigned int s_id=0; s_id<n_sys; s_id++)
    {
      const std::string& s_name = equation_systems.get_system(s_id).name();
      if (s_name.rfind("NP:", 0)==0)
      {
        np_sys_names.push_back(s_name);
      }
    }
    equation_systems.parameters.set<Real>("dt") =
      this->get_min_dt(equation_systems, np_sys_names);
    PMToolBox::output_message(std::string("dt required by NP systems = ")
    +std::to_string(equation_systems.parameters.set<Real>("dt")), *comm_in);
    equation_systems.parameters.set<std::vector<unsigned int> >(
            "boundary_id_dirichlet_np") = boundary_id_dirichlet_np;
    equation_systems.parameters.set<std::vector<std::vector<Real>> >(
            "boundary_value_dirichlet_np") = boundary_value_dirichlet_np;
  }
  equation_systems.parameters.set<int> ("o_precision") = o_precision;
}

void CopssPointParticleSystem::update_object()
{
  if (point_particle_model == "polymer_chain")
  {
    chain_broken = polymer_chain->check_chain(max_spring_len);

    if (chain_broken) {
      ss <<"   ********** warning: Polymer chain is broken ---> bead position is corrected by scaling the chain length and moving the particle according to periodicity";
      PMToolBox::output_message(ss, *comm_in);
    }
  }
}

void CopssPointParticleSystem::write_object(unsigned int step_id)
{
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Allgather the distributed vector ROUT to local vector lvec on all
       processors
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       */
  std::vector<Real> lvec;

  brownian_sys->vector_transform(lvec, &ROUT, "backward"); // ROUT -> lvec

  // write output file
  if (point_particle_model == "polymer_chain") {
    polymer_chain->write_polymer_chain(step_id,
                                       o_step,
                                       real_time,
                                       center0,
                                       lvec,
                                       output_file,
                                       comm_in->rank());
  }
  else {
    polymer_chain->write_bead(step_id,
                              o_step,
                              real_time,
                              center0,
                              lvec,
                              output_file,
                              comm_in->rank());
  } // end else
}

// ============================================================================
void CopssPointParticleSystem::run(EquationSystems& equation_systems) {
  PerfLog perf_log("Copss-Hydrodynamics-PointParticleSystem");

  // get stokes system from equation systems
  PMSystemStokes& system = equation_systems.get_system<PMSystemStokes>("Stokes");

  // validate Stokes against analytic solution
  if (simulation_name == "ggem_validation")
  {
    perf_log.push("GGEM validation");
    system.reinit_system(neighbor_list_update_flag, "disturbed");
    system.test_velocity_profile();
    system.test_nodal_error();
    perf_log.pop("GGEM validation");
    return;
  }

  // validate Poisson against analytic solution
  if (simulation_name=="ggem_validation_poisson")
  {
    perf_log.push("GGEMPoisson validation");
    system.reinit_system(neighbor_list_update_flag, "disturbed");
    PMSystemPoisson& system_poisson =
      equation_systems.get_system<PMSystemPoisson>("Poisson");
    system_poisson.test_potential_profile();
    system_poisson.test_nodal_error();
    perf_log.pop("GGEMPoisson validation");
    return;
  }

  // validate NP system (diffusion only) with analytic solution
  if (simulation_name=="np_validation_analytic")
  {
    perf_log.push("NP (diffusion only) validation");

    // test concentration profile of all NP systems
    unsigned int n_sys = equation_systems.n_systems();
    for (unsigned int s_id=0; s_id<n_sys; s_id++)
    {
      const std::string& s_name = equation_systems.get_system(s_id).name();
      if (s_name.rfind("NP:", 0)==0)
      {
        PMSystemNP& np_system = equation_systems.get_system<PMSystemNP>(s_name);
        np_system.test_concentration_profile();
      }
    }

    perf_log.pop("NP (diffusion only) validation");
    return;
  }


  ss << "============================4. Start moving particles ====================";
  PMToolBox::output_message(ss, *comm_in);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     Parameters for dynamic process
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  NP    = point_mesh->num_particles();
  n_vec = dim * NP;
  hmin  = hminf;
  hmax  = hmaxf;

  if (!with_brownian and with_hi) {
    for (int i = 1; i < max_dr_coeff.size(); i++) max_dr_coeff[i] *= hmin;  // only
                                                                            // magnify
                                                                            // max_dr_coeff[1]
                                                                            // and
                                                                            // max_dr_coeff[2]
  }

  if (update_neighbor_list_everyStep) {
    ss << "====> neighbor_list is updated at every time step (including half step of fixman if available)\n";
    PMToolBox::output_message(ss, *comm_in);
  }
  else {
    neighbor_list_update_interval =
      (real_time <
       max_dr_coeff[0]) ? int(search_radius_p / 2. /
                              max_dr_coeff[1]) : int(search_radius_p / 2. /
                                                     max_dr_coeff[2]);
    ss << "====> neighbor_list is updated every " << neighbor_list_update_interval << " steps\n"
       << "Warning: be careful of using this option. Although the difference between results from updating neighborList every some steps and from"
       << "updating neighborList at each step seems tiny, but we have not fully validated it.\n";
    PMToolBox::output_message(ss, *comm_in);
  }

  // Get a better conformation of polymer chains before simulation.
  this->update_object();

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     Compute undisturbed velocity field without particles.
     NOTE: We MUST re-init particle-mesh before solving Stokes
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  ss << "==>(1/3) Solve the undisturbed system";
  PMToolBox::output_message(ss, *comm_in);
  perf_log.push("solve_undisturbed_system");
  this->solve_undisturbed_system(equation_systems);
  perf_log.pop("solve_undisturbed_system");

  // create Brownian system for simulation
  ss << "==>(2/3) Prepare RIN & ROUT and Brownian_system in binary format at step 0";
  PMToolBox::output_message(ss, *comm_in);
  this->create_brownian_system(equation_systems);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     Advancing in time. Fixman Mid-Point algorithm
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  ss << "==>(3/3) Start calculating dynamics and advancing time steps";
  PMToolBox::output_message(ss, *comm_in);
  vel0.resize(n_vec);
  vel1.resize(n_vec);

  if (adaptive_dt and point_particle_model == "polymer_chain") {
    for (int i = 0; i < max_dr_coeff.size();
         i++) max_dr_coeff[i] *= Ss2 / Rb / Rb;
  }
  // attach max_dr_coeff to equation systems for future use
  equation_systems.parameters.set<std::vector<Real> >("max_dr_coeff") = max_dr_coeff;

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
  } // end step integration
  perf_log.pop("integration");

  // destroy objects after integration
  this->destroy();
}
} // end of namespace
