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

CopssPointParticleSystem::~CopssPointParticleSystem() {
  delete polymer_chain;
  delete mesh;
  delete point_mesh;
  delete pm_periodic_boundary;
  delete brownian_sys;
  delete fix_factory;
  polymer_chain        = NULL;
  point_mesh           = NULL;
  pm_periodic_boundary = NULL;
  mesh                 = NULL;
  pm_periodic_boundary = NULL;
  brownian_sys         = NULL;
  fix_factory          = NULL;

  for (int i = 0; i < fixes.size(); i++)
  {
    delete fixes[i];
    fixes[i] = NULL;
  }
  fixes.clear();
}

// ==========================================================================
void CopssPointParticleSystem::read_ggem_info() {
  alpha = input_file("alpha", 0.1);
  std::ostringstream ss;
  ss << "##########################################################\n" 
     << "#                 GGEM information                        \n" 
     << "##########################################################\n"
     << "-----------> the smoothing parameter in GGEM alpha = " << alpha << "\n" 
     << "-----------> recommend meshsize <= " << 1. / (std::sqrt(2)* alpha);
  PMToolBox::output_message(ss.str(), *comm_in); 
}

// ==========================================================================
void CopssPointParticleSystem::read_particle_info() {
  particle_type = input_file("particle_type", "other");

  if (particle_type != "point_particle") {
    error_msg = "invalid particle type (" + particle_type + ") defined\n";
    PMToolBox::output_message(error_msg, *comm_in);
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
    error_msg = "	Invalid point_particle_model !!!";
    PMToolBox::output_message(error_msg, *comm_in);
    libmesh_error();
  }
} // end read_particle_parameter()

// ==========================================================================
void CopssPointParticleSystem::create_object() {
  const unsigned int chain_id = 0;

  polymer_chain = new PolymerChain(chain_id, *pm_periodic_boundary);
  std::ostringstream pfilename;
  pfilename << "point_particle_data.in";
  PMToolBox::output_message(
    "--------------> skip generating datafile, will read in existed data file: "
    + pfilename.str(), *comm_in);
  polymer_chain->read_particles_data(pfilename.str());
  if (restart)
  {
    pfilename.str("");
    pfilename.clear();
    PMToolBox::output_message(
      "--------------> in restart mode, load particle positions from saved restart files",
      *comm_in);
    if (point_particle_model == "polymer_chain") {
      pfilename << "output_polymer_" << o_step << ".vtk";
      polymer_chain->read_particles_data_restart_vtk(pfilename.str());
    }
    else if (point_particle_model == "bead") {
      pfilename << "output_bead_" << o_step << ".csv";
      polymer_chain->read_particles_data_restart_csv(pfilename.str());
    }
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
  std::ostringstream ss;
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
  PMToolBox::output_message(ss.str(), *comm_in);
  comm_in->barrier();
} // end function

// =====================================================================
void CopssPointParticleSystem::create_object_mesh() {
  // prepare domain and objects
  PMToolBox::output_message("==>(1/4) Generate/Create domain Mesh", *comm_in);
  this->create_domain_mesh();
  PMToolBox::output_message("==>(2/4) Create periodic box", *comm_in);
  this->create_periodic_boundary();
  PMToolBox::output_message("==>(3/4) Create polymer chain object \
(for beads or polymer_chain)", *comm_in);
  this->create_object();
  PMToolBox::output_message("==>(4/4) Create point_mesh object", *comm_in);
  point_mesh = new PointMesh<3>(*mesh,
                                *polymer_chain,
                                search_radius_p,
                                search_radius_e);
  point_mesh->add_periodic_boundary(*pm_periodic_boundary);

  // reinit point mesh (including particles and neighbor list)
  point_mesh->reinit(neighbor_list_update_flag, true);
  std::ostringstream ss;
  ss << "-------------> Reinit point mesh object, finished! \n"
     << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n"
     << "### The point-mesh info:\n"
     << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n"
     << "Total number of point particles:" << point_mesh->num_particles() <<"\n"
     << "search_radius_p = " << search_radius_p << "\n"
     << "search_radius_e = " << search_radius_e << "\n";
  PMToolBox::output_message(ss.str(), *comm_in);
} // end function

// ==================================================================================
void CopssPointParticleSystem::attach_object_mesh(PMLinearImplicitSystem& system)
{
  system.attach_point_mesh(point_mesh);
}

// ======================================================================================
void CopssPointParticleSystem::set_parameters(EquationSystems& equation_systems) {
  equation_systems.parameters.set<unsigned int>(
                                        "linear solver maximum iterations") =
    max_linear_iterations;
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
  equation_systems.parameters.set<SystemSolverType>(
    "solver_type_poisson") = solver_type_poisson;
  equation_systems.parameters.set<bool>("module_poisson")
    = module_poisson;
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
  equation_systems.parameters.set<Real>(  "phi0")
    = phi0;
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
  equation_systems.parameters.set<std::vector<unsigned int> >(
    "boundary_id_dirichlet_poisson") = boundary_id_dirichlet_poisson;
  equation_systems.parameters.set<std::vector<unsigned int> >(
    "boundary_id_neumann_poisson") = boundary_id_neumann_poisson;
  equation_systems.parameters.set<std::vector<Real> >(
    "boundary_value_dirichlet_poisson") = boundary_value_dirichlet_poisson;
  equation_systems.parameters.set<std::vector<Real> >(
    "boundary_value_neumann_poisson") = boundary_value_neumann_poisson;
}

void CopssPointParticleSystem::update_object()
{
  if (point_particle_model == "polymer_chain")
  {
    chain_broken = polymer_chain->check_chain(max_spring_len);

    if (chain_broken) {
      output_msg =
        "   ********** warning: Polymer chain is broken ---> bead position is corrected by scaling the chain length and moving the particle according to periodicity";
      PMToolBox::output_message(output_msg, *comm_in);
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

  // validate GGEMStokes if simulation_name = ggem_validation
  if (simulation_name == "ggem_validation") {
    perf_log.push("GGEM validation");
    system.reinit_system(neighbor_list_update_flag, build_elem_neighbor_list);
    system.test_velocity_profile();
    perf_log.pop("GGEM validation");
    return;
  }

  // validate GGEMPoisson if simulation_name = ggem_validation_poisson
  if (simulation_name == "ggem_validation_poisson") {
    perf_log.push("GGEMPoisson validation");
    PMSystemPoisson& system_poisson =
      equation_systems.get_system<PMSystemPoisson>("Poisson");
    // Build neighbor list, will this update point_mesh in PMSystemPoisson?
    system.reinit_system(neighbor_list_update_flag, build_elem_neighbor_list);
    system_poisson.test_potential_profile();
    perf_log.pop("GGEMPoisson validation");
    return;
  }

  PMToolBox::output_message(
    "============================4. Start moving particles ====================",
    *comm_in);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     Parameters for dynamic process
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  NP    = point_mesh->num_particles();
  n_vec = dim * NP;
  hmin  = hminf;
  hmax  = hmaxf;

  if (with_brownian == false and with_hi == true) {
    for (int i = 1; i < max_dr_coeff.size(); i++) max_dr_coeff[i] *= hmin;  // only
                                                                            // magnify
                                                                            // max_dr_coeff[1]
                                                                            // and
                                                                            // max_dr_coeff[2]
  }

  if (update_neighbor_list_everyStep) {
    PMToolBox::output_message(
      "====> neighbor_list is updated at every time step (including half step of fixman if available)\n",
      *comm_in);
  }
  else {
    neighbor_list_update_interval =
      (real_time <
       max_dr_coeff[0]) ? int(search_radius_p / 2. /
                              max_dr_coeff[1]) : int(search_radius_p / 2. /
                                                     max_dr_coeff[2]);
    std::ostringstream ss;
    ss << "====> neighbor_list is updated every " << neighbor_list_update_interval << " steps\n"
       << "Warning: be careful of using this option. Although the difference between results from updating neighborList every some steps and from"
       << "updating neighborList at each step seems tiny, but we have not fully validated it.\n";
    PMToolBox::output_message(ss.str(), *comm_in);
  }

  // Get a better conformation of polymer chains before simulation.
  this->update_object();

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     Compute undisturbed velocity field without particles.
     NOTE: We MUST re-init particle-mesh before solving Stokes
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  PMToolBox::output_message("==>(1/3) Solve the undisturbed system", *comm_in);
  perf_log.push("solve_undisturbed_system");
  this->solve_undisturbed_system(equation_systems);
  perf_log.pop("solve_undisturbed_system");

  // create Brownian system for simulation
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

  if (adaptive_dt == true and point_particle_model == "polymer_chain") {
    for (int i = 0; i < max_dr_coeff.size();
         i++) max_dr_coeff[i] *= Ss2 / Rb / Rb;
  }

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
