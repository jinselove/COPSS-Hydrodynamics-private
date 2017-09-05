#include "copss_point_particle_system.h"

using std::cout;
using std::endl;
using std::string;


namespace libMesh{
//==========================================================================
CopssPointParticleSystem::CopssPointParticleSystem(CopssInit& init)
:Copss(init)
{

}

CopssPointParticleSystem::~CopssPointParticleSystem(){
	delete polymer_chain;
	delete point_mesh;
	polymer_chain = NULL;
	point_mesh = NULL;
}


//==========================================================================
void CopssPointParticleSystem::read_particle_info(){
	if (particle_type != "point_particle"){
		error_msg = "invalid particle type ("+particle_type+") defined\n";
		PMToolBox::output_message(error_msg,comm_in);
		libmesh_error();
	}
	point_particle_model	= input_file("point_particle_model", "other");
	if(point_particle_model == "bead"){
	}
	else if(point_particle_model == "polymer_chain"){
		bk = input_file("bk", 1E-6);//Kuhn length(um)
		Nks = input_file("Nks", 1E-6);// # of Kuhn length per spring
		Ss2 = Nks*bk*bk/6.; // (um^2)
		q0 = Nks * bk;// maximum spring length (um)
    max_spring_len = q0/Rb;     // non-dimensional max spring length
	}
	else{
		error_msg = "	Invalid point_particle_model !!!";
		PMToolBox::output_message(error_msg, comm_in);
		libmesh_error();
	}
}// end read_particle_parameter()


//==========================================================================
void CopssPointParticleSystem::create_object(){
  const unsigned int chain_id = 0;
  polymer_chain = new PolymerChain(chain_id, *pm_periodic_boundary);
  std::ostringstream pfilename;
  if(restart)
  {
    if(point_particle_model == "polymer_chain"){
    	pfilename << "output_polymer_"<< restart_step << ".vtk";
      polymer_chain->read_data_vtk(pfilename.str());
    }
    else if (point_particle_model == "bead"){
      pfilename << "output_bead_" << restart_step << ".csv";
      polymer_chain->read_data_csv(pfilename.str());
    }
    cout <<"-------------> read "<< point_particle_model << "data from " << pfilename.str() << " in restart mode" << endl;

  } 
  else
  {
  	pfilename << "point_particle_data.in";
    cout<<"--------------> skip generating datafile, will read in existed pizza file: "<<pfilename.str()<<endl;
    polymer_chain->read_data_pizza(pfilename.str());
    cout<<"--------------> Polymer_chain class is built!\n";
  }
  // output data read from file
  if(point_particle_model == "bead"){
    Nb = polymer_chain -> n_beads();
    Ns = Nb - 1;
  }
  else if (point_particle_model == "polymer_chain"){
    Nb = polymer_chain -> n_beads();
    nChains = polymer_chain -> n_chains();
    nBonds = polymer_chain -> n_bonds();
    Ns = nBonds / nChains;
    chain_length = Ns * q0; // contour length of the spring (um)
    Dc = Db / Real(Nb); // Diffusivity of the chain (um^2/s)
  }
  // for particular models
  cout<<"##########################################################\n"
           <<"#                  Particle Parameters                    \n"
           <<"##########################################################\n\n"
           <<"   particle type             : " << particle_type.c_str() << endl
           <<"   point Particle model      : " << point_particle_model.c_str() << endl
           <<"   number of point particles Nb = " << Nb << endl;
  if(point_particle_model == "polymer_chain"){
  cout<<"   number of springs per Chain       Ns  = " << Ns << endl
           <<"   number of Chains              nChains = " << nChains << endl
           <<"   Kuhn length                      bk  = " <<bk <<" (um)\n"
           <<"   # of Kuhn segment per spring      Nks = " << Nks << "\n"
           <<"   second moment of polymer chain    Ss2 = " << Ss2 << " (um^2)\n"
           <<"   maximum spring length             q0  = " << q0  << " (um)\n"
           <<"   chain length of polymer           Lc  = " << chain_length <<" (um)\n"
           <<"   chain diffusivity                 Dc  = " << Dc <<" (um^2/s)\n";
  }

  cout << "------------> The non-dimensional variables:\n";


  cout<<"   non-dimensional bead radius      a0     = " << 1.0 << "\n"
           <<"   non-dimensional ksi = sqrt(PI)/(3a0)    = " << std::sqrt(PI) / 3. <<"\n";
  if(point_particle_model == "polymer_chain"){
    cout <<"   non-dimensional Kuhn length    bk/a     = " <<bk/Rb <<"\n"
              <<"   non-dimensional spring length  q0/a     = " <<q0/Rb <<"\n"
              <<"   non-dimensional contour length Lc/a     = " <<chain_length/Rb <<"\n"
              <<"   non-dimensional Ss/a = sqrt(Ss2/a^2)    = " <<std::sqrt(Ss2/Rb/Rb) <<"\n"
              <<"   non-dimensional ksi = sqrt(PI)/(3a0)    = " <<std::sqrt(PI)/(3.) <<"\n";
  }
  pfilename.str(""); pfilename.clear();
  comm_in.barrier();
}//end function

//=====================================================================
void CopssPointParticleSystem::create_object_mesh(){
  // prepare domain and objects
  cout << "\n==>(1/4) Generate/Create domain Mesh\n";
  this -> create_domain_mesh();
  cout << "\n==>(2/4) Create periodic box \n ";
  this -> create_periodic_boundary();
  cout << "\n==>(3/4) Create polymer chain object (for beads or polymer_chain) \n";
  this -> create_object();

  cout << "\n==>(4/4) Create point_mesh object \n";
  // create object mesh
  search_radius_p = 4.0/alpha;
  search_radius_e = 0.5*max_mesh_size + search_radius_p;

  point_mesh = new PointMesh<3> (*mesh, *polymer_chain, search_radius_p, search_radius_e);
  //point_mesh = std::unique_ptr<PointMesh<3> >(new PointMesh<3> (*mesh, *polymer_chain, search_radius_p, search_radius_e));

  point_mesh->add_periodic_boundary(*pm_periodic_boundary);

  point_mesh->reinit();
  cout <<"-------------> Reinit point mesh object, finished! \n"
            <<"- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
		        <<"### The point-mesh info:\n"
		        <<"- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
		        <<"Total number of point particles:"<<point_mesh->num_particles() << endl
		        <<"search_radius_p = "<< search_radius_p << endl
            <<"search_radius_e = "<< search_radius_e << endl;
  if(print_info == true){
    if (comm_in.rank() == 0) point_mesh->print_point_info();
  } 
} // end function


//==================================================================================
void CopssPointParticleSystem::attach_object_mesh(PMLinearImplicitSystem& system)
{
	system.attach_point_mesh(point_mesh);
}

//======================================================================================
void CopssPointParticleSystem::set_parameters(EquationSystems& equation_systems){
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
  equation_systems.parameters.set<Real>       ("viscosity_0")  = muc;
  equation_systems.parameters.set<Real>               ("br0")  = 1.0;
  equation_systems.parameters.set<Real>               ("bk")   = bk;
  equation_systems.parameters.set<Real>               ("q0")   = q0;
  equation_systems.parameters.set<Real>       ("bead radius")  = Rb;
  equation_systems.parameters.set<Real>              ("drag")  = drag_c;
  equation_systems.parameters.set<Real>                ("tc")  = tc;
  equation_systems.parameters.set<Real>               ("Nks")  = Nks;
  equation_systems.parameters.set<Real>               ("Ss2")  = Ss2;
  equation_systems.parameters.set<string> ("particle_type")  = particle_type;
  equation_systems.parameters.set<string> ("point_particle_model") = point_particle_model;
  // Attach force fields
  equation_systems.parameters.set<std::vector<string>> ("pp_force_types") = pp_force_type;
  for (int i=0; i<num_pp_force; i++) equation_systems.parameters.set<std::vector<Real>> (pp_force[i].first) = pp_force[i].second;
  equation_systems.parameters.set<std::vector<string>> ("pw_force_types") = pw_force_type;
  for (int i=0; i<num_pw_force; i++) equation_systems.parameters.set<std::vector<Real>> (pw_force[i].first) = pw_force[i].second;
  equation_systems.parameters.set<string> ("test_name") = test_name;
  equation_systems.parameters.set<string> ("wall_type") = wall_type;
  equation_systems.parameters.set<std::vector<Real>> (wall_type) = wall_params;
}

void CopssPointParticleSystem::update_object(std::string stage)
{
  if(point_particle_model == "polymer_chain")
  {
    chain_broken = polymer_chain->check_chain(max_spring_len);
    if(chain_broken) {
      output_msg = "   ********** warning: Polymer chain is broken " + stage +"\n"+
                   "   ********** warning: bead position is corrected by scaling the chain length and moving the particle according to periodicity";
      PMToolBox::output_message(output_msg, comm_in);
    }
  }  
}


void CopssPointParticleSystem::write_object(std::size_t step_id)
{
  if(comm_in.rank()==0){
  /*---------------------------------------------------------------------------------------
   * output polymer chain / bead data in the VTK format at step i
   -----------------------------------------------------------------------------------------*/
   if(point_particle_model == "polymer_chain"){
     oss << "output_polymer_" << o_step << ".vtk";
     polymer_chain->write_polymer_chain( oss.str() );
   }
   else{
     oss << "output_bead_" << o_step <<".csv";
     polymer_chain->write_bead(oss.str());
   } // end else
   /*----------------------------------------------------------------------------------------------------
   * Output mean square displacement, radius of gyration, chain stretch, and center of mass at step i
   --------------------------------------------------------------------------------------------------- */
  } // end if comm_in.rank() == 0
  brownian_sys->output_statistics_stepi(out_msd_flag, out_stretch_flag, out_gyration_flag, out_com_flag,
                                             step_id, real_time, center0, ROUT);
  oss.str(""); oss.clear();
  /*----------------------------------------------------------------------------------------------------
   * Write out ROUT for restart mode at step i
   ----------------------------------------------------------------------------------------------------*/
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vector_ROUT.dat",FILE_MODE_WRITE,&viewer);
  VecView(ROUT,viewer);
}

void CopssPointParticleSystem::run(EquationSystems& equation_systems){
  PerfLog perf_log("Copss-Hydrodynamics-PointParticleSystem");
  cout<<endl<<"============================4. Start moving particles ============================"<<endl<<endl;
  // get stokes system from equation systems
  PMLinearImplicitSystem& system = equation_systems.get_system<PMLinearImplicitSystem> ("Stokes");
   
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Parameters for dynamic process
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  NP     = point_mesh->num_particles();
  n_vec  = dim*NP;


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
  if(restart==false)
  {
    if(comm_in.rank()==0){
	     if(point_particle_model == "polymer_chain") polymer_chain->write_polymer_chain("output_polymer_0.vtk");
	     else polymer_chain->write_bead("output_bead_0.csv");
    }
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
  if (adaptive_dt == true and point_particle_model == "polymer_chain") max_dr_coeff = 0.1 * Ss2 / Rb / Rb;


  //start integration

  perf_log.push ("integration");
  for(unsigned int i=istart; i<=iend; ++i)
  {
    cout << "\nStarting Fixman Mid-Point algorithm at step "<< i << endl;
    // integrate particle movement using fixman's mid point scheme
    this -> fixman_integrate(equation_systems, i);  
  } // end step integration
  
  perf_log.pop ("integration");
}

} // end of namespace
