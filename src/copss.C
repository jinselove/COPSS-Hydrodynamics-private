// Parallel Finite Element-General Geometry Ewald-like Method.
// Copyright (C) 2015-2016 Xujun Zhao, Jiyuan Li, Xikai Jiang

// This code is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.


// This code is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.


// You should have received a copy of the GNU General Public
// License along with this code; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#include "copss.h"
using std::cout;
using std::endl;

namespace libMesh
{

//==========================================================================
Copss::Copss(const CopssInit& init)
{
  comm_in = init.comm();
  cout << endl <<"============================0.  Initialize libMesh.===========================" << endl;
  this -> check_libmesh();
}

Copss::~Copss()
{
}


int Copss::check_libmesh(){

  // This example NaNs with the Eigen sparse linear solvers and Trilinos solvers,
  // but should work OK with either PETSc or Laspack.
  libmesh_example_requires(libMesh::default_solver_package() != EIGEN_SOLVERS,
                           "--enable-petsc or --enable-laspack");
  libmesh_example_requires(libMesh::default_solver_package() != TRILINOS_SOLVERS,
                           "--enable-petsc or --enable-laspack");
  
  // check libmesh enables
  #ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
  #endif
  
  #ifndef LIBMESH_HAVE_SLEPC
  libmesh_example_requires(false, "--enable-slepc");
  #endif
  
  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(3 == LIBMESH_DIM, "--3D support");

  return 0;
}


//==========================================================================
void Copss::start_time(struct tm * timeinfo){
  cout <<"\n---------------------------------------------------------------------\n"
         <<"The program starts. \n"
         <<"The current date/time is: " << asctime (timeinfo)  << "\n"
         <<"---------------------------------------------------------------------\n\n";
}


//==========================================================================
void Copss::end_time(struct tm * timeinfo){
    cout << "\n---------------------------------------------------------------------\n"
              <<"The current date/time is: " <<asctime (timeinfo) << "\n"
              <<"The program ends. \n"
              <<"---------------------------------------------------------------------\n\n";
}

//====================================================================
EquationSystems Copss::init_system(std::string _control_fileName){
   control_fileName = _control_fileName;
   cout <<"\n============================1. read input parameters ============================\n";
   this -> read_input();
   cout <<"\n============================2. Create Point-mesh object =========================\n";
   this -> create_object_mesh();
   cout <<"\n==========3. Create equation_systems object (of type EquationSystems) ===========\n";   
   return this -> create_equation_systems();
}

//====================================================================
void Copss::read_input()
{
  const GetPot tmp(control_fileName);
  input_file = tmp;
  this -> read_system_info();
  this -> read_physical_info();
  this -> read_particle_info();
  this -> read_domain_info();
  this -> read_force_info();
  this -> read_ggem_info();
  this -> read_stokes_solver_info();
  this -> read_chebyshev_info();
  this -> read_run_info();
} // end read_data function

//====================================================================
void Copss::read_system_info()
{
    test_name = input_file("test_name", "validation");
    cout <<"\n##########################################################\n"
         <<  "#                       system_name                       \n"
         <<  "##########################################################\n\n"
         <<  "-----------> system_name: "<< test_name.c_str() << endl;
    print_info = input_file("print_info", false);
}

  /*
   * Read physical parameters 
   */
void Copss::read_physical_info()
{
  T = input_file("temperature", 297);// K
  kBT    = kB * T; //(N*um)
  viscosity            = input_file("viscosity", 1.0); // viscosity (cP = N*s/um^2)
  Rb                   = input_file("radius", 0.10); // radius of the bead (um)
  drag_c      = 6.*PI*viscosity*Rb;    // Drag coefficient (N*s/um)
  Db          = kBT/drag_c;     // diffusivity of a bead (um^2/s)  

  tc   = drag_c*Rb*Rb/kBT;       // diffusion time (s)
  uc   = kBT/(drag_c*Rb);        // characteristic velocity (um/s)
  fc   = kBT/Rb;                 // characteristic force (N)
  muc  = 1./(6.*PI);             // non-dimensional viscosity

  // print out physical parameters information
  cout<<"\n ##########################################################\n"
           <<" #                  System Physical Parameters             \n"
           <<" ##########################################################\n\n"
           <<"   temperature           T   = " << T <<"(K)\n"
           <<"   viscosity             mu  = " << viscosity << " (cP = N*s/um^2)\n"
           <<"   Energy unit           kBT = " << kBT << " (N*um = N*um)\n"
           <<"   Radius of the bead     a  = " << Rb <<" (um)\n"
           <<"   bead diffusivity      Db  = " << Db <<"(um^2/s)\n"
           <<"   HI Drag coefficient  zeta = 6*PI*mu*a = "<<drag_c <<" (N*s/um)\n"
           <<"   ksi = sqrt(PI)/(3a)       =  " <<std::sqrt(PI)/(3. * Rb) << " (1/um)\n"
           <<"   ------------> The characteristic variables:\n"
           <<"   characteristic time          = " << tc <<" (s)\n"
           <<"   characteristic velocity      = " << uc <<" (um/s)\n"
           <<"   characteristic force         = " << fc <<" (N)\n"; 
}// end read_physical_parameter()


  /*
   * Read Geometry infomation
   */

void Copss::read_domain_info()
{
  dim = input_file("dimension", 3);
  //=============== wall type and wall params
  wall_type = input_file("wall_type", "not_defined");
  if (wall_type == "slit" or wall_type == "sphere"){
    wall_params.resize(input_file.vector_variable_size(wall_type));
    for (unsigned int j = 0; j < wall_params.size(); j++) wall_params[j] = input_file(wall_type,0.0,j);
  }
  else{
    std::cout << "*** Error: wall_type ("<< wall_type <<") is not supported !!! (in check_wall)" << std::endl;
    libmesh_error();    
  }
  //=============== periodicity
  periodicity.resize(input_file.vector_variable_size("periodicity"));
  if(periodicity.size() != dim){
    cout <<"Warning: periodicity has to be defined for all dimensions" << endl;
    libmesh_error();
  }
  for (unsigned int i=0; i < periodicity.size(); i++){ periodicity[i] = input_file("periodicity", false, i); }
  if(periodicity[0]==true and periodicity[1] == true and periodicity[2]==true){
    error_msg = "warning: The box cannot be periodic on all directions at the same time. (required by FEM)";
    PMToolBox::output_message(error_msg, comm_in);
    libmesh_error();
  }
  //============== inlet 
  inlet.resize(input_file.vector_variable_size("inlet"));
  if(inlet.size() != dim){
    cout <<"Warning: inlet has to be defined for all dimensions" << endl;
    libmesh_error();
  }
  for (unsigned int i=0; i < inlet.size(); i++){ 
     inlet[i] = input_file("inlet", false, i);
     if(inlet[i]==true and periodicity[i]==true) {
      error_msg = "warning: A inlet direction has to be non-periodic";
      PMToolBox::output_message(error_msg,comm_in);
      libmesh_error();
     }
  }
  //============== inlet pressure
  inlet_pressure.resize(input_file.vector_variable_size("inlet_pressure"));
  for (unsigned int i=0; i < inlet_pressure.size(); i++){ inlet_pressure[i] = input_file("inlet_pressure", 0, i); }

  //============== shear
  std::cout << "\n\nWarning: shear flow has not been implemented \n\n";
  shear.resize(input_file.vector_variable_size("shear"));
  if(shear.size() != dim*2){
    cout <<"Warning: shear has to be defined for all dimensions and both upper wall and bottom wall" << endl;
    libmesh_error();
  }
  for (unsigned int i=0; i < shear.size(); i++){ 
     shear[i] = input_file("shear", false, i);
     if(shear[i]==true and periodicity[int(i/2)]==true) {
      error_msg = "warning: A shear direction has to be non-periodic";
      PMToolBox::output_message(error_msg,comm_in);
      libmesh_error();
     }
  }
  // =============== shear rate
  shear_rate.resize(input_file.vector_variable_size("shear_rate"));
  for (unsigned int i=0; i < shear_rate.size(); i++){ shear_rate[i] = input_file("shear_rate", 0, i); }

  cout <<endl<< "##########################################################"<<endl
       << "#                  Geometry information                   " <<endl
       << "##########################################################"<<endl<<endl
       << "-----------> Dimension: " << dim << endl
       << "-----------> Wall type: " << wall_type << endl;
  cout << "-----------> Wall size parameters: ";
       for (int i = 0; i < wall_params.size(); ++i)
       {
         cout << wall_params[i] <<"   ";
       }
       cout << endl;
  cout << "-----------> Periodicity of the box: ";
  for (int i=0; i<dim; i++) cout << std::boolalpha <<periodicity[i] <<", ";
  cout <<endl;
  cout << "-----------> Inlet/Outlet of the box: ";
  for (int i=0; i<dim; i++) cout << std::boolalpha << inlet[i] <<"(pressure = " <<inlet_pressure[i] <<" ), ";
  cout <<endl;
  cout << "-----------> shear of the box: ";
  for (int i=0; i<dim*2; i++) cout << std::boolalpha << shear[i] <<"(shear_rate = " <<shear_rate[i] <<" ), ";
  cout <<endl;
  cout <<endl<< "##########################################################"<<endl
       << "#                  Domain Mesh information                     " <<endl
       << "##########################################################"<<endl<<endl;
  generate_mesh = input_file("generate_mesh", true); 
  if (generate_mesh){
    n_mesh.resize(input_file.vector_variable_size("n_mesh"));
    for (unsigned int i=0; i < n_mesh.size(); i++){ n_mesh[i] = input_file("n_mesh", 1, i); }
    cout << "------------> Generate Mesh using COPSS: n_mesh = " << n_mesh[0] <<";" << n_mesh[1] <<";" << n_mesh[2] << endl;
  }
  else{
    domain_mesh_file = input_file("domain_mesh_file" , "nothing");
    cout <<"------------> Load mesh file from "<< domain_mesh_file.c_str() << endl;
  } // end else
} // end read_domain_info()

  /*
   * read force types
   */
void Copss::read_force_info(){
  // read particle-particle force types
  numForceTypes = input_file.vector_variable_size("force_field");
  forceTypes.resize(numForceTypes);
  forces.resize(numForceTypes);
  for (unsigned int i=0; i < numForceTypes; i++){    
    forceTypes[i] = input_file("force_field", "nothing", i);
    std::vector<Real> params(input_file.vector_variable_size(forceTypes[i]));
    if(forceTypes[i] != "nothing"){
      for (unsigned int j = 0; j < params.size(); j++){
        params[j] = input_file(forceTypes[i],0.0,j);
      }
    }
    forces[i].first = forceTypes[i];
    forces[i].second = params;
  }

  cout <<endl<< "##########################################################"<<endl
       << "#    Force field (fixes)                " <<endl
       << "##########################################################"<<endl<<endl;
  for (int i = 0; i < numForceTypes; i++){
    cout << "-----------> ";
    cout << forces[i].first <<"= '";
    for (int j = 0; j < forces[i].second.size(); j++){
          cout << forces[i].second[j] <<"  ";       
    }
    cout << "'" <<endl;
  }

} // end read_force_info()

/*
 * read Stokes Solver  
 */
void Copss::read_stokes_solver_info(){
  max_linear_iterations = input_file("max_linear_iterations", 100);
  linear_solver_rtol   = input_file("linear_solver_rtol", 1E-6);
  linear_solver_atol   = input_file("linear_solver_atol", 1E-6);
  user_defined_pc            = input_file("user_defined_pc", true);
  schur_user_ksp       = input_file("schur_user_ksp", false);
  schur_user_ksp_rtol  = input_file("schur_user_ksp_rtol", 1E-6);
  schur_user_ksp_atol  = input_file("schur_user_ksp_atol", 1E-6);
  schur_pc_type = input_file("schur_pc_type", "SMp");
  stokes_solver_type = input_file("stokes_solver", "superLU_dist");
  if(stokes_solver_type=="superLU_dist") {
    solver_type = superLU_dist;
    user_defined_pc = false;
  }
  else if(stokes_solver_type=="field_split") {
    solver_type = field_split;
    user_defined_pc = true;
  }
  else {
    solver_type = user_define;
  }  

  cout <<endl<< "##########################################################"<<endl
       << "#                 Solver information                      " <<endl
       << "##########################################################"<<endl<<endl;
  cout << "-----------> Stokes solver type = " << stokes_solver_type <<endl;
  if (stokes_solver_type=="field_split"){
    cout<<"-----------> FieldSplit Schur Complement Reduction Solver"<<endl;
    cout<<"-----------> schur_pc_type = " << schur_pc_type << endl;
      if(schur_user_ksp){
            cout<<"----------->  user defined KSP is used for Schur Complement!"<< endl;
            cout<<"----------->  KSP rel tolerance for Schur Complement solver is = " << schur_user_ksp_rtol <<endl;
            cout<<"----------->  KSP abs tolerance for Schur Complement solver is = " << schur_user_ksp_atol <<endl;
        }// end if(schur_user_ksp)
    }// end if (stokes_solver_type == "field_split")
}// end read_stokes_solver_info()

/*
 * read Chebyshev info
 */
void Copss::read_chebyshev_info(){
  with_hi        = input_file("with_hi", false);
  with_brownian  = input_file("with_brownian", false);
  restart       = input_file("restart", false); 
  compute_eigen = input_file("compute_eigen", false);
  eig_min = 0.;
  eig_max = 0.;
  if(with_hi and with_brownian and restart and !compute_eigen){
    read_eigen = true;  
    this -> read_restart_eigenvalue();
  }
  max_n_cheb = input_file("max_n_cheb", 10);
  tol_cheb = input_file("tol_cheb", 0.1);
  eig_factor = input_file("eig_factor", 1.05);
  tol_eigen = input_file("tol_eigen", 0.01);
  if(with_brownian and with_hi){
    cout <<endl<< "##########################################################"<<endl
         << "#   Chebyshev information (only for brownian System)            " <<endl
         << "##########################################################"<<endl<<endl;  
    cout << "-----------> compute eigen values  = " <<std::boolalpha << compute_eigen <<endl;
    cout << "-----------> max number of chebyshev polynomial = " <<max_n_cheb <<endl;
    cout << "-----------> tolerance of chebyshev polynomial = " <<tol_cheb <<endl;
    cout << "-----------> factor of eigenvalues range = " <<eig_factor <<endl;
    cout << "-----------> tolerance of eigenvalues convergence = " <<tol_eigen <<endl;
    if (read_eigen){
      cout << "-----------> read eigen value range from 'out.eigenvalue = '"<<std::boolalpha<<read_eigen << endl
           << "-----------> initial eigen value range = ( " << eig_min << ", " << eig_max << " )" << endl;
    }
  }
} // end read_chebyshev_info()


/*
 * read run time info 
 */
void Copss::read_run_info(){
  //############## Without Brownian ###############################
  // For polymer_chain and bead: maximum displacement (non dimensional) of one step = 0.1 * fluid mesh size minimum (hmin)
  //############## With Brownian ##################################
  // For polymer_chain: maximum displacement (non dimensional) of one step = 0.1 * Ss2/Rb/Rb
  // For bead: maximum displacement (non_dimensional) of one step = 0.1
  if(with_brownian){
    random_seed   = input_file("random_seed",111);
  }
  max_dr_coeff.resize(input_file.vector_variable_size("max_dr_coeff"));
  for (unsigned int i=0; i<max_dr_coeff.size();i++){
    max_dr_coeff[i] = input_file("max_dr_coeff", 0.1, i);
  }
  adaptive_dt    = input_file("adaptive_dt", true);  
  if (restart){
    this->read_restart_time();
    istart = restart_step;
  }
  else{
    restart_step = -1;
    istart = 0;
    o_step = 0;
    real_time = 0.;
  }
  // neighbor list update
  neighbor_list_update_flag = true; // set this to be ture before simulation
  update_neighbor_list_everyStep = input_file("update_neighbor_list_everyStep", true);
  // total number of steps 
  nstep = input_file("nstep", 1);
  // free-draining steps when Chebysheve cannot converge even when eigenvalue is recomputed
  n_relax_step = input_file("n_relax_step", 10);
  // write interval
  write_interval = input_file("write_interval", 1);
  //output file flags
  output_file.resize(input_file.vector_variable_size("output_file"));
  for (unsigned int i=0; i < output_file.size(); i++){
    output_file[i] = input_file("output_file","not defined", i);
  }
  // debug info
  debug_info   = input_file("debug_info", false);
  // write to screen
  cout <<"\n##########################################################\n"
       << "#                 Run information                      \n"
       << "##########################################################\n\n"
       << "-----------> adaptive_dt: " << std::boolalpha << adaptive_dt << endl
       << "-----------> max_dr_coeff at step 0(this value will be multiplied by hmin if with_brownian == false and with_hi == true): " << max_dr_coeff[0] << endl
       << "-----------> max_dr_coeff at step i > 0(this value will be multiplied by hmin if with_brownian == false and with_hi == true): " << max_dr_coeff.back() << endl
       << "-----------> debug_info: " << std::boolalpha << debug_info << endl; 
  if (with_hi) cout << "-----------> with_hi: " <<std::boolalpha <<with_hi <<endl;
  if (with_brownian){
	 cout << "-----------> with_brownian: " <<std::boolalpha<<with_brownian <<endl
        << "-----------> random seed: " <<random_seed <<endl;
  }
  cout << "-----------> Restart mode: "<<std::boolalpha << restart <<endl;
  if(restart){
   cout <<"-----------> Restart step: "<<restart_step << endl;
   cout <<"-----------> Restart from real_time: "<<real_time <<endl;
   cout <<"-----------> Restart particle data read from o_step: " << o_step << endl;
  }
  cout << "-----------> nstep: " <<nstep <<endl
       << "-----------> write_interval: " << write_interval << endl
       << "-----------> write output file: " << endl;
  for(int i = 0; i < output_file.size(); i++){
    cout << "                               " << output_file[i] << endl;
  }
} // end read_run_info()

//============================================================================
void Copss::read_restart_time()
{
  // Open the local file and check the existance
  std::cout <<"\n###Read time output: out.time"<<std::endl;
  std::ifstream fin;
  fin.open ("out.time", std::ios_base::ate);
  char c = '\0';
  if( fin.is_open() ){
    const unsigned int length = fin.tellg();
    cout << "length = " << length <<endl;
    for (int i = length-2; i>0;i--){
      fin.seekg(i);
      c = fin.get();
      if(c=='\r' || c=='\n') // new line?
        break;
    }
    fin >> restart_step >> o_step >> real_time;
    fin.close();
  } // end if
  else{
    printf("***warning: read_restart_time can NOT read the time output: out.time");
    libmesh_error();
  }
}

//============================================================================
void Copss::read_restart_eigenvalue()
{
  //Open the local file and check the existance
  std::cout <<"\n###Read eigenvalue from previous simulation to restart: out.eigenvalue" << std::endl;
  std::ifstream fin;
  fin.open ("out.eigenvalue", std::ios_base::ate);
  std::string tmp;
  char c = '\0';
  if( fin.is_open() ){
    const unsigned int length = fin.tellg();
    cout << "length = " << length <<endl;
    for (int i = length-2; i>0;i--){
      fin.seekg(i);
      c = fin.get();
      if(c=='\r' || c=='\n') // new line?
        break;
    }
    fin >> eig_min >> eig_max;
    fin.close();
  } // end if
  else{
    printf("***warning: read_restart_eigenvalue() can NOT read the time output: out.eigenvalue");
    libmesh_error();
  }  
}

//============================================================================
void Copss::create_domain_mesh()
{
  if(dim == 2) {
    error_msg = "Copss::create_mesh() only works for 3D systems; 2D simulation needs extra implementation";
    PMToolBox::output_message(error_msg, comm_in);
    libmesh_error();
  }
  mesh = new SerialMesh(comm_in);
  //mesh = std::unique_ptr<SerialMesh> (new SerialMesh (comm_in));
  if(generate_mesh){
    if (wall_type == "slit"){      
        const Real meshsize_x   = (wall_params[1] - wall_params[0])/Real( n_mesh[0] );
        const Real meshsize_y   = (wall_params[3] - wall_params[2])/Real( n_mesh[1] );
        const Real meshsize_z   = (wall_params[5] - wall_params[4])/Real( n_mesh[2] );
//        cout << "mesh_size = " <<meshsize_x  << "; "<<meshsize_y << "; "<<meshsize_z << endl;      
        hminf  = std::min(meshsize_x, meshsize_y);
        hminf  = std::min(hminf, meshsize_z);
        hmaxf  = std::max(meshsize_x, meshsize_y);
        hmaxf  = std::max(hmaxf, meshsize_z);
        cout<<"\n##########################################################\n"
            <<"#                 The created mesh information              \n"
            << "########################################################## \n\n"
            << "   nx_mesh = " << n_mesh[0] <<", Lx = " << wall_params[1]-wall_params[0] <<", hx = "<< meshsize_x <<endl
            << "   ny_mesh = " << n_mesh[1] <<", Ly = " << wall_params[3]-wall_params[2] <<", hy = "<< meshsize_y <<endl
            << "   nz_mesh = " << n_mesh[2] <<", Lz = " << wall_params[5]-wall_params[4] <<", hz = "<< meshsize_z <<endl
            << "   minimum mesh size of fluid: hminf = " << hminf << endl
            << "   maximum mesh size of fluid: hmaxf = " << hmaxf <<endl;
     /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     * Create a mesh, distributed across the default MPI communicator.
     * We build a mesh with Quad9(8) elements for 2D and HEX27(20) element for 3D
    / - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        if(dim==2)
        {
        MeshTools::Generation::build_square (*mesh, n_mesh[0], n_mesh[1],
                                           wall_params[0], wall_params[1], wall_params[2], wall_params[3], QUAD8);        // QUAD8/9
        }else{
          MeshTools::Generation::build_cube (*mesh, n_mesh[0], n_mesh[1], n_mesh[2],
                                            wall_params[0], wall_params[1], wall_params[2], wall_params[3], wall_params[4], wall_params[5], HEX20);  // HEX20/27
        }
    } // end wall_type = "slit"
    else if (wall_type == "sphere"){
      cout << "create_domain_mesh () does not support sphere mesh so far ! " << endl;
      libmesh_error (); 
    }
  }// end if (generate_mesh)
  else{
    if(domain_mesh_file != "nothing"){
          mesh->read(domain_mesh_file);    
          mesh->all_second_order();
          mesh->prepare_for_use();
          const std::vector<Real> mesh_size = PMToolBox::mesh_size(*mesh);
          hminf = mesh_size[0];
          hmaxf = mesh_size[1];
          cout <<endl<< "##########################################################"<<endl
           << "#                 The Read-in mesh information                      " <<endl
           << "##########################################################"<<endl<<endl;
            cout << "   minimum mesh size of fluid: hminf = " << hminf << endl;
            cout << "   maximum mesh size of fliud: hmaxf = " << hmaxf << endl;
          }
    else{
        cout <<"**************************************warning***********************************"<<endl;
        cout << "domain_mesh_file has to be specified" << endl;
        cout <<"********************************************************************************"<<endl;    
        libmesh_error();
    } 
  } // end else (generate mesh)

 // initialize search radius
  // if(with_hi == false){
  //   search_radius_p = input_file("search_radius_p", 0.0);
  //   if (search_radius_p <= 0.) {
  //     cout <<"  Warning: for Free draining systems, a search_radius_p(> 0) needs to be given to build particle-particle neighborlist" << endl;
  //     libmesh_error();
  //   }
  // }
  // else{
  //   search_radius_p = 4.0/alpha;
  // }
  search_radius_p = 4./alpha; 
  search_radius_e = 0.5*hmaxf + 4./alpha;
  // print mesh info
  mesh -> print_info();
} // end function

//============================================================================
void Copss::create_periodic_boundary(){
  if (wall_type == "slit"){
    const Point bbox_pmin(wall_params[0], wall_params[2], wall_params[4]);
    const Point bbox_pmax(wall_params[1], wall_params[3], wall_params[5]);
    // construct PMPeriodicBoundary class using info above  
    pm_periodic_boundary = new PMPeriodicBoundary(bbox_pmin, bbox_pmax, periodicity, inlet, inlet_pressure);
  }
  else if (wall_type == "sphere"){
    // check: No PBC, No inlet/outlet
    bool error_flag = false;
    for (int i=0; i<dim; i++){
      if(periodicity[i] or inlet[i] or shear[2*i] or shear[2*i+1]) error_flag = true;
    }
    if (error_flag){
      cout << "spherical domain cannot have PBC or inlet/outlet, check control file" << endl;
      libmesh_error();
    }
    cout << "spherical domain cannot have PBC, but we need to create a PBC object using a non-existed cubic box to keep COPSS running !" << endl; 
    const Point bbox_pmin(-Real(wall_params[0]/2.),-Real(wall_params[0]/2.),-Real(wall_params[0]/2.));
    const Point bbox_pmax(Real(wall_params[0]/2.),Real(wall_params[0]/2.),Real(wall_params[0]/2.));
    // construct PMPeriodicBoundary class using info above  
    pm_periodic_boundary = new PMPeriodicBoundary(bbox_pmin, bbox_pmax, periodicity, inlet, inlet_pressure);    
  }
  else{
    cout << "COPSS::create_periodic_boundary() only support 'slit' or 'sphere' wall_type for now !" << endl;
  }
} // end function


//=============================================================================
EquationSystems Copss::create_equation_systems()
{
  // Initialize equation_systems object using the 'mesh' we created before
  //equation_systems = new EquationSystems(*mesh);
  cout << "==>(1/8) Initialize equation_systems object using the 'mesh' we created before" <<endl;
  EquationSystems equation_systems(*mesh);
  // Add 'Stokes' system (of PMLinearImplicitSystem) to the 'equation_systems'
  cout << "==>(2/8) Add 'Stokes' system (of PMLinearImplicitSystem) to the 'equation_systems'" <<endl;
  PMLinearImplicitSystem& system = equation_systems.add_system<PMLinearImplicitSystem> ("Stokes");
  cout << "==>(3/8) Add variables to 'Stokes' system" <<endl;
  //Add variables to 'Stokes' system"
  u_var = system.add_variable ("u", SECOND);
  v_var = system.add_variable ("v", SECOND);
  if(dim==3)  w_var  = system.add_variable ("w", SECOND);
  const unsigned int p_var = system.add_variable ("p", FIRST);

  // attach object_mesh to pm_linear_implicit_system
  cout << "==>(4/8) Attach object_mesh to the system" <<endl;  
  this->attach_object_mesh(system);

  // attach period boudary to pm_linear_implicit_system
  cout<<"==>(5/8) Add period boundary conditions to 'Stokes' system"<<endl;
  this->attach_period_boundary(system);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initialize the Preconditioning matrix for saddle point problems if required.
   Initialize the equation system and zero the preconditioning matrix
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if( user_defined_pc ){
    system.add_matrix("Preconditioner");
  } 
  /* Initialize the data structures for the equation system. */
  cout<<"==>(6/8) Init equation_systems (libmesh function, to init all systems in equation_systems)"<<endl;
  equation_systems.init();
  // zero the PC matrix, which MUST be done after es.init()
  if( user_defined_pc ) {
    cout<<"==> (If user_defined_pc) Zero preconditioner matrix"<<endl;
    system.get_matrix("Preconditioner").zero();
  }
  cout<<"--------------> Equation systems are initialized:\n"<<std::endl;
  // set parameters for equation systems
  cout<<"==>(7/8) Set parameters of equation_systems"<<endl;
  this -> set_parameters(equation_systems);
  // initialized force field
  cout<<"==>(8/8) Attach fixes to 'stokes' system"<<endl;
  this -> attach_fixes(system); 
  /* Print information about the mesh and system to the screen. */
  cout << endl <<"--------------> Print equation systems info" <<endl;
    equation_systems.print_info();
    cout <<"  System has: "<< mesh->n_elem()<<" elements,\n"
             <<"              "<< mesh->n_nodes()<<" nodes,\n"
             <<"              "<< equation_systems.n_dofs()<<" degrees of freedom.\n"
             <<"              "<< equation_systems.n_active_dofs()<<" active degrees of freedom.\n"
             <<"              "<< point_mesh->num_particles()<<" particles.\n" << std::endl;
  return equation_systems;
}

//=======================================================================================
void Copss::attach_fixes(PMLinearImplicitSystem& pm_system)
{

  // Fix *f;
  fix_factory = new FixFactory();
  fixes.resize(numForceTypes);
  for (int i=0; i < numForceTypes; i++){
    fixes[i] = fix_factory -> buildFix(forceTypes[i], pm_system);
  }
  pm_system.attach_fixes(fixes);
}


//===============================================================================
void Copss::attach_period_boundary(PMLinearImplicitSystem& system)
{
  cout << "--------------> Get dof_map of 'Stokes' system"<<endl;  
  DofMap& dof_map = system.get_dof_map();
  cout << "--------------> Add periodicBoundary object to 'dof_map'"<<endl;
  /*** set PBC in x-direction ***/
  if (periodicity[0])
  {
    PeriodicBoundary pbcx(RealVectorValue(wall_params[1]-wall_params[0], 0., 0.));
    pbcx.set_variable(u_var);
    pbcx.set_variable(v_var);
    if(dim==3) pbcx.set_variable(w_var);
    //pbcx.set_variable(p_var); //*** NOT include p!
    
    // is this boundary number still true for 3D?
    if(dim==2)
    {
      pbcx.myboundary = 3;
      pbcx.pairedboundary = 1;
    }
    else if(dim==3)
    {
      pbcx.myboundary = 4; // left face (viewing from X negative to X positive)
      pbcx.pairedboundary = 2; // right face (viewing from X positive to X negative)

    } // end if
    
    dof_map.add_periodic_boundary(pbcx);    
    // check
    if (search_radius_p>=(wall_params[1]-wall_params[0])/2.)
    {
      output_msg = std::string("****************************** warning: ********************************\n")+
                   "**** The search radius is larger than half domain length in x direction\n!"+
                   "**** search radius = "+ std::to_string(search_radius_p) + 
                   ", half domain size Lx/2 =" + std::to_string((wall_params[1]-wall_params[0])/2.)+"\n"+
                   "************************************************************************\n\n";
      PMToolBox::output_message(output_msg, comm_in);
      libmesh_error();
    }
  }
  /*** set PBC in y-direction ***/
  if (periodicity[1])
  {
    PeriodicBoundary pbcy(RealVectorValue(0., wall_params[3]-wall_params[2], 0.));
    pbcy.set_variable(u_var);
    pbcy.set_variable(v_var);
    if(dim==3) pbcy.set_variable(w_var);
    //pbcy.set_variable(p_var); //*** NOT include p!
    
    if(dim==2)
    {
      pbcy.myboundary = 0;
      pbcy.pairedboundary = 2;
    }
    else if(dim==3)
    {
      pbcy.myboundary = 1; // bottom face, viewing from Y negative to Y positive
      pbcy.pairedboundary = 3; // top face, viewing from Y positive to Y negative
    } // end if
    
    dof_map.add_periodic_boundary(pbcy);
    // check
    if (search_radius_p>=(wall_params[3]-wall_params[2])/2.)
    {
      output_msg = std::string("****************************** warning: ********************************\n")+
                   "**** The search radius is larger than half domain length in y direction\n!"+
                   "**** search radius = "+ std::to_string(search_radius_p) + 
                   ", half domain size Ly/2 =" + std::to_string((wall_params[3]-wall_params[2])/2.)+"\n"+
                   "************************************************************************\n\n";
      PMToolBox::output_message(output_msg, comm_in);
      libmesh_error();
    }
  }
  /*** set PBC in z-direction ***/
  if (periodicity[2])
  {
    PeriodicBoundary pbcz(RealVectorValue(0., 0., wall_params[5]-wall_params[4]));
    pbcz.set_variable(u_var);
    pbcz.set_variable(v_var);
    if(dim==3) pbcz.set_variable(w_var);
    //pbcz.set_variable(p_var); //*** NOT include p!
    
    if(dim==3)
    {
      pbcz.myboundary = 0; // bottom face, viewing from Z negative to Z positive
      pbcz.pairedboundary = 5; // top face, viewing from Z positive to Z negative 
    } // end if
    
    dof_map.add_periodic_boundary(pbcz);
    // check
    if (search_radius_p>=(wall_params[5]-wall_params[4])/2.)
    {
      output_msg = std::string("****************************** warning: ********************************\n")+
                   "**** The search radius is larger than half domain length in z direction\n!"+
                   "**** search radius = "+ std::to_string(search_radius_p) + 
                   ", half domain size Lx/2 =" + std::to_string((wall_params[5]-wall_params[4])/2.)+"\n"+
                   "************************************************************************\n\n";
      PMToolBox::output_message(output_msg, comm_in);
      libmesh_error();
    }
  }  
}

//============================================================================================
void Copss::solve_undisturbed_system(EquationSystems& equation_systems)
{
   // get stokes system from equation systems
  PMLinearImplicitSystem& system = equation_systems.get_system<PMLinearImplicitSystem> ("Stokes");
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Compute undisturbed velocity field without particles.
   NOTE: We MUST re-init particle-mesh before solving Stokes
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (update_neighbor_list_everyStep) neighbor_list_update_flag = true;
  system.reinit_hi_system(neighbor_list_update_flag); // neighbor_list_update is ture here
  if(print_info){
    cout << "print particle information after first reinit\n";
    if(comm_in.rank()==0)
    point_mesh->print_point_info();
  }
  reinit_stokes = true;
  system.solve_stokes("undisturbed",reinit_stokes);
  v0_ptr = system.solution->clone(); // backup v0
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write out the equation systems at Step 0 (undisturbed field)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  exodus_ptr = new ExodusII_IO(*mesh);
  if(std::find(output_file.begin(), output_file.end(), "equation_systems") != output_file.end() && restart==false)
  {
    //system.add_local_solution(); // Don't add local solution for undisturbed system!
#ifdef LIBMESH_HAVE_EXODUS_API
    exodus_ptr->write_equation_systems(out_system_filename, equation_systems);
#endif
  }
}

//============================================================================================
void Copss::create_brownian_system(EquationSystems& equation_systems)
{
 /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create vectors and Shell Mat for use:
   U0:          particle velocity vector;
   R0/R_mid:    particle position vector;
   dw/dw_mid:   random vector;
   RIN/ROUT:    the initial and intermediate particle postion vector for msd output
   RIN will not change, and ROUT excludes pbc
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  brownian_sys = new BrownianSystem(equation_systems);
  brownian_sys->init_petsc_random(&rand_ctx);
  brownian_sys->_create_shell_mat(n_vec, &M);
  brownian_sys->_create_petsc_vec(n_vec,&R0);
  VecDuplicate(R0,&U0);
  VecDuplicate(R0,&R_mid);
  VecDuplicate(R0,&dw_mid);
  brownian_sys->extract_particle_vector(&ROUT,"coordinate","extract");
  VecDuplicate(ROUT,&RIN);
  VecCopy(ROUT,RIN);  // RIN = ROUT = the initial position vector
  brownian_sys->set_std_random_seed(random_seed);  // random seed
  if(restart)
  {
    // read RIN & ROUT from local file output during the previous simulation
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vector_RIN.dat",FILE_MODE_READ,&viewer);
    VecLoad(RIN,viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vector_ROUT.dat",FILE_MODE_READ,&viewer);
    VecLoad(ROUT,viewer);
  }
  else
  {
    // write out binary file of RIN, which may be used at restart mode.
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vector_RIN.dat",FILE_MODE_WRITE,&viewer);
    VecView(RIN,viewer);
  }
  PetscViewerDestroy(&viewer);
  comm_in.barrier();
}

//==============================================================================================
void Copss::fixman_integrate(EquationSystems& equation_systems, unsigned int& i)
{
//  PerfLog perf_log("integration");
  // get stokes system from equation systems
  PMLinearImplicitSystem& system = equation_systems.get_system<PMLinearImplicitSystem> ("Stokes");
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Compute the "disturbed" particle velocity + "undisturbed" velocity = U0
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //cout <<"Compute the disturbed particle velocity at step "<<i+1<<endl;
  // perf_log.push("reinit sytem");
  if(i>0){
    *(system.solution) = *v0_ptr; // re-assign the undisturbed solution
    // Update the local values to reflect the solution on neighboring processors
    system.update();
    if(update_neighbor_list_everyStep){
      neighbor_list_update_flag = true;
    }
    else if (timestep_duration > neighbor_list_update_interval){
       cout << "======> update neighbor list at timestep " << i << endl;
       neighbor_list_update_flag = true;
       timestep_duration = 0;
    }
    system.reinit_hi_system(neighbor_list_update_flag);
  }
  if(print_info){
  if(comm_in.rank()==0)
    point_mesh->print_point_info();
  }
  // perf_log.pop("reinit sytem");
//  perf_log.push("compute_point_velocity undisturbed");
  // compute undisturbed velocity of points 
  system.compute_point_velocity("undisturbed", vel0);
//  perf_log.pop("compute_point_velocity undisturbed");

  reinit_stokes = false;
//  perf_log.push("solve_stokes disturbed");
  system.solve_stokes("disturbed",reinit_stokes); // Using StokesSolver
//  perf_log.pop("solve_stokes disturbed");

  // compute distrubed velocity of points
//  perf_log.push("compute_point_velocity disturbed");
  system.compute_point_velocity("disturbed", vel1);
//  perf_log.pop("compute_point_velocity disturbed");
  // add up undistrubed and disturbed velocity of points
  for(std::size_t j=0; j<vel1.size();++j) vel1[j] += vel0[j];
 // transform total point velocity to U0 in Brownian_system
  brownian_sys->vector_transform(vel1, &U0, "forward");
  // assign vel1 to particle velocity
  point_mesh->set_bead_velocity(vel1);
  if (debug_info){
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     * ---> test: output the particle velocity
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    const Real v0_min = v0_ptr->min();
    const Real v0_max = v0_ptr->max();
    const Real v0_sum = v0_ptr->sum();
    if(comm_in.rank()==0){
      for(unsigned int j=0; j<NP;++j)
      {    
        std::vector<Real> vtest0(dim), vtest1(dim);
        for(std::size_t k=0; k<dim;++k){
          vtest0[k] = vel0[j*dim+k];
          vtest1[k] = vel1[j*dim+k];
        }
        Point vtest2 = point_mesh->particles()[j]->particle_velocity();
        printf("--->test in test_move_particles(): velocity on the %u-th point:\n",j);
        printf("            U0 = (%.12e,%.12e,%.12e)\n",   vtest0[0],vtest0[1],vtest0[2]);
        printf("       U0 + U1 = (%.12e,%.12e,%.12e)\n", vtest1[0],vtest1[1],vtest1[2]);
        printf("     bead velocity(should be equal to U0+U1) = (%.12e,%.12e,%.12e)\n\n", vtest2(0), vtest2(1), vtest2(2));
      }    
      printf("fluid field v0_min = %.12e, v0_max = %.12e, v0_sum = %.12e)\n",v0_min,v0_max,v0_sum);
    }    
  }
  /*---------------------------------------------------------------------------------------
   * write equation system at step i
   * print out information at step 0
   * do not print out information at the first step when restart since it is identical to the 
   * last step before restart.
  -----------------------------------------------------------------------------------------*/
  if(i%write_interval == 0){
    cout << "\nStarting Fixman integration at step "<< i << endl;
    if (i != restart_step){
      /*
       * write equation system to output file
       */
      if (i != 0){
        if(std::find(output_file.begin(), output_file.end(), "equation_systems") != output_file.end())
        {
          system.add_local_solution();  // add local solution for the disturbed system
          system.solution->add(*v0_ptr);// add the undisturbed solution
  #ifdef LIBMESH_HAVE_EXODUS_API
          exodus_ptr->append(true);
          exodus_ptr->write_timestep(out_system_filename,equation_systems,o_step,o_step);
  #endif
        } // end if (write es)  
      }
      /*
       * write particle to output file
       */ 
      this -> write_object(i);
      /*----------------------------------------------------------------------------------------------------
       * Write out ROUT for restart mode at step i
       ----------------------------------------------------------------------------------------------------*/
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vector_ROUT.dat",FILE_MODE_WRITE,&viewer);
      VecView(ROUT,viewer);
      PetscViewerDestroy(&viewer);        
    }
    // update o_step
    o_step++;
  } // end if (i % write_interval == 0 )
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Adaptive time step.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Real dt = 0;
  if(adaptive_dt){
    Real vp_max = point_mesh->maximum_bead_velocity();
    Real vp_min = point_mesh->minimum_bead_velocity();
    if (i==0) dt = (vp_max <= 1.) ? (max_dr_coeff[0]) : (max_dr_coeff[0] * 1. / vp_max);  
    else dt = (vp_max <= 1.) ? (max_dr_coeff.back()) : (max_dr_coeff.back() * 1. / vp_max);  
    if(i % write_interval == 0){
      cout << "       ##############################################################################################################" << endl
           << "       # Max velocity magnitude is " <<std::setprecision(o_precision) <<std::fixed << vp_max << endl
           << "       # Min velocity magnitude is " <<std::setprecision(o_precision) <<std::fixed << vp_min << endl
           << "       # minimum mesh size = " << hmin << endl
           << "       # The adaptive time increment at step "<< i << " is dt = " <<std::setprecision(o_precision) <<std::fixed << dt<<endl
           << "       # (with Brownian) adapting_time_step = max_dr_coeff * bead_radius / (max_bead_velocity at t_i)" << endl
           << "       # (without Brownian) adapting_time_step = max_dr_coeff * mesh_size_min / (max_bead_velocity at t_i) "<<endl      
           << "       # Chebyshev failure steps = " << n_chebyshev_failure << endl
           << "       ##############################################################################################################" << endl;
    } // end if (i% write_interval)  
  }
  else{
    if (i==0) dt = max_dr_coeff[0] * 1.;
    else dt = max_dr_coeff.back() * 1.; 
    if(i % write_interval == 0){
      cout << "       ##############################################################################################################" << endl
           << "       # The fixed time increment at step "<< i << " is dt = " <<std::setprecision(o_precision)<<std::fixed << dt<<endl
     	     << "       # (With Brownian) fixed_time_step = max_dr_coeff * bead_radius / 1.0"<<endl
           << "       # (Without Brownian) fixed_time_step = max_dr_coeff * mesh_size_min / 1.0 "<<endl
           << "       # Chebyshev failure steps = " << n_chebyshev_failure << endl
           << "       ##############################################################################################################" << endl;
    } // end if (i % write_interval == 0)
  } // end if (adaptive_dt)
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    If with Brownian motion, we use midpoint scheme
    If without Brownian motion, we use normal stepping: dR = Utotal*dt
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */  
  if (with_brownian)
  {
    /*-------------------------------------------
    * For restarted Brownian system, to make sure the trajectory after is exactly 
    * the same as a continus simulation, we need to let the random number generator
    * run "restart_step" steps
    * ------------------------------------------*/
    if (i==restart_step){
      cout << "       ##############################################################################################################" << endl
           << "       # Run RandomGenerator for "<<restart_step<<" steps to make sure the trajectroy after restart is the same as a " <<endl
	   << "       # continuous simulation (this is only done once at the first step of restarted simulation)" <<endl
           << "       ##############################################################################################################" << endl;
      Vec dw_prep; // this Vec is never used
      for (unsigned int j=0; j<restart_step; j++){
        brownian_sys->std_random_vector(0.0, 1.0,"gaussian", &dw_prep);
	VecDestroy(&dw_prep);
      }
    }
    // PerfLog perf_log("BD step");
    // perf_log.push("bd");
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Generate random vector dw whose mean = 0, variance = sqrt(2*dt)
     petsc_random_vector generates a uniform distribution [0 1] whose
     mean = 0.5 and variance = 1/12, so we need a shift and scale operation.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    Real mean_dw = 0.0, variance_dw = 0.0;
    // A more precise way is to construct a random vector with gaussian distribution
    const Real std_dev  = std::sqrt(dt);
    brownian_sys->std_random_vector(0.0,std_dev,"gaussian",&dw);
    brownian_sys->_vector_mean_variance(dw, mean_dw, variance_dw);
    VecScale(dw,std::sqrt(2.0));
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Print out the mean and variance or view the generated vector.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//      PetscPrintf(PETSC_COMM_WORLD,
//                 "Generated random_vector: mean = %f, variance = %f\n",
//                 mean_dw, variance_dw);
//      PetscPrintf(PETSC_COMM_WORLD,
//                 "Predicted random_vector: mean = %f, variance = %f\n",
//                 0., std::sqrt(2.*dt));
    // Compute dw = B^-1 * dw using Chebyshev polynomial, dw will be changed!
    VecCopy (dw,dw_mid);  // save dw to dw_mid, which will be used for Chebyshev
    // perf_log.push("cheb_converge");
    for(std::size_t j=0; j<2; j++)
    {
    //  cout << "j = "<<j<<endl;
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Compute the max/min eigenvalues if needed. Otherwise, magnify the interval.
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      if(compute_eigen){
        cout <<"compute_eigen = "<<compute_eigen << endl;
        //cout << "Compute the max & min eigenvalues for Chebyshev polynomial at step "<<i+1<<endl;
        brownian_sys->compute_eigenvalues(eig_min,eig_max,tol_eigen);
        // Magnify the spectral range by a factor (1.05 by default).
	      eig_max *= eig_factor; eig_min /= eig_factor;
        if(eig_min < 0 or eig_max < 0 or eig_min > eig_max){
          cout <<"--->Invalid eigenvalue range: eig_min = "<<eig_min <<"; eig_max = "<<eig_max << endl;
          cheb_converge = false;
          break;  
        } 
        else{
          cout <<"--->Valid eigenvalue range: recomputed eigen values and magnify the range by a factor eig_factor = " << eig_factor 
               <<"; eig_min = " << eig_min
               <<"; eig_max = " << eig_max
               <<"; tol_cheb = "<< tol_cheb
               <<"; max_n_cheb = " << max_n_cheb << endl
               <<"--->Write this eigenvalue range to out.eigenvalue for restart purpose" << endl;
          std::ofstream out_file;  
          if(comm_in.rank() == 0){
            if(i == 0){
              out_file.open("out.eigenvalue", std::ios_base::out);
              out_file << "eigen_min" <<"  " <<"eigen_max" <<"\n";      
            }
            else{
              out_file.open("out.eigenvalue", std::ios_base::app);
            }
            out_file.precision(o_precision);
            out_file <<eig_min <<"  " <<eig_max <<"\n";
            out_file.close();  
          }
        }
      }
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Compute the Brownian displacement B^-1 * dw using Chebyshev approximation.
       Here dw is both input and output variables, so it will be changed.
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      cheb_converge = brownian_sys->chebyshev_polynomial_approximation(max_n_cheb,
                                                                      eig_min,eig_max,tol_cheb,&dw);  
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       If converged, dw returns the Brownian displacement, then break the j-loop;
       Otherwise, recompute eigenvalues
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      if(cheb_converge){
        compute_eigen = false; break;
      }
      else{
        compute_eigen = true;
        VecCopy(dw_mid,dw); /*copy back, recompute eigenvalues*/
        //cout << "It is necessry to re-compute the eigenvalues at step " <<i+1<<endl;
      }
    } // end for j-loop
    // perf_log.pop("cheb_converge");
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Double-check the convergence of Chebyshev polynomial approximation
     If cheb_converge = true continue fixman integration
     otherwise, relax the system for serval free-draining steps
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    // perf_log.push("fixman step");
    if(cheb_converge){
      // Compute dw_mid = D*B^-1*dw, which can be obtained by solving the Stokes
      brownian_sys->hi_ewald(M,dw,dw_mid);  // dw_mid = D * dw
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Particle coordinate vector R0.
       Move the particle R_mid = R0 + 0.5*(U0+U1)*dt (deterministic)
       and R_mid = R_mid + 0.5*sqrt(2)*D*B^-1*dw     (stochastic)
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      brownian_sys->extract_particle_vector(&R0,"coordinate","extract");
      VecWAXPY(R_mid,0.5*dt,U0,R0);  // R_mid = R0 + 0.5*dt*(U0+U1)  (R0 and U0 do NOT change)
      coef = 0.5;                    // coefficient. sqrt(2) is introduced when generating dw
      VecAXPY(R_mid,coef,dw_mid);    // R_mid = R_mid + 0.5*sqrt(2)*D*B^-1*dw
      brownian_sys->extract_particle_vector(&R_mid,"coordinate","assign"); // Update mid-point coords && apply pbc to particle position
      this -> update_object(); 
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Update the particle mesh for the mid-point step,
       and recompute U0 + U1_mid, D_mid*(B^-1*dw)
       NOTE: the FEM solution of undisturbed field doesn't change, but particles
       move, so U0 needs to be re-evaluated at the new position.
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      *(system.solution) = *v0_ptr;       // re-assign the undisturbed solution
      system.update();
      // comment the line below if not update neighbor list at each time step
      if(update_neighbor_list_everyStep) neighbor_list_update_flag = true;
      system.reinit_hi_system(neighbor_list_update_flag);
      system.compute_point_velocity("undisturbed", vel0);
      reinit_stokes = false;
      system.solve_stokes("disturbed",reinit_stokes);   // solve the disturbed solution
      system.compute_point_velocity("disturbed", vel1);
      for(std::size_t j=0; j<vel1.size();++j) vel1[j] += vel0[j];
      brownian_sys->vector_transform(vel1, &U0, "forward"); // (U0+U1)_mid
      brownian_sys->hi_ewald(M,dw,dw_mid);  // dw_mid = D_mid*dw, where dw=B^-1*dw computed above
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       the mid-point to the NEW point, and update the particle coordinates
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      //cout << "Update from mid-point to the NEW particle coordinates at step " <<i+1<<endl;
      VecWAXPY(R_mid,dt,U0,R0);         // R_mid = R0 + dt*U0_mid
      VecAXPY(R_mid,2.0*coef,dw_mid); // R_mid = R_mid + sqrt(2)*D_mid*B^-1*dw
      brownian_sys->extract_particle_vector(&R_mid,"coordinate","assign"); // Update mid-point coords && apply pbc to particle position
      this -> update_object();
      // Update ROUT (position vector excluding pbc) at the i-th step
      VecAXPY(ROUT,dt,U0);            // ROUT = ROUT + dt*U0_mid
      VecAXPY(ROUT,2.0*coef,dw_mid);  // ROUT = ROUT + sqrt(2)*D_mid*B^-1*dw
    }// end if cheb_converge
    else {
      n_chebyshev_failure += 1;
      PMToolBox::output_message("****** Warning: After recomputing eigenvalues, Chebysheve failed to converge at step " +std::to_string(i)+"; relax the system for "+std::to_string(n_relax_step)+" steps", comm_in);
      for (int relax_step=0; relax_step<n_relax_step; relax_step++){
        langevin_integrate(equation_systems, i);
        i++;
      }
    }// end else cheb_converge
    // perf_log.pop("fixman step");
    // perf_log.pop("bd");
  } // end if Brownian
  else{ // if without Brownian
    // Move the particle R_mid = R0 + (U0+U1)*dt (deterministic)
    brownian_sys->extract_particle_vector(&R0,"coordinate","extract");
    VecWAXPY(R_mid,dt,U0,R0);  // R_mid = R0 + dt*Utotal (U0 is actually Utotal)
    brownian_sys->extract_particle_vector(&R_mid,"coordinate","assign"); // Update mid-point coords 
    this -> update_object();
    // Update ROUT (position vector excluding pbc) at the i-th step
    VecAXPY(ROUT,dt,U0); // ROUT = ROUT + dt*U0_mid   
  } // end else (without_brownian)
  real_time += dt;
  timestep_duration += 1;
}

void Copss::langevin_integrate(EquationSystems& equation_systems, unsigned int& i)
{
  // get stokes system from equation systems
  PMLinearImplicitSystem& system = equation_systems.get_system<PMLinearImplicitSystem> ("Stokes");
  if(i>0){
    system.update();
    if(update_neighbor_list_everyStep){
      neighbor_list_update_flag = true;
    }
    else if (timestep_duration > neighbor_list_update_interval){
       neighbor_list_update_flag = true;
       timestep_duration = 0;
    }
    // whether or not reinit neighbor list depends on the neighbor_list_update_flag
    system.reinit_fd_system(neighbor_list_update_flag);
  }
  Point p_velocity(0.);
  for (std::size_t p_id = 0; p_id < NP; p_id++) {
    for (int _dim = 0; _dim < dim; _dim++){
      p_velocity(_dim) = point_mesh->particles()[p_id]->particle_force()(_dim);        
      vel1[dim*p_id+_dim] = p_velocity(_dim);
    } 
  }
  // set vel1 to particle velocity
  point_mesh->set_bead_velocity(vel1);
  // transform total point velocity to U0 in Brownian_system
  brownian_sys->vector_transform(vel1, &U0, "forward");
  /*---------------------------------------------------------------------------------------
   * write equation system at step i
   * print out information at step 0
   * do not print out information at the first step when restart since it is identical to the 
   * last step before restart.
  -----------------------------------------------------------------------------------------*/
  if(i%write_interval == 0){
    cout << "\nStarting Langevin integration at step "<< i << endl;
    if (i != restart_step){  
      /*
       * write particle to output file
       */ 
      this -> write_object(i);
      /*----------------------------------------------------------------------------------------------------
       * Write out ROUT for restart mode at step i
       ----------------------------------------------------------------------------------------------------*/
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vector_ROUT.dat",FILE_MODE_WRITE,&viewer);
      VecView(ROUT,viewer);
      PetscViewerDestroy(&viewer);        
    }
    // update o_step
    o_step++;
  } // end if (i % write_interval == 0 )
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Adaptive time step.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Real dt = 0;
  if(adaptive_dt){
    Real vp_max = point_mesh->maximum_bead_velocity();
    Real vp_min = point_mesh->minimum_bead_velocity();
    if (i==0) dt = (vp_max <= 1.) ? (max_dr_coeff[0]) : (max_dr_coeff[0] * 1. / vp_max);  
    else dt = (vp_max <= 1.) ? (max_dr_coeff.back()) : (max_dr_coeff.back() * 1. / vp_max);  
    if(i % write_interval == 0){
      cout << "       ##############################################################################################################" << endl
           << "       # Max velocity magnitude is " << std::setprecision(o_precision) <<std::fixed <<vp_max << endl
           << "       # Min velocity magnitude is " << std::setprecision(o_precision) <<std::fixed <<vp_min << endl
           << "       # minimum mesh size = " << hmin << endl
           << "       # The adaptive time increment at step "<< i << " is dt = " << std::setprecision(o_precision) <<std::fixed <<dt<<endl
           << "       # adapting_time_step = max_dr_coeff * bead_radius / (max_bead_velocity at t_i)" << endl
           << "       ##############################################################################################################" << endl;
    } // end if (i% write_interval)  
  }
  else{
    if (i==0) dt = max_dr_coeff[0] * 1.;
    else dt = max_dr_coeff.back() * 1.; 
    if(i % write_interval == 0){
      cout << "       ##############################################################################################################" << endl
           << "       # The fixed time increment at step "<< i << " is dt = " << std::setprecision(o_precision) <<std::fixed <<dt<<endl
           << "       # fixed_time_step = max_dr_coeff * bead_radius / 1.0"<<endl
           << "       ##############################################################################################################" << endl;
    } // end if (i % write_interval == 0)
  } // end if (adaptive_dt)
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    If with Brownian motion, we use midpoint scheme
    If without Brownian motion, we use normal stepping: dR = Utotal*dt
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */  
  if (with_brownian)
  {
    /*-------------------------------------------
    * For restarted Brownian system, to make sure the trajectory after is exactly 
    * the same as a continus simulation, we need to let the random number generator
    * run "restart_step" steps
    * ------------------------------------------*/
    if (i==restart_step){
      Vec dw_prep; // this Vec is never used
      for (unsigned int j=0; j<restart_step; j++){
        brownian_sys->std_random_vector(0.0, 1.0,"gaussian", &dw_prep);
	VecDestroy(&dw_prep);
      }
    }
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Generate random vector dw whose mean = 0, variance = sqrt(2*dt)
     petsc_random_vector generates a uniform distribution [0 1] whose
     mean = 0.5 and variance = 1/12, so we need a shift and scale operation.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    Real mean_dw = 0.0, variance_dw = 0.0;
    // A more precise way is to construct a random vector with gaussian distribution
    const Real std_dev  = std::sqrt(dt);
    brownian_sys->std_random_vector(0.0,std_dev,"gaussian",&dw);
    brownian_sys->_vector_mean_variance(dw, mean_dw, variance_dw);
    VecScale(dw,std::sqrt(2.0));
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Print out the mean and variance or view the generated vector.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//    PetscPrintf(PETSC_COMM_WORLD,
//               "Generated random_vector: mean = %f, variance = %f\n",
//               mean_dw, variance_dw);
//    PetscPrintf(PETSC_COMM_WORLD,
//               "Predicted random_vector: mean = %f, variance = %f\n",
//               0., std::sqrt(2.*dt));
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Particle coordinate vector R0.
     Move the particle R_mid = R0 + U0*dt (deterministic)
     and R_mid = R_mid + dw, where dw = sqrt(2*dt)*R, where <R>=0, <R^2>=1     (stochastic)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    brownian_sys->extract_particle_vector(&R0,"coordinate","extract");
    VecWAXPY(R_mid,dt,U0,R0);  // R_mid = R0 + dt * U0  (R0 and U0 do NOT change)
    coef = 1.;                    // coefficient. 
    VecAXPY(R_mid,coef,dw);    // R_mid = R_mid + dw, where dw = sqrt(2*dt)*R, where <R>=0, <R^2>=1
    brownian_sys->extract_particle_vector(&R_mid,"coordinate","assign"); // Update mid-point coords && apply pbc to particle position
    this -> update_object(); 
    // Update ROUT (position vector excluding pbc) at the i-th step
    VecAXPY(ROUT,dt,U0);            // ROUT = ROUT + dt*U0_mid
    VecAXPY(ROUT,coef,dw);   // ROUT = ROUT + dw, where dw = sqrt(2*dt)*R, where <R>=0, <R^2>=1
  } // end if Brownian
  
  else{ // if without Brownian
    // Move the particle R_mid = R0 + (U0+U1)*dt (deterministic)
    brownian_sys->extract_particle_vector(&R0,"coordinate","extract");
    VecWAXPY(R_mid,dt,U0,R0);  // R_mid = R0 + dt*Utotal (U0 is actually Utotal)
    brownian_sys->extract_particle_vector(&R_mid,"coordinate","assign"); // Update mid-point coords
    this -> update_object();
    // Update ROUT (position vector excluding pbc) at the i-th step
    VecAXPY(ROUT,dt,U0); // ROUT = ROUT + dt*U0_mid   
  } // end else (without_brownian)
  real_time += dt;
  timestep_duration += 1;
}


void Copss::destroy()
{
  MatDestroy(&M);
  VecDestroy(&U0);
  VecDestroy(&R0);
  VecDestroy(&R_mid);
  VecDestroy(&dw_mid);
  PetscRandomDestroy(&rand_ctx);
  if(with_brownian){
    VecDestroy(&dw);
  }
  if(exodus_ptr and with_hi) {
    delete exodus_ptr;
  }
  PetscViewerDestroy(&viewer);
}


} // end namespace





