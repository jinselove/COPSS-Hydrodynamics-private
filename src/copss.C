
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

namespace libMesh {
// ==========================================================================
Copss::Copss(CopssInit& init)
{
  comm_in = &init.comm();
  ss << "============================0.  Initialize libMesh.===========================";
  PMToolBox::output_message(ss, *comm_in);
  this->check_libmesh();
}

Copss::~Copss()
{ 
  delete mesh; mesh = nullptr;
  delete point_mesh; point_mesh = nullptr;
  delete pm_periodic_boundary; pm_periodic_boundary = nullptr;
  delete brownian_sys; brownian_sys = nullptr;
  delete fix_factory; fix_factory = nullptr;
  for (int i = 0; i < fixes.size(); i++) 
    {
      delete fixes[i]; fixes[i] = nullptr;
    }
  fixes.clear(); 
}

int Copss::check_libmesh() {
  // This example NaNs with the Eigen sparse linear solvers and Trilinos
  // solvers,
  // but should work OK with either PETSc or Laspack.
  libmesh_example_requires(libMesh::default_solver_package() != EIGEN_SOLVERS,
                           "--enable-petsc or --enable-laspack");
  libmesh_example_requires(libMesh::default_solver_package() != TRILINOS_SOLVERS,
                           "--enable-petsc or --enable-laspack");

  // check libmesh enables
  #ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
  #endif // ifndef LIBMESH_ENABLE_AMR

  #ifndef LIBMESH_HAVE_SLEPC
  libmesh_example_requires(false, "--enable-slepc");
  #endif // ifndef LIBMESH_HAVE_SLEPC

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(3 == LIBMESH_DIM, "--3D support");

  return 0;
}

// ==========================================================================
void Copss::start_time(struct tm *timeinfo) {
  ss << "---------------------------------------------------------------------\n"
     << "The program starts. \n"
     << "The current date/time is: " << asctime(timeinfo)  << "\n"
     << "---------------------------------------------------------------------\n";
  PMToolBox::output_message(ss, *comm_in);
}

// ==========================================================================
void Copss::end_time(struct tm *timeinfo) {
  ss << "\n---------------------------------------------------------------------\n"
     << "The current date/time is: " << asctime(timeinfo) << "\n"
     << "The program ends. \n"
     << "---------------------------------------------------------------------\n";
  PMToolBox::output_message(ss, *comm_in);
}

// ====================================================================
EquationSystems Copss::init_system(std::string _control_fileName) {
  control_fileName = _control_fileName;
  ss << "============================1. read input parameters ============================";
  PMToolBox::output_message(ss, *comm_in);
  this->read_input();
  ss << "============================2. Create Point-mesh object =========================";
  PMToolBox::output_message(ss, *comm_in);
  this->create_object_mesh();
  ss << "==========3. Create equation_systems object (of type EquationSystems) ===========";
  PMToolBox::output_message(ss, *comm_in);
  return this->create_equation_systems();
}

// ====================================================================
void Copss::read_input()
{
  const GetPot tmp(control_fileName);

  input_file = tmp;
  this->read_system_info();
  this->read_module_info();
  this->read_physical_info();
  this->read_particle_info();
  this->read_domain_info();
  this->read_force_info();
  this->read_ggem_info();
  this->read_solver_info();
  this->read_chebyshev_info();
  this->read_run_info();
} // end read_data function

// ====================================================================
void Copss::read_system_info()
{
  simulation_name = input_file("simulation_name", "validation");
  ss << "##########################################################\n"
     <<  "#                       simulation_name                    \n"
     <<  "###########################################################\n"
     <<  "-----------> simulation_name: " << simulation_name.c_str();
  PMToolBox::output_message(ss, *comm_in);
  print_info = input_file("print_info", false);
}

// ====================================================================
void Copss::read_module_info()
{
  module_poisson = input_file("module_poisson", false);
  module_np = input_file("module_np", false);
  ss << "##########################################################\n"
     << "#                       module info                    \n"
     << "###########################################################\n"
     << "-----------> Stokes module: true (HI can be turned off by setting "
        "with_hi = false)\n";
  if (module_poisson) ss << "-----------> Poisson module: true\n";
  if (module_np) ss << "-----------> Nernst-Planck module: true\n";
  PMToolBox::output_message(ss, *comm_in);
  print_info = input_file("print_info", false);
}

// ====================================================================
void Copss::read_physical_info()
{
  // temperature T, [T] = K
  T = input_file("temperature", 297);
  // characteristic energy, [kBT] = (N*um)
  kBT = kB * T;
  // viscosity, [viscosity] = N*s/um^2
  viscosity = input_file("viscosity", 1.0);
  // characteristic length, [Rb] = um
  Rb = input_file("radius", 0.10);
  // drag coefficient, [drag_c] = N*s/um
  drag_c    = 6. * PI * viscosity * Rb;
  // Characteristic diffusion coefficient, i.e., diffusion of a bead of
  // radius Rb in infinite dilution, [Db] = um^2/s
  Db        = kBT / drag_c;
  // characteristic time, i.e., the time that a bead of radius Rb diffuses a
  // bead radius away, [tc] = s
  tc  = drag_c * Rb * Rb / kBT;
  // characteristic velocity, [uc] = um/s
  uc  = kBT / (drag_c * Rb);
  // characteristic force, [fc] = N
  fc  = kBT / Rb;
  // non-dimensional viscosity of fluid, [muc] = 1
  muc = 1. / (6. * PI);
  // read parameters related to module_poisson
  if (module_poisson)
  {
    // Relative permittivity of the fluid
    epsilon = input_file("epsilon", 1.);
    // Characteristic electrostatic potential (V)
    phi0    = elementary_charge / (4. * PI * epsilon * epsilon_0 * Rb * 1E-6);
    // characteristic electrostatic field, [efield_c] = V/um
    efield0 = phi0 / (Rb);
    // characteristic volume charge density, [charge_rho0] = C/um^3
    charge_rho0 = elementary_charge / (Rb * Rb * Rb);
    // characteristic surface charge density, [charge_sigma0] = C/um^2
    charge_sigma0 = elementary_charge / (Rb * Rb);
  } // end if module_poisson

  // read parameters related to
  if (module_np)
  {
    // name of all ion species
    ion_name.resize(input_file.vector_variable_size("ion_name"));
    // diffusivity of all ion species (unit=um^2/s)
    ion_diffusivity.resize(input_file.vector_variable_size("ion_diffusivity"));
    // valence of all ion species (1)
    ion_valence.resize(input_file.vector_variable_size("ion_valence"));
    // real in data for individual ions
    if (ion_name.size() ==ion_diffusivity.size()
      and ion_name.size() == ion_valence.size())
    {
      for (unsigned int j = 0; j < ion_name.size(); j++)
      {
        ion_name[j] = input_file("ion_name", "", j);
        // get ion_diffusivity and normalize it using characteristic
        // diffusion coefficient
        ion_diffusivity[j] = input_file("ion_diffusivity", 0.0, j) / Db;
        ion_valence[j] = input_file("ion_valence", 0, j);
      }
    }
    else
    {
      ss << "Error: ion parameters ('ion_name', 'ion_diffusivity', "
            "'ion_valence') are not correctly set. "
            "Exiting ...\n";
      PMToolBox::output_message(ss, *comm_in);
      libmesh_error();
    }
    // np system relaxation time for initialization (unit = tc)
    np_system_relaxation_time = input_file("np_system_relaxation_time",
      2. / (*std::min_element(ion_diffusivity.begin(), ion_diffusivity.end())));
  } // end if module_np

  // print out physical parameters information related to Stokes System
  ss << "\n ##########################################################\n"
     << " #                  System Physical Parameters             \n"
     << " ##########################################################\n\n"
     << "   temperature           T   = " << T << "(K)\n"
     << "   viscosity             mu  = " << viscosity << " (cP = N*s/um^2)\n"
     << "   Energy unit           kBT = " << kBT << " (N*um = N*um)\n"
     << "   Radius of the bead     a  = " << Rb << " (um)\n"
     << "   bead diffusivity      Db  = " << Db << "(um^2/s)\n"
     << "   HI Drag coefficient  zeta = 6*PI*mu*a = " << drag_c << " (N*s/um)\n"
     << "   ksi = sqrt(PI)/(3a)       =  " << std::sqrt(PI) / (3. * Rb) << " (1/um)\n"
     << "   ------------> The characteristic variables:\n"
     << "   characteristic time          = " << tc << " (s)\n"
     << "   characteristic velocity      = " << uc << " (um/s)\n"
     << "   characteristic force         = " << fc << " (N)\n";
  // print out physical parameters information related to Poisson System
  if (module_poisson)
  {
    ss << "   characteristic electrostatic potential = " << phi0 <<  " (uV)\n"
       << "   relative dielectric permittivity = " << epsilon << " (1)\n"
       << "   characteristic electrostatic potential: phi0 = " << phi0 << " V\n"
       << "   characteristic electrostatic potential: phi0 = " << phi0 << " V\n"
       << "   characteristic electrostatic field: efield0 = " << efield0 <<" V/um\n"
       << "   characteristic volume charge density: charge_rho0 = " <<
       charge_rho0 << " C/um^3\n"
       << "   characteristic surface charge density: charge_sigma0 = " <<
       charge_sigma0 << " C/um^2\n";
  } // end if module_poisson

  // print out physical parameters information related to NP System
  if (module_np)
  {
    ss << "   characteristic ion concentration = " << c0 << " (M)\n"
       << "   relaxation time for system initialization, "
          "np_system_relaxation_time = " << np_system_relaxation_time
          << " (unit=characteristic time tc ~= O(R_ion/R_bead) ~= O"
             "(D_bead/D_ion))\n"
       << "   ------------> parameters of all ion species:\n";
    for (int i=0; i<ion_name.size(); i++)
    {
      ss << "   " << ion_name[i] << " : "
         << "ion diffusivity = " << ion_diffusivity[i]
         << "(unit = bead diffusivity Db ~= O(R_bead/R_ion)), "
         << "valence = " << ion_valence[i] << " (unit = 1)";
    }

  } // end if module_np
  // output message
  PMToolBox::output_message(ss, *comm_in);
} // end read_physical_parameter()

// ====================================================================
void Copss::read_domain_info()
{
  dim = input_file("dimension", 3);

  // =============== wall type and wall params
  wall_type = input_file("wall_type", "not_defined");

  if (wall_type == "slit" or wall_type == "sphere") {
    wall_params.resize(input_file.vector_variable_size(wall_type));

    for (unsigned int j = 0; j < wall_params.size();
         j++) wall_params[j] = input_file(wall_type, 0.0, j);
  }
  else {
    ss << "Error: wall_type (" << wall_type <<") is not supported !!! (in "
                                               "check_wall)";
    PMToolBox::output_message(ss, *comm_in);
    libmesh_error();
  }

  // =============== periodicity
  periodicity.resize(input_file.vector_variable_size("periodicity"));

  if (periodicity.size() != dim) {
    ss << "Error: periodicity has to be defined for all dimensions";
    PMToolBox::output_message(ss, *comm_in);
    libmesh_error();
  }

  for (unsigned int i = 0; i < periodicity.size();
       i++) periodicity[i] = input_file("periodicity", false, i);

  if (periodicity[0] == true and periodicity[1] == true and periodicity[2] ==
      true) {
    ss << "Error: The box cannot be periodic on all directions at the same "
         "time. (required by FEM). Exiting...";
    PMToolBox::output_message(ss, *comm_in);
    libmesh_error();
  }

  // ============== inlet
  inlet.resize(input_file.vector_variable_size("inlet"));

  if (inlet.size() != dim) {
    ss << "Error: inlet has to be defined for all dimensions. Exiting...",
    PMToolBox::output_message(ss, *comm_in);
    libmesh_error();
  }

  for (unsigned int i = 0; i < inlet.size(); i++) {
    inlet[i] = input_file("inlet", false, i);

    if (inlet[i] == true and periodicity[i] == true) {
      ss << "Error: A inlet direction has to be non-periodic. Exiting...";
      PMToolBox::output_message(ss, *comm_in);
      libmesh_error();
    }
  }

  // ============== inlet pressure
  inlet_pressure.resize(input_file.vector_variable_size("inlet_pressure"));

  for (unsigned int i = 0; i < inlet_pressure.size(); i++)
  {
    inlet_pressure[i] = input_file("inlet_pressure", 0., i);
  }
  // ============== shear or not on each boundary pair
  shear.resize(input_file.vector_variable_size("shear"));
  if (shear.size() != dim) {
    ss << "Error: shear has to be defined for all boundary pairs. Exiting...";
    PMToolBox::output_message(ss, *comm_in);
    libmesh_error();
  }
  for (unsigned int i = 0; i < shear.size(); i++) {
    shear[i] = input_file("shear", false, i);

    if (shear[i] == true and periodicity[i] == true) {
      ss << "Error: A shear direction has to be non-periodic. Exiting...";
      PMToolBox::output_message(ss, *comm_in);
      libmesh_error();
    }
  }

  // =============== shear stress on each boundary
  shear_rate.resize(input_file.vector_variable_size("shear_rate"));

  for (unsigned int i = 0; i < shear_rate.size(); i++) shear_rate[i] = input_file(
      "shear_rate",
      0.0,
      i);

  // =============== shear direction on each boundary
  shear_direction.resize(input_file.vector_variable_size("shear_direction"));

  for (unsigned int i = 0; i < shear_direction.size();
       i++) shear_direction[i] = input_file("shear_direction", 0, i);

  // boundary condition for Poisson Module
  if (module_poisson)
  {
    // Read Dirichlet BC: boundary_ids and boundary_values (fixed
    // electrostatic potential at boundaries)
    boundary_id_dirichlet_poisson.resize(input_file.vector_variable_size(
            "boundary_id_dirichlet_poisson"));
    boundary_value_dirichlet_poisson.resize(input_file.vector_variable_size(
            "boundary_value_dirichlet_poisson"));
    // validation
    if (boundary_id_dirichlet_poisson.size() !=
        boundary_value_dirichlet_poisson.size()) {
      ss << "Error: size of 'boundary_id_dirichlet_poisson' does not match "
            "with boundary_value_dirichlet_poisson'. Exiting ..." << "\n";
      PMToolBox::output_message(ss, *comm_in);
      libmesh_error();
    }
    else{
      for (unsigned int i = 0; i < boundary_id_dirichlet_poisson.size(); i++) {
        boundary_id_dirichlet_poisson[i] = input_file(
                "boundary_id_dirichlet_poisson", 0, i);
        boundary_value_dirichlet_poisson[i] = input_file(
                "boundary_value_dirichlet_poisson", 0.0, i);
      }
    } // end else

    // Read Neumann BC: boundary_ids and boundary values (fixed surface
    // charge density at boundaries)
    boundary_id_neumann_poisson.resize(input_file.vector_variable_size(
            "boundary_id_neumann_poisson"));
    boundary_value_neumann_poisson.resize(input_file.vector_variable_size(
            "boundary_value_neumann_poisson"));
    // validation
    if (boundary_id_neumann_poisson.size() !=
        boundary_value_neumann_poisson.size()) {
      ss << "Error: size of 'boundary_id_neumann_poisson' does not match with "
          << "'boundary_value_neumann_poisson'. Exiting ..." << "\n";
      PMToolBox::output_message(ss, *comm_in);
      libmesh_error();
    }
    else{
      for (unsigned int i = 0; i < boundary_id_neumann_poisson.size(); i++) {
        boundary_id_neumann_poisson[i] = input_file(
          "boundary_id_neumann_poisson", 0, i);
        boundary_value_neumann_poisson[i] = input_file(
          "boundary_value_neumann_poisson", 0.0, i);
      }
    } // end else
  }

  // Boundary condition for Nernst-Planck system
  if (module_np)
  {
    // if check charge neutrality for BC
    check_charge_neutrality = input_file("check_charge_neutrality", true);

    // read Dirichlet Boundary ids of NP system
    if (input_file.have_variable("boundary_id_dirichlet_np"))
    {
      boundary_id_dirichlet_np.resize(input_file.vector_variable_size(
        "boundary_id_dirichlet_np"));
      for (int i = 0; i < boundary_id_dirichlet_np.size(); i++) {
        boundary_id_dirichlet_np[i] = input_file("boundary_id_dirichlet_np",
                                                 0, i);
      }
    }
    else
    {
      ss << "Error: boundary_id_dirichlet_np is not given. Exiting...";
      PMToolBox::output_message(ss, *comm_in);
      libmesh_error();
    }

    // get Dirichlet boundary value for all ion species
    boundary_value_dirichlet_np.resize(ion_name.size());
    for (int ion_id=0; ion_id<ion_name.size(); ion_id++)
    {
      // check the boundary value is given
      if (input_file.have_variable
      ("boundary_value_dirichlet_np_"+ion_name[ion_id]))
      {
        // temporary placeholder of boundary value for this ion
        std::vector<Real> tmp_c;
        tmp_c.resize(input_file.vector_variable_size(
          "boundary_value_dirichlet_np_"+ion_name[ion_id]));
        // check the size is equal to the size of boundary_id_dirichlet_np
        if (tmp_c.size() != boundary_id_dirichlet_np.size())
        {
          ss << "Error: boundary_value_dirichlet_np_" << ion_name[ion_id]
             << " has the wrong size. Exiting...";
          PMToolBox::output_message(ss, *comm_in);
          libmesh_error();
        }
        else
        {
          // read boundary value to tmp_c
          for (int i=0; i<tmp_c.size(); i++)
          {
            tmp_c[i] = input_file
              ("boundary_value_dirichlet_np_"+ion_name[ion_id], 0.0, i);
          }
          // put tmp_c to boundary_value_dirichlet_np[ion_id]
          boundary_value_dirichlet_np[ion_id] = tmp_c;
        }
      }
      else
      {
        ss << "Error: boundary_value_dirichlet_np_" << ion_name[ion_id]
           << " is not given. Exiting...";
        PMToolBox::output_message(ss, *comm_in);
        libmesh_error();
      }
    } // end for loop over ion_id

    // check charge neutrality for Boundary Values
    if (check_charge_neutrality)
    {
      for (int i=0; i<boundary_id_dirichlet_np.size(); i++)
      {
        Real total_c = 0;
        for (int ion_id=0; ion_id<boundary_value_dirichlet_np.size(); ion_id++)
        {
          total_c += boundary_value_dirichlet_np[ion_id][i] *
                     ion_valence[ion_id];
        }
        if (abs(total_c) > 1.e-6)
        {
          ss << "Error: total ion charge (concentration*valence) is not neutral"
                " at boundary id = "<< boundary_id_dirichlet_np[i]
             << ". Exiting...";
          PMToolBox::output_message(ss, *comm_in);
          libmesh_error();
        }
      }
    }
  } // end if molule_np

  // print out Geometry information
  ss << "##########################################################\n"
     << "#                  Geometry information                  \n"
     << "##########################################################\n"
     << "-----------> Dimension: " << dim << "\n"
     << "-----------> Wall type: " << wall_type << "\n"
     << "-----------> Wall size parameters: ";

  for (int i = 0; i < wall_params.size(); ++i)
    ss << wall_params[i] << "   ";
  ss << "\n";

  ss << "-----------> Periodicity of the box: ";
  for (int i = 0; i < dim; i++)
    ss << std::boolalpha << periodicity[i] << ", ";
  ss << "\n";

  ss << "-----------> Inlet/Outlet of the box: ";
  for (int i = 0; i < dim; i++)
    ss << inlet[i] << "(pressure = " << inlet_pressure[i] << " ), ";
  ss << "\n";

  ss << "-----------> shear on boundary pairs: ";
  for (int i = 0; i < dim; i++)
    ss << shear[i] << "(shear_rate = " << shear_rate[i] << " ), ";
  ss << "\n";
  // Poisson System BC
  if (module_poisson)
  {
    ss << "-----------> Dirichlet boundary condition for Poisson (Potential)"
          ":\n";
    for (int i = 0; i < boundary_id_dirichlet_poisson.size(); i++)
      ss << "Boundary ID " << boundary_id_dirichlet_poisson[i]
         << " : value = " << boundary_value_dirichlet_poisson[i] << "\n";
    ss << "\n";
    ss << "-----------> Neumann boundary condition for Poisson (Surface Charge Density):\n";
    for (int i = 0; i < boundary_id_neumann_poisson.size(); i++)
      ss << "Boundary ID " << boundary_id_neumann_poisson[i]
         << " : value = " << boundary_value_neumann_poisson[i] << "\n";
    ss << "\n";
  } // end if module_poisson

  // NP system BC
  if (module_np)
  {
    ss << "-----------> Dirichlet boundary condition for Nernst-Planck "
          "(ion concentration):\n";
    for (int i=0; i<boundary_id_dirichlet_np.size(); i++)
    {
      // print concentration of all ion species at this boundary
      ss << "At boundary ID " << boundary_id_dirichlet_np[i] << " : ";
      for (int ion_id=0; ion_id<ion_name.size(); ion_id++)
      {
        ss << ion_name[ion_id] << "(" <<
        boundary_value_dirichlet_np[ion_id][i] <<"[M]); ";
      }
      ss << "\n";
    }
  } // end if module_np

  // Domain mesh information
  ss << "##########################################################" << "\n"
     << "#                  Domain Mesh information                " << "\n"
     << "##########################################################" << "\n";
  generate_mesh = input_file("generate_mesh", true);
  if (generate_mesh)
  {
    n_mesh.resize(input_file.vector_variable_size("n_mesh"));
    for (unsigned int i = 0; i < n_mesh.size(); i++) n_mesh[i] = input_file(
        "n_mesh", 1, i);
    ss << "------------> Generate Mesh using COPSS: n_mesh = " << n_mesh[0] <<
      ";" << n_mesh[1] << ";" << n_mesh[2] << "\n";
  }
  else {
    domain_mesh_file = input_file("domain_mesh_file", "nothing");
    ss << "------------> Load mesh file from " << domain_mesh_file.c_str() << "\n";
  } // end else

  // print out information
  PMToolBox::output_message(ss, *comm_in);
}   // end read_domain_info()

// ====================================================================
void Copss::read_force_info() {
  // read particle-particle force types
  numForceTypes = input_file.vector_variable_size("force_field");
  forceTypes.resize(numForceTypes);
  forces.resize(numForceTypes);

  for (unsigned int i = 0; i < numForceTypes; i++) {
    forceTypes[i] = input_file("force_field", "nothing", i);
    std::vector<Real> params(input_file.vector_variable_size(forceTypes[i]));

    if (forceTypes[i] != "nothing") {
      for (unsigned int j = 0; j < params.size(); j++) {
        params[j] = input_file(forceTypes[i], 0.0, j);
      }
    }
    forces[i].first  = forceTypes[i];
    forces[i].second = params;
  }
  ss << "##########################################################\n" 
     << "#    Force field (fixes)                \n"
     << "##########################################################\n";

  for (int i = 0; i < numForceTypes; i++) {
    ss << "-----------> " << forces[i].first << "= '";
    for (int j = 0; j < forces[i].second.size(); j++) {
      ss << forces[i].second[j] << "  ";
    }
    ss << "'" << "\n";
  }
  PMToolBox::output_message(ss, *comm_in);
} // end read_force_info()

// ====================================================================
void Copss::read_solver_info() {
  max_linear_iterations = input_file("max_linear_iterations", 100);
  linear_solver_rtol    = input_file("linear_solver_rtol", 1E-6);
  linear_solver_atol    = input_file("linear_solver_atol", 1E-6);
  user_defined_pc       = input_file("user_defined_pc", true);
  schur_user_ksp        = input_file("schur_user_ksp", false);
  schur_user_ksp_rtol   = input_file("schur_user_ksp_rtol", 1E-6);
  schur_user_ksp_atol   = input_file("schur_user_ksp_atol", 1E-6);
  schur_pc_type         = input_file("schur_pc_type", "SMp");

  solver_stokes = input_file("solver_stokes", "superLU_dist");

  if (solver_stokes == "superLU_dist") {
    solver_type_stokes = superLU_dist;
    user_defined_pc    = false;
  }
  else if (solver_stokes == "field_split") {
    solver_type_stokes = field_split;
    user_defined_pc    = true;
  }
  else {
    solver_type_stokes = user_define;
  }

  // solver info for Poisson system
  if (module_poisson)
  {
    solver_poisson = input_file("solver_poisson", "superLU_dist");
    if (solver_poisson == "superLU_dist") {
      solver_type_poisson = superLU_dist;
      user_defined_pc_poisson = false;
    }
    else if (solver_poisson == "field_split") {
      solver_type_poisson = field_split;
      user_defined_pc_poisson = false;
      ss << "Warning: solver_poisson is set to be "
            "'field_split' but user-defined preconditioner for Poisson System "
            "is not yet supported. Setting user_defined_pc_poisson to be "
            "false...\n";
      PMToolBox::output_message(ss, *comm_in);
    }
    else {
      solver_type_poisson = user_define;
    }
  } // end if module_poisson

  // solver info for Nernst Planck System
  if (module_np)
  {
    solver_np = input_file("solver_np", "superLU_dist");

    if (solver_np == "superLU_dist") {
      solver_type_np = superLU_dist;
      user_defined_pc_np     = false;
    }
    else if (solver_np == "field_split") {
      solver_type_np = field_split;
      user_defined_pc_np     = false;
      ss << "Warning: solver_np is set to be "
            "'field_split' but user-defined preconditioner for Nernst Planck "
            "System is not yet supported. Setting user_defined_pc_np to be "
            "false...\n";
      PMToolBox::output_message(ss, *comm_in);
    }
    else {
      solver_type_np = user_define;
    }
  }

  ss << "##########################################################\n"
     << "#                 Solver information                      \n"
     << "##########################################################\n"
     << "-----------> Stokes solver type = " << solver_stokes << "\n";

  if (solver_stokes == "field_split") {
    ss << "-----------> FieldSplit Schur Complement Reduction Solver\n"
       << "-----------> schur_pc_type = " << schur_pc_type << "\n";

    if (schur_user_ksp) {
      ss << "----------->  user defined KSP is used for Schur Complement!\n" 
         << "----------->  KSP rel tolerance for Schur Complement solver is = "
         << schur_user_ksp_rtol << "\n"
         << "----------->  KSP abs tolerance for Schur Complement solver is = "
         << schur_user_ksp_atol << "\n";
    }
  }

  if (module_poisson)
  {
    ss << "-----------> Poisson solver type = " << solver_poisson << "\n";
  }

  if (module_np)
  {
    ss << "-----------> Nernst Planck solver type = " << solver_np << "\n";
  }

  PMToolBox::output_message(ss, *comm_in);
} // end read_solver_info()

// ====================================================================
void Copss::read_chebyshev_info() {
  with_hi       = input_file("with_hi", false);
  with_brownian = input_file("with_brownian", false);
  restart       = input_file("restart", false);
  compute_eigen = input_file("compute_eigen", false);
  eig_min       = 0.;
  eig_max       = 0.;

  if (with_hi and with_brownian and restart and !compute_eigen) {
    read_eigen = true;
    this->read_restart_eigenvalue();
  }
  max_n_cheb = input_file("max_n_cheb", 10);
  tol_cheb   = input_file("tol_cheb", 0.1);
  eig_factor = input_file("eig_factor", 1.05);
  tol_eigen  = input_file("tol_eigen", 0.01);
  if (with_brownian and with_hi) {
    ss << "##########################################################\n"
       << "#   Chebyshev information (only for brownian System)      \n"
       << "##########################################################\n"
       << "-----------> compute eigen values  = " << compute_eigen << "\n"
       << "-----------> max number of chebyshev polynomial = " << max_n_cheb <<"\n"
       << "-----------> tolerance of chebyshev polynomial = " << tol_cheb << "\n"
       << "-----------> factor of eigenvalues range = " << eig_factor << "\n"
       << "-----------> tolerance of eigenvalues convergence = " << tol_eigen << "\n";

    if (read_eigen) {
      ss << "-----------> read eigen value range from 'out.eigenvalue = '" 
         << read_eigen << "\n"
         << "-----------> initial eigen value range = ( " 
         << eig_min << ", " << eig_max << " )" << "\n";
    }
  }
  PMToolBox::output_message(ss, *comm_in);
  // build elem-particle neighbor list if either with_hi or module_poisson or
  // module_np is true
  build_elem_neighbor_list = with_hi || module_poisson || module_np;
} // end read_chebyshev_info()

// ====================================================================
void Copss::read_run_info() {
  // ############## Without Brownian ###############################
  // For polymer_chain and bead: maximum displacement (non dimensional) of one
  // step = 0.1 * fluid mesh size minimum (hmin)
  // ############## With Brownian ##################################
  // For polymer_chain: maximum displacement (non dimensional) of one step = 0.1
  // * Ss2/Rb/Rb
  // For bead: maximum displacement (non_dimensional) of one step = 0.1
  if (with_brownian) {
    random_seed = input_file("random_seed", 111);
  }
  max_dr_coeff.resize(input_file.vector_variable_size("max_dr_coeff"));
  if (max_dr_coeff.size() == 3) {
    for (unsigned int i = 0; i < max_dr_coeff.size(); i++) {
      max_dr_coeff[i] = input_file("max_dr_coeff", 0.1, i);
    }
  }
  else {
    ss << "Error: max_dr_coeff needs three parameters: t_milestone, "
           "max_dr_coeff before milestone, max_dr_coeff after milestone. "
           "Exiting...";
    PMToolBox::output_message(ss, *comm_in);
    libmesh_error();
  }
  adaptive_dt = input_file("adaptive_dt", true);
  if (restart) {
    this->read_restart_time();
    istart = restart_step;
  }
  else {
    restart_step = -1;
    istart       = 0;
    o_step       = 0;
    real_time    = 0.;
  }

  // neighbor list update
  neighbor_list_update_flag      = true; // set this to be ture before
                                         // simulation
  update_neighbor_list_everyStep = input_file("update_neighbor_list_everyStep",
                                              true);
  // total number of steps
  nstep = input_file("nstep", 1);
  // free-draining steps when Chebysheve cannot converge even when eigenvalue is
  // recomputed
  n_relax_step = input_file("n_relax_step", 10);
  // write interval
  write_interval = input_file("write_interval", 1);
  // output file flags
  output_file.resize(input_file.vector_variable_size("output_file"));
  for (unsigned int i = 0; i < output_file.size(); i++) {
    output_file[i] = input_file("output_file", "not defined", i);
  }
  // debug info
  debug_info = input_file("debug_info", false);
  // output precision
  o_precision = input_file("output_precision", 6);
  // write to screen
  ss << "\n##########################################################\n"
     << "#                 Run information                      \n"
     << "##########################################################\n\n"
     << "-----------> adaptive_dt: " << adaptive_dt << "\n"
     << "-----------> simulation milestone t_milestone = " << max_dr_coeff[0] << "\n"
     << "-----------> before t_milestone, max_dr_coeff = " << max_dr_coeff[1] 
     << "; after milestone, max_dr_coeff = " << max_dr_coeff[2] << "\n"
     << "-----------> debug_info: " << debug_info << "\n"
     << "-----------> with_hi: "  << with_hi << "\n";

  if (with_brownian) {
    ss << "-----------> with_brownian: " << with_brownian << "\n"
       << "-----------> random seed: " << random_seed << "\n";
  }
  ss << "-----------> Restart mode: " << restart << "\n";
  if (restart) {
    ss << "-----------> Restart step: " << restart_step << "\n"
       << "-----------> Restart from real_time: " << real_time << "\n"
       << "-----------> Restart particle data read from o_step: " << o_step << "\n";
  }
  ss << "-----------> nstep: " << nstep << "\n"
     << "-----------> write_interval: " << write_interval << "\n"
     << "-----------> write output file: " << "\n";

  for (int i = 0; i < output_file.size(); i++) {
    ss << "                              " << output_file[i] << "\n";
  }
  ss << "-----------> output precision: " << o_precision << "\n";
  PMToolBox::output_message(ss, *comm_in);
} // end read_run_info()

// ============================================================================
void Copss::read_restart_time()
{
  // Open the local file and check the existance
  ss << "###Read time output: out.time" << "\n";
  std::ifstream fin;
  fin.open("out.time", std::ios_base::ate);
  char c = '\0';
  if (fin.is_open()) {
    const unsigned int length = fin.tellg();
    for (int i = length - 2; i > 0; i--) {
      fin.seekg(i);
      c = fin.get();
      if ((c == '\r') || (c == '\n')) // new line?
        break;
    }
    fin >> restart_step >> o_step >> real_time;
    fin.close();
    ss << "restart time_step: " << restart_step
       << "restart o_step: " << o_step
       << "restart real_time: "<< real_time;
    PMToolBox::output_message(ss, *comm_in); 
  } // end if
  else {
    ss << "Error: read_restart_time can NOT read the time output: out.time. "
          "Exiting...";
    PMToolBox::output_message(ss, *comm_in);
    libmesh_error();
  }
}

// ============================================================================
void Copss::read_restart_eigenvalue()
{
  // Open the local file and check the existance
  ss << "###Read eigenvalue from previous simulation to restart: out.eigenvalue\n";
  std::ifstream fin;
  fin.open("out.eigenvalue", std::ios_base::ate);
  std::string tmp;
  char c = '\0';
  if (fin.is_open()) {
    const unsigned int length = fin.tellg();
    for (int i = length - 2; i > 0; i--) {
      fin.seekg(i);
      c = fin.get();
      if ((c == '\r') || (c == '\n')) // new line?
        break;
    }
    fin >> eig_min >> eig_max;
    fin.close();
    ss << "restart eig_min: " << eig_min
       << "restart eig_max: " << eig_max;
    PMToolBox::output_message(ss, *comm_in);
  } // end if
  else {
    ss << "Error: read_restart_eigenvalue() can NOT read the time output: out"
         ".eigenvalue. Exiting...";
    PMToolBox::output_message(ss, *comm_in);
    libmesh_error();
  }
}

// ============================================================================
void Copss::create_domain_mesh()
{
  if (dim == 2) 
  {
    ss << "Error::Copss::create_domain_mesh() only works for 3D systems.\
      2D simulation needs extra implementation. Exiting...";
    PMToolBox::output_message(ss, *comm_in);
    libmesh_error();
  }
  mesh = new SerialMesh(*comm_in);
  if (generate_mesh) 
  {
    if (wall_type == "slit") 
    {
      const std::vector<Real> mesh_size = PMToolBox::mesh_size(*mesh);
      const Real meshsize_x = (wall_params[1] - wall_params[0]) / Real(n_mesh[0]);
      const Real meshsize_y = (wall_params[3] - wall_params[2]) / Real(n_mesh[1]);
      const Real meshsize_z = (wall_params[5] - wall_params[4]) / Real(n_mesh[2]);
      hminf = std::min(meshsize_x, meshsize_y);
      hminf = std::min(hminf, meshsize_z);
      hmaxf = std::max(meshsize_x, meshsize_y);
      hmaxf = std::max(hmaxf, meshsize_z);
      // Build HEX20 3D mesh
      MeshTools::Generation::build_cube(*mesh,
                                        n_mesh[0],
                                        n_mesh[1],
                                        n_mesh[2],
                                        wall_params[0],
                                        wall_params[1],
                                        wall_params[2],
                                        wall_params[3],
                                        wall_params[4],
                                        wall_params[5],
                                        HEX20); // HEX20/27
     ss << "\n##########################################################\n"
        << "#                 The created mesh information              \n"
        << "########################################################## \n\n"
        << "   minimum mesh size of fluid: hminf = " << hminf << "\n"
        << "   maximum mesh size of fluid: hmaxf = " << hmaxf << "\n";
      // need to first modify the boundary_id to a temporary code
      // otherwise we might make mistake like change 1 -> 3 -> 5 -> 1
      for (int i = 0; i < copss_slitMesh_boundary_id.size(); i++) {
        libMesh::MeshTools::Modification::change_boundary_id(*mesh,
          copss_slitMesh_boundary_id[i], -copss_slitMesh_boundary_id[i]);
      }
      for (int i = 0; i < copss_slitMesh_boundary_id.size(); i++) {
        libMesh::MeshTools::Modification::change_boundary_id(*mesh,
          -copss_slitMesh_boundary_id[i], slitMesh_boundary_id[i]);              
      }
    } 
    else 
    {
      ss << "Error: COPSS only supports generating domain mesh for 'slit' "
             "wall. Please load the domain mesh file for other wall types. "
             "Exiting ...";
      PMToolBox::output_message(ss, *comm_in);
      libmesh_error();
    }
  } 
  else 
  {
    if (domain_mesh_file != "nothing") {
      mesh->read(domain_mesh_file);
      mesh->all_second_order();
      mesh->prepare_for_use();
      const std::vector<Real> mesh_size = PMToolBox::mesh_size(*mesh);
      hminf = mesh_size[0];
      hmaxf = mesh_size[1];
      ss << "##########################################################" << "\n"
         << "#              The Read-in mesh information               " << "\n"
         << "##########################################################" << "\n"
         << "   minimum mesh size of fluid: hminf = " << hminf << "\n"
         << "   maximum mesh size of fliud: hmaxf = " << hmaxf << "\n";
    }
    else 
    {
      PMToolBox::output_message("Error: 'domain_mesh_file' needs to be specified. Exiting ...",
        *comm_in);
      libmesh_error();
    }
  }
  // MeshBase::const_element_iterator el =
  //   mesh->active_local_elements_begin();
  // const MeshBase::const_element_iterator end_el =
  //   mesh->active_local_elements_end();
  // for (; el != end_el; ++el)
  // {
  //   // Store a pointer to the element we are currently working on.
  //   const Elem *elem           = *el;
  //   const unsigned int elem_id = elem->id();
  //   std::cout<<"elem_id = "<<elem_id <<"\n";
  //   for (unsigned int s = 0; s < elem->n_sides(); s++)
  //   {
  //     std::cout<<"side id = " << s <<", boundary_ids = ";
  //     std::vector<boundary_id_type> boundary_ids;
  //     mesh->get_boundary_info().boundary_ids(elem, s, boundary_ids);
  //     for (int i = 0; i<boundary_ids.size(); i++){
  //       std::cout<<boundary_ids[i]<<"; ";
  //     }
  //     std::cout<<"\n";
  //     // If this side is on the boundary
  //   }
  // } // end for elem-loop
  search_radius_p = 4. / alpha;
  search_radius_e = 0.5 * hmaxf + 4. / alpha;
  // print mesh info
  mesh->print_info();
  mesh->get_boundary_info().print_summary();
  PMToolBox::output_message(ss, *comm_in);
} // end function

// ============================================================================
void Copss::create_periodic_boundary() {
  if (wall_type == "slit") {
    const Point bbox_pmin(wall_params[0], wall_params[2], wall_params[4]);
    const Point bbox_pmax(wall_params[1], wall_params[3], wall_params[5]);

    // construct PMPeriodicBoundary class using info above
    pm_periodic_boundary = new PMPeriodicBoundary(bbox_pmin,
                                                  bbox_pmax,
                                                  periodicity,
                                                  inlet,
                                                  inlet_pressure);
  }
  else if (wall_type == "sphere") {
    // check: No PBC, No inlet/outlet
    bool error_flag = false;
    for (int i = 0; i < dim; i++) {
      if (periodicity[i] or inlet[i] or shear[i]) error_flag = true;
    }
    if (error_flag) {
      ss << "Error: spherical domain cannot have PBC or inlet/outlet, Please "
             "check control file. Exiting...";
      PMToolBox::output_message(ss, *comm_in);
      libmesh_error();
    }
    // spherical domain cannot have PBC, but we need to create a PBC object 
    // using a non-existed cubic box to keep COPSS running !"
    const Point bbox_pmin(-Real(wall_params[0] / 2.),
                          -Real(wall_params[0] / 2.),
                          -Real(wall_params[0] / 2.));
    const Point bbox_pmax(Real(wall_params[0] / 2.),
                          Real(wall_params[0] / 2.),
                          Real(wall_params[0] / 2.));

    // construct PMPeriodicBoundary class using info above
    pm_periodic_boundary = new PMPeriodicBoundary(bbox_pmin,
                                                  bbox_pmax,
                                                  periodicity,
                                                  inlet,
                                                  inlet_pressure);
  }
  else {
    ss << "Error: COPSS::create_periodic_boundary() only support 'slit' or 'sphere' wall_type for now !";
    PMToolBox::output_message(ss, *comm_in);
  }
} // end function

// =============================================================================
EquationSystems Copss::create_equation_systems()
{
  // Initialize equation_systems object using the 'mesh' we created before
  // equation_systems = new EquationSystems(*mesh);
  ss << "==> Initialize equation_systems object using the 'mesh' we created "
        "before";
  PMToolBox::output_message(ss, *comm_in);
  EquationSystems equation_systems(*mesh);

  // Add 'Stokes' system (of PMSystemStokes) to the 'equation_systems'
  ss << "==> Add 'Stokes' system (of PMSystemStokes) to the 'equation_systems'";
  PMToolBox::output_message(ss, *comm_in);
  PMSystemStokes& system = equation_systems.add_system<PMSystemStokes>("Stokes");
  
  // Add variables to 'Stokes' system"
  ss << "----> add variables to the system";
  PMToolBox::output_message(ss, *comm_in);
  u_var = system.add_variable("u", SECOND);
  v_var = system.add_variable("v", SECOND);
  if (dim == 3) w_var = system.add_variable("w", SECOND);
  const unsigned int p_var = system.add_variable("p", FIRST);

  ss << "----> attach object_mesh to the system";
  PMToolBox::output_message(ss, *comm_in);
  this->attach_object_mesh(system);

  // attach period boundary to the system
  ss << "----> add period boundary conditions to the system";
  PMToolBox::output_message(ss, *comm_in);
  // don't include "p" in the periodic variables
  std::vector<std::string> periodic_vars {"u", "v", "w"};
  this->attach_period_boundary(system, periodic_vars);

  // Initialize the Preconditioning matrix for saddle point problems if
  // required. Initialize the equation system and zero the preconditioning
  // matrix.
  if (user_defined_pc) {
    ss << "----> add Preconditioner when user_defined_pc is true";
    PMToolBox::output_message(ss, *comm_in);
    system.add_matrix("Preconditioner");
  }

  // Poisson Module
  if (module_poisson)
  {
    // Add 'Poisson' system to the 'equation_systems'
    ss << "==> Add 'Poisson' system to the 'equation_systems'";
    PMToolBox::output_message(ss, *comm_in);
    PMSystemPoisson& system_poisson =
      equation_systems.add_system<PMSystemPoisson>("Poisson");

    // Add variables to the system
    ss << "----> add variables to the system";
    PMToolBox::output_message(ss, *comm_in);
    phi_var = system_poisson.add_variable("phi", SECOND);

    // attach object_mesh to the system
    ss << "----> attach object_mesh to the system";
    PMToolBox::output_message(ss, *comm_in);
    this->attach_object_mesh(system_poisson);

    // attach period boundary to the system
    ss << "----> add period boundary conditions to the system";
    PMToolBox::output_message(ss, *comm_in);
    std::vector<std::string> periodic_vars {"phi"};
    this->attach_period_boundary(system_poisson, periodic_vars);
  }

  // Nernst-Planck Module
  if (module_np)
  {
    // fixme: remove the following if condition after final implementation
    if(with_brownian & with_hi){
      PMToolBox::output_message("Error: cannot support module_np with both "
                                "Brownian and HI on for now because of the "
                                "half-point scheme in Fixman's integration. "
                                "", *comm_in);
      libmesh_error();
    }
    // add a NP system to equation systems for each ion
    for (int ion_id=0 ; ion_id<ion_name.size(); ion_id++)
    {
      std::string np_sys_name = std::string("NP") + ":" + ion_name[ion_id];
      std::string np_var_name = std::string("c") + ":" + ion_name[ion_id];

      // Add This NP system to equation_systems and get a reference
      ss << "==> Add " << np_sys_name << " system to the 'equation_systems'";
      PMToolBox::output_message(ss, *comm_in);
      PMSystemNP& system_np = equation_systems.add_system<PMSystemNP>
        (np_sys_name);

      // Add variables to the system
      ss << "----> add variables to the system";
      PMToolBox::output_message(ss, *comm_in);
      c_var = system_np.add_variable(np_var_name, SECOND);

      // attach object_mesh to the system
      ss << "----> attach object_mesh to the system";
      PMToolBox::output_message(ss, *comm_in);
      this->attach_object_mesh(system_np);

      // attach period boundary to the system
      ss << "----> add period boundary conditions to the system";
      PMToolBox::output_message(ss, *comm_in);
      std::vector<std::string> periodic_vars {np_var_name};
      this->attach_period_boundary(system_np, periodic_vars);

      // attach ion type (id and name) to the system
      ss << "----> attach ion type to the system";
      PMToolBox::output_message(ss, *comm_in);
      system_np.attach_ion_type(ion_id, ion_name[ion_id]);
    }
  }

  /* Initialize the data structures for the equation system. */
  ss << "==> Init equation_systems (libmesh function, to init all systems in "
        "equation_systems)";
  PMToolBox::output_message(ss, *comm_in);
  equation_systems.init();
  // zero the PC matrix, which MUST be done after es.init()
  if (user_defined_pc) {
    ss << "==> (If user_defined_pc) Zero preconditioner matrix";
    PMToolBox::output_message(ss, *comm_in);
    system.get_matrix("Preconditioner").zero();
  }

  // set parameters for equation systems
  ss << "==> Set parameters of equation_systems";
  PMToolBox::output_message(ss, *comm_in);
  this->set_parameters(equation_systems);

  // attach fixes to Stokes system for force field calculation
  ss << "==> attach fixes to 'stokes' system";
  PMToolBox::output_message(ss, *comm_in);
  this->attach_fixes(system);

  // Print information about the mesh and system to the screen.
  ss << "--------------> Print equation systems info" << "\n"
     << "  System has: " << mesh->n_elem() << " elements,\n"
     << "              " << mesh->n_nodes() << " nodes,\n"
     << "              " << equation_systems.n_dofs() << " degrees of freedom.\n"
     << "              " << equation_systems.n_active_dofs() << " active degrees of freedom.\n"
     << "              " << point_mesh->num_particles() << " particles.\n";
  PMToolBox::output_message(ss.str(), *comm_in);
  equation_systems.print_info();

  return equation_systems;
}

// =======================================================================================
void Copss::attach_fixes(PMLinearImplicitSystem& pm_system)
{
  // Fix *f;
  fix_factory = new FixFactory();
  fixes.resize(numForceTypes);

  for (int i = 0; i < numForceTypes; i++) {
    fixes[i] = fix_factory->buildFix(forceTypes[i], pm_system);
  }
  pm_system.attach_fixes(fixes);
}

// ===============================================================================
void Copss::attach_period_boundary(PMLinearImplicitSystem& system,
  std::vector<std::string>& periodic_vars)
{
  ss << "--------> Get dof_map of " << system.name() << " system.";
  PMToolBox::output_message(ss, *comm_in);
  DofMap& dof_map = system.get_dof_map();
  ss << "--------> Create periodicBoundary object and add to dof_map";
  PMToolBox::output_message(ss, *comm_in);
  // loop over each direction
  for (int dim_i = 0; dim_i < dim; dim_i++)
  {
    // proceed if this direction is periodic
    if (periodicity[dim_i])
    {
      ss << "---------------> direction: " << dim_i;
      PMToolBox::output_message(ss, *comm_in);
      // create a PeriodicBoundary object
      PeriodicBoundary pbc(RealVectorValue(
        (wall_params[1] - wall_params[0]) * (dim_i==0),
        (wall_params[3] - wall_params[2]) * (dim_i==1),
        (wall_params[5] - wall_params[4]) * (dim_i==2)));
      // add all variables to the pbc object
      for (unsigned int var_id=0; var_id < periodic_vars.size(); var_id++)
      {
        pbc.set_variable(system.variable_number(periodic_vars[var_id]));
      }
      // set myboundary and pairedboundary
      pbc.myboundary     = slitMesh_boundary_id[2 * dim_i];
      pbc.pairedboundary = slitMesh_boundary_id[2 * dim_i + 1];
      // check if the PBC is valid (particle search radius = 4/alpha <=
      // half of the box length in this direction)
      if (search_radius_p >= (wall_params[2 * dim_i + 1] - wall_params[2 *
      dim_i]) / 2.)
      {
        ss << "Error: GGEM particle search radius is larger than half domain "
              "length in " << dim_i << "direction! "
              << "(search radius = 4./alpha" << search_radius_p
              << ", half domain size Lx/2 = " << (wall_params[dim_i * 2 + 1]
              - wall_params[dim_i * 2]) / 2. << "). Exiting ...";
        PMToolBox::output_message(ss, *comm_in);
        libmesh_error();
      }
      // add pbc object to dof_map
      dof_map.add_periodic_boundary(pbc);
    } // end if (periodicity[dim_i])
  } // end for loop over dim_i
}

// ============================================================================================
void Copss::solve_undisturbed_system(EquationSystems& equation_systems)
{
  // get stokes system from equation systems
  PMSystemStokes& system = equation_systems.get_system<PMSystemStokes>("Stokes");
  // if either with_hi or module_poisson is true, build_elem_neighbor_list will
  // be true
  if (update_neighbor_list_everyStep) neighbor_list_update_flag = true;
  system.reinit_system(neighbor_list_update_flag, build_elem_neighbor_list, "undisturbed");
  if (print_info) 
  {
    PMToolBox::output_message("--------After reinit undisturbed system---------", *comm_in);
    if (comm_in->rank() == 0) point_mesh->print_point_info();
  }
  // if with_hi is true, solve the undisturbed Stokes equation and backup the solution
  if (with_hi) {
    system.solve("undisturbed");
  }
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     write out the equation systems at Step 0 (undisturbed field)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  PMToolBox::output_message("start writing undisturbed solution to file", *comm_in);
  if ((std::find(output_file.begin(), output_file.end(),
                 "equation_systems") != output_file.end()) && !restart)
  {
    ExodusII_IO(*mesh).write_equation_systems("output_equation_systems_undisturbed.e",
      equation_systems);
    // system.write_equation_systems(0, 0., "undisturbed", 
      // "output_equation_systems_undisturbed");
  }
  PMToolBox::output_message("end writing undisturbed solution to file", *comm_in);
}

// ============================================================================================
void Copss::create_brownian_system(EquationSystems& equation_systems)
{
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     Create vectors and Shell Mat for use:
     U0:          particle velocity vector;
     R0/R_mid:    particle position vector;
     dw/dw_mid:   random vector;
     RIN/ROUT:    the initial and intermediate particle postion vector for msd
       output
     RIN will not change, and ROUT excludes pbc
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  brownian_sys = new BrownianSystem(equation_systems);
  brownian_sys->init_petsc_random(&rand_ctx);
  brownian_sys->_create_shell_mat(n_vec, &M);
  brownian_sys->_create_petsc_vec(n_vec, &R0);
  VecDuplicate(R0, &U0);
  VecDuplicate(R0, &R_mid);
  VecDuplicate(R0, &dw_mid);
  brownian_sys->extract_particle_vector(&ROUT, "coordinate", "extract");
  VecDuplicate(ROUT, &RIN);
  VecCopy(ROUT, RIN);                             // RIN = ROUT = the initial
                                                  // position vector
  brownian_sys->set_std_random_seed(random_seed); // random seed

  if (restart)
  {
    // read RIN & ROUT from local file output during the previous simulation
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,
                          "vector_RIN.dat",
                          FILE_MODE_READ,
                          &viewer);
    VecLoad(RIN, viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,
                          "vector_ROUT.dat",
                          FILE_MODE_READ,
                          &viewer);
    VecLoad(ROUT, viewer);
  }
  else
  {
    // write out binary file of RIN, which may be used at restart mode.
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,
                          "vector_RIN.dat",
                          FILE_MODE_WRITE,
                          &viewer);
    VecView(RIN, viewer);
  }
  PetscViewerDestroy(&viewer);
  comm_in->barrier();
}

// ===================================================================
const Real Copss::get_min_dt(EquationSystems &es,
                             std::vector<std::string> sys_names)
{
  // Initialize min_dt as a large FLOAT number
  Real min_dt = std::numeric_limits<double>::max();
  // If sys_names is not given, get all system names from es
  if (sys_names.size()==0)
  {
    unsigned int n_sys = es.n_systems();
    sys_names.resize(n_sys);
    for (unsigned int s_id=0; s_id<n_sys; s_id++)
      sys_names[s_id] = es.get_system(s_id).name();
  }
  // loop over all systems in sys_names and find the minimum dt
  for (unsigned int s_id=0; s_id<sys_names.size(); s_id++)
  {
    // get a reference to the sub system
    PMLinearImplicitSystem &sys = es.get_system<PMLinearImplicitSystem>
      (sys_names[s_id]);
    // get dt from this system
    const Real &dt = sys.get_dt();
    min_dt = (dt < min_dt) ? dt : min_dt;
//    ss << "dt(" << sys_names[s_id] << "): " << dt;
//    PMToolBox::output_message(ss, *comm_in);
  }
  // return minimum dt
  return min_dt;
}

// =========================================================================
void Copss::fixman_integrate(EquationSystems& equation_systems, unsigned int& i)
{
  // PerfLog perf_log("integration");
  PMSystemStokes& system = equation_systems.get_system<PMSystemStokes>("Stokes");
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the "disturbed" particle velocity + "undisturbed" velocity = U0
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (i > 0) {
    *(system.solution) = *(system.undisturbed_solution); // re-assign the undisturbed solution
    // Update the local values to reflect the solution on neighboring processors
    system.update();
    if (update_neighbor_list_everyStep) {
      neighbor_list_update_flag = true;
    }
    else if (timestep_duration > neighbor_list_update_interval) {
      ss << "======> update neighbor list at timestep " + std::to_string(i);
      PMToolBox::output_message(ss, *comm_in);
      neighbor_list_update_flag = true;
      timestep_duration         = 0;
    }
  }
  system.reinit_system(neighbor_list_update_flag, build_elem_neighbor_list, "disturbed");
  if (print_info) 
  {
    PMToolBox::output_message("--------After reinit system at integration step " + std::to_string(i) + "---------", *comm_in);
    if (comm_in->rank() == 0) point_mesh->print_point_info();
  }
  // compute undisturbed velocity of points
  system.compute_point_velocity("undisturbed", vel0);
  // compute disturbed velocity of points
  system.solve("disturbed"); // Using SolverStokes
  system.compute_point_velocity("disturbed", vel1);
  // add up undisturbed and disturbed velocity of points
  for (std::size_t j = 0; j < vel1.size(); ++j) vel1[j] += vel0[j];
  // transform total point velocity to U0 in Brownian_system
  brownian_sys->vector_transform(vel1, &U0, "forward");
  // assign vel1 to particle velocity
  point_mesh->set_bead_velocity(vel1);
  /*---------------------------------------------------------------------------------------
   * write equation system at step i
   * print out information at step 0
   * do not print out information at the first step when restart since it is
   * identical to the
   * last step before restart.
     -----------------------------------------------------------------------------------------*/
  if (i % write_interval == 0) 
  {
    ss << "Starting Fixman integration at step " + std::to_string(i);
    PMToolBox::output_message(ss, *comm_in);
    if (i != restart_step) 
    {
      if (std::find(output_file.begin(), output_file.end(), "equation_systems")
        != output_file.end())
      {
        system.write_equation_systems(o_step, real_time, "total");
      }
      /*
       * write particle to output file
       */
      this->write_object(i);
      /*----------------------------------------------------------------------------------------------------
      * Write out ROUT for restart mode at step i
         ----------------------------------------------------------------------------------------------------*/
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,
                            "vector_ROUT.dat",
                            FILE_MODE_WRITE,
                            &viewer);
      VecView(ROUT, viewer);
      PetscViewerDestroy(&viewer);
    }

    // update o_step
    o_step++;
  } // end if (i % write_interval == 0 )

  // get time step dt
  const Real& dt = this->get_min_dt(equation_systems);
//  ss << "min{dt} of all systems: " << dt;
//  PMToolBox::output_message(ss, *comm_in);
  equation_systems.parameters.set<Real>("dt") = dt;

  // use MidPoint scheme if Brownian is on
  if (with_brownian)
  {
    /*-------------------------------------------
     * For restarted Brownian system, to make sure the trajectory after is
     *exactly
     * the same as a continus simulation, we need to let the random number
     *generator
     * run "restart_step" steps
     * ------------------------------------------*/
    if (i == restart_step) {
      ss << "       ##########################################################\n"
         << "       # Run RandomGenerator for " << restart_step << " steps to make sure the trajectroy after restart is the same as a \n" 
         << "       # continuous simulation (this is only done once at the first step of restarted simulation)\n"
         << "       ###########################################################\n";
      PMToolBox::output_message(ss, *comm_in);
      Vec dw_prep; // this Vec is never used
      for (unsigned int j = 0; j < restart_step; j++) {
        brownian_sys->std_random_vector(0.0, 1.0, "gaussian", &dw_prep);
        VecDestroy(&dw_prep);
      }
    }

    // PerfLog perf_log("BD step");
    // perf_log.push("bd");

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Generate random vector dw whose mean = 0, variance = sqrt(2*dt)
       petsc_random_vector generates a uniform distribution [0 1] whose
       mean = 0.5 and variance = 1/12, so we need a shift and scale operation.
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    Real mean_dw = 0.0, variance_dw = 0.0;

    // A more precise way is to construct a random vector with gaussian
    // distribution
    const Real std_dev = std::sqrt(dt);
    brownian_sys->std_random_vector(0.0, std_dev, "gaussian", &dw);
    brownian_sys->_vector_mean_variance(dw, mean_dw, variance_dw);
    VecScale(dw, std::sqrt(2.0));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
       Print out the mean and variance or view the generated vector.
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    //      PetscPrintf(PETSC_COMM_WORLD,
    //                 "Generated random_vector: mean = %f, variance = %f\n",
    //                 mean_dw, variance_dw);
    //      PetscPrintf(PETSC_COMM_WORLD,
    //                 "Predicted random_vector: mean = %f, variance = %f\n",
    //                 0., std::sqrt(2.*dt));
    // Compute dw = B^-1 * dw using Chebyshev polynomial, dw will be changed!
    VecCopy(dw, dw_mid); // save dw to dw_mid, which will be used for Chebyshev

    // perf_log.push("cheb_converge");
    for (std::size_t j = 0; j < 2; j++)
    {
      //  ss << "j = "<<j<<"\n";

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        -
         Compute the max/min eigenvalues if needed. Otherwise, magnify the
           interval.
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           - */
      if (compute_eigen) {
        ss << "compute_eigen = " + std::to_string(compute_eigen);
        PMToolBox::output_message(ss, *comm_in);

        // ss << "Compute the max & min eigenvalues for Chebyshev polynomial
        // at step "<<i+1<<"\n";
        brownian_sys->compute_eigenvalues(eig_min, eig_max, tol_eigen);

        // Magnify the spectral range by a factor (1.05 by default).
        eig_max *= eig_factor; eig_min /= eig_factor;

        if (eig_min < 0 or eig_max<0 or eig_min>eig_max) {
          ss << "--->Invalid eigenvalue range: eig_min = " + std::to_string(eig_min)
             << "; eig_max = " + std::to_string(eig_max) + "\n";
          PMToolBox::output_message(ss, *comm_in);
          cheb_converge = false;
          break;
        }
        else {
          ss << "--->Valid eigenvalue range: recomputed eigen values and magnify the range by a factor eig_factor = "<< eig_factor
             << "; eig_min = " << eig_min
             << "; eig_max = " << eig_max
             << "; tol_cheb = " << tol_cheb
             << "; max_n_cheb = " << max_n_cheb << "\n"
             << "--->Write this eigenvalue range to out.eigenvalue for restart purpose";
          PMToolBox::output_message(ss, *comm_in);
          std::ofstream out_file;

          if (comm_in->rank() == 0) {
            if (i == 0) {
              out_file.open("out.eigenvalue", std::ios_base::out);
              out_file << "eigen_min" << "  " << "eigen_max" << "\n";
            }
            else {
              out_file.open("out.eigenvalue", std::ios_base::app);
            }
            out_file.precision(o_precision);
            out_file << eig_min << "  " << eig_max << "\n";
            out_file.close();
          }
        }
      }

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        -
         Compute the Brownian displacement B^-1 * dw using Chebyshev
           approximation.
         Here dw is both input and output variables, so it will be changed.
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           - */
      cheb_converge = brownian_sys->chebyshev_polynomial_approximation(max_n_cheb,
                                                                       eig_min,
                                                                       eig_max,
                                                                       tol_cheb,
                                                                       &dw);

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        -
         If converged, dw returns the Brownian displacement, then break the
           j-loop;
         Otherwise, recompute eigenvalues
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           - */
      if (cheb_converge) {
        compute_eigen = false; break;
      }
      else {
        compute_eigen = true;
        VecCopy(dw_mid, dw); /*copy back, recompute eigenvalues*/

        // ss << "It is necessry to re-compute the eigenvalues at step "
        // <<i+1<<"\n";
      }
    } // end for j-loop

    // perf_log.pop("cheb_converge");

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Double-check the convergence of Chebyshev polynomial approximation
       If cheb_converge = true continue fixman integration
       otherwise, relax the system for serval free-draining steps
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         */

    // perf_log.push("fixman step");
    if (cheb_converge) 
    {
      // Compute dw_mid = D*B^-1*dw, which can be obtained by solving the Stokes
      brownian_sys->hi_ewald(M, dw, dw_mid); // dw_mid = D * dw

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Particle coordinate vector R0.
         Move the particle R_mid = R0 + 0.5*(U0+U1)*dt (deterministic)
         and R_mid = R_mid + 0.5*sqrt(2)*D*B^-1*dw     (stochastic)
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            - */
      brownian_sys->extract_particle_vector(&R0, "coordinate", "extract");
      // R_mid = R0 + 0.5*dt*(U0+U1) (R0 and U0 do NOT change)
      VecWAXPY(R_mid, 0.5 * dt, U0, R0);                  
      // coefficient.                   
      coef = 0.5;                                                            
      // R_mid = R_mid + 0.5*sqrt(2)*D*B^-1*dw; 
      // sqrt(2) is introduced when generating dw
      VecAXPY(R_mid, coef, dw_mid);
      // Update mid-point coords && apply pbc to particle position                                          
      brownian_sys->extract_particle_vector(&R_mid, "coordinate", "assign"); 
      this->update_object();
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        -
         Update the particle mesh for the mid-point step,
         and recompute U0 + U1_mid, D_mid*(B^-1*dw)
         NOTE: the FEM solution of undisturbed field doesn't change, but
           particles
         move, so U0 needs to be re-evaluated at the new position.
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           - */
      *(system.solution) = *(system.undisturbed_solution); // re-assign the undisturbed solution
      system.update();

      // comment the line below if not update neighbor list at each time step
      if (update_neighbor_list_everyStep) neighbor_list_update_flag = true;
      system.reinit_system(neighbor_list_update_flag, build_elem_neighbor_list, "disturbed");
      system.compute_point_velocity("undisturbed", vel0);
      system.solve("disturbed"); // solve the disturbed solution
      system.compute_point_velocity("disturbed", vel1);

      for (std::size_t j = 0; j < vel1.size(); ++j) vel1[j] += vel0[j];
      brownian_sys->vector_transform(vel1, &U0, "forward"); 
      // dw_mid = D_mid*dw, where dw=B^-1*dw 
      brownian_sys->hi_ewald(M, dw, dw_mid);                
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        -
         the mid-point to the NEW point, and update the particle coordinates
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           - */

      // ss << "Update from mid-point to the NEW particle coordinates at step
      // " <<i+1<<"\n";
      // R_mid = R0 + dt*U0_mid
      VecWAXPY(R_mid, dt, U0, R0);                           
      // R_mid = R_mid + sqrt(2)*D_mid*B^-1*dw                
      VecAXPY(R_mid, 2.0 * coef, dw_mid);
      // Update mid-point coords && apply pbc to particle position                                    
      brownian_sys->extract_particle_vector(&R_mid, "coordinate", "assign"); 
      this->update_object();

      // Update ROUT (position vector excluding pbc) at the i-th step
      VecAXPY(ROUT, dt,         U0);     // ROUT = ROUT + dt*U0_mid
      VecAXPY(ROUT, 2.0 * coef, dw_mid); // ROUT = ROUT + sqrt(2)*D_mid*B^-1*dw
    }// end if cheb_converge
    // Use dR = Utotal*dt if Brownian is off
    else 
    {
      n_chebyshev_failure += 1;
      ss << "****** Warning: After recomputing eigenvalues, Chebysheve failed to converge at step "
         << std::to_string(i) << "; relax the system for " 
         << std::to_string(n_relax_step) + " steps";
      PMToolBox::output_message(ss, *comm_in);

      for (int relax_step = 0; relax_step < n_relax_step; relax_step++) {
        langevin_integrate(equation_systems, i);
        i++;
      }
    } // end else cheb_converge
     // perf_log.pop("fixman step");
     // perf_log.pop("bd");
  }   // end if Brownian
  else 
  { // if without Brownian
    // Move the particle R_mid = R0 + (U0+U1)*dt (deterministic)
    brownian_sys->extract_particle_vector(&R0, "coordinate", "extract");
    // R_mid = R0 + dt*Utotal (U0 is actually Utotal)
    VecWAXPY(R_mid, dt, U0, R0);
    // Update mid-point coords                                           
    brownian_sys->extract_particle_vector(&R_mid, "coordinate", "assign"); 
    this->update_object();
    // Update ROUT (position vector excluding pbc) at the i-th step
    VecAXPY(ROUT, dt, U0); // ROUT = ROUT + dt*U0_mid
  } // end else (without_brownian)
  real_time += dt;
  timestep_duration += 1;
}

void Copss::langevin_integrate(EquationSystems& equation_systems, unsigned int& i)
{
  // get stokes system from equation systems
  PMSystemStokes& system = equation_systems.get_system<PMSystemStokes>("Stokes");
  if (i > 0) {
    system.update();

    if (update_neighbor_list_everyStep) {
      neighbor_list_update_flag = true;
    }
    else if (timestep_duration > neighbor_list_update_interval) {
      neighbor_list_update_flag = true;
      timestep_duration         = 0;
    }

    // whether or not reinit neighbor list depends on the
    // neighbor_list_update_flag
  }
  system.reinit_system(neighbor_list_update_flag, build_elem_neighbor_list, "disturbed");
  if (print_info) 
  {
    PMToolBox::output_message("--------After reinit system at integration step " + std::to_string(i) + "---------", *comm_in);
    if (comm_in->rank() == 0) point_mesh->print_point_info();
  }

  Point p_velocity(0.);

  for (std::size_t p_id = 0; p_id < NP; p_id++) {
    for (int _dim = 0; _dim < dim; _dim++) {
      p_velocity(_dim) =
        point_mesh->particles()[p_id]->particle_force()(_dim);
      vel1[dim * p_id + _dim] = p_velocity(_dim);
    }
  }

  // set vel1 to particle velocity
  point_mesh->set_bead_velocity(vel1);

  // transform total point velocity to U0 in Brownian_system
  brownian_sys->vector_transform(vel1, &U0, "forward");

  /*---------------------------------------------------------------------------------------
   * write equation system at step i
   * print out information at step 0
   * do not print out information at the first step when restart since it is
   *identical to the
   * last step before restart.
     -----------------------------------------------------------------------------------------*/
  if (i % write_interval == 0) {
    ss <<"Starting Langevin integration at step " + std::to_string(i);
    PMToolBox::output_message(ss, *comm_in);
    if (i != restart_step) {
        if (std::find(output_file.begin(), output_file.end(), "equation_systems") 
          != output_file.end())
        {
          system.write_equation_systems(o_step, real_time, "total");
        } // end if (write es)
      // write particle to output file
      this->write_object(i);
      /*----------------------------------------------------------------------------------------------------
      * Write out ROUT for restart mode at step i
         ----------------------------------------------------------------------------------------------------*/
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,
                            "vector_ROUT.dat",
                            FILE_MODE_WRITE,
                            &viewer);
      VecView(ROUT, viewer);
      PetscViewerDestroy(&viewer);
    }

    // update o_step
    o_step++;
  } // end if (i % write_interval == 0 )

  // get time step dt
  const Real& dt = this->get_min_dt(equation_systems);
//  ss << "min{dt} of all systems: " << dt;
//  PMToolBox::output_message(ss, *comm_in);
  equation_systems.parameters.set<Real>("dt") = dt;


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     If with Brownian motion, we use midpoint scheme
     If without Brownian motion, we use normal stepping: dR = Utotal*dt
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       */
  if (with_brownian)
  {
    /*-------------------------------------------
     * For restarted Brownian system, to make sure the trajectory after is
     *exactly
     * the same as a continus simulation, we need to let the random number
     *generator
     * run "restart_step" steps
     * ------------------------------------------*/
    if (i == restart_step) {
      Vec dw_prep; // this Vec is never used

      for (unsigned int j = 0; j < restart_step; j++) {
        brownian_sys->std_random_vector(0.0, 1.0, "gaussian", &dw_prep);
        VecDestroy(&dw_prep);
      }
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Generate random vector dw whose mean = 0, variance = sqrt(2*dt)
       petsc_random_vector generates a uniform distribution [0 1] whose
       mean = 0.5 and variance = 1/12, so we need a shift and scale operation.
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         */
    Real mean_dw = 0.0, variance_dw = 0.0;

    // A more precise way is to construct a random vector with gaussian
    // distribution
    const Real std_dev = std::sqrt(dt);
    brownian_sys->std_random_vector(0.0, std_dev, "gaussian", &dw);
    brownian_sys->_vector_mean_variance(dw, mean_dw, variance_dw);
    VecScale(dw, std::sqrt(2.0));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Print out the mean and variance or view the generated vector.
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         */

    //    PetscPrintf(PETSC_COMM_WORLD,
    //               "Generated random_vector: mean = %f, variance = %f\n",
    //               mean_dw, variance_dw);
    //    PetscPrintf(PETSC_COMM_WORLD,
    //               "Predicted random_vector: mean = %f, variance = %f\n",
    //               0., std::sqrt(2.*dt));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Particle coordinate vector R0.
       Move the particle R_mid = R0 + U0*dt (deterministic)
       and R_mid = R_mid + dw, where dw = sqrt(2*dt)*R, where <R>=0, <R^2>=1
              (stochastic)
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          */
    brownian_sys->extract_particle_vector(&R0, "coordinate", "extract");
    // R_mid = R0 + dt * U0 (R0 and U0 do NOT change)
    VecWAXPY(R_mid, dt, U0, R0);                     
    // coefficient.                      
    coef = 1.;                
    // R_mid = R_mid + dw, where dw = sqrt(2*dt)*R, where <R>=0, <R^2>=1                                             
    VecAXPY(R_mid, coef, dw);                                    
    // Update mid-point coords && apply pbc to particle position          
    brownian_sys->extract_particle_vector(&R_mid, "coordinate", "assign"); 
    this->update_object();

    // Update ROUT (position vector excluding pbc) at the i-th step
    // ROUT = ROUT + dt*U0_mid
    VecAXPY(ROUT, dt,   U0);                                
    // ROUT = ROUT + dw, where dw = sqrt(2*dt)*R, where <R>=0, <R^2>=1              
    VecAXPY(ROUT, coef, dw);                                               
  } // end if Brownian
  else {                                                                   
    // Move the particle R_mid = R0 + (U0+U1)*dt (deterministic)
    brownian_sys->extract_particle_vector(&R0, "coordinate", "extract");
    // R_mid = R0 + dt*Utotal (U0 is actually Utotal)
    VecWAXPY(R_mid, dt, U0, R0);                                
    // Update mid-point coords           
    brownian_sys->extract_particle_vector(&R_mid, "coordinate", "assign"); 
    this->update_object();
    // Update ROUT (position vector excluding pbc) at the i-th step
    VecAXPY(ROUT, dt, U0); // ROUT = ROUT + dt*U0_mid
  } // end else (without_brownian)
  real_time         += dt;
  timestep_duration += 1;
}

void Copss::destroy()
{
  VecDestroy(&U0);
  VecDestroy(&R0);
  VecDestroy(&R_mid);
  VecDestroy(&dw_mid);
  PetscRandomDestroy(&rand_ctx);
  if (with_brownian) {
    VecDestroy(&dw);
  }
  if (&viewer) {
    PetscViewerDestroy(&viewer);
  }
  PetscViewerDestroy(&viewer);
}
} // end namespace
