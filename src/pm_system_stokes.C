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


// std C++
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <time.h>

// libmesh head files
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmv_io.h"
#include "libmesh/vtk_io.h"

// local header files
#include "pm_toolbox.h"
#include "ggem_stokes.h"
#include "brownian_system.h"
#include "pm_system_stokes.h"
#include "pm_system_poisson.h"

namespace libMesh {
// ======================================================================================
PMSystemStokes::PMSystemStokes(EquationSystems  & es,
                               const std::string& name,
                               const unsigned int number)
  : PMLinearImplicitSystem(es, name, number),
  _solver_stokes(es)
{
  if (name != "Stokes") libmesh_error();
  assemble_stokes     = (new AssembleStokes(es, name));
  analytical_solution = assemble_stokes->get_analytical_solution();
  ggem_stokes         = assemble_stokes->get_ggem_stokes();
}

// ==================================================================================
PMSystemStokes::~PMSystemStokes()
{
  // Clear data
  this->clear();
}

// ==================================================================================
void PMSystemStokes::clear()
{
  // delete the pointer
  if (assemble_stokes) {
    delete assemble_stokes;
  }

  // analytical_solution and ggem_stokes pointer will be destroyed in
  // AssembleStokes
}

// ==================================================================================
void PMSystemStokes::reinit_system(bool      & neighbor_list_update_flag,
                                   const bool& build_elem_neighbor_list)
{
  START_LOG("reinit_system()", "PMSystemStokes");
  this->comm().barrier(); // Is this at the beginning or the end necessary?
  // reinit point-mesh system, including
  // (1) build the point-point neighbor list according to search radius;
  // (2) build the element-point neighbor list according to search radius;
  // (3) evaluate forces
  // perf_log.push("reinit point_mesh");
  _point_mesh->reinit(neighbor_list_update_flag, build_elem_neighbor_list);

  // perf_log.pop("reinit point_mesh");
  // perf_log.push("fix compute");
  // update the tracking points position on the mesh if needed
  if (_particle_mesh != NULL) {
    // update node positions on particle mesh using updated point mesh
    // information
    _point_mesh->update_particle_mesh(_particle_mesh);

    // zero force density of each rigid particle
    _particle_mesh->zero_node_force();

    // need to check if particles are on the pbc, if so, rebuild the particle
    // mesh
    _fixes[0]->check_pbc_pre_fix();

    // compute forces on surface nodes
    for (std::size_t i = 0; i < _fixes.size(); i++) {
      // _fixes[i]->print_fix();
      _fixes[i]->compute();
    }

    // sync forces on nodes to PointParticles in point mesh
    _fixes[0]->sync_node_to_pointmesh();

    // need to restore particle mesh after applying all fixes if particles are
    // on pbc
    _fixes[0]->check_pbc_post_fix();
  }
  else
  {
    for (std::size_t i = 0; i < _fixes.size(); i++) {
      _fixes[i]->compute();
    }
  }
  // multi-Physics coupling
  if (this->get_equation_systems().parameters.get<bool>("module_poisson")) {
    this->couple_poisson(false);  
  }

  // perf_log.pop("fix compute");
  STOP_LOG("reinit_system()", "PMSystemStokes");
}

// ===========================================================
void PMSystemStokes::assemble_matrix(const std::string& system_name,
                                     const std::string& option)
{
  libmesh_assert(this->matrix);
  libmesh_assert(this->matrix->initialized());

  START_LOG("assemble_matrix()", "PMSystemStokes");

  // init the matrices: global stiffness and PC matrix (if required)
  this->matrix->zero();
  const bool user_defined_pc = this->get_equation_systems().parameters.get<bool>(
    "user_defined_pc");

  if (user_defined_pc) {
    std::cout << "--->test in PMSystemStokes::assemble_matrix(): "
              << "Initialize the preconditioning matrix (all zeros) \n";
    this->get_matrix("Preconditioner").zero();
  }

  // Call assemble function to assemble matrix
  // assemble_matrix_sedimentation_ex1(this->get_equation_systems(),
  // system_name, option);
  assemble_stokes->assemble_global_K(system_name, option);

  // close the matrices
  this->matrix->close(); // close the matrix

  if (user_defined_pc) this->get_matrix("Preconditioner").close();

  STOP_LOG("assemble_matrix()", "PMSystemStokes");
}

// ==================================================================================
void PMSystemStokes::assemble_rhs(const std::string& system_name,
                                  const std::string& option)
{
  libmesh_assert(this->rhs);
  libmesh_assert(this->rhs->initialized());

  START_LOG("assemble_rhs()", "PMSystemStokes");

  // PerfLog perf_log("assemble_rhs");
  // init, assemble, and close the rhs vector
  // perf_log.push("zero");
  this->rhs->zero();
  // perf_log.pop("zero");
  // assemble_rhs_sedimentation_ex1 (this->get_equation_systems(), system_name,
  // option);
  // perf_log.push("assemble");
  assemble_stokes->assemble_global_F(system_name, option);

  // perf_log.pop("assemble");
  // perf_log.push("close");
  this->rhs->close();

  // perf_log.pop("close");

  STOP_LOG("assemble_rhs()", "PMSystemStokes");
}

// ==================================================================================
void PMSystemStokes::solve(const std::string& option)
{
  START_LOG("solve()", "PMSystemStokes");

  // PerfLog perf_log("solve_stokes");
  // Real t1, t2;
  // std::string msg = "---> solve Stokes";
  // PMToolBox::output_message(msg, this->comm());

  // Assemble the global matrix and pc matrix, and record the CPU wall time.
  if (_re_init)
  {
    // perf_log.push("assemble_matrix (undisturbed)");
    // t1 = MPI_Wtime();

    // set the solver type for the Stokes equation
    const SystemSolverType solver_type =
      this->get_equation_systems().parameters.get<SystemSolverType>(
        "solver_type_stokes");
    _solver_stokes.set_solver_type(solver_type);

    // Assemble the global matrix, and init the KSP solver
    this->assemble_matrix("Stokes", option);
    _solver_stokes.init_ksp_solver();
    
    // set re_init to false once K matrix is built
    _re_init = false;

    // perf_log.pop("assemble_matrix (undisturbed)");

    // t2 = MPI_Wtime();
    // std::cout << "Time used to assemble the global matrix and reinit KSP is "
    // <<t2-t1<<" s\n\n";
  }

  // assemble the rhs vector, and record the CPU wall time.
  // t1 = MPI_Wtime();
  // perf_log.push("assemble_rhs");
  this->assemble_rhs("Stokes", option);

  // perf_log.pop("assemble_rhs");

  // t2 = MPI_Wtime();
  // std::cout << "Time used to assemble the right-hand-side vector is "
  // <<t2-t1<<" s\n";

  // solve the problem
  // perf_log.push("solve()");
  _solver_stokes.solve();

  // perf_log.pop("solve()");

  STOP_LOG("solve()", "PMSystemStokes");
}

// ==================================================================================
void PMSystemStokes::add_local_solution()
{
  START_LOG("add_local_solution()", "PMSystemStokes");

  // Check if the system solution vector is closed or not
  if ((this->solution->closed()) == false) this->solution->close();

  // this->update();

  // Get the parameters and Initialize the quantities
  const bool test_output            = false;
  MeshBase & mesh                   = this->get_mesh();
  const std::size_t   dim           = mesh.mesh_dimension();
  const std::size_t   n_local_nodes = mesh.n_local_nodes();
  std::vector<Number> local_solution(dim * n_local_nodes);
  std::vector<numeric_index_type> dof_indices(dim * n_local_nodes);

  // printf("--->test in add_local_solution() n_local_nodes = %lu on the
  // processor %u\n",
  //       n_local_nodes,this->comm().rank());

  // Update the system solution by adding the local solution (from Green's
  // function)
  MeshBase::node_iterator nd           = mesh.local_nodes_begin();
  const MeshBase::node_iterator end_nd = mesh.local_nodes_end();
  std::size_t local_count              = 0;

  for (; nd != end_nd; ++nd)
  {
    // Store a pointer to the current node, and extract a point
    Node *node = *nd;
    Point pt;

    for (unsigned int i = 0; i < dim; ++i) pt(i) =  (*node)(i);

    // this is a test for dof_number at each node
    if (test_output)
    {
      const unsigned int node_id = node->id();
      std::ostringstream oss;
      oss << "          NODE " << node_id;
      PMToolBox::output_message(oss, this->comm());
      node->print_info();

      if (this->comm().rank() == 0) printf("--->test: nodal dof number :");

      for (unsigned int i = 0; i < dim; ++i)
      {
        dof_id_type dof_num = node->dof_number(this->number(), i, 0);

        if (this->comm().rank() == 0) printf(" %u", dof_num);
      }

      if (this->comm().rank() == 0) printf(" \n");
    }

    // get the dof numbers at this node (only for velocity)
    std::vector<dof_id_type> dof_nums(dim);

    for (unsigned int i = 0; i < dim; ++i) { // var = 0, 1, 2 = i
      dof_nums[i] = node->dof_number(this->number(), i, 0);
    }

    // compute the local velocity of fluid at the current node
    const std::vector<Real> Ulocal =
      this->local_velocity_fluid(pt, "regularized");

    // store the local velocity and dof indices
    for (unsigned int i = 0; i < dim; ++i)
    {
      local_solution[local_count * dim + i] =  Ulocal[i];
      dof_indices[local_count * dim + i]    =  dof_nums[i];
    }
    local_count++;
  } // end for

  // printf("--->test in add_local_solution() local_count = %lu on the processor
  // %u\n",
  //       n_local_nodes,this->comm().rank());

  // add the local to the global
  // this->solution->zero();
  this->solution->add_vector(local_solution, dof_indices);
  this->solution->close();
  this->update();
  
  // multi-Physics coupling
  if (this->get_equation_systems().parameters.get<bool>("module_poisson")) {
    this->couple_poisson(true);  
  }
  
  STOP_LOG("add_local_solution()", "PMSystemStokes");
}

// ==================================================================================
void PMSystemStokes::test_l2_norm(bool& neighbor_list_update_flag)
{
  START_LOG("test_l2_norm()", "PMSystemStokes");
  std::ostringstream ss;
  ss << "--->test in PMSystemStokes::test_l2_norm(): \n";
  PMToolBox::output_message(ss, this->comm());
  bool build_elem_neighbor_list = true;

  // Numerical solution: Global(FEM) + Local(Analytical)
  this->reinit_system(neighbor_list_update_flag, build_elem_neighbor_list);
  _re_init = true;
  this->solve("disturbed");
  this->add_local_solution();

  // Theoretical solution(only velocity)
  const unsigned int dim = 3;
  MeshBase& mesh = this->get_mesh();
  const unsigned int n_nodes = mesh.n_nodes();
  Real val0_norm = 0., val1_norm = 0.;
  Real val2_norm = 0., val3_norm = 0.;

  // AnalyticalSolution analytical_solution(*this);

  // Loop over each node and compute the nodal velocity value
  MeshBase::node_iterator nd           = mesh.local_nodes_begin();
  const MeshBase::node_iterator end_nd = mesh.local_nodes_end();

  for (; nd != end_nd; ++nd)
  {
    // Store a pointer to the current node, and extract a point
    Node *node = *nd;
    Point pt;

    for (unsigned int i = 0; i < dim; ++i) pt(i) =  (*node)(i);

    // get the dof numbers at this node (only for velocity)
    std::vector<dof_id_type> dof_nums(dim);

    for (unsigned int i = 0; i < dim; ++i) { // var = 0, 1, 2 = i
      dof_nums[i] = node->dof_number(this->number(), i, 0);
    }

    // Get the numerical solution
    std::vector<Real> Unum;
    this->solution->get(dof_nums, Unum);

    // compute the local velocity of fluid at the current node
    const std::vector<Real> Uexact =
      analytical_solution->exact_solution_infinite_domain(*ggem_stokes, pt);

    // compute the errors
    for (unsigned int i = 0; i < dim; ++i) {
      Real tmpt = std::abs(Unum[i] - Uexact[i]);
      val0_norm += tmpt;
      val1_norm += std::abs(Uexact[i]);

      val2_norm += tmpt * tmpt;
      val3_norm += Uexact[i] * Uexact[i];
    }
  } // end for nd-loop

  // Compute the error: l1 and l2-norm of errors
  this->comm().sum(val0_norm);
  this->comm().sum(val1_norm);
  this->comm().sum(val2_norm);
  this->comm().sum(val3_norm);

  const Real l1_norm = val0_norm / val1_norm;
  const Real l2_norm = std::sqrt(val2_norm) / std::sqrt(val3_norm);
  ss << "--->test in test_l1_norm: l1_norm = " << l1_norm
     << "; l2_norm = " << l2_norm << "\n";
  PMToolBox::output_message(ss, this->comm());
  STOP_LOG("test_l2_norm()", "PMSystemStokes");
}

// ==================================================================================
void PMSystemStokes::write_equation_systems(const std::size_t  time_step,
                                            const std::string& output_filename,
                                            const std::string& output_format)
{
  START_LOG("write_equation_systems()", "PMSystemStokes");

  // Write out the FEM results: global solution
  MeshBase & mesh           = this->get_mesh();
  const bool fem_sol_output = true;

  if (fem_sol_output)
  {
    std::ostringstream file_name_fem;
    file_name_fem << output_filename + "_fem";

    if (output_format == "EXODUS")
    {
#ifdef LIBMESH_HAVE_EXODUS_API

      if (time_step == 0)
      {
        file_name_fem << ".e";
        ExodusII_IO(mesh).write_equation_systems(file_name_fem.str(),
                                                 this->get_equation_systems());
      }
      else
      {
        // file_name_fem << ".e-s." << std::setw(8) << std::setfill('0') <<
        // std::right << time_step;
        // ExodusII_IO(mesh).write_equation_systems(file_name_fem.str(),this->get_equation_systems());

        file_name_fem << ".e";
        ExodusII_IO exodus_IO(mesh);
        exodus_IO.append(true);
        exodus_IO.write_timestep(file_name_fem.str(),
                                 this->get_equation_systems(),
                                 time_step + 1,
                                 time_step + 1);
      } // end if-else
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
    }
    else if (output_format == "VTK")
    {
#ifdef LIBMESH_HAVE_VTK
      file_name_fem << "_" << std::setw(8) << std::setfill('0') << std::right <<
        time_step << ".vtu";
      VTKIO(mesh).write_equation_systems(file_name_fem.str(),
                                         this->get_equation_systems());
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
    }
    else
    {
      file_name_fem << "_" << std::setw(8) << std::setfill('0') << std::right <<
        time_step << ".gmv.";
      GMVIO(mesh).write_equation_systems(file_name_fem.str(),
                                         this->get_equation_systems());
    }
  }

  // Update the system solution by adding the local solution (from Green's
  // function)
  this->add_local_solution();

  // Write out the FEM results: global solution
  std::ostringstream file_name;
  file_name << output_filename + "_total";

  if (output_format == "EXODUS")
  {
#ifdef LIBMESH_HAVE_EXODUS_API

    if (time_step == 0)
    {
      file_name << ".e";
      ExodusII_IO(mesh).write_equation_systems(file_name.str(),
                                               this->get_equation_systems());
    }
    else
    {
      // file_name << ".e-s." << std::setw(8) << std::setfill('0') << std::right
      // << time_step;
      // ExodusII_IO(mesh).write_equation_systems(file_name.str(),this->get_equation_systems());

      file_name << ".e";
      ExodusII_IO exodus_IO(mesh);
      exodus_IO.append(true);
      exodus_IO.write_timestep(file_name.str(), this->get_equation_systems(),
                               time_step + 1, time_step + 1);
    } // end if-else
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
  }
  else if (output_format == "VTK")
  {
#ifdef LIBMESH_HAVE_VTK
    file_name << "_" << std::setw(8) << std::setfill('0') << std::right <<
      time_step << ".vtu";
    VTKIO(mesh).write_equation_systems(file_name.str(),
                                       this->get_equation_systems());
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
  }
  else
  {
    file_name << "_" << std::setw(8) << std::setfill('0') << std::right <<
      time_step << ".gmv";
    GMVIO(mesh).write_equation_systems(file_name.str(),
                                       this->get_equation_systems());
  }

  STOP_LOG("write_equation_systems()", "PMSystemStokes");
}

// ==================================================================================
void PMSystemStokes::compute_point_velocity(const std::string& option,
                                            std::vector<Real>& pv)
{
  START_LOG("compute_point_velocity()", "PMSystemStokes");

  // dim*NP: size of the velocity vector
  const MeshBase  & mesh   = this->get_mesh();
  const std::size_t NP     =  _point_mesh->num_particles();
  const std::size_t dim    = mesh.mesh_dimension();
  const dof_id_type n_elem = mesh.n_elem();

  // std::vector<Real> pvlocal(dim*NP,0.);  // declared on each processor
  std::vector<Real> _pv_send_list;                            // point velocity
                                                              // send list
  std::vector<Real> _pglobal_send_list;
  std::vector<Real> _plocal_send_list;
  std::vector<std::size_t> _pid_send_list;                    // point id send
                                                              // list
  const unsigned int u_var      = this->variable_number("u"); // u_var = 0
  const unsigned int v_var      = this->variable_number("v"); // v_var = 1
  const unsigned int w_var      = this->variable_number("w"); // w_var = 2
  const std::string  force_type = "regularized";

  // Loop over each point, and compute its velocity.
  // Collect the global velocity from FEM through Allgather operation
  for (std::size_t i = 0; i < NP; ++i)
  {
    // 0. point coordinates & residing element
    const Point pt               = _point_mesh->particles()[i]->point();
    const Point fv               = _point_mesh->particles()[i]->particle_force();
    const PointType   point_type = _point_mesh->particles()[i]->point_type();
    const dof_id_type elem_id    = _point_mesh->particles()[i]->elem_id();

    if (elem_id > n_elem)
    {
      std::cout << "---> Error in PMSystemStokes::compute_point_velocity():\n"
                << "      " << i << "-th particle (position = " << pt(0) <<
        ", " << pt(1) << ", " << pt(2) << ") is out of domain"
                << std::endl << std::endl;
    }
    const Elem *elem = mesh.elem(elem_id);

    // 1. Global(FEM) solution at the current point. This is done on local
    // processors
    Point Uglobal;

    // 2. Local velocity of particle i, This is also done on local processors
    Point Ulocal;

    // 3. total velocity
    Point U;

    if (elem->processor_id() == this->processor_id())
    {
      // get global_velocity
      Uglobal(0) = this->point_value(u_var, pt, *elem);
      Uglobal(1) = this->point_value(v_var, pt, *elem);
      Uglobal(2) = this->point_value(w_var, pt, *elem);

      // get local velocity
      if (option == "disturbed") {
        Ulocal = this->local_velocity_bead(i, force_type);

        // exclude self exclusioon term for point_type = "POLYMER_BEAD"
        if (point_type == POLYMER_BEAD) {
          const Point Uex = this->global_self_exclusion(i);
          U = Uglobal + Ulocal - Uex + fv; // drag0 by default = 1.
        }
        else if (point_type == LAGRANGIAN_POINT) U = Uglobal + Ulocal;
      }
      else if (option == "undisturbed") {
        U = Uglobal;
      }
      else {
        libmesh_error();
      }

      // pack particle id and its velocity
      _pid_send_list.push_back(i);

      for (int j = 0; j < dim; j++) {
        _pv_send_list.push_back(U(j));

        // _pglobal_send_list.push_back(Uglobal(j));
        // _plocal_send_list.push_back(Ulocal(j));
      }

      // printf("i = %i, processor_id = %i, Ulocal = %f, %f, %f\n", i,
      // this->processor_id(), Ulocal(0), Ulocal(1), Ulocal(2));
    } // end if (elem->processor_id() == this->processor_id)
  }   // end for i-loop

  // Check the size of local_pv and the size of list on each process after
  // allgather
  this->comm().allgather(_pid_send_list); // allgather the particle id
  this->comm().allgather(_pv_send_list); // allgather the particle velocity

  // this->comm().allgather(_pglobal_send_list);
  // this->comm().allgather(_plocal_send_list);
  if (_pid_send_list.size() != NP)
  {
    std::cout << "---> Error in PMSystemStokes::compute_point_velocity: \n"
              << "       _pid_send_list.size() != NP"
              << std::endl << std::endl;
    libmesh_error();
  }

  for (std::size_t i = 0; i < NP; ++i) {
    const std::size_t p_id = _pid_send_list[i];

    for (int j = 0; j < dim; ++j) {
      pv[dim * p_id + j] = _pv_send_list[i * dim + j];
    }

    // // ---------------------------- output for test
    // -----------------------------
    // if (this->comm().rank()==0 && option=="disturbed")
    // {
    //   printf("\n--->test in compute_point_velocity(): output point
    // velocity:\n");
    //   printf("point %lu: Uglobal (FEM) = (%E, %E, %E)\n",
    //          p_id, _pglobal_send_list[dim*i], _pglobal_send_list[dim*i+1],
    // _pglobal_send_list[dim*i+2] );
    //   printf("              Ulocal (Green Function) = (%E, %E, %E)\n",
    //  
    //   
    //   
    //  _plocal_send_list[dim*i],_plocal_send_list[dim*i+1],_plocal_send_list[dim*i+2]
    // );
    // //  printf("              Uexc    = (%E, %E, %E)\n", Uex[0], Uex[1],
    // Uex[2] );
    // //  printf("       Stokes drag    = (%E, %E, %E)\n", fv[0], fv[1], fv[2]
    // );
    // //  printf("           ---Utotal  = (%E, %E, %E)\n\n",
    //   //       pv[dim*pid], pv[dim*pid+1], pv[dim*pid+2] );
    // }
  }

  // std::cout <<"pv = ";
  // for(std::size_t i=0; i<pv.size(); i++){
  //   std::cout << pv[i] <<";";
  // }
  // std::cout<<std::endl;

  STOP_LOG("compute_point_velocity()", "PMSystemStokes");
}

// ==================================================================================
std::vector<Real>PMSystemStokes::compute_unperturbed_point_velocity()
{
  START_LOG("compute_unperturbed_point_velocity()", "PMSystemStokes");

  const std::size_t NP  =  _point_mesh->num_particles();
  const std::size_t dim =  this->get_mesh().mesh_dimension();

  // first solve the undisturbed flow field.
  _re_init = true;
  this->solve("undisturbed");
  std::vector<Real> pv_unperturbed(NP * dim, 0);

  // only FEM solution without particles! (undistrubed solution)
  this->compute_point_velocity("undisturbed", pv_unperturbed);

  STOP_LOG("compute_unperturbed_point_velocity()", "PMSystemStokes");
  return pv_unperturbed;
}

// ==================================================================================
std::vector<Real>PMSystemStokes::point_velocity(
  const std::vector<Real>& vel_beads,
  const std::size_t        i) const
{
  START_LOG("point_velocity()", "PMSystemStokes");

  const std::size_t dim =  this->get_mesh().mesh_dimension();
  std::vector<Real> point_v(dim);

  for (std::size_t k = 0; k < dim; ++k) {
    point_v[k] = vel_beads[i * dim + k];
  }

  STOP_LOG("point_velocity()", "PMSystemStokes");
  return point_v;
}

// ==================================================================================
std::vector<Number>PMSystemStokes::local_velocity_fluid(const Point      & p,
                                                        const std::string& force_type)
const
{
  START_LOG("local_velocity_fluid()", "PMSystemStokes");

  std::vector<Real> Ulocal = ggem_stokes->local_velocity_fluid(_point_mesh,
                                                               p,
                                                               force_type);

  STOP_LOG("local_velocity_fluid()", "PMSystemStokes");

  return Ulocal;
}

// ==================================================================================
std::vector<Number>PMSystemStokes::local_velocity_fluid(const Elem        *elem,
                                                        const Point      & p,
                                                        const std::string& force_type)
const
{
  START_LOG("local_velocity_fluid()", "PMSystemStokes");

  std::vector<Real> Ulocal = ggem_stokes->local_velocity_fluid(_point_mesh,
                                                               elem,
                                                               p,
                                                               force_type);

  STOP_LOG("local_velocity_fluid()", "PMSystemStokes");

  return Ulocal;
}

// ==================================================================================
Point PMSystemStokes::local_velocity_bead(const std::size_t& bead_id,
                                          const std::string& force_type) const
{
  START_LOG("local_velocity_bead()", "PMSystemStokes");

  Point Ulocal =
    ggem_stokes->local_velocity_bead(_point_mesh, bead_id, force_type);

  STOP_LOG("local_velocity_bead()", "PMSystemStokes");
  return Ulocal;
}

// ==================================================================================
Point PMSystemStokes::global_self_exclusion(const std::size_t p_id) const
{
  START_LOG("global_self_exclusion()", "PMSystemStokes");

  Point self_v = ggem_stokes->global_self_exclusion(_point_mesh, p_id);

  STOP_LOG("global_self_exclusion()", "PMSystemStokes");
  return self_v;
}

// ==================================================================================
void PMSystemStokes::test_velocity_profile()
{
  START_LOG("test_velocity_profile()", "PMSystemStokes");
  bool neighbor_list_update_flag = true;
  std::ostringstream ss;
  // solve the disturbed velocity field: global solution(FEM solution)
  // In this test, we assume that undisturbed velocity is zero!
  ss <<"========>2. Test in PMSystemStokes:: test_velocity_profile()";
  PMToolBox::output_message(ss, this->comm());
  // build both particle-particle and particle-elem neighbor list
  _re_init = true;
  this->solve("disturbed");
  // output the velocity profiles along xyz directions, global + local
  // solutions.
  const Point& box_min = _point_mesh->pm_periodic_boundary()->box_min();
  const Point& box_len = _point_mesh->pm_periodic_boundary()->box_length();
  const Real   xn = 200, yn = 80, zn = 80;
  const Real   dx = box_len(0) / xn, dy = box_len(1) / yn, dz = box_len(2) / zn;
  const unsigned int NP = _point_mesh->num_particles();

  // AnalyticalSolution analytical_solution(*this);

  std::ofstream outfile;
  int o_width = 12, o_precision = 9;

  // 1. write out the velocity profile along x-direction.
  std::ostringstream filenamex;
  filenamex << "output_velocity_profile_x_" << NP << "P.txt";
  outfile.open(filenamex.str(), std::ios_base::out);

  for (std::size_t i = 0; i < xn + 1; ++i)
  {
    Point pt(box_min(0) + Real(i) * dx, 0., 0.);

    // global velocity from FEM
    std::vector<Real>  Uglobal(3), Utotal(3);
    const unsigned int u_var = this->variable_number("u"); // u_var = 0
    const unsigned int v_var = this->variable_number("v"); // v_var = 1
    const unsigned int w_var = this->variable_number("w"); // w_var = 2
    Uglobal[0] = this->point_value(u_var, pt);             // this is slow, but
                                                           // only for test.
    Uglobal[1] = this->point_value(v_var, pt);
    Uglobal[2] = this->point_value(w_var, pt);

    // local velocity from analytical function
    const std::vector<Real> Ulocal =
      this->local_velocity_fluid(pt, "regularized");

    for (std::size_t j = 0; j < 3; ++j) Utotal[j] = Uglobal[j] + Ulocal[j];

    // Exact solution for an unbounded domain
    const std::vector<Real> Uexact =
      analytical_solution->exact_solution_infinite_domain(*ggem_stokes, pt);

    // write the velocity, x vx vy vz
    outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
    outfile.precision(o_precision);   outfile.width(o_width);

    if (this->comm().rank() == 0)
      outfile << pt(0) << "  " << Utotal[0] << "  " << Utotal[1] << "  " <<
        Utotal[2]
              << "  " << Uglobal[0] << "  " << Uglobal[1] << "  " << Uglobal[2]
              << "  " << Ulocal[0] << "  " << Ulocal[1] << "  " << Ulocal[2]
              << "  " << Uexact[0] << "  " << Uexact[1] << "  " << Uexact[2] <<
        "\n";
  } // end for i-loop
  outfile.close();

  // 2. write out the velocity profile along y-direction.
  std::ostringstream filenamey;
  filenamey << "output_velocity_profile_y_" << NP << "P.txt";
  outfile.open(filenamey.str(), std::ios_base::out);

  for (std::size_t i = 0; i < yn + 1; ++i)
  {
    Point pt(0., box_min(1) + Real(i) * dy, 0.);

    // global velocity from FEM
    std::vector<Real>  Uglobal(3), Utotal(3);
    const unsigned int u_var = this->variable_number("u"); // u_var = 0
    const unsigned int v_var = this->variable_number("v"); // v_var = 1
    const unsigned int w_var = this->variable_number("w"); // w_var = 2
    Uglobal[0] = this->point_value(u_var, pt);             // this is slow, but
                                                           // only for test.
    Uglobal[1] = this->point_value(v_var, pt);
    Uglobal[2] = this->point_value(w_var, pt);

    // local velocity from analytical function
    const std::vector<Real> Ulocal =
      this->local_velocity_fluid(pt, "regularized");

    for (std::size_t j = 0; j < 3; ++j) Utotal[j] = Uglobal[j] + Ulocal[j];

    // Exact solution for an unbounded domain
    const std::vector<Real> Uexact =
      analytical_solution->exact_solution_infinite_domain(*ggem_stokes, pt);

    // write the velocity, y vx vy vz
    outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
    outfile.precision(o_precision);   outfile.width(o_width);

    if (this->comm().rank() == 0)
      outfile << pt(1) << "  " << Utotal[0] << "  " << Utotal[1] << "  " <<
        Utotal[2]
              << "  " << Uglobal[0] << "  " << Uglobal[1] << "  " << Uglobal[2]
              << "  " << Ulocal[0] << "  " << Ulocal[1] << "  " << Ulocal[2]
              << "  " << Uexact[0] << "  " << Uexact[1] << "  " << Uexact[2] <<
        "\n";
  } // end for i-loop
  outfile.close();

  // 3. write out the velocity profile along z-direction.
  std::ostringstream filenamez;
  filenamez << "output_velocity_profile_z_" << NP << "P.txt";
  outfile.open(filenamez.str(), std::ios_base::out);

  for (std::size_t i = 0; i < zn + 1; ++i)
  {
    Point pt(0., 0., box_min(2) + Real(i) * dz);

    // global velocity from FEM
    std::vector<Real>  Uglobal(3), Utotal(3);
    const unsigned int u_var = this->variable_number("u"); // u_var = 0
    const unsigned int v_var = this->variable_number("v"); // v_var = 1
    const unsigned int w_var = this->variable_number("w"); // w_var = 2
    Uglobal[0] = this->point_value(u_var, pt);             // this is slow, but
                                                           // only for test.
    Uglobal[1] = this->point_value(v_var, pt);
    Uglobal[2] = this->point_value(w_var, pt);

    // local velocity from analytical function
    const std::vector<Real> Ulocal =
      this->local_velocity_fluid(pt, "regularized");

    for (std::size_t j = 0; j < 3; ++j) Utotal[j] = Uglobal[j] + Ulocal[j];

    // Exact solution for an unbounded domain
    const std::vector<Real> Uexact =
      analytical_solution->exact_solution_infinite_domain(*ggem_stokes, pt);

    // write the velocity, z vx vy vz
    outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
    outfile.precision(o_precision);   outfile.width(o_width);

    if (this->comm().rank() == 0)
      outfile << pt(2) << "  " << Utotal[0] << "  " << Utotal[1] << "  " <<
        Utotal[2]
              << "  " << Uglobal[0] << "  " << Uglobal[1] << "  " << Uglobal[2]
              << "  " << Ulocal[0] << "  " << Ulocal[1] << "  " <<  Ulocal[2]
              << "  " << Uexact[0] << "  " << Uexact[1] << "  " << Uexact[2] <<
        "\n";
  } // end for i-loop
  outfile.close();
  // done and write out the results
  std::ostringstream output_filename;
  output_filename << "output_velocity_profile_" << NP << "P.e";
#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO(this->get_mesh()).write_equation_systems(output_filename.str(),
                                                       this->get_equation_systems());
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
  STOP_LOG("test_velocity_profile()", "PMSystemStokes");
}

// ==================================================================================
// const std::vector<Real> PMSystemStokes::exact_solution(const Point& pt0)
// const
// {
//  START_LOG("exact_solution()", "PMSystemStokes");
//
//  // ksi's value should be consistent with that in
// GGEMStokes::regularization_parameter()
//  const Real ksi = std::sqrt(libMesh::pi)/3.0;  // = 0.591
//  const Real muc = 1.0/(6*libMesh::pi);
//  const unsigned int dim = 3;
//  std::vector<Real> UA(dim,0.);
//  DenseMatrix<Number> GT;
//
//  // GGEM object and number of points in the system
//  GGEMStokes ggem_stokes;
//  const std::size_t n_points = _point_mesh->num_particles();
//
//  // loop over each point
//  for(std::size_t i=0; i<n_points; ++i)
//  {
//    const Point pti = _point_mesh->particles()[i]->point();
//    const Point x   = pt0 - pti;
//
//    bool  zero_limit  = false;
//    if(x.size()<1E-6) zero_limit  = true;
//
//    // use ksi instead of alpha
//    GT = ggem_stokes.green_tensor_exp(x,ksi,muc,dim,zero_limit);
//    const std::vector<Real> fv =
// _point_mesh->particles()[i]->particle_force();
//    //printf("--->test in exact_solution(): i = %lu, fv = (%f,%f,%f)\n",
// i,fv[0],fv[1],fv[2]);
//
//    // 3. compute u due to this particle
//    for (std::size_t k=0; k<dim; ++k){
//      for (std::size_t l=0; l<dim; ++l){
//        UA[k] += GT(k,l)*fv[l];
//      } // end for l
//    } // end for k
//  } // end for i
//
//  STOP_LOG("exact_solution()", "PMSystemStokes");
//  return UA;
// }


// ==================================================================================
void PMSystemStokes::write_fluid_velocity_data(const std::string& filename)
{
  START_LOG("write_fluid_velocity_data()", "PMSystemStokes");

  // Check if the system solution vector is closed or not
  if ((this->solution->closed()) == false) this->solution->close();
  std::vector<Real> global_solution;
  this->solution->localize(global_solution);

  // Get the parameters
  MeshBase& mesh        = this->get_mesh();
  const std::size_t dim = mesh.mesh_dimension();

  // Use only one processor
  if (this->comm().rank() == 0) {
    // Output to data file
    std::ofstream out_file;
    out_file.open(filename, std::ios_base::out);
    out_file.setf(std::ios::right); out_file.setf(std::ios::fixed);
    out_file.precision(8); out_file.width(10);
    out_file << "% " << "x y z vx vy vz" << std::endl;

    // Loop through all nodes
    MeshBase::node_iterator nd           = mesh.nodes_begin();
    const MeshBase::node_iterator end_nd = mesh.nodes_end();

    for (; nd != end_nd; ++nd)
    {
      // Output nodal coordinates
      Node *node = *nd;
      out_file << (*node)(0) << " " << (*node)(1) << " " << (*node)(2) << " ";

      // Get the dof numbers at this node (only for velocity)
      std::vector<dof_id_type> dof_nums(dim);

      for (unsigned int i = 0; i < dim; ++i) { // var = 0, 1, 2 = i
        dof_nums[i] = node->dof_number(this->number(), i, 0);
        out_file << global_solution[dof_nums[i]] << " ";
      }

      out_file << std::endl;
    } // End looping nodes

    out_file.close();
  }
  this->comm().barrier();

  STOP_LOG("write_fluid_velocity_data()", "PMSystemStokes");
}

// ============================================================================
void PMSystemStokes::couple_poisson(const bool &add_local_solution_to_output)
{
  START_LOG("couple_poisson()", "PMSystemStokes");
  
  if (add_local_solution_to_output) {
    this->get_equation_systems().get_system<PMSystemPoisson>("Poisson").
      add_local_solution();
  }
  else {
    this->get_equation_systems().get_system<PMSystemPoisson>("Poisson").solve(
      "unused");
    this->get_equation_systems().get_system<PMSystemPoisson>("Poisson").
      add_electrostatic_forces();
  }    

  STOP_LOG("couple_poisson()", "PMSystemStokes");
}
} // end namespace libMesh
