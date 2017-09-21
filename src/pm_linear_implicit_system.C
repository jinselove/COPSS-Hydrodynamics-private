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

// libmesh headers
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmv_io.h"
#include "libmesh/vtk_io.h"


// particle-mesh header files
#include "pm_toolbox.h"
#include "ggem_system.h"
#include "brownian_system.h"
#include "pm_linear_implicit_system.h"
#include "analytical_solution.h"

// assemble functions
#include "assemble_navier_stokes.h"

namespace libMesh
{



// ======================================================================================
PMLinearImplicitSystem::PMLinearImplicitSystem(EquationSystems& es,
                                               const std::string& name,
                                               const unsigned int number)
: Parent (es, name, number),
  _stokes_solver(es),
  _point_mesh(NULL),
  _particle_mesh(NULL)
{
  // Assemble Navier Stokes
  _assemble_ns = ( new AssembleNS(es) );
}



// ==================================================================================
PMLinearImplicitSystem::~PMLinearImplicitSystem()
{
  // Clear data
  this->clear();
}




// ==================================================================================
void PMLinearImplicitSystem::clear ()
{
  // clear the parent data
  Parent::clear();
  
  // delete the pointer
  if(_assemble_ns) {
    delete _assemble_ns;
  }
}

  

// ===========================================================
void PMLinearImplicitSystem::assemble_matrix(const std::string& system_name,
                			     const std::string& option)
{
  libmesh_assert (this->matrix);
  libmesh_assert (this->matrix->initialized());
  
  START_LOG("assemble_matrix()", "PMLinearImplicitSystem");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   init the matrices: global stiffness and PC matrix (if required)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  this->matrix->zero();
  const bool user_defined_pc = this->get_equation_systems().parameters.get<bool>("user_defined_pc");
  if( user_defined_pc ) {
    std::cout << "--->test in PMLinearImplicitSystem::assemble_matrix(): "
              << "Initialize the preconditioning matrix (all zeros) \n";
    this->get_matrix("Preconditioner").zero();
  }
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Call assemble function to assemble matrix
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //assemble_matrix_sedimentation_ex1(this->get_equation_systems(), system_name, option);
  _assemble_ns->assemble_global_K(system_name, option);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   close the matrices
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  this->matrix->close();  // close the matrix
  if( user_defined_pc ) this->get_matrix("Preconditioner").close();
  
  
  STOP_LOG("assemble_matrix()", "PMLinearImplicitSystem");
}



// ==================================================================================
void PMLinearImplicitSystem::assemble_rhs(const std::string& system_name,
                                          const std::string& option)
{
  libmesh_assert (this->rhs);
  libmesh_assert (this->rhs->initialized());
  
  START_LOG("assemble_rhs()", "PMLinearImplicitSystem");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   init, assemble, and close the rhs vector
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  this->rhs->zero();
  //assemble_rhs_sedimentation_ex1 (this->get_equation_systems(), system_name, option);
  _assemble_ns->assemble_global_F(system_name, option);
  this->rhs->close();
  
  STOP_LOG("assemble_rhs()", "PMLinearImplicitSystem");
}
  

    
  
// ==================================================================================
void PMLinearImplicitSystem::solve_stokes (const std::string& option,
                                           const bool& re_init)
{
  START_LOG("solve_stokes()", "PMLinearImplicitSystem");
//  PerfLog perf_log("solve_stokes");
//  Real t1, t2;
  //std::string msg = "---> solve Stokes";
  //PMToolBox::output_message(msg, this->comm());
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Assemble the global matrix and pc matrix, and record the CPU wall time.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if(re_init)
  {
//    perf_log.push("assemble_matrix (undisturbed)");
  //  t1 = MPI_Wtime();
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     set the solver type for the Stokes equation
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    const StokesSolverType solver_type  = this->get_equation_systems().parameters.get<StokesSolverType> ("solver_type");
    _stokes_solver.set_solver_type(solver_type);
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Assemble the global matrix, and init the KSP solver
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    this->assemble_matrix("Stokes",option);
    _stokes_solver.init_ksp_solver();
//    perf_log.pop("assemble_matrix (undisturbed)");
    
  //  t2 = MPI_Wtime();
    //std::cout << "Time used to assemble the global matrix and reinit KSP is " <<t2-t1<<" s\n\n";
  }
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   assemble the rhs vector, and record the CPU wall time.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//  t1 = MPI_Wtime();
//  perf_log.push("assemble_rhs");
  this->assemble_rhs ("Stokes",option);
//  perf_log.pop("assemble_rhs");

//  t2 = MPI_Wtime();
  //std::cout << "Time used to assemble the right-hand-side vector is " <<t2-t1<<" s\n";
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   solve the problem
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//  perf_log.push("solve()");

  _stokes_solver.solve();
//  perf_log.pop("solve()");

  
  STOP_LOG("solve_stokes()", "PMLinearImplicitSystem");
}


  
  
// ==================================================================================
void PMLinearImplicitSystem::compute_point_velocity(const std::string& option,
                                                    std::vector<Real>& pv)
{
  //START_LOG("compute_point_velocity()", "PMLinearImplicitSystem");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   dim*NP: size of the velocity vector
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const std::size_t NP     =  _point_mesh->num_particles();
  const std::size_t dim    =  this->get_mesh().mesh_dimension();
  const dof_id_type n_elem =  this->get_mesh().n_elem();
  std::vector<Real> pvlocal(dim*NP,0.);  // declared on each processor
  std::vector<Real> _pv_send_list;                      // point velocity send list
  std::vector<std::size_t> _pid_send_list;              // point id send list
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over each point, and compute its velocity.
   Collect the global velocity from FEM through Allgather operation
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for(std::size_t i=0; i<NP; ++i)
  {
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     0. point coordinates & residing element
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    const Point pt            = _point_mesh->particles()[i]->point();
    const dof_id_type elem_id = _point_mesh->particles()[i]->elem_id();
    if(elem_id>n_elem)
    {
      printf("--->error in PMLinearImplicitSystem::compute_point_velocity() :\n");
      printf("    Point id = %u, elem id = %u, is out of simulation domain!\n",
             _point_mesh->particles()[i]->id(),elem_id);
      printf("    Total number of element in this mesh is %u!\n", n_elem);
    }
    const Elem* elem          = this->get_mesh().elem(elem_id);
    
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     1. Global(FEM) solution at the current point. This is done on local processors
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    std::vector<Real> Uglobal(dim);
    if(elem && (elem->processor_id()==this->processor_id()) )
    {
//      printf("--->test in compute_particle_velocity() particle id = %lu: (%f, %f, %f) rank = %i\n",
//             i,pt(0),pt(1),pt(2),this->processor_id() );
      
      // Numeric ids corresponding to each variable in the system
      const unsigned int u_var = this->variable_number ("u");      // u_var = 0
      const unsigned int v_var = this->variable_number ("v");      // v_var = 1
      Uglobal[0] = this->point_value(u_var, pt, *elem);
      Uglobal[1] = this->point_value(v_var, pt, *elem);
      if(dim==3)
      {
        const unsigned int w_var = this->variable_number ("w");      // w_var = 2
        Uglobal[2] = this->point_value(w_var, pt, *elem);
      }
      
      // Pack the particle id and its velocity
      _pid_send_list.push_back(i);
      for(std::size_t j=0; j<dim; ++j)
        _pv_send_list.push_back(Uglobal[j]);
      
      //printf("particle id = %lu  on processor %i\n", i, this->processor_id() );
      //printf("      Uglobal = (%f, %f, %f)\n", Uglobal[0], Uglobal[1], Uglobal[2]);
    } // end if
    
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     2. get the Free space solution (local solution) on all processors!
     No need to exclude the "self-term" for the local solution. => set "false"
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    if (option=="disturbed")
    {
      const std::string force_type = "regularized";
      const std::vector<Real> Ulocal = this->local_velocity_bead(i,force_type);

      // pass Ulacal to pvlocal
      for(std::size_t j=0; j<dim; ++j)
        pvlocal[dim*i+j] = Ulocal[j];
    } // end if(option)
    
  } // end for i-loop

  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Check the size of local_pv and the size of list on each process after allgather
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  this->comm().allgather(_pid_send_list); // allgather the particle id
  this->comm().allgather(_pv_send_list);  // allgather the particle velocity
  if (_pid_send_list.size() != NP)
  {
    libmesh_assert("*** error in PMLinearImplicitSystem::compute_particle_velocity()!");
    libmesh_error();
  }
//  printf("size of _pv_send_list is %lu, _pid_send_list.size() = %lu on the processor %i\n",
//         _pv_send_list.size(), _pid_send_list.size(), this->comm().rank() );
  
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Obtain the particle velocity vector excluding the global self exclusion.
   FIXME: Polymer beads and tracking points have the same radius value (=1)???
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //const Real mu0    = this->get_equation_systems().parameters.get<Real>("viscosity_0");
  const Real drag0  = 1.0;  // nondimensional drag coef = 6*PI*mu0*a where mu0 = 1/(6*PI), a=1
  for(std::size_t i=0; i<NP; ++i)
  {
    // id and position for the current bead
    const std::size_t pid      = _pid_send_list[i];
    const std::vector<Real> fv = _point_mesh->particles()[pid]->particle_force();
    const PointType point_type = _point_mesh->particles()[pid]->point_type();
    
    // exclusion of "self-term" for the global (FEM) solution
    const std::vector<Real> Uex = this->global_self_exclusion(pid);
    
    // Ug + Ul - U_self + Stokes drag
    for(std::size_t j=0; j<dim; ++j)
    {
      if (option=="disturbed") {
        if(point_type==POLYMER_BEAD) {          // the point is "bead"
          pv[dim*pid+j] = _pv_send_list[dim*i+j] + pvlocal[dim*pid+j]  - Uex[j] + fv[j]/drag0;
        }
        else if(point_type==LAGRANGIAN_POINT) { // the point is a tracking point
//          const Real hmin   = this->get_equation_systems().parameters.get<Real> ("solid mesh size");
//          const Real bead_r = std::sqrt(libMesh::pi)*hmin/4.0;
          pv[dim*pid+j] = _pv_send_list[dim*i+j] + pvlocal[dim*pid+j];// + fv[j]/(drag0*bead_r);
        }
      }
      else if (option=="undisturbed") {
        pv[dim*pid+j] = _pv_send_list[dim*i+j];
      }
      else {
        libmesh_error();
      } // end if-else
    } // end for j-loop

    // ---------------------------- output for test -----------------------------
//    if (this->comm().rank()==0 && option=="disturbed")
//    {
//      printf("\n--->test in compute_point_velocity(): output point velocity:\n");
//      printf("i = %lu, point %lu: Uglobal (FEM) = (%E, %E, %E)\n",
//             i, pid, _pv_send_list[dim*i], _pv_send_list[dim*i+1], _pv_send_list[dim*i+2] );
//      printf("              Ulocal (Green Function) = (%E, %E, %E)\n",
//             pvlocal[dim*i],pvlocal[dim*i+1],pvlocal[dim*i+2] );
//      printf("              Uexc    = (%E, %E, %E)\n", Uex[0], Uex[1], Uex[2] );
//      printf("       Stokes drag    = (%E, %E, %E)\n", fv[0], fv[1], fv[2] );
//      printf("           ---Utotal  = (%E, %E, %E)\n\n",
//             pv[dim*pid], pv[dim*pid+1], pv[dim*pid+2] );
//    }
    // --------------------------------------------------------------------------
  } // end for i-loop
  //STOP_LOG("compute_point_velocity()", "PMLinearImplicitSystem");
}
  
  
  
// ==================================================================================
std::vector<Real> PMLinearImplicitSystem::compute_unperturbed_point_velocity()
{
  START_LOG("compute_unperturbed_point_velocity()", "PMLinearImplicitSystem");
  const std::size_t NP     =  _point_mesh->num_particles();
  const std::size_t dim    =  this->get_mesh().mesh_dimension(); 
  // first solve the undisturbed flow field.
  const bool re_init = true;
  this->solve_stokes("undisturbed",re_init);
  std::vector<Real> pv_unperturbed(NP*dim, 0);
  // only FEM solution without particles! (undistrubed solution)
  this->compute_point_velocity("undisturbed", pv_unperturbed);
  
  STOP_LOG("compute_unperturbed_point_velocity()", "PMLinearImplicitSystem");
  return pv_unperturbed;
}

  
  
  
// ==================================================================================
std::vector<Real> PMLinearImplicitSystem::point_velocity(const std::vector<Real>& vel_beads,
                                                         const std::size_t i) const
{
  START_LOG("point_velocity()", "PMLinearImplicitSystem");
  
  const std::size_t dim  =  this->get_mesh().mesh_dimension();
  std::vector<Real> point_v(dim);
  for(std::size_t k=0; k<dim; ++k){
    point_v[k] = vel_beads[i*dim + k];
  }
  
  STOP_LOG("point_velocity()", "PMLinearImplicitSystem");
  return point_v;
}

// ==================================================================================
void PMLinearImplicitSystem::reinit_system()
{
  START_LOG("reinit_system()", "PMLinearImplicitSystem");
  this->comm().barrier(); // Is this at the beginning or the end necessary?
  
  // reinit point-mesh system, including
  // (1) build the point-point neighbor list according to search radius;
  // (2) build the element-point neighbor list according to search radius;
  // (3) evaluate forces
  const bool with_hi = true;
  _point_mesh->reinit(with_hi);  
  // update the tracking points position on the mesh if needed
  if(_particle_mesh != NULL){
    _point_mesh->update_particle_mesh(_particle_mesh);
  }   
  for (std::size_t i = 0; i < _fixes.size(); i++){
    // _fixes[i]->print_fix();
    _fixes[i]->compute();
  }

  STOP_LOG("reinit_system()", "PMLinearImplicitSystem");
}

// ==================================================================================
void PMLinearImplicitSystem::reinit_fd_system()
{
  START_LOG("reinit_fd_system()", "PMLinearImplicitSystem");
  this->comm().barrier(); // Is this at the beginning or the end necessary?
  
  // reinit point-mesh system, including
  // (1) build the point-point neighbor list according to search radius;
  // (2) set the elem_id and proc_id for points
  const bool with_hi = false;
  _point_mesh->reinit(with_hi);  
  for (std::size_t i = 0; i < _fixes.size(); i++){
    // _fixes[i]->print_fix();
    _fixes[i]->compute();
  }

  STOP_LOG("reinit_fd_system()", "PMLinearImplicitSystem");
}

  

// ==================================================================================
std::vector<Number> PMLinearImplicitSystem::local_velocity_fluid(const Point &p,
                                                                 const std::string& force_type) const
{
  START_LOG("local_velocity_fluid()", "PMLinearImplicitSystem");
  
  /* ----------------------------------------------------------------------------------------
   * get parameters of equation systems:
   * mu = 1/(6*PI) and bead_r0 = 1 are normalized viscosity and bead radius, respectively!
   ---------------------------------------------------------------------------------------- */
  const std::size_t dim  =  this->get_mesh().mesh_dimension();
  const Real         mu  =  this->get_equation_systems().parameters.get<Real> ("viscosity_0");
  const Real      alpha  =  this->get_equation_systems().parameters.get<Real> ("alpha");
  const Real    bead_r0  =  this->get_equation_systems().parameters.get<Real> ("br0");
  
  
  /* ----------------------------------------------------------------------------------------
   FIXME: When system contains different "particle_type", for example, interaction of
   particle and polymer chain, this becomes invalid!
   ---------------------------------------------------------------------------------------- */
  const std::string particle_type
           = this->get_equation_systems().parameters.get<std::string> ("particle_type");
  
  
  /* ----------------------------------------------------------------------------------------
   Mesh size for fluid and solid, which will be used to get the local solution.
   For point particle simulation, there is no solid mesh. "hmin" actually is NOT used.
   For non-point particle cases, hmin will be used to evaluate the screening parameter ksi
   ---------------------------------------------------------------------------------------- */
  Real  hmin =  this->get_equation_systems().parameters.get<Real> ("minimum fluid mesh size");
  if(particle_type=="rigid_particle"){
    hmin  =  this->get_equation_systems().parameters.get<Real> ("minimum solid mesh size");
  }
  
  
  /* ----------------------------------------------------------------------------------------
   compute the local velocity corresponding to the unbounded domain
   ---------------------------------------------------------------------------------------- */
  GGEMSystem ggem_sys;
  std::vector<Real> Ulocal = ggem_sys.local_velocity_fluid(_point_mesh,p,alpha,mu,bead_r0,
                                                           hmin,dim,force_type);
 
  
  STOP_LOG("local_velocity_fluid()", "PMLinearImplicitSystem");
  return Ulocal;
}
  
  

// ==================================================================================
std::vector<Number> PMLinearImplicitSystem::local_velocity_bead(const std::size_t& bead_id,
                                                                const std::string& force_type) const
{
  START_LOG("local_velocity_bead()", "PMLinearImplicitSystem");
  
  /* ----------------------------------------------------------------------------------------
   * get parameters of equation systems:
   * mu = 1/(6*PI) and bead_r0 = 1 are normalized viscosity and bead radius, respectively!
   ---------------------------------------------------------------------------------------- */
  const std::size_t dim  =  this->get_mesh().mesh_dimension();
  const Real         mu  =  this->get_equation_systems().parameters.get<Real> ("viscosity_0");
  const Real      alpha  =  this->get_equation_systems().parameters.get<Real> ("alpha");
  const Real    bead_r0  =  this->get_equation_systems().parameters.get<Real> ("br0");
  
  
  /* ----------------------------------------------------------------------------------------
   FIXME: When system contains different "particle_type", for example, interaction of
   // particle and polymer chain, this becomes invalid!
   ---------------------------------------------------------------------------------------- */
  const std::string particle_type
        = this->get_equation_systems().parameters.get<std::string> ("particle_type");
  
  
  /* ----------------------------------------------------------------------------------------
   Mesh size for fluid and solid, which will be used to get the local solution.
   For point particle simulation, there is no solid mesh. "hmin" actually is NOT used.
   For non-point particle cases, hmin will be used to evaluate the screening parameter ksi
   ---------------------------------------------------------------------------------------- */
  Real  hmin =  this->get_equation_systems().parameters.get<Real> ("minimum fluid mesh size");
  if(particle_type == "rigid_particle"){
    hmin  =  this->get_equation_systems().parameters.get<Real> ("minimum solid mesh size");
  }
  
  
  /* ----------------------------------------------------------------------------------------
   compute the local velocity corresponding to the unbounded domain
   ---------------------------------------------------------------------------------------- */
  GGEMSystem ggem_sys;
  std::vector<Real> Ulocal = ggem_sys.local_velocity_bead(_point_mesh,bead_id,alpha,mu,bead_r0,
                                                          hmin,dim,force_type);
  
  
  STOP_LOG("local_velocity_bead()", "PMLinearImplicitSystem");
  return Ulocal;
}

  

// ==================================================================================
std::vector<Real> PMLinearImplicitSystem::global_self_exclusion(const std::size_t p_id) const
{
  START_LOG("global_self_exclusion()", "PMLinearImplicitSystem");
  
  /* ----------------------------------------------------------------------------------------
   get parameters of equation systems.
   ---------------------------------------------------------------------------------------- */
  const std::size_t dim  =  this->get_mesh().mesh_dimension();
  const Real mu      =  this->get_equation_systems().parameters.get<Real> ("viscosity_0");
  const Real alpha   =  this->get_equation_systems().parameters.get<Real> ("alpha");
  
  
  /* ----------------------------------------------------------------------------------------
   compute the global self-exclusion velocity
   ---------------------------------------------------------------------------------------- */
  GGEMSystem ggem_sys;
  std::vector<Real> self_v = ggem_sys.global_self_exclusion(_point_mesh,p_id,alpha,mu,dim);
  
  STOP_LOG("global_self_exclusion()", "PMLinearImplicitSystem");
  return self_v;
}

  
  
// ==================================================================================
void PMLinearImplicitSystem::test_l2_norm()
{
  START_LOG("test_l2_norm()", "PMLinearImplicitSystem");
  std::string msg = "--->test in PMLinearImplicitSystem::test_l2_norm(): \n";
  PMToolBox::output_message(msg, this->comm());
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Numerical solution: Global(FEM) + Local(Analytical)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  this->reinit_system();       // re-init particle-mesh before start
  const bool re_init = true;
  this->solve_stokes("disturbed",re_init);
  this->add_local_solution();
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Theoretical solution(only velocity)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const unsigned int dim        = 3;
  MeshBase&         mesh        = this->get_mesh();
  const unsigned int n_nodes    = mesh.n_nodes();
  Real val0_norm = 0., val1_norm = 0.;
  Real val2_norm = 0., val3_norm = 0.;
  AnalyticalSolution analytical_solution(*this);
  
  // Loop over each node and compute the nodal velocity value
  MeshBase::node_iterator       nd     = mesh.local_nodes_begin();
  const MeshBase::node_iterator end_nd = mesh.local_nodes_end();
  for ( ; nd != end_nd; ++nd)
  {
    // Store a pointer to the current node, and extract a point
    Node* node = *nd;
    Point pt;
    for(unsigned int i=0; i<dim; ++i)pt(i) =  (*node)(i) ;
    
    // get the dof numbers at this node (only for velocity)
    std::vector<dof_id_type> dof_nums(dim);
    for(unsigned int i=0; i<dim; ++i){ // var = 0, 1, 2 = i
      dof_nums[i] = node->dof_number(this->number(), i, 0);
    }
    
    // Get the numerical solution
    std::vector<Real> Unum;
    this->solution->get(dof_nums, Unum);
    
    // compute the local velocity of fluid at the current node
    //const std::vector<Real> Uexact = this->exact_solution(pt);
    const std::vector<Real> Uexact = analytical_solution.exact_solution_infinite_domain(pt);
    
    // compute the errors
    for(unsigned int i=0; i<dim; ++i){
      Real tmpt  = std::abs(Unum[i] - Uexact[i]);
      val0_norm += tmpt;
      val1_norm += std::abs(Uexact[i]);
      
      val2_norm += tmpt*tmpt;
      val3_norm += Uexact[i]*Uexact[i];
    }
    
  } // end for nd-loop

  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Compute the error: l1 and l2-norm of errors
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  this->comm().sum(val0_norm);
  this->comm().sum(val1_norm);
  this->comm().sum(val2_norm);
  this->comm().sum(val3_norm);
  
  const Real l1_norm = val0_norm/val1_norm;
  const Real l2_norm = std::sqrt(val2_norm)/std::sqrt(val3_norm);
  printf("--->test in test_l1_norm: l1_norm = %E, l2_norm = %E\n", l1_norm, l2_norm);
  
  STOP_LOG("test_l2_norm()", "PMLinearImplicitSystem");
}
  

  
// ==================================================================================
void PMLinearImplicitSystem::test_velocity_profile()
{
  START_LOG("test_velocity_profile()", "PMLinearImplicitSystem");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   solve the disturbed velocity field: global solution(FEM solution)
   In this test, we assume that undisturbed velocity is zero!
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  std::cout<< "========>2. Test in PMLinearImplicitSystem::test_velocity_profile(): \n";
  
  this->reinit_system();       // re-init particle-mesh before start
  const bool re_init = true;
  this->solve_stokes("disturbed",re_init);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   output the velocity profiles along xyz directions, global + local solutions.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const Point& box_min = _point_mesh->pm_periodic_boundary()->box_min();
  const Point& box_len = _point_mesh->pm_periodic_boundary()->box_length();
  const Real xn = 200, yn = 80, zn = 80;
  const Real dx = box_len(0)/xn, dy = box_len(1)/yn, dz = box_len(2)/zn;
  const unsigned int NP = _point_mesh->num_particles();
  AnalyticalSolution analytical_solution(*this);
  
  std::ofstream outfile;
  int o_width = 12, o_precision = 9;
  
  // 1. write out the velocity profile along x-direction.
  std::ostringstream filenamex;
  filenamex << "output_velocity_profile_x_" << NP << "P.txt";
  outfile.open(filenamex.str(), std::ios_base::out);
  for(std::size_t i=0; i<xn+1; ++i)
  {
    Point pt (box_min(0)+Real(i)*dx, 0. ,0.);
    
    // global velocity from FEM
    std::vector<Real> Uglobal(3), Utotal(3);
    const unsigned int u_var = this->variable_number ("u");      // u_var = 0
    const unsigned int v_var = this->variable_number ("v");      // v_var = 1
    const unsigned int w_var = this->variable_number ("w");      // w_var = 2
    Uglobal[0] = this->point_value(u_var, pt);  // this is slow, but only for test.
    Uglobal[1] = this->point_value(v_var, pt);
    Uglobal[2] = this->point_value(w_var, pt);
    
    // local velocity from analytical function
    const std::vector<Real> Ulocal = this->local_velocity_fluid(pt,"regularized");
    for(std::size_t j=0; j<3; ++j) Utotal[j] = Uglobal[j] + Ulocal[j];
    
    // Exact solution for an unbounded domain
    const std::vector<Real> Uexact = analytical_solution.exact_solution_infinite_domain(pt);
    
    // write the velocity, x vx vy vz
    outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
    outfile.precision(o_precision);   outfile.width(o_width);
    if( this->comm().rank()==0 )
      outfile << pt(0) << "  " << Utotal[0] << "  " << Utotal[1] << "  " << Utotal[2]
                       << "  " << Uglobal[0]<< "  " << Uglobal[1]<< "  " << Uglobal[2]
                       << "  " << Ulocal[0] << "  " << Ulocal[1] << "  " << Ulocal[2]
                       << "  " << Uexact[0] << "  " << Uexact[1] << "  " << Uexact[2] <<"\n";
  } // end for i-loop
  outfile.close();
  
  // 2. write out the velocity profile along y-direction.
  std::ostringstream filenamey;
  filenamey << "output_velocity_profile_y_" << NP << "P.txt";
  outfile.open(filenamey.str(), std::ios_base::out);
  for(std::size_t i=0; i<yn+1; ++i)
  {
    Point pt(0., box_min(1)+Real(i)*dy, 0.);
    
    // global velocity from FEM
    std::vector<Real> Uglobal(3), Utotal(3);
    const unsigned int u_var = this->variable_number ("u");      // u_var = 0
    const unsigned int v_var = this->variable_number ("v");      // v_var = 1
    const unsigned int w_var = this->variable_number ("w");      // w_var = 2
    Uglobal[0] = this->point_value(u_var, pt);  // this is slow, but only for test.
    Uglobal[1] = this->point_value(v_var, pt);
    Uglobal[2] = this->point_value(w_var, pt);
    
    // local velocity from analytical function
    const std::vector<Real> Ulocal = this->local_velocity_fluid(pt,"regularized");
    for(std::size_t j=0; j<3; ++j) Utotal[j] = Uglobal[j] + Ulocal[j];
    
    // Exact solution for an unbounded domain
    const std::vector<Real> Uexact = analytical_solution.exact_solution_infinite_domain(pt);
    
    // write the velocity, y vx vy vz
    outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
    outfile.precision(o_precision);   outfile.width(o_width);
    if( this->comm().rank()==0 )
      outfile << pt(1) << "  " << Utotal[0] << "  " << Utotal[1] << "  " << Utotal[2]
                       << "  " << Uglobal[0]<< "  " << Uglobal[1]<< "  " << Uglobal[2]
                       << "  " << Ulocal[0] << "  " << Ulocal[1] << "  " << Ulocal[2]
                       << "  " << Uexact[0] << "  " << Uexact[1] << "  " << Uexact[2] <<"\n";
  } // end for i-loop
  outfile.close();
  
  // 3. write out the velocity profile along z-direction.
  std::ostringstream filenamez;
  filenamez << "output_velocity_profile_z_" << NP << "P.txt";
  outfile.open(filenamez.str(), std::ios_base::out);
  for(std::size_t i=0; i<zn+1; ++i)
  {
    Point pt(0., 0., box_min(2)+Real(i)*dz);
    
    // global velocity from FEM
    std::vector<Real> Uglobal(3), Utotal(3);
    const unsigned int u_var = this->variable_number ("u");      // u_var = 0
    const unsigned int v_var = this->variable_number ("v");      // v_var = 1
    const unsigned int w_var = this->variable_number ("w");      // w_var = 2
    Uglobal[0] = this->point_value(u_var, pt);  // this is slow, but only for test.
    Uglobal[1] = this->point_value(v_var, pt);
    Uglobal[2] = this->point_value(w_var, pt);
    
    // local velocity from analytical function
    const std::vector<Real> Ulocal = this->local_velocity_fluid(pt,"regularized");
    for(std::size_t j=0; j<3; ++j) Utotal[j] = Uglobal[j] + Ulocal[j];
    
    // Exact solution for an unbounded domain
    const std::vector<Real> Uexact = analytical_solution.exact_solution_infinite_domain(pt);
    
    // write the velocity, z vx vy vz
    outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
    outfile.precision(o_precision);   outfile.width(o_width);
    if( this->comm().rank()==0 )
      outfile << pt(2) << "  " << Utotal[0] << "  " << Utotal[1]<< "  " <<  Utotal[2]
                       << "  " << Uglobal[0]<< "  " << Uglobal[1]<< "  " << Uglobal[2]
                       << "  " << Ulocal[0] << "  " << Ulocal[1]<< "  " <<  Ulocal[2]
                       << "  " << Uexact[0] << "  " << Uexact[1] << "  " << Uexact[2] <<"\n";
  } // end for i-loop
  outfile.close();
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   done and write out the results
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  std::ostringstream output_filename;
  output_filename << "output_velocity_profile_" << NP << "P.e";
#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO(this->get_mesh()).write_equation_systems(output_filename.str(),
                                                   this->get_equation_systems() );
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
  
  
  STOP_LOG("test_velocity_profile()", "PMLinearImplicitSystem");
}

  

// ==================================================================================
//const std::vector<Real> PMLinearImplicitSystem::exact_solution(const Point& pt0) const
//{
//  START_LOG("exact_solution()", "PMLinearImplicitSystem");
//  
//  // ksi's value should be consistent with that in GGEMSystem::regularization_parameter()
//  const Real ksi = std::sqrt(libMesh::pi)/3.0;  // = 0.591
//  const Real muc = 1.0/(6*libMesh::pi);
//  const unsigned int dim = 3;
//  std::vector<Real> UA(dim,0.);
//  DenseMatrix<Number> GT;
//  
//  // GGEM object and number of points in the system
//  GGEMSystem ggem_system;
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
//    GT = ggem_system.green_tensor_exp(x,ksi,muc,dim,zero_limit);
//    const std::vector<Real> fv = _point_mesh->particles()[i]->particle_force();
//    //printf("--->test in exact_solution(): i = %lu, fv = (%f,%f,%f)\n", i,fv[0],fv[1],fv[2]);
//    
//    // 3. compute u due to this particle
//    for (std::size_t k=0; k<dim; ++k){
//      for (std::size_t l=0; l<dim; ++l){
//        UA[k] += GT(k,l)*fv[l];
//      } // end for l
//    } // end for k
//  } // end for i
//  
//  STOP_LOG("exact_solution()", "PMLinearImplicitSystem");
//  return UA;
//}
  
  

// ==================================================================================
void PMLinearImplicitSystem::write_out_single_particle(const Point& coords,
                                                       const std::vector<Real>& vel,
                                                       const int i_step,
                                                       const Real time) const
{
  START_LOG("write_out_single_particle()", "PMLinearImplicitSystem");

  this->comm().barrier();
  std::size_t dim = this->get_mesh().mesh_dimension();
  std::string filename = "single_particle_history.txt";
  std::ofstream outfile;
  int o_width = 15, o_precision = 9;
  
  // the first step, ios_base::out, File open for writing
  if(i_step<=0)
  {
    outfile.open(filename,std::ios_base::out);
    outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
    outfile.precision(o_precision);   outfile.width(o_width/3);
    if( this->comm().rank()==0 ) outfile << i_step << "  " << time;
    
    for ( std::size_t i=0; i<dim; ++i )
    {
      outfile.setf(std::ios::right);  outfile.setf(std::ios::fixed);
      outfile.precision(o_precision); outfile.width(o_width);
      if( this->comm().rank()==0 )    outfile << coords(i);
    }
    for ( std::size_t i=0; i<dim; ++i )
    {
      outfile.setf(std::ios::right);  outfile.setf(std::ios::fixed);
      outfile.precision(o_precision); outfile.width(o_width);
      if( this->comm().rank()==0 )    outfile << vel[i];
    }
    if( this->comm().rank()==0 ) outfile << "\n";
  }
  
  // the following steps, ios_base::app, output operations happen at the end of the file,
  // appending to its existing contents.
  if(i_step>0)
  {
    outfile.open(filename,std::ios_base::app);
    outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
    outfile.precision(o_precision);   outfile.width(o_width);
    if( this->comm().rank()==0 ) outfile << i_step << "  " << time;
    
    for ( std::size_t i=0; i<dim; ++i )
    {
      outfile.setf(std::ios::right);  outfile.setf(std::ios::fixed);
      outfile.precision(o_precision); outfile.width(o_width);
      if( this->comm().rank()==0 )    outfile << coords(i);
    }
    for ( std::size_t i=0; i<dim; ++i )
    {
      outfile.setf(std::ios::right);  outfile.setf(std::ios::fixed);
      outfile.precision(o_precision); outfile.width(o_width);
      if( this->comm().rank()==0 )    outfile << vel[i];
    }
    if( this->comm().rank()==0 ) outfile << "\n";
  }
  
  // close the file
  outfile.close();
  this->comm().barrier();
  
  STOP_LOG("write_out_single_particle()", "PMLinearImplicitSystem");
}



// ==================================================================================
void PMLinearImplicitSystem::write_out_point_coordinate(Vec* vin,
                                                        const std::size_t istep,
                                                        const Real& time,
                                                        const std::string& filename,
                                                        const std::string& openmode) const
{
  START_LOG("write_out_particle_coordinate()", "PMLinearImplicitSystem");
  
  VecScatter ctx;
  Vec vout;
  const PetscScalar     *px;
  
  //  PetscPrintf(PETSC_COMM_WORLD,"--->test in write_out_particle_coordinate vin = \n");
  //  VecView(*vin,PETSC_VIEWER_STDOUT_WORLD);  // View the random vector
  
  VecScatterCreateToZero(*vin,&ctx,&vout);
  VecScatterBegin(ctx,*vin,vout,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,*vin,vout,INSERT_VALUES,SCATTER_FORWARD);
  VecGetArrayRead(vout,&px);
  
  //  PetscPrintf(PETSC_COMM_WORLD,"--->test in write_out_particle_coordinate vout = \n");
  //  VecView(vout,PETSC_VIEWER_STDOUT_WORLD);  // View the random vector
  
  //
  const std::size_t NP = _point_mesh->num_particles();
  const std::size_t dim = this->get_mesh().mesh_dimension();
  
  std::ofstream outfile;
  const int o_width = 5, o_precision = 9;
  if(this->comm().rank()==0)
  {
    if(openmode=="out"){
      outfile.open(filename,std::ios_base::out);
    }
    else if(openmode=="app"){
      outfile.open(filename,std::ios_base::app);
    }
    else{
      outfile.open(filename,std::ios_base::app);
    }
    // end of if-else
    
    if(istep==0) outfile << NP << "\n";
    outfile << "Step " << istep << " time = " <<time<< "\n";
    for (std::size_t i=0; i<NP; ++i)
    {
      outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
      outfile.precision(o_precision);   outfile.width(o_width);
      outfile << i << " ";
      for(std::size_t j=0; j<dim; ++j)
      {
        Real xyz = Real( px[i*dim+j] );
        outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
        outfile.precision(o_precision);   outfile.width(o_width);
        outfile << xyz << "  ";
      }
      outfile << "\n";
    }
    
    outfile << "End Step " << istep << "\n\n";
    outfile.close();
  }
  
  // restore vector and destroy ctx
  VecRestoreArrayRead(vout,&px);
  VecScatterDestroy(&ctx);
  VecDestroy(&vout);
  
  STOP_LOG("write_out_particle_coordinate()", "PMLinearImplicitSystem");
}



  
// ==================================================================================
void PMLinearImplicitSystem::write_point_csv(const std::string& filename,
                                             const std::vector<Real>& pv,
                                             const bool write_velocity) const
{
  START_LOG("write_point_csv()", "PMLinearImplicitSystem");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Write out the CSV file on the 0-th processor
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const std::size_t NP  = _point_mesh->num_particles();
  const std::size_t dim = this->get_mesh().mesh_dimension();
  std::ofstream outfile;
  const int o_width = 5, o_precision = 9;
  if(this->comm().rank()==0)
  {
    outfile.open(filename,std::ios_base::out);
    outfile << "Point #, "<< "X Coord, "<< "Y Coord, "<< "Z Coord";
    if(write_velocity)
      outfile << ",  Vx,      "<< "Vy,       "<< "Vz,      "<< "Vmag ";
    outfile << "\n";
    
    for (std::size_t i=0; i<NP; ++i)
    {
      // write the point id and coordinates
      outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
      outfile.precision(o_precision);   outfile.width(o_width);
      outfile << i << " ";
      for(std::size_t j=0; j<dim; ++j)
      {
        const Real xyz = _point_mesh->particles()[i]->point()(j);
        outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
        outfile.precision(o_precision);   outfile.width(o_width);
        outfile << ",  "<< xyz;
      } // end for j-loop
      
      // write the point velocity if needed
      if(write_velocity)
      {
        Real Vmag = 0.0;
        for(std::size_t j=0; j<dim; ++j)
        {
          const Real vxyz = pv[i*dim + j];
          outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
          outfile.precision(o_precision);   outfile.width(o_width);
          outfile << ",  "<< vxyz;
          Vmag += vxyz*vxyz;
        }
        Vmag = std::sqrt(Vmag);
        outfile << ",  "<< Vmag;
      }
      
      outfile << "\n";
    }
    outfile.close();
  }
  
  
  STOP_LOG("write_point_csv()", "PMLinearImplicitSystem");
}



// ==================================================================================
PetscErrorCode PMLinearImplicitSystem::write_point_csv(const std::string& filename,
                                                       Vec * petsc_vec,
                                                       const bool write_velocity) const
{
  PetscInt          low, high, nlocal, vsize;
  PetscScalar       *pv;
  PetscErrorCode    ierr;
  PetscFunctionBeginUser;
  START_LOG ("write_point_csv()", "PMLinearImplicitSystem");
  
  std::vector<Real> std_vec;
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Get the Ownership range and local components, then copy to a std::vector!
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (petsc_vec && write_velocity)
  {
    ierr = VecGetSize(*petsc_vec, &vsize);CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(*petsc_vec,&low,&high); CHKERRQ(ierr);
    ierr = VecGetLocalSize(*petsc_vec,&nlocal);         CHKERRQ(ierr);
    ierr = VecGetArray(*petsc_vec,&pv);                 CHKERRQ(ierr);
    
    std_vec.resize( (std::size_t)nlocal );
    for(int i=0; i<nlocal; ++i) std_vec[i] = pv[i];
    this->comm().allgather(std_vec);
    ierr = VecRestoreArray(*petsc_vec,&pv);             CHKERRQ(ierr);
  }
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Call the member function to write out point CSV file
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  this->write_point_csv(filename, std_vec, write_velocity);
  
  STOP_LOG("write_point_csv()", "PMLinearImplicitSystem");
  PetscFunctionReturn(0);
}


  
// ==================================================================================
void PMLinearImplicitSystem::write_equation_systems(const std::size_t time_step,
                                                    const std::string& output_filename,
                                                    const std::string& output_format)
{
  START_LOG("write_equation_systems()", "PMLinearImplicitSystem");
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Write out the FEM results: global solution
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  MeshBase&         mesh    = this->get_mesh();
  const bool fem_sol_output = true;
  if(fem_sol_output)
  {
    std::ostringstream file_name_fem;
    file_name_fem << output_filename + "_fem";
    
    if (output_format=="EXODUS")
    {
#ifdef LIBMESH_HAVE_EXODUS_API
      
      if(time_step==0)
      {
        file_name_fem << ".e";
        ExodusII_IO(mesh).write_equation_systems(file_name_fem.str(),this->get_equation_systems());
      }
      else
      {
//        file_name_fem << ".e-s." << std::setw(8) << std::setfill('0') << std::right << time_step;
//        ExodusII_IO(mesh).write_equation_systems(file_name_fem.str(),this->get_equation_systems());
        
        file_name_fem << ".e";
        ExodusII_IO exodus_IO(mesh);
        exodus_IO.append(true);
        exodus_IO.write_timestep(file_name_fem.str(),this->get_equation_systems(),
                                 time_step+1,time_step+1);
      } // end if-else
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
    }
    else if (output_format=="VTK")
    {
#ifdef LIBMESH_HAVE_VTK
      file_name_fem <<"_"<< std::setw(8) << std::setfill('0') << std::right << time_step<< ".vtu";
      VTKIO(mesh).write_equation_systems (file_name_fem.str(), this->get_equation_systems());
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
    }
    else
    {
      file_name_fem <<"_"<< std::setw(8) << std::setfill('0') << std::right << time_step<< ".gmv.";
      GMVIO(mesh).write_equation_systems (file_name_fem.str(), this->get_equation_systems());
    }
  }
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Update the system solution by adding the local solution (from Green's function)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  this->add_local_solution();
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Write out the FEM results: global solution
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  std::ostringstream file_name;
  file_name << output_filename + "_total";
  
  if (output_format=="EXODUS")
  {
#ifdef LIBMESH_HAVE_EXODUS_API
    if(time_step==0)
    {
      file_name << ".e";
      ExodusII_IO(mesh).write_equation_systems(file_name.str(),this->get_equation_systems());
    }
    else
    {
//      file_name << ".e-s." << std::setw(8) << std::setfill('0') << std::right << time_step;
//      ExodusII_IO(mesh).write_equation_systems(file_name.str(),this->get_equation_systems());
      
      file_name << ".e";
      ExodusII_IO exodus_IO(mesh);
      exodus_IO.append(true);
      exodus_IO.write_timestep(file_name.str(),this->get_equation_systems(),
                               time_step+1,time_step+1);
    } // end if-else
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
  }
  else if (output_format=="VTK")
  {
#ifdef LIBMESH_HAVE_VTK
    file_name <<"_"<< std::setw(8) << std::setfill('0') << std::right << time_step<< ".vtu";
    VTKIO(mesh).write_equation_systems (file_name.str(), this->get_equation_systems());
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
  }
  else
  {
    file_name <<"_"<< std::setw(8) << std::setfill('0') << std::right << time_step<< ".gmv";
    GMVIO(mesh).write_equation_systems (file_name.str(), this->get_equation_systems());
  }
  
    
  STOP_LOG("write_equation_systems()", "PMLinearImplicitSystem");
}



  
// ==================================================================================
void PMLinearImplicitSystem::add_local_solution()
{
  START_LOG("add_local_solution()", "PMLinearImplicitSystem");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Check if the system solution vector is closed or not
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  if( (this->solution->closed())==false ) this->solution->close();
  //this->update();
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Get the parameters and Initialize the quantities
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  const bool test_output    = false;
  MeshBase&         mesh    = this->get_mesh();
  const std::size_t  dim    = mesh.mesh_dimension();
  const std::size_t n_local_nodes = mesh.n_local_nodes();
  std::vector<Number>    local_solution(dim*n_local_nodes);
  std::vector<numeric_index_type> dof_indices(dim*n_local_nodes);
  //printf("--->test in add_local_solution() n_local_nodes = %lu on the processor %u\n",
  //       n_local_nodes,this->comm().rank());
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Update the system solution by adding the local solution (from Green's function)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  MeshBase::node_iterator       nd     = mesh.local_nodes_begin();
  const MeshBase::node_iterator end_nd = mesh.local_nodes_end();
  std::size_t local_count = 0;
  for ( ; nd != end_nd; ++nd)
  {
    // Store a pointer to the current node, and extract a point
    Node* node = *nd;
    Point pt;
    for(unsigned int i=0; i<dim; ++i)pt(i) =  (*node)(i) ;
     
    // ---------- this is a test for dof_number at each node -----------
    if (test_output)
    {
      const unsigned int node_id = node->id();
      std::ostringstream oss;
      oss << "          NODE " << node_id;
      PMToolBox::output_message(oss.str(), this->comm());
      node->print_info();
      
      if(this->comm().rank()==0) printf("--->test: nodal dof number :");
      for(unsigned int i=0; i<dim; ++i)
      {
        dof_id_type dof_num = node->dof_number(this->number(), i, 0);
        if(this->comm().rank()==0) printf(" %u", dof_num);
      }
      if(this->comm().rank()==0) printf(" \n");
    }
     
    // get the dof numbers at this node (only for velocity)
    std::vector<dof_id_type> dof_nums(dim);
    for(unsigned int i=0; i<dim; ++i){ // var = 0, 1, 2 = i
      dof_nums[i] = node->dof_number(this->number(), i, 0);
    }
    
    // compute the local velocity of fluid at the current node
    const std::vector<Real> Ulocal = this->local_velocity_fluid(pt,"regularized");

    // store the local velocity and dof indices
    for(unsigned int i=0; i<dim; ++i)
    {
      local_solution[local_count*dim+i]  =  Ulocal[i];
      dof_indices   [local_count*dim+i]  =  dof_nums[i];
    }
    local_count++;
  } // end for

  //printf("--->test in add_local_solution() local_count = %lu on the processor %u\n",
  //       n_local_nodes,this->comm().rank());
  
  // add the local to the global
  //this->solution->zero();
  this->solution->add_vector(local_solution, dof_indices);
  this->solution->close();
  this->update();
 
  STOP_LOG("add_local_solution()", "PMLinearImplicitSystem");
}
 

 
// ==================================================================================
void PMLinearImplicitSystem::write_fluid_velocity_data(const std::string& filename)
{
  START_LOG("write_fluid_velocity_data()", "PMLinearImplicitSystem");
 
  // Check if the system solution vector is closed or not
  if( (this->solution->closed())==false ) this->solution->close();
  std::vector<Real> global_solution;
  this->solution->localize(global_solution);

  // Get the parameters
  MeshBase&         mesh = this->get_mesh();
  const std::size_t  dim = mesh.mesh_dimension();

  // Use only one processor
  if(this->comm().rank()==0){
    // Output to data file
    std::ofstream out_file;
    out_file.open(filename, std::ios_base::out);
    out_file.setf(std::ios::right); out_file.setf(std::ios::fixed);
    out_file.precision(8); out_file.width(10);
    out_file << "% " << "x y z vx vy vz" << std::endl;
 
    // Loop through all nodes
    MeshBase::node_iterator       nd     = mesh.nodes_begin();
    const MeshBase::node_iterator end_nd = mesh.nodes_end();
    for ( ; nd != end_nd; ++nd)
    {
      // Output nodal coordinates
      Node* node = *nd;
      out_file << (*node)(0) << " " << (*node)(1) << " " << (*node)(2) << " ";
  
      // Get the dof numbers at this node (only for velocity)
      std::vector<dof_id_type> dof_nums(dim);
      for(unsigned int i=0; i<dim; ++i){ // var = 0, 1, 2 = i
        dof_nums[i] = node->dof_number(this->number(), i, 0);
        out_file << global_solution[dof_nums[i]] << " ";
      }
  
      out_file << std::endl;
  
    } // End looping nodes
  
    out_file.close();
  }
  this->comm().barrier();

  STOP_LOG("write_fluid_velocity_data()", "PMLinearImplicitSystem");
}


} // end namespace libMesh

