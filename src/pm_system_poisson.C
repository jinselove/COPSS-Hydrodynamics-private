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
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"

// local header files
#include "pm_toolbox.h"
#include "ggem_poisson.h"
#include "pm_system_poisson.h"


namespace libMesh {
// ======================================================================================
PMSystemPoisson::PMSystemPoisson(EquationSystems  & es,
                                 const std::string& name,
                                 const unsigned int number)
  : PMLinearImplicitSystem(es, name, number),
  _solver_poisson(es)
{
  if (name != "Poisson") 
  {
    PMToolBox::output_message("Error: system name in PMSystemPoisson initializer is not 'Poisson'. Exiting ..."
      , this->comm());
    libmesh_error();
  }
  // Poisson equation assembly
  _assemble_poisson   = new AssemblePoisson(es, name);
  analytical_solution = _assemble_poisson->get_analytical_solution();
  ggem_poisson        = _assemble_poisson->get_ggem_poisson();
  o_precision = 6;
}

// ==================================================================================
PMSystemPoisson::~PMSystemPoisson()
{
  // Clear data
  this->clear();
}

// ==================================================================================
void PMSystemPoisson::clear()
{
  // delete the pointer
  if (_assemble_poisson) {
    delete _assemble_poisson;
  }
}

// ===========================================================
void PMSystemPoisson::assemble_matrix(const std::string& system_name,
                                      const std::string& option)
{
  libmesh_assert(this->matrix);
  libmesh_assert(this->matrix->initialized());

  START_LOG("assemble_matrix()", "PMSystemPoisson");

  // init the matrices: global stiffness and PC matrix (if required)
  this->matrix->zero();
  // Call assemble function to assemble matrix
  _assemble_poisson->assemble_global_K(system_name, option);
  // close the matrices
  this->matrix->close(); // close the matrix
  STOP_LOG("assemble_matrix()", "PMSystemPoisson");
}

// ==================================================================================
void PMSystemPoisson::assemble_rhs(const std::string& system_name,
                                   const std::string& option)
{
  libmesh_assert(this->rhs);
  libmesh_assert(this->rhs->initialized());

  START_LOG("assemble_rhs()", "PMSystemPoisson");

  // init, assemble, and close the rhs vector
  this->rhs->zero();
  // assemble rhs vector
  _assemble_poisson->assemble_global_F(system_name, option);
  // close rhs vector
  this->rhs->close();

  STOP_LOG("assemble_rhs()", "PMSystemPoisson");
}

// ==================================================================================
void PMSystemPoisson::solve(const std::string& option)
{
  START_LOG("solve()", "PMSystemPoisson");

  // Assemble the global matrix and pc matrix at the first step, when
  // re_init=true.
  if (_re_init)
  {
    // Set the solver type for the Poisson equation
    const SystemSolverType solver_type =
      this->get_equation_systems().parameters.get<SystemSolverType>(
        "solver_type_poisson");
    _solver_poisson.set_solver_type(solver_type);
    // Assemble the global matrix, and init the KSP solver
    this->assemble_matrix(this->name(), option);
    _solver_poisson.init_ksp_solver(this->name());
    //set _re_init to false once K matrix is built
    _re_init = false;
  }
  
  // assemble the rhs vector, and record the CPU wall time.
  // t1 = MPI_Wtime();
  this->assemble_rhs(this->name(), option);
  // t2 = MPI_Wtime();
  // std::cout << "For Poisson equation, time used to assemble the
  // right-hand-side vector is " <<t2-t1<<" s\n";

  // solve the problem
  _solver_poisson.solve(this->name());
  

  STOP_LOG("solve()", "PMSystemPoisson");
}

// ==================================================================================
//void PMSystemPoisson::compute_point_potential(std::vector<Real>& pv)
//{
//  START_LOG("compute_point_potential()", "PMSystemPoisson");
//
//  // NP: size of the potential vector
//  const MeshBase  & mesh   = this->get_mesh();
//  const std::size_t NP     =  _point_mesh->num_particles();
//  const std::size_t dim    = mesh.mesh_dimension();
//  const dof_id_type& n_elem = mesh.n_elem();
//
//  std::vector<Real> _pv_send_list;                               // point
//                                                                 // potential
//                                                                 // send list
//  std::vector<Real> _pglobal_send_list;
//  std::vector<Real> _plocal_send_list;
//  std::vector<std::size_t> _pid_send_list;                       // point id
//                                                                 // send list
//  const unsigned int phi_var     = this->variable_number("phi"); // phi_var = 0
//  const std::string  charge_type = "regularized";
//
//  // Loop over each point, and compute electrical potential at its location.
//  // Collect the global potential from FEM through Allgather operation
//  for (std::size_t i = 0; i < NP; ++i)
//  {
//    // 0. point coordinates & residing element
//    const Point pt               = _point_mesh->particles()[i]->point();
//    const PointType   point_type = _point_mesh->particles()[i]->point_type();
//    const dof_id_type elem_id    = _point_mesh->particles()[i]->elem_id();
//
//    if (elem_id > n_elem)
//    {
//      std::cout << "---> Error in PMSystemPoisson::compute_point_potential():\n"
//                << "      " << i << "-th particle (position = " << pt(0) <<
//        ", " << pt(1) << ", " << pt(2) << ") is out of domain"
//                << std::endl << std::endl;
//    }
//    const Elem *elem = mesh.elem(elem_id);
//
//    // 1. Global(FEM) solution at the current point. This is done on local
//    // processors
//    Real phi_global;
//
//    // 2. Local potential at particle i, This is also done on local processors
//    Real phi_local;
//
//    // 3. total potential
//    Real phi_total;
//
//    if (elem->processor_id() == this->processor_id())
//    {
//      // get global potential
//      phi_global = this->point_value(phi_var, pt, *elem);
//
//      // get local velocity
//      phi_local = this->local_potential_bead(i, charge_type);
//
//      // FIXME: do we need exclude self-exclusion term for point_type =
//      // "POLYMER_BEAD"?
//      if (point_type == POLYMER_BEAD) {
//        phi_total = phi_global + phi_local;
//      }
//      else if (point_type == LAGRANGIAN_POINT) {
//        std::cout << "---> Error in PMSystemPoisson::compute_point_potential: \n"
//                  <<
//          "     point_type LAGRANGIAN_POINT is not supported yet for electrostatics."
//                  << std::endl << std::endl;
//
//        libmesh_error();
//      }
//
//      // pack particle id and its velocity
//      _pid_send_list.push_back(i);
//      _pv_send_list.push_back(phi_total);
//
//      // printf("i = %i, processor_id = %i, phi_local = %f\n", i,
//      // this->processor_id(), phi_local);
//    } // end if (elem->processor_id() == this->processor_id)
//  }   // end for i-loop
//
//  // Check the size of local_pv and the size of list on each process after
//  // allgather
//  this->comm().allgather(_pid_send_list); // allgather the particle id
//  this->comm().allgather(_pv_send_list); // allgather the electrical potential
//
//  if (_pid_send_list.size() != NP)
//  {
//    std::cout << "---> Error in PMSystemPoisson::compute_point_potential: \n"
//              << "       _pid_send_list.size() != NP"
//              << std::endl << std::endl;
//    libmesh_error();
//  }
//
//  for (std::size_t i = 0; i < NP; ++i) {
//    const std::size_t p_id = _pid_send_list[i];
//    pv[p_id] = _pv_send_list[i];
//
//    // ---------------------------- output for debug
//    // -----------------------------
//    // if (this->comm().rank()==0)
//    // {
//    //   printf("\n--->test in compute_point_potential(): output electrical
//    // potential at the point:\n");
//    //   printf("point %lu: phi_global (FEM) = (%E)\n", p_id,
//    // _pglobal_send_list[i] );
//    //   printf("              phi_local (Green Function) = (%E)\n",
//    // _plocal_send_list[i] );
//    //   printf("           ---phi_total  = (%E)\n\n", pv[pid] );
//    // }
//  }
//
//  STOP_LOG("compute_point_potential()", "PMSystemPoisson");
//}

// ==================================================================================
void PMSystemPoisson::compute_point_efield(std::vector<Point>& pv)
{
  START_LOG("compute_point_efield()", "PMSystemPoisson");

  // dim*NP: size of the electric field vector
  const MeshBase  & mesh   = this->get_mesh();
  const std::size_t NP     =  _point_mesh->num_particles();
  const std::size_t dim    = mesh.mesh_dimension();
  const dof_id_type& n_elem = mesh.n_elem();
  const unsigned int phi_var     = this->variable_number("phi"); // phi_var = 0
  const std::string  charge_type = "regularized";

  // create local buffer of total efield and
  // particle id. These vectors will be evaluated in each processor locally
  // and then gathered together using this->comm().allgather function
  std::vector<Point> efield_total_send_list;
  std::vector<size_t> pid_send_list;

  // Loop over each point, and compute electric field in its location
  // Collect the global electrical potential from FEM through Allgather
  // operation
  for (std::size_t i = 0; i < NP; ++i)
  {
    // 0. point coordinates & its residing element
    const Point pt               = _point_mesh->particles()[i]->point();
    const PointType   point_type = _point_mesh->particles()[i]->point_type();
    const dof_id_type elem_id    = _point_mesh->particles()[i]->elem_id();
    if (elem_id > n_elem)
    {
      std::cout << "---> Error in PMSystemPoisson::compute_point_efield():\n"
                << "      " << i << "-th particle (position = " << pt(0) <<
        ", " << pt(1) << ", " << pt(2) << ") is out of domain"
                << std::endl << std::endl;
    }
    const Elem *elem = mesh.elem(elem_id);

    // global potential gradient
    Gradient efield_global;
    // local solution, but we only need to compute the gradient
    std::pair<Real, Point> sol_local;
    const std::string sol_type = "grad";
    // 3. total electric field
    Point efield_total;

    if (elem->processor_id() == this->processor_id())
    {
      // Global electric field
      efield_global = this->point_gradient(phi_var, pt, *elem);
      // evaluate local solution
      this->local_potential_field(pt, charge_type, sol_type, sol_local,
        elem_id);
      // FIXME: do we need to exclude self exclusion term for point_type =
      // "POLYMER_BEAD"
      if (point_type == POLYMER_BEAD) {
        efield_total(0) = efield_global(0) + sol_local.second(0);
        efield_total(1) = efield_global(1) + sol_local.second(1);
        efield_total(2) = efield_global(2) + sol_local.second(2);
      }
      else if (point_type == LAGRANGIAN_POINT) {
        std::cout << "---> Error in PMSystemPoisson::compute_point_potential: \n"
                  <<"     point_type LAGRANGIAN_POINT is not supported yet for electrostatics."
                  << std::endl << std::endl;
        libmesh_error();
      }
      // pack particle id
      pid_send_list.push_back(i);
      efield_total_send_list.push_back(efield_total);
    }// end if (elem->processor_id() == this->processor_id)
  } // end for i-loop

  // Check the size of local_pv and the size of list on each process after
  // allgather
  this->comm().allgather(pid_send_list);
  this->comm().allgather(efield_total_send_list);
  if (pid_send_list.size() != NP)
  {
    std::cout << "---> Error in PMSystemPoisson::compute_point_efield: \n"
              << "       pid_send_list.size() != NP"
              << std::endl << std::endl;
    libmesh_error();
  }

  // Fill pv
  for (std::size_t i = 0; i < NP; ++i)
    pv[pid_send_list[i]] = efield_total_send_list[pid_send_list[i]];

  STOP_LOG("compute_point_efield()", "PMSystemPoisson");
}

// ==================================================================================
void PMSystemPoisson::add_electrostatic_forces()
{
  START_LOG("add_electrostatic_forces()", "PMSystemPoisson");

  const MeshBase  & mesh = this->get_mesh();
  const std::size_t NP   = _point_mesh->num_particles();
  const std::size_t dim  = mesh.mesh_dimension();

  // Compute electric field at every bead's location
  std::vector<Point> pv(NP, 0.);
  this->compute_point_efield(pv);

  Point pforce;
  Real  charge;
  for (std::size_t i = 0; i < NP; ++i) {
    charge = _point_mesh->particles()[i]->charge();
    pforce = pv[i] * charge;
    _point_mesh->particles()[i]->add_particle_force(pforce);
  }

  STOP_LOG("add_electrostatic_forces()", "PMSystemPoisson");
}

// ==========================================================================
Real PMSystemPoisson::local_potential_field(const Point      & p,
                                            const std::string& charge_type,
                                            dof_id_type p_elem_id)
                                            const
{
  START_LOG("local_potential_field()", "PMSystemPoisson");

  // locate point element id if it's not given
  if (p_elem_id==-1)
  {
    const MeshBase& mesh = this->get_mesh();
    p_elem_id = mesh.point_locator().operator()(p)->id();
  }

  Real phi_local =
    ggem_poisson->local_solution_field(_point_mesh, p, charge_type, p_elem_id);

  STOP_LOG("local_potential_field()", "PMSystemPoisson");
  return phi_local;
}

// ===========================================================================
void PMSystemPoisson::local_potential_field(const Point &p,
                                            const std::string &charge_type,
                                            const std::string& sol_option,
                                            std::pair<Real, Point>& local_sol,
                                            dof_id_type p_elem_id)const
{
  START_LOG("local_potential_field()", "PMSystemPoisson");

  // locate point element id if it's not given
  if (p_elem_id==-1)
  {
    const MeshBase& mesh = this->get_mesh();
    p_elem_id = mesh.point_locator().operator()(p)->id();
  }

  ggem_poisson->local_solution_field(_point_mesh, p, charge_type, p_elem_id,
    sol_option, local_sol);

  STOP_LOG("local_potential_field()", "PMSystemPoisson");
}


// ===========================================================================
Real PMSystemPoisson::total_potential_laplacian_field(const Point &pt,
                                                      const std::string&charge_type,
                                                      dof_id_type elem_id) const
{
  START_LOG("total_potential_laplacian_field()", "PMSystemPoisson");

  // locate point element id if it's not given
  if (elem_id==-1)
  {
    const MeshBase& mesh = this->get_mesh();
    elem_id = mesh.point_locator().operator()(pt)->id();
  }

  const std::vector<dof_id_type>& point_nb_list =
    _point_mesh->get_elem_point_neighbor_list(elem_id);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Loop over all the neighbor list beads, and
     and compute the laplacian, which is identical to -4*pi*charge_density
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  // Initialize variables
  Real laplacian = 0., rho=0.;

  // loop over all neighbor points of this field point
  for (std::size_t v = 0; v < point_nb_list.size(); ++v)
  {
    // 0. particle id and position, vector x = ptx - pt0
    const dof_id_type& p_id = point_nb_list[v];
    const Point& ptj = _point_mesh->particles()[p_id]->point();
    const Point& x = _point_mesh->pm_periodic_boundary()->point_vector(pt, ptj);

    // 1. compute the Green function of particle-v
    if (charge_type == "regularized") {
        rho += ggem_poisson->regularized_charge_exp(x) *
          (_point_mesh->particles()[p_id]->charge());
    }
    else {
      std::cout<<"Error: invalid charge_type. Exiting..."<<std::endl;
      libmesh_error();
    } // end if-else
  }

  laplacian = -4 * pi * rho;

  STOP_LOG("total_potential_laplacian_field()", "PMSystemPoisson");

  return laplacian;
}

// ===========================================================================
Real PMSystemPoisson::local_potential_bead(const std::size_t& bead_id,
                                           const std::string& charge_type) const
{
  START_LOG("local_potential_bead()", "PMSystemPoisson");

  Real phi_local = ggem_poisson->local_solution_bead(_point_mesh, bead_id,
    charge_type);

  STOP_LOG("local_potential_bead()", "PMSystemPoisson");
  return phi_local;
}

// ===========================================================================
void PMSystemPoisson::local_potential_bead(const std::size_t& bead_id,
                                           const std::string& charge_type,
                                           const std::string& sol_type,
                                           std::pair<Real, Point>& local_sol)
                                           const
{
  START_LOG("local_potential_bead()", "PMSystemPoisson");

  ggem_poisson->local_solution_bead(_point_mesh, bead_id, charge_type,
    sol_type, local_sol);

  STOP_LOG("local_potential_bead()", "PMSystemPoisson");
}

// ===========================================================================
void PMSystemPoisson::test_potential_profile()
{
  START_LOG("test_potential_profile()", "PMSystemPoisson");

  PMToolBox::output_message("========>2. Test in PMSystemPoisson::"
                            "test_potential_profile(): \n", this->comm());

  // Solve the electrical potential field: global solution by FEM
  _re_init = true; // assemble global matrix and init ksp_solver
  this->solve("unused");

  // get reference to system essentials
  MeshBase& mesh = this->get_mesh();
  const std::size_t& dim = mesh.mesh_dimension();
  const unsigned int& phi_var = this->variable_number("phi"); // phi_var = 0

  // Output electrical potential profiles along xyz directions, global + local
  // solutions
  const Point& box_min = _point_mesh->pm_periodic_boundary()->box_min();
  const Point& box_len = _point_mesh->pm_periodic_boundary()->box_length();
  // nn: number of spacial points to check
  const unsigned int ns = 200;
  // NP: number of point particles (beads) in the system
  const unsigned int NP = _point_mesh->num_particles();

  // define output file format parameters
  std::ofstream outfile;
  o_precision = this->get_equation_systems().parameters.get<int>
          ("o_precision");
  // define output filename for x, y, z directions
  std::vector<std::string> filenames;
  filenames.push_back("x");
  filenames.push_back("y");
  if (dim==3) filenames.push_back("z");

  // loop over all [dim_i] directions and write the output
  for (int dim_i=0; dim_i<dim; dim_i++)
  {
    // create output file
    std::ostringstream oss;
    oss << "output_potential_profile_" << filenames[dim_i] << "_" << NP <<
    "P.csv";
    // open output file
    PMToolBox::output_message("writing output to " + oss.str() + "\n",
            this->comm());
    outfile.open(oss.str(), std::ios_base::out);
    // set output file format (fixed + precision)
    outfile.setf(std::ios::fixed);
    outfile.precision(o_precision);
    // write column header
    if (this->comm().rank() == 0)
      outfile <<filenames[dim_i]<<","
        <<"phi_total,phi_global,phi_local,phi_exact,"
        <<"phi_grad_total(x),phi_grad_total(y),phi_grad_total(z),"
        <<"phi_grad_global(x),phi_grad_global(y),phi_grad_global(z),"
        <<"phi_grad_local(x),phi_grad_local(y),phi_grad_local(z),"
        <<"phi_grad_exact(x),phi_grad_exact(y),phi_grad_exact(z),"
        <<"phi_laplacian_total,phi_laplacian_exact"<<"\n";

    // separation between spacial points in this direction
    const Real ds = box_len(dim_i) / Real(ns);
    // create helper vector
    std::vector<int> helper_vec(dim, 0);
    helper_vec[dim_i] = 1;
    // loop over all [ns] space points along this direction, the coordinate of
    // other two directions are set to be zero
    for (std::size_t i = 0; i < ns + 1; ++i)
    {
      // create the i-th spacial point
      Point pt((box_min(0) + Real(i) * ds) * helper_vec[0],
               (box_min(1) + Real(i) * ds) * helper_vec[1],
               (box_min(2) + Real(i) * ds) * helper_vec[2]);
      // Initialize solution
      Real phi_global=0., phi_total=0., phi_exact;
      // Initialize solution gradient
      Point phi_grad_global, phi_grad_total, phi_grad_exact;

      // evaluate global solution from FEM solution, this is slow, but we
      // tolerate it only for test
      phi_global = this->point_value(phi_var, pt);
      phi_grad_global = this->point_gradient(phi_var, pt);

      // evaluate local potential from GGEM
      std::pair<Real, Point> local_sol;
      this->local_potential_field(pt, "regularized", "phi&grad", local_sol);

      // get total solution at this point
      phi_total = phi_global + local_sol.first;
      phi_grad_total = phi_grad_global + local_sol.second;

      // get exact solution for an unbounded domain from analytical solution
      phi_exact = analytical_solution->exact_solution_infinite_domain(
        *ggem_poisson, pt);
      phi_grad_exact =
        analytical_solution->exact_solution_infinite_domain_grad
        (*ggem_poisson, pt);

      // get laplacian
      Real phi_exact_laplacian =
        analytical_solution->exact_solution_infinite_domain_laplacian
        (*ggem_poisson, pt);
      Real phi_total_laplacian = this->total_potential_laplacian_field(pt,
        "regularized");

      // write the result to output file
      if (this->comm().rank() == 0)
        outfile <<pt(dim_i)<<","
        <<phi_total<<","
        <<phi_global<<","
        <<local_sol.first<<","
        <<phi_exact<<","
        <<phi_grad_total(0)<<","<<phi_grad_total(1)<<","<<phi_grad_total(2)<<","
        <<phi_grad_global(0)<<","<<phi_grad_global(1)<<","<<phi_grad_global(2)<<","
        <<local_sol.second(0)<<","<<local_sol.second(1)<<","<<local_sol.second(2)<<","
        <<phi_grad_exact(0)<<","<<phi_grad_exact(1)<<","<<phi_grad_exact(2)<<","
        <<phi_total_laplacian<<","<<phi_exact_laplacian
        <<"\n";
    } // end for i-loop
    outfile.close();
  } // end loop dim_i

  // write equation system to output file

#ifdef LIBMESH_HAVE_EXODUS_API
  PMToolBox::output_message("writing global solution to global_solution.e",
    this->comm());
  std::string output_filename = "global_solution.e";
  ExodusII_IO(this->get_mesh()).write_equation_systems(output_filename,
                                                       this->get_equation_systems());
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // write points value (for test only)
//  const int n_pts = 200;
//  std::vector<Point> pts(n_pts);
//  for (int i=0; i<pts.size(); i++)
//    pts[i]=Point(box_min(0) + i*box_len(0)/Real(n_pts), 0., 0.);
//  this->output_point_solution(pts, "test_point_solution.csv");

//  PMToolBox::output_message("writing nodal global solution to "
//                            "poisson_global_solution.out", this->comm());
//  this->output_nodal_solution("poisson_global_solution.out");

  STOP_LOG("test_potential_profile()", "PMSystemPoisson");
}

// =========================================================================
void PMSystemPoisson::test_nodal_error()
{
  START_LOG("test_nodal_error()", "PMSystemPoisson");
  std::ostringstream ss;
  ss << "===> test in PMSystemPoisson::test_nodal_error(): \n";
  PMToolBox::output_message(ss, this->comm());

  // neighbor list is already updated before calling this function

  // solve the disturbed test system to get global solution
  _re_init = true;
  this->solve("disturbed");

  //update system solution to total
  this->update_solution_to_total();

  // Theoretical solution
  const unsigned int dim = 3;
  MeshBase& mesh = this->get_mesh();
  const unsigned int n_nodes = mesh.n_nodes();
  Real max_abs_error=0.;
//  Real max_abs_error_grad=0.;
  Real val0_norm = 0.;
  Real val1_norm = 0.;
  Real val2_norm = 0.;
  Real val3_norm = 0.;

  const unsigned int phi_var = this->variable_number("phi");

  // Loop over each node and compute the nodal solution
  MeshBase::node_iterator nd           = mesh.local_nodes_begin();
  const MeshBase::node_iterator end_nd = mesh.local_nodes_end();

  for (; nd != end_nd; ++nd)
  {
    // Store a pointer to the current node, and extract a point
    Node *node = *nd;
    Point pt;

    for (unsigned int i = 0; i < dim; ++i) pt(i) =  (*node)(i);

    // get the dof index at this node
    const dof_id_type& dof_index = node->dof_number(this->number(), phi_var, 0);

    // get the element id of this node from stored mapping
    const dof_id_type &elem_id = _point_mesh->get_node_elem_id(node->id());
    const Elem* elem = mesh.elem_ptr(elem_id);

    const Real& phi_total = this->solution->operator()(dof_index);
//    const RealGradient& phi_total_grad = this->point_gradient(phi_var, pt,
//      *elem);

    // compute the local velocity of fluid at the current node
    Real phi_exact =
      analytical_solution->exact_solution_infinite_domain(*ggem_poisson, pt);
//    Point phi_exact_grad =
//      analytical_solution->exact_solution_infinite_domain_grad(*ggem_poisson,
//        pt);

    // compute the errors
    max_abs_error = std::max(std::abs(phi_total-phi_exact), max_abs_error);
//    max_abs_error_grad =
//      std::max({std::abs(phi_total_grad(0)-phi_exact_grad(0)),std::abs
//      (phi_total_grad(1)-phi_exact_grad(1)), std::abs(phi_total_grad(2)
//      -phi_exact_grad(2)), max_abs_error_grad});
    Real tmpt = std::abs(phi_total-phi_exact);
    val0_norm += tmpt;
    val1_norm += std::abs(phi_exact);
    val2_norm += tmpt * tmpt;
    val3_norm += phi_exact * phi_exact;
  } // end for nd-loop

  // Compute the error: l1 and l2-norm of errors
  this->comm().sum(val0_norm);
  this->comm().sum(val1_norm);
  this->comm().sum(val2_norm);
  this->comm().sum(val3_norm);
  this->comm().max(max_abs_error);
//  this->comm().max(max_abs_error_grad);

  const Real l1_norm = val0_norm / val1_norm;
  const Real l2_norm = std::sqrt(val2_norm) / std::sqrt(val3_norm);
  ss << ">>> sum(|phi_total-phi_exact|) / sum(|phi_exact|) = " << l1_norm
     << "; sum(|phi_total-phi_exact|^2) / sum(|phi_exact|^2) = " << l2_norm
     << "; max(|phi_total-phi_exact|) = "<<max_abs_error
//     <<"; max(phi_grad_total-phi_grad_exact|) = "<<max_abs_error_grad
     <<"\n";
  PMToolBox::output_message(ss, this->comm());

  // write equation system to output file
#ifdef LIBMESH_HAVE_EXODUS_API
  PMToolBox::output_message("writing total solution to total_solution.e",
this->comm());
  std::string output_filename = "total_solution.e";
ExodusII_IO(this->get_mesh()).write_equation_systems(output_filename,
                                               this->get_equation_systems());
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // output nodal solution
//  PMToolBox::output_message("writing nodal total solution to "
//                            "poisson_total_solution.out", this->comm());
//  this->output_nodal_solution("poisson_total_solution.out");

  // resume solution to global
  this->resume_solution_to_global();

  STOP_LOG("test_nodal_error()", "PMSystemPoisson");
}

// ===========================================================================
void PMSystemPoisson::output_nodal_solution(const std::string& output_filename)
{
  START_LOG("output_nodal_solution()", "PMSystemPoisson");
  std::filebuf fb;
  fb.open (output_filename,std::ios::out);
  std::ostream os(&fb);

  std::vector<Real> v(this->solution->size());
  this->solution->localize(v);
  // right now we only want one copy of the output
  if (this->solution->processor_id())
    return;

  // get a reference to mesh
  const MeshBase& mesh = this->get_mesh();
  // loop over all nodes
  os <<"node_dof_id,x,y,z,phi\n";
  for(dof_id_type i=0; i!=v.size(); i++){
    const Node* node = mesh.node_ptr(i);
    os<<i<<","
    <<node->operator()(0)<<","<<node->operator()(1)<<","<<node->operator()(2)
    <<","<<v[i]<<"\n";
  }

  fb.close();
  STOP_LOG("output_nodal_solution()", "PMSystemPoisson");
}

// ===========================================================================

// ===========================================================================
void PMSystemPoisson::output_point_solution(const std::vector<Point>& pts,
                          const std::string& filename)
{
  START_LOG("output_point_solution()", "PMSystemPoisson");

  // get system dimension
  MeshBase& mesh = this->get_mesh();
  const std::size_t dim = mesh.mesh_dimension();

  // define output file format parameters
  o_precision = this->get_equation_systems().parameters.get<int>
    ("o_precision");
  std::ofstream outfile;
  outfile.open(filename, std::ios_base::out);
  // set output file format (fixed + precision)
  outfile.setf(std::ios::fixed);
  outfile.precision(o_precision);

  if (this->comm().rank() == 0)
    outfile << "x,y,z,phi_local,phi_global,phi_total\n";

  for (std::size_t i = 0; i < pts.size(); i++)
  {
    // global potential from FEM
    Real phi_global, phi_total;
    const unsigned int phi_var = this->variable_number("phi"); // phi_var = 0
    // get global phi from FEM solution, this is slow, but we tolerate it
    // only for test
    phi_global = this->point_value(phi_var, pts[i]);

    // get local potential from analytical function
    const Real phi_local = this->local_potential_field(pts[i], "regularized");

    // get total solution at this point
    phi_total = phi_global + phi_local;

    // write the result to output file
    if (this->comm().rank() == 0)
      outfile << pts[i](0)<<","<<pts[i](1)<<","<<pts[i](2)<<","
              << phi_local<<","<<phi_global<<","<<phi_total<<"\n";
  } // end for i-loop
  outfile.close();

  STOP_LOG("output_point_solution()", "PMSystemPoisson");
}

// ===========================================================================
void PMSystemPoisson::update_solution_to_total()
{
  START_LOG("update_solution_to_total()", "PMSystemPoisson");
  // Check if the system solution vector is closed or not
  if ((this->solution->closed()) == false) this->solution->close();

  // clone solution to global solution before make changes
  this->_global_solution = this->solution->clone();

  // Add local_disturbed_solution for each nodes

  // Get the parameters and Initialize the quantities
  MeshBase& mesh = this->get_mesh();
  const unsigned int& dim = mesh.mesh_dimension();

  // loop over local nodes and update local solution vector according to the
  // dof_indices of each node
  MeshBase::node_iterator nd           = mesh.local_nodes_begin();
  const MeshBase::node_iterator end_nd = mesh.local_nodes_end();

  for (; nd != end_nd; ++nd) {
    // Store a pointer to the current node, and extract a point
    Node *node = *nd;
    Point pt;
    for (unsigned int i = 0; i < dim; ++i) pt(i) = (*node)(i);

    // get the element id of this node from stored mapping
    const dof_id_type &elem_id = _point_mesh->get_node_elem_id(node->id());

    // get the dof index at this node
    const unsigned int phi_var = this->variable_number("phi");
    const dof_id_type& dof_index = node->dof_number(this->number(), phi_var, 0);

    // get the local solution on this node calling GGEM
    std::pair<Real, Point> sol_local;
    this->local_potential_field(pt, "regularized", "phi&grad", sol_local, elem_id);

    // add local solution to system solution
    this->solution->add(dof_index, sol_local.first);
  } // end loop over local nodes

  this->solution->close();

  // update the local values in current_local_solution to reflect the
  // solution on neighboring processors since we make changes to this->solution
  this->update();

  STOP_LOG("update_solution_to_total()", "PMSystemPoisson");
}

// ===========================================================================
void PMSystemPoisson::resume_solution_to_global()
{
  START_LOG("resume_solution_to_global()", "PMSystemPoisson");

  *(this->solution) = *(this->_global_solution);
  this->update();

  STOP_LOG("resume_solution_to_global()", "PMSystemPoisson");
}


} // end of namespace
