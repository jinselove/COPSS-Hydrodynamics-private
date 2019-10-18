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
  _solver_poisson.solve();
  

  STOP_LOG("solve()", "PMSystemPoisson");
}

// ==================================================================================
void PMSystemPoisson::add_local_solution()
{
  START_LOG("add_local_solution()", "PMSystemPoisson");

  // Check if the system solution vector is closed or not
  if ((this->solution->closed()) == false) this->solution->close();

  // Get parameters and initialize the quantities
  MeshBase& mesh                    = this->get_mesh();
  const std::size_t   dim           = mesh.mesh_dimension();
  const std::size_t   n_local_nodes = mesh.n_local_nodes();
  std::vector<Number> local_solution(n_local_nodes);
  std::vector<numeric_index_type> dof_indices(n_local_nodes);

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
    // Store a pointer to the current node, and extract a point coordinate
    Node *node = *nd;
    Point pt;

    for (unsigned int i = 0; i < dim; ++i) pt(i) = (*node)(i);

    // this is a test for dof_number at each node
    // const unsigned int node_id = node->id();
    // std::ostringstream oss;
    // oss << "          NODE " << node_id;
    // PMToolBox::output_message(oss.str(), this->comm());
    // node->print_info();
    // if(this->comm().rank()==0) printf("--->test: nodal dof number :");
    // for(unsigned int i=0; i<dim; ++i)
    // {
    //  dof_id_type dof_num = node->dof_number(this->number(), i, 0);
    //  if(this->comm().rank()==0) printf(" %u", dof_num);
    // }
    // if(this->comm().rank()==0) printf(" \n");

    // Store local electrical potential and dof indices
    local_solution[local_count] = this->local_potential_field(pt, "regularized");
    dof_indices[local_count]    = node->dof_number(this->number(), 0, 0);

    local_count++;
  }

  // printf("--->test in add_local_solution() local_count = %lu on the processor
  // %u\n",
  //       n_local_nodes,this->comm().rank());

  // add local potential to the global potential
  this->solution->add_vector(local_solution, dof_indices);
  this->solution->close();
  this->update();

  STOP_LOG("add_local_solution()", "PMSystemPoisson");
}

// ==================================================================================
void PMSystemPoisson::test_l2_norm(bool& neighbor_list_update_flag)
{
  START_LOG("test_l2_norm()", "PMSystemPoisson");
  std::string msg = "--->test in PMSystemPoisson::test_l2_norm(): \n";
  PMToolBox::output_message(msg, this->comm());

  // FIXME: to be implemented

  STOP_LOG("test_l2_norm()", "PMSystemPoisson");
}


// ==================================================================================
void PMSystemPoisson::compute_point_potential(std::vector<Real>& pv)
{
  START_LOG("compute_point_potential()", "PMSystemPoisson");

  // NP: size of the potential vector
  const MeshBase  & mesh   = this->get_mesh();
  const std::size_t NP     =  _point_mesh->num_particles();
  const std::size_t dim    = mesh.mesh_dimension();
  const dof_id_type n_elem = mesh.n_elem();

  std::vector<Real> _pv_send_list;                               // point
                                                                 // potential
                                                                 // send list
  std::vector<Real> _pglobal_send_list;
  std::vector<Real> _plocal_send_list;
  std::vector<std::size_t> _pid_send_list;                       // point id
                                                                 // send list
  const unsigned int phi_var     = this->variable_number("phi"); // phi_var = 0
  const std::string  charge_type = "regularized";

  // Loop over each point, and compute electrical potential at its location.
  // Collect the global potential from FEM through Allgather operation
  for (std::size_t i = 0; i < NP; ++i)
  {
    // 0. point coordinates & residing element
    const Point pt               = _point_mesh->particles()[i]->point();
    const PointType   point_type = _point_mesh->particles()[i]->point_type();
    const dof_id_type elem_id    = _point_mesh->particles()[i]->elem_id();

    if (elem_id > n_elem)
    {
      std::cout << "---> Error in PMSystemPoisson::compute_point_potential():\n"
                << "      " << i << "-th particle (position = " << pt(0) <<
        ", " << pt(1) << ", " << pt(2) << ") is out of domain"
                << std::endl << std::endl;
    }
    const Elem *elem = mesh.elem(elem_id);

    // 1. Global(FEM) solution at the current point. This is done on local
    // processors
    Real phi_global;

    // 2. Local potential at particle i, This is also done on local processors
    Real phi_local;

    // 3. total potential
    Real phi_total;

    if (elem->processor_id() == this->processor_id())
    {
      // get global potential
      phi_global = this->point_value(phi_var, pt, *elem);

      // get local velocity
      phi_local = this->local_potential_bead(i, charge_type);

      // FIXME: do we need exclude self-exclusion term for point_type =
      // "POLYMER_BEAD"?
      if (point_type == POLYMER_BEAD) {
        phi_total = phi_global + phi_local;
      }
      else if (point_type == LAGRANGIAN_POINT) {
        std::cout << "---> Error in PMSystemPoisson::compute_point_potential: \n"
                  <<
          "     point_type LAGRANGIAN_POINT is not supported yet for electrostatics."
                  << std::endl << std::endl;

        libmesh_error();
      }

      // pack particle id and its velocity
      _pid_send_list.push_back(i);
      _pv_send_list.push_back(phi_total);

      // printf("i = %i, processor_id = %i, phi_local = %f\n", i,
      // this->processor_id(), phi_local);
    } // end if (elem->processor_id() == this->processor_id)
  }   // end for i-loop

  // Check the size of local_pv and the size of list on each process after
  // allgather
  this->comm().allgather(_pid_send_list); // allgather the particle id
  this->comm().allgather(_pv_send_list); // allgather the electrical potential

  if (_pid_send_list.size() != NP)
  {
    std::cout << "---> Error in PMSystemPoisson::compute_point_potential: \n"
              << "       _pid_send_list.size() != NP"
              << std::endl << std::endl;
    libmesh_error();
  }

  for (std::size_t i = 0; i < NP; ++i) {
    const std::size_t p_id = _pid_send_list[i];
    pv[p_id] = _pv_send_list[i];

    // ---------------------------- output for debug
    // -----------------------------
    // if (this->comm().rank()==0)
    // {
    //   printf("\n--->test in compute_point_potential(): output electrical
    // potential at the point:\n");
    //   printf("point %lu: phi_global (FEM) = (%E)\n", p_id,
    // _pglobal_send_list[i] );
    //   printf("              phi_local (Green Function) = (%E)\n",
    // _plocal_send_list[i] );
    //   printf("           ---phi_total  = (%E)\n\n", pv[pid] );
    // }
  }

  STOP_LOG("compute_point_potential()", "PMSystemPoisson");
}

// ==================================================================================
void PMSystemPoisson::compute_point_efield(std::vector<Real>& pv)
{
  START_LOG("compute_point_efield()", "PMSystemPoisson");

  // dim*NP: size of the electric field vector
  const MeshBase  & mesh   = this->get_mesh();
  const std::size_t NP     =  _point_mesh->num_particles();
  const std::size_t dim    = mesh.mesh_dimension();
  const dof_id_type n_elem = mesh.n_elem();

  // std::vector<Real> pvlocal(dim*NP,0.);  // declared on each processor
  std::vector<Real> _pv_send_list;                               // point
                                                                 // electric
                                                                 // field send
                                                                 // list
  std::vector<Real> _pglobal_send_list;
  std::vector<Real> _plocal_send_list;
  std::vector<std::size_t> _pid_send_list;                       // point id
                                                                 // send list
  const unsigned int phi_var     = this->variable_number("phi"); // phi_var = 0
  const std::string  charge_type = "regularized";

  // Get Finite Element type for "phi"
  FEType fe_phi_type = this->variable_type(phi_var);

  // Build a Finite Element object of the specified type for the potential
  // variable
  UniquePtr<FEBase> fe_phi(FEBase::build(dim, fe_phi_type));

  // Gauss quadrature rule for numerical integration
  QGauss qrule(dim, SECOND);

  // Tell finite element objects to use the quadrature rule
  fe_phi->attach_quadrature_rule(&qrule);

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

    // 1. Global(FEM) solution at the current point. This is done on local
    // processors
    Gradient efield_global;

    // 2. Local electric field at particle i, This is also done on local
    // processors
    RealGradient efield_local;

    // 3. total electric field
    Point efield_total;

    if (elem->processor_id() == this->processor_id())
    {
      // Global electric field
      efield_global = this->point_gradient(phi_var, pt, *elem);

      // Extract shape function derivatives to be evaluated at the bead location
      const std::vector<std::vector<RealGradient> >& dphi = fe_phi->get_dphi();

      // Extract bead coordinates in the reference frame
      std::vector<Point> bead_ptx;
      bead_ptx.resize(1);
      bead_ptx[0] = FE<3, LAGRANGE>::inverse_map(elem, pt);

      // Evaluate the shape functions derivatives at this bead location
      fe_phi->reinit(elem, &bead_ptx);

      // Local electrical potential on all nodes in this element
      std::vector<Real> phi_local;
      phi_local.resize(elem->n_nodes());

      // Loop through nodes on this element
      for (unsigned int nn = 0; nn < elem->n_nodes(); nn++)
      {
        // Coordinate of the node
        const Point ptx = elem->point(nn);

        // Evaluate local electrical potential (phi) on each node
        phi_local[nn] = this->local_potential_field(elem, ptx, "regularized");

        // phi_local[nn] = this->local_potential_field(ptx, "regularized");
      }

      // Interpolate local electric field on this bead location
      for (unsigned int j = 0; j < phi_local.size();
           j++) efield_local += phi_local[j] * dphi[j][0];

      // FIXME: do we need to exclude self exclusion term for point_type =
      // "POLYMER_BEAD"
      if (point_type == POLYMER_BEAD) {
        efield_total(0) = efield_global(0) + efield_local(0);
        efield_total(1) = efield_global(1) + efield_local(1);
        efield_total(2) = efield_global(2) + efield_local(2);
      }
      else if (point_type == LAGRANGIAN_POINT) {
        std::cout << "---> Error in PMSystemPoisson::compute_point_potential: \n"
                  <<
          "     point_type LAGRANGIAN_POINT is not supported yet for electrostatics."
                  << std::endl << std::endl;

        libmesh_error();
      }

      // pack particle id and its velocity
      _pid_send_list.push_back(i);

      for (unsigned int j = 0; j < dim; j++) {
        _pv_send_list.push_back(efield_total(j));

        // _pglobal_send_list.push_back(efield_global(j));
        // _plocal_send_list.push_back(efield_local(j));
      }

      // printf("i = %i, processor_id = %i, efield_local = %f, %f, %f\n", i,
      // this->processor_id(), efield_local(0), efield_local(1),
      // efield_local(2));
    } // end if (elem->processor_id() == this->processor_id)
  }   // end for i-loop

  // Check the size of local_pv and the size of list on each process after
  // allgather
  this->comm().allgather(_pid_send_list); // allgather the particle id
  this->comm().allgather(_pv_send_list); // allgather the total electric field

  // this->comm().allgather(_pglobal_send_list);
  // this->comm().allgather(_plocal_send_list);
  if (_pid_send_list.size() != NP)
  {
    std::cout << "---> Error in PMSystemPoisson::compute_point_efield: \n"
              << "       _pid_send_list.size() != NP"
              << std::endl << std::endl;
    libmesh_error();
  }

  for (std::size_t i = 0; i < NP; ++i) {
    const std::size_t p_id = _pid_send_list[i];

    for (unsigned int j = 0; j < dim; ++j) {
      pv[dim * p_id + j] = _pv_send_list[i * dim + j];
    }
  }

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
  std::vector<Real> pv(NP * dim, 0.);
  this->compute_point_efield(pv);

  Point pforce;
  Real  charge;

  for (std::size_t i = 0; i < NP; ++i) {
    charge = _point_mesh->particles()[i]->charge();

    pforce(0) = charge * pv[dim * i];
    pforce(1) = charge * pv[dim * i + 1];
    pforce(2) = charge * pv[dim * i + 2];

    _point_mesh->particles()[i]->add_particle_force(pforce);
  }

  STOP_LOG("add_electrostatic_forces()", "PMSystemPoisson");
}

// ==================================================================================
Real PMSystemPoisson::local_potential_field(const Point      & p,
                                            const std::string& charge_type) const
{
  START_LOG("local_potential_field()", "PMSystemPoisson");

  Real phi_local =
    ggem_poisson->local_potential_field(_point_mesh, p, charge_type);

  STOP_LOG("local_potential_field()", "PMSystemPoisson");
  return phi_local;
}

// ==================================================================================
Real PMSystemPoisson::local_potential_field(const Elem        *elem,
                                            const Point      & p,
                                            const std::string& charge_type) const
{
  START_LOG("local_potential_field()", "PMSystemPoisson");

  Real phi_local = ggem_poisson->local_potential_field(_point_mesh,
                                                       elem,
                                                       p,
                                                       charge_type);

  STOP_LOG("local_potential_field()", "PMSystemPoisson");
  return phi_local;
}

// ==================================================================================
Real PMSystemPoisson::local_potential_bead(const std::size_t& bead_id,
                                           const std::string& charge_type) const
{
  START_LOG("local_potential_bead()", "PMSystemPoisson");

  Real phi_local = ggem_poisson->local_potential_bead(_point_mesh,
                                                      bead_id,
                                                      charge_type);

  STOP_LOG("local_potential_bead()", "PMSystemPoisson");
  return phi_local;
}

// ==================================================================================
void PMSystemPoisson::test_potential_profile()
{
  START_LOG("test_potential_profile()", "PMSystemPoisson");

  PMToolBox::output_message("========>2. Test in PMSystemPoisson::"
                            "test_potential_profile(): \n", this->comm());
  // Solve the electrical potential field: global solution by FEM
  _re_init = true; // assemble global matrix and init ksp_solver
  this->solve("unused");
  // get system dimension
  MeshBase& mesh = this->get_mesh();
  const std::size_t dim = mesh.mesh_dimension();
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
      outfile << filenames[dim_i] << ",phi_total,phi_global,phi_local,"
                                     "phi_exact\n";
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
      // global potential from FEM
      Real phi_global, phi_total;
      const unsigned int phi_var = this->variable_number("phi"); // phi_var = 0
      // get global phi from FEM solution, this is slow, but we tolerate it
      // only for test
      phi_global = this->point_value(phi_var, pt);
      // get local potential from analytical function
      const Real phi_local = this->local_potential_field(pt, "regularized");
      // get total solution at this point
      phi_total = phi_global + phi_local;
      // get exact solution for an unbounded domain from analytical solution
      const Real phi_exact = analytical_solution->exact_solution_infinite_domain(
              *ggem_poisson, pt);
      // write the result to output file
      if (this->comm().rank() == 0)
        outfile << pt(dim_i) << "," << phi_total
                << "," << phi_global
                << "," << phi_local
                << "," << phi_exact << "\n";
    } // end for i-loop
    outfile.close();
  } // end loop dim_i

  // write equation system to output file
  std::ostringstream output_filename;
  output_filename << "output_potential_profile_" << NP << "P.e";
#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO(this->get_mesh()).write_equation_systems(output_filename.str(),
                                                       this->get_equation_systems());
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  STOP_LOG("test_potential_profile()", "PMSystemPoisson");
}



// =============================================================================
void PMSystemPoisson::update_solution_for_output(const std::string&
  solution_name)

{
  START_LOG("update_solution_for_output()", "PMSystemPoisson");

  // clone Poisson FEM System Solution to solution_backup
  this->solution_backup = this->solution->clone();
  // update system solution based on the "solution_name"
  if (solution_name == "disturbed_global")
  {
    // do nothing
  }
  else if (solution_name == "disturbed_total" or solution_name == "total")
  {
    this->add_local_solution();
  }
  else
  {
    std::ostringstream ss;
    ss << "Error: invalid solution_name: " << solution_name
       << "; Suggested options are: "
       << "'disturbed_global', 'disturbed_total'"
       << ", 'undisturbed', 'total'. Exiting ...";
    PMToolBox::output_message(ss.str(), this->comm());
  }

  STOP_LOG("update_solution_for_output()", "PMSystemPoisson");
}

} // end of namespace
