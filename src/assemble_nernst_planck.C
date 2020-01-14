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


// C++ includes
#include <algorithm>
#include <math.h>

// Libmesh includes
#include "libmesh/exodusII_io.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/auto_ptr.h"

// For systems of equations the DenseSubMatrix and DenseSubVector provide
// convenient ways
// for assembling the element matrix and vector on a component-by-component
// basis.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/mesh.h"
#include <fstream>
// User defined header includes
#include "pm_toolbox.h"
#include "assemble_nernst_planck.h"
#include "pm_system_np.h"
#include "pm_system_poisson.h"
#include "pm_system_stokes.h"

// ==================================================================================
AssembleNP::AssembleNP(EquationSystems& es,
                       const std::string& name)
  : AssembleSystem(es)
{
  libmesh_assert_equal_to(name, "NP");
  analytical_solution = new AnalyticalSolutionNP(name);
  // Initialize boundary sides for BC
  _boundary_sides_dirichlet_np.resize(_mesh.n_elem());
}

// ==================================================================================
AssembleNP::~AssembleNP()
{
  delete analytical_solution; analytical_solution = nullptr;
}

// ==================================================================================
void AssembleNP::assemble_global_K(const std::string& system_name,
                                   const std::string& option)
{
  START_LOG("assemble_global_K()", "AssembleNP");

  (void) option;

  // Get a Reference to the LinearImplicitSystem we are solving
  PMSystemNP& np_system = _eqn_sys.get_system<PMSystemNP> (system_name);

  // Get a const reference to the Finite element type for the first (and only
  // for NP system) variable in the system
  FEType fe_type = np_system.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe (FEBase::build(_dim, fe_type));
  UniquePtr<FEBase> fe_face (FEBase::build(_dim, fe_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  QGauss qrule (_dim, fe_type.default_quadrature_order());
  QGauss qface (_dim-1, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);
  fe_face->attach_quadrature_rule(&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.  We will start
  // with the element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW      = fe->get_JxW();
  const std::vector<Real> & JxW_face = fe_face->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();
  const std::vector<std::vector<Real> > & psi = fe_face->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // The XY locations of the quadrature points used for face integration
  const std::vector<Point> & qface_points = fe_face->get_xyz();

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  // in future examples.
  const DofMap & dof_map = np_system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // initialize node penalty
  _node_penalty.resize(_mesh.n_elem());

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.
  MeshBase::const_element_iterator el = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh
    .active_local_elements_end();

  for(; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem * elem = *el;
    const dof_id_type&  elem_id = elem->id();

    // Fill _boundary_sides_dirichlet_np
    this->select_boundary_side(elem, system_name);

    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices (elem, dof_indices);

    // Compute the element-specific data for the current
    // element.  This involves computing the location of the
    // quadrature points (q_point) and the shape functions
    // (phi, dphi) for the current element.
    fe->reinit (elem);

    // Zero the element matrix and right-hand side before
    // summing them.  We use the resize member here because
    // the number of degrees of freedom might have changed from
    // the last element.  Note that this will be the case if the
    // element type is different (i.e. the last element was a
    // triangle, now we are on a quadrilateral).
    Ke.resize (dof_indices.size(),
               dof_indices.size());
//    Fe.resize (dof_indices.size());

    // Now we will build the element matrix.
    // Constructing the RHS requires the solution and its
    // gradient from the previous timestep.  This myst be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    // loop over Gauss quadrature points
    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {
      // loop over rows of Ke
      for (unsigned int i=0; i<phi.size(); i++)
      {
        // loop over cols of Ke
        for (unsigned int j=0; j<phi.size(); j++)
        {
          Ke(i,j) += JxW[qp] * (
                                  // Mass term
                                  phi[i][qp]*phi[j][qp]
                                  // Diffusion term when using semi implicit
                                  // Euler for Diffusion
                                  + 0.5 * np_system.dt * np_system
                                  .ion_diffusivity * dphi[i][qp] *
                                  dphi[j][qp]
                               );
        } // end loop over cols of Ke
      } // end loop over rows of Ke
    } // end loop for qp

    // apply bc by penalty on K matrix
    {
      // The following loops over the sides of the element.
      // If the element has no neighbor on a side then that
      // side MUST live on a boundary of the domain.
      for (unsigned int s=0; s<elem->n_sides(); s++)
      {
        if (elem->neighbor_ptr(s) == libmesh_nullptr) {
          fe_face->reinit(elem, s);
          for (unsigned int qp = 0; qp < qface.n_points(); qp++)
          {
            // Matrix contribution
            for (unsigned int i = 0; i < psi.size(); i++)
              for (unsigned int j = 0; j < psi.size(); j++)
                Ke(i, j) += penalty * JxW_face[qp] * psi[i][qp] * psi[j][qp];
          } // end for loop over qp
        } // end if condition
      } // end for loop over all sides

      // resize the penalty vector of this element. This has to be done after
      // fe_face->reinit(elem, s) so that psi.size() is not zero. And notice
      // that psi.size() should not change for different s so we can do it
      // after the last s. For elements without sides on the boundaries, psi
      // .size() will be just zero.
      _node_penalty[elem_id].resize(psi.size(), 0.);
    }
    // If this assembly program were to be used on an adaptive mesh,
    // we would have to apply any hanging node constraint equations.
    dof_map.constrain_element_matrix(Ke, dof_indices);

    // Add the element matrix to the global system.
    np_system.matrix->add_matrix(Ke, dof_indices);
  } // end for loop over local elements

  STOP_LOG("assemble_global_K()", "AssembleNP");
}

// ==================================================================================
void AssembleNP::assemble_global_F(const std::string& system_name,
                                   const std::string& option) {
  START_LOG("assemble_global_F()", "AssembleNP");

  // Get a Reference to the LinearImplicitSystem we are solving
  PMSystemNP &np_system = _eqn_sys.get_system<PMSystemNP>(system_name);

  // Poisson system
  PMSystemPoisson* poisson_system = nullptr;
  DofMap* poisson_dof_map = nullptr;
  unsigned short int phi_var;


  // Stokes system
  PMSystemStokes* stokes_system = nullptr;
  std::vector<unsigned short int> vel_vars(_dim);
  DofMap* stokes_dof_map = nullptr;

  // initialize a pointer to all PMSystemNP systems; no need to clean these
  // pointers afterwards since we will only attach the reference here and the
  // actual object will be destroyed somewhere else.
  std::vector<PMSystemNP*> np_systems;

  // loop over all systems and push the references of all NP systems and
  // their dof maps to np_systems and np_dof_maps vector respectively
  unsigned int n_sys = _eqn_sys.n_systems();
  for (unsigned int s_id=0; s_id<n_sys; s_id++){
    const std::string& s_name = _eqn_sys.get_system(s_id).name();
    if (s_name.rfind("NP:", 0)==0){
      // np system
      np_systems.push_back(&(_eqn_sys.get_system<PMSystemNP>(s_name)));
//      // np system dof map
//      np_dof_maps.push_back(&(_eqn_sys.get_system<PMSystemNP>(s_name)
//        .get_dof_map()));
    }
  }

  // Reassign pointers depending on 'option'; notice that these pointers will
  // be destroyed when the systems they point to are destroyed at the end of
  // the simulation so that there will not be data leak; cannot use smarter
  // pointer here since when smarter pointer at the end of this function,
  // the system they point to will be destroyed too, thus the simulation
  // won't continue
  if (option == "diffusion")
  {
    // no need to do anything
  }
  else if (option == "diffusion&convection"){
    stokes_system = &(_eqn_sys.get_system<PMSystemStokes>("Stokes"));
  }
  else if (option == "diffusion&electrostatics"){
    poisson_system = &(_eqn_sys.get_system<PMSystemPoisson>("Poisson"));
  }
  else if (option == "diffusion&electrostatics&convection"){
    poisson_system = &(_eqn_sys.get_system<PMSystemPoisson>("Poisson"));
    stokes_system = &(_eqn_sys.get_system<PMSystemStokes>("Stokes"));
  }
  else{
    std::stringstream oss;
    std::cout << "Error: invalid option in AssembleNP::assemble_global_F"
           "(system_name, option), option specified is :" << option
            << ". Exiting...";
    libmesh_error();
  }
  // prepare FEM objects for Stokes system if existed
  if (stokes_system != nullptr){
    stokes_dof_map = &(stokes_system->get_dof_map());
    vel_vars[0] = stokes_system->variable_number("u");
    vel_vars[1] = stokes_system->variable_number("v");
    vel_vars[2] = stokes_system->variable_number("w");
  }
  if (poisson_system != nullptr){
    poisson_dof_map = &(poisson_system->get_dof_map());
    phi_var = poisson_system->variable_number("phi");
  }

  // Get a const reference to the Finite element type for the first (and only
  // for NP system) variable in the system
  FEType fe_type = np_system.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe (FEBase::build(_dim, fe_type));
  UniquePtr<FEBase> fe_face (FEBase::build(_dim, fe_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  QGauss qrule (_dim, fe_type.default_quadrature_order());
  QGauss qface (_dim-1, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);
  fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.  We will start
  // with the element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW      = fe->get_JxW();
  const std::vector<Real> & JxW_face = fe_face->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();
  const std::vector<std::vector<Real> > & psi = fe_face->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // The XY locations of the quadrature points of 3D elements, this will give
  // the coordinates of quadrature points in real space instead of reference
  // space
  const std::vector<Point> & q_points = fe->get_xyz();

  // The XY locations of the quadrature points used for face integration
  const std::vector<Point> & qface_points = fe_face->get_xyz();

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  // in future examples.
  const DofMap & dof_map = np_system.get_dof_map();

  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.

  const dof_id_type& n_elem = _mesh.n_elem();
//  std::ofstream outfile;
//  outfile.open("debug.csv", std::ios_base::out);
//  outfile<<"np_coeff,dt,2nd(ion+particle),2nd(ion),1st_x,1st_y,1st_z,coeff_2\n";

  MeshBase::const_element_iterator el = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh
    .active_local_elements_end();

  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem * elem = *el;
    const dof_id_type& elem_id = elem->id();

    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices (elem, dof_indices);

    // Compute the element-specific data for the current
    // element.  This involves computing the location of the
    // quadrature points (q_point) and the shape functions
    // (phi, dphi) for the current element.
    fe->reinit (elem);

    // Zero the element matrix and right-hand side before
    // summing them.  We use the resize member here because
    // the number of degrees of freedom might have changed from
    // the last element.  Note that this will be the case if the
    // element type is different (i.e. the last element was a
    // triangle, now we are on a quadrilateral).
    Fe.resize (dof_indices.size());

    // Now we will build the element matrix and right-hand-side.
    // Constructing the RHS requires the solution and its
    // gradient from the previous timestep.  This myst be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    // for Diffusion only system
    if ((stokes_system == nullptr) and (poisson_dof_map == nullptr))
    {
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {
        // Values to hold the old system solution & its gradient at this qp point
        Number c_old = 0.;
        Gradient gradient_c_old;

        // compute the old solution & its gradient at this qp point
        for (unsigned int l = 0; l < phi.size(); l++)
        {
          // get old solution at this qp point
          // notice that we are not using Libmesh::TransientSystem, thus the
          // current_solution of our system is actually solution from previous
          // time step
          c_old +=
            phi[l][qp] * (*np_system.current_local_solution)(dof_indices[l]);

          // get solution gradient at this qp_point
          gradient_c_old.add_scaled(dphi[l][qp],
                                    (*np_system.current_local_solution)(
                                      dof_indices[l]));
        }
        // Now compute the element rhs
        const Gradient coeff_1 = -0.5 * np_system.dt * np_system
          .ion_diffusivity * gradient_c_old;

        for (unsigned int i = 0; i < phi.size(); i++)
        {
          Fe(i) += JxW[qp] * (
            // Mass term
            c_old * phi[i][qp]
            // diffusion term when using semi-implicit Euler for diffusion
            + coeff_1 * dphi[i][qp]
            );
        }
      } // end loop over qp
    }
    // for diffusion+convection system
    else if((stokes_system != nullptr) and (poisson_dof_map == nullptr))
    {
      // get dof_indices associated with u, v, w of
      // the stokes system for this element
      std::vector<std::vector<dof_id_type>> stokes_dof_indices(_dim,
        std::vector<dof_id_type>(phi.size()));
      for (int dim_i=0; dim_i<_dim; dim_i++)
        stokes_dof_map->dof_indices(elem,
          stokes_dof_indices[dim_i], vel_vars[dim_i]);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {
        // Values to hold the old system solution & its gradient at this qp point
        Number c_old = 0.;
        Gradient gradient_c_old;

        // Initialize global, local and undisturbed velocity at this qp point
        Point total_vel_old;

        // compute the old solution & its gradient at this qp point
        for (unsigned int l = 0; l < phi.size(); l++)
        {
          // get old solution at this qp point
          // notice that we are not using Libmesh::TransientSystem, thus the
          // current_solution of our system is actually solution from previous
          // time step
          c_old +=
            phi[l][qp] * (*np_system.current_local_solution)(dof_indices[l]);

          // get solution gradient at this qp_point
          gradient_c_old.add_scaled(dphi[l][qp],
                                    (*np_system.current_local_solution)(
                                      dof_indices[l]));

          // get velocity vector (size=_dim) from stokes solution at this
          // qp point, this only gives us the global + undisturbed
          // contribution of velocity
          for (int dim_i=0; dim_i<_dim; dim_i++) {
            total_vel_old(dim_i) += phi[l][qp] *
              (
                stokes_system->current_local_solution->operator()
                  (stokes_dof_indices[dim_i][l])
                + stokes_system->local_undisturbed_solution
                [stokes_dof_indices[dim_i][l]]
              );
          }
        }

        // add local velocity of this qp point
        total_vel_old += stokes_system->local_velocity_fluid(q_points[qp],
          "regularized", elem_id);

        // Now compute the element rhs
        const Gradient coeff_1 = -0.5 * np_system.dt * np_system
          .ion_diffusivity * gradient_c_old;
        const Real coeff_2 = -np_system.dt * (total_vel_old * gradient_c_old);
        for (unsigned int i = 0; i < phi.size(); i++)
        {
          Fe(i) += JxW[qp] * (
            // Mass term
            c_old * phi[i][qp]
            // diffusion term when using semi-implicit Euler for diffusion
            + coeff_1 * dphi[i][qp]
            // convection term when using fully explicit Euler for convection
            + coeff_2 * phi[i][qp]
            );
        }
      } // end loop over qp
    }
      // for diffusion&electrostatics&convection system
    else if((stokes_system == nullptr) and (poisson_system != nullptr))
    {
      // get dof_indices associated with "phi" of the Poisson system
      std::vector<dof_id_type> poisson_dof_indices(phi.size());
      poisson_dof_map->dof_indices(elem, poisson_dof_indices);

      // loop over all qp points to compute rhs
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {
        // Values to hold the old system solution & its gradient at this qp point
        Number c_old = 0.;
        Gradient gradient_c_old;

        // evaluate laplacian of the total potential on this qp point
        // (contribution from discrete charges)
        Real total_potential_laplacian_old =
          poisson_system->total_potential_laplacian_field(q_points[qp],
                                                          "regularized", elem_id);



        // evaluate gradient of the local potential on this qp point
        std::pair<Real, Point> local_potential_old;
        poisson_system->local_potential_field(q_points[qp], "regularized",
                                              "grad", local_potential_old, elem_id);

        // initialize the gradient of the global potential on this qp point,
        // the value of this gradient will be calculated using interpolation
        Gradient global_potential_gradient_old;

        // charge density contributed by all ions at this qp point
        Real charge_density_ion_old = 0.;

        // compute the old solution & its gradient at this qp point
        for (unsigned int l = 0; l < phi.size(); l++)
        {
          // get old solution at this qp point
          // notice that we are not using Libmesh::TransientSystem, thus the
          // current_solution of our system is actually solution from previous
          // time step
          c_old +=
            phi[l][qp] * (*np_system.current_local_solution)(dof_indices[l]);

          // get solution gradient at this qp_point
          gradient_c_old.add_scaled(dphi[l][qp],
                                    (*np_system.current_local_solution)(
                                      dof_indices[l]));

          // interpolate the global potential gradient at this qp point
          global_potential_gradient_old.add_scaled(dphi[l][qp],
            poisson_system->current_local_solution->operator()
            (poisson_dof_indices[l]));

          // interpolate all np systems to get the total ion concentration *
          // valence at this qp point
          Real c_z = 0.;
          for (short int s_id=0; s_id<np_systems.size(); s_id++)
          {
            c_z += (phi[l][qp] *
              np_systems[s_id]->current_local_solution->operator()
              (dof_indices[l])) * np_systems[s_id]->ion_valence;
          }
          charge_density_ion_old += c_z;
        }
        charge_density_ion_old *= _eqn_sys.parameters.get<Real>
          ("coeff_ion_charge_density");
        total_potential_laplacian_old += -4. * pi * charge_density_ion_old;

        // coefficient for diffusion term
        const Gradient coeff_1 = -0.5 * np_system.dt * np_system
          .ion_diffusivity * gradient_c_old;
        // coefficient for electrostatics term
//        std::cout<<"np_system.np_coeff="<<np_system.np_coeff<<";"
//          <<"dt="<<np_system.dt<<";first_order="<<gradient_c_old *
//                                                  (local_potential_old.second + global_potential_gradient_old)
//          <<";second_order="<<gradient_c_old*total_potential_laplacian_old;
        const Real coeff_2 = np_system.np_coeff * np_system.dt * (c_old *
          total_potential_laplacian_old + gradient_c_old *
          (local_potential_old.second + global_potential_gradient_old));
//        std::cout<<"coeff_2="<<coeff_2;
//        outfile<<np_system.np_coeff<<","<<np_system.dt<<","
//          <<total_potential_laplacian_old<<","
//          <<-4.*pi*charge_density_ion_old<<","
//          << (local_potential_old.second + global_potential_gradient_old)(0)<<","
//          << (local_potential_old.second + global_potential_gradient_old)(1)
//          <<","
//          << (local_potential_old.second + global_potential_gradient_old)(2)
//          <<","
//          <<coeff_2<<"\n";
        for (unsigned int i = 0; i < phi.size(); i++)
        {
          Fe(i) += JxW[qp] * (
            // Mass term
            c_old * phi[i][qp]
            // diffusion term when using semi-implicit Euler for diffusion
            + coeff_1 * dphi[i][qp]
            // electrostatics term using full explicit Euler for electrostatics
            + coeff_2 * phi[i][qp]
          );
        }
      } // end loop over qp
    }
    // for diffusion&electrostatics&convection system
    else
    {
      // get dof_indices associated with u, v, w of
      // the stokes system for this element
      std::vector<std::vector<dof_id_type>> stokes_dof_indices(_dim,
        std::vector<dof_id_type>(phi.size()));
      for (int dim_i=0; dim_i<_dim; dim_i++)
        stokes_dof_map->dof_indices(elem,
          stokes_dof_indices[dim_i], vel_vars[dim_i]);

      // get dof_indices associated with "phi" of the Poisson system
      std::vector<dof_id_type> poisson_dof_indices(phi.size());
      poisson_dof_map->dof_indices(elem, poisson_dof_indices);

      // loop over all qp points to compute rhs
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {
        // Values to hold the old system solution & its gradient at this qp point
        Number c_old = 0.;
        Gradient gradient_c_old;

        // Initialize global, local and undisturbed velocity at this qp point
        Point total_vel_old;

        // evaluate laplacian of the total potential on this qp point
        Real total_potential_laplacian_old =
          poisson_system->total_potential_laplacian_field(q_points[qp],
            "regularized", elem_id);

        // evaluate gradient of the local potential on this qp point
        std::pair<Real, Point> local_potential_old;
        poisson_system->local_potential_field(q_points[qp], "regularized",
          "grad", local_potential_old, elem_id);
        // initialize the gradient of the global potential on this qp point,
        // the value of this gradient will be calculated using interpolation
        Gradient global_potential_gradient_old;

        // compute the old solution & its gradient at this qp point
        for (unsigned int l = 0; l < phi.size(); l++)
        {
          // get old solution at this qp point
          // notice that we are not using Libmesh::TransientSystem, thus the
          // current_solution of our system is actually solution from previous
          // time step
          c_old +=
            phi[l][qp] * (*np_system.current_local_solution)(dof_indices[l]);

          // get solution gradient at this qp_point
          gradient_c_old.add_scaled(dphi[l][qp],
                                    (*np_system.current_local_solution)(
                                      dof_indices[l]));

          // get velocity vector (size=_dim) from stokes solution at this
          // qp point, this only gives us the global + undisturbed
          // contribution of velocity
          for (int dim_i=0; dim_i<_dim; dim_i++) {
            total_vel_old(dim_i) += phi[l][qp] *
                                    (
                                      stokes_system->current_local_solution->operator()
                                        (stokes_dof_indices[dim_i][l])
                                      + stokes_system->local_undisturbed_solution
                                      [stokes_dof_indices[dim_i][l]]
                                    );
          }

          // interpolate the global potential gradient at this qp point
          global_potential_gradient_old.add_scaled(dphi[l][qp],
            poisson_system->current_local_solution->operator()
            (poisson_dof_indices[l]));
        }

        // add local velocity of this qp point
        total_vel_old += stokes_system->local_velocity_fluid(q_points[qp],
          "regularized", elem_id);

        // coefficient for diffusion term
        const Gradient coeff_1 = -0.5 * np_system.dt * np_system
          .ion_diffusivity * gradient_c_old;
        // coefficient for convection term
        const Real coeff_2 = -np_system.dt * (total_vel_old * gradient_c_old);
        // coefficient for electrostatics term
        const Real coeff_3 = np_system.np_coeff * np_system.dt * (c_old *
          total_potential_laplacian_old + gradient_c_old *
          (local_potential_old.second + global_potential_gradient_old));

        for (unsigned int i = 0; i < phi.size(); i++)
        {
          Fe(i) += JxW[qp] * (
            // Mass term
            c_old * phi[i][qp]
            // diffusion term when using semi-implicit Euler for diffusion
            + coeff_1 * dphi[i][qp]
            // convection term when using fully explicit Euler for convection
            + coeff_2 * phi[i][qp]
            // electrostatics term using full explicit Euler for electrostatics
            + coeff_3 * phi[i][qp]
          );
        }
      } // end loop over qp
    }

    // apply bc by penalty
    {
      // if _reinit_node_penalty is true, we recalculate penalty on the boundary
      // nodes of this element
      if (_reinit_node_penalty)
      {
        // reset the penalty value of this element to 0.
        std::fill(_node_penalty[elem_id].begin(),
          _node_penalty[elem_id].end(), 0.);
        // loop over all boundary sides of this element
        for (dof_id_type s_id=0;
          s_id<_boundary_sides_dirichlet_np[elem_id].size(); s_id++)
        {
          // extract the tuple for this side of this element, use referece to avoid data copying
          // get<0>(t) : side id
          // get<1>(t) : the dirichlet boundary id associated with this side
          // get<2>(t) : the dirichlet value associated with this side
          const auto &t = _boundary_sides_dirichlet_np[elem_id][s_id];
          // get side
          const dof_id_type& s = std::get<0>(t);
          // get the Dirichlet Boundary Value on this side
          Real boundary_value = std::get<2>(t);

          // reinit fe_face on this surface
          fe_face->reinit(elem, s);

          // loop over all qp points on this surface to integrate the penalty
          for (unsigned int qp=0; qp<qface.n_points(); qp++)
          {
            // reset boundary value to be the exact solution on boundaries
            // for the test system
            if (_eqn_sys.parameters.get<std::string>("simulation_name")=="np_validation_analytic")
              boundary_value =
                analytical_solution->exact_solution_infinite_domain(
                  qface_points[qp],
                  _eqn_sys.parameters.get<Real>("real_time"),
                  _eqn_sys.parameters.get<std::vector<Real>>
                  ("ion_diffusivity")[0]);

            // RHS contribution
            for (dof_id_type i=0; i<psi.size(); i++)
              _node_penalty[elem_id][i] += penalty * JxW_face[qp] *
                boundary_value * psi[i][qp];
          } // end loop over qp
        } // end loop over s_id (boundary side id of this element)
      } // end if _reinit_node_penalty

      // add penalty to each node of the Fe vector; notice that we cannot use
      // i<psi.size() as the for loop condition because when
      // _reinit_node_penalty is false, then psi will not be initialized, i.e
      // ., psi.size() = 0.
      for (unsigned int i=0; i<_node_penalty[elem_id].size(); i++){
//        std::cout<<"elem_id = " << elem_id
//          <<"; _node_penalty[elem_id].size()=" << _node_penalty[elem_id].size()
//          <<"; i=" << i
//          <<"; node_penalty = " << _node_penalty[elem_id][i] << std::endl;
        Fe(i) += _node_penalty[elem_id][i];
      }
    }// end apply penalty

    // If this assembly program were to be used on an adaptive mesh,
    // we would have to apply any hanging node constraint equations.
    // perf_log.push("finish");
    dof_map.constrain_element_vector(Fe, dof_indices);

    // Add the element matrix and rhs vector to the global system.
    np_system.rhs->add_vector(Fe, dof_indices);
  } // end loop over local elements
//  outfile.close();

  // set _reinit_node_penalty to false to avoid duplicate calculation
  // of _node_penalty unless for the test case where Dirichlet Boundary
  // value is not constant; Notice that this can only be done outside of the
  // elem loop, otherwise one cpu might change the variable status and affect
  // the operations on other elements, for example, cpu 0 set
  // _reinit_node_penalty to false and cpu 1 might stop calculating penalty
  // because the global status has changed.
  if (_eqn_sys.parameters.get<std::string>("simulation_name")
      =="np_validation_analytic")
    _reinit_node_penalty = true;
  else
    _reinit_node_penalty = false;
  STOP_LOG("assemble_global_F()", "AssembleNP");
}

// ==================================================================================
void AssembleNP::select_boundary_side(const Elem *elem,
                                      const std::string& system_name)
{
  START_LOG("select_boundary_side()", "AssembleNP");
  // Get a reference to the Particle-Mesh System.
  PMSystemNP & pm_system = _eqn_sys.get_system<PMSystemNP>(system_name);
  const std::size_t elem_id   = elem->id();
  _boundary_sides_dirichlet_np[elem_id].resize(0);

  // The following loops over the sides of the element. If the element has NO
  // neighbors on a side then that side MUST live on a boundary of the domain.
  for (dof_id_type s = 0; s < elem->n_sides(); s++)
  {
    // If this side is on the boundary
    if (elem->neighbor(s) == NULL)
    {
      // Check if the s-th face of this element is associated with any Dirichlet boundary id
      // if any, stop looking
      for (dof_id_type i = 0; i < pm_system.boundary_id_dirichlet_np.size();
        i++)
      {
        if (_mesh.get_boundary_info().has_boundary_id(elem, s,
                                                      pm_system.boundary_id_dirichlet_np[i]))
        {
          _boundary_sides_dirichlet_np[elem_id].push_back(
            std::make_tuple(s,
                            pm_system.boundary_id_dirichlet_np[i],
                            pm_system.boundary_value_dirichlet_np[i]));
          break;
        }
      }
    }
  }

  STOP_LOG("select_boundary_side()", "AssembleNP");
}