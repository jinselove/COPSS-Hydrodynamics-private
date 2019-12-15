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


#include "libmesh/libmesh_logging.h"
#include "ggem_poisson.h"


namespace libMesh {
// ======================================================================
GGEMPoisson::GGEMPoisson()
  : GGEMSystem() {}


// ======================================================================
GGEMPoisson::~GGEMPoisson()
{
  // do nothing
}

// ======================================================================
void GGEMPoisson::set_alpha(const Real& _alpha)
{
  START_LOG("set_alpha()", "GGEMPoisson");

  alpha             = _alpha;
  alpha2            = alpha * alpha;
  alpha3            = alpha2 * alpha;
  alpha3_pi_23      = alpha3 / pi_23;
  two_alpha_sqrt_pi = 2. * alpha / sqrt_pi;

  STOP_LOG("set_alpha()", "GGEMPoisson");
}

// ======================================================================
void GGEMPoisson::set_ksi()
{
  START_LOG("set_ksi()", "GGEMPoisson");

  if ((point_type == NOT_DEFINED_POINTTYPE) || (br0 == NAN))
  {
    std::cout <<
      "Error: error in GGEMPoisson::set_ksi(), required parameters\
                      are not initialized yet." << std::endl;
    libmesh_error();
  }

  switch (point_type)
  {
  case POLYMER_BEAD: // Case I: bead type point
  {
    ksi = 3. / br0;
    break;
  }

  case POINT_PARTICLE: // Case III: point particle
  {
    ksi = 3. / br0;
    break;
  }

  case LAGRANGIAN_POINT: // Case II: Lagrangian point of immersed body
  {
    std::cout << "Error: ksi is not defined for LAGRANGIAN_POINT yet" <<
        std::endl;
    libmesh_error();
  }

  default:
  {
    std::cout << "Error: error in GGEMPoisson::set_ksi(): Undefined PointType!" <<
        std::endl;
    libmesh_error();
  }
  }

  // init related helper parameter
  ksi2            = ksi * ksi;
  ksi3            = ksi2 * ksi;
  ksi3_pi_23      = ksi3 / pi_23;
  two_ksi_sqrt_pi = 2. * ksi / sqrt_pi;
  coeff_1         = 2. * (ksi - alpha) / sqrt_pi;

  STOP_LOG("set_ksi", "GGEMPoisson");
}

// ======================================================================
Real GGEMPoisson::smoothed_charge_exp(const Real& r) const
{
  START_LOG("smoothed_force_exp()", "GGEMSystem");

  const Real r2   = r * r;
  const Real a2r2 = alpha2 * r2;

  Real g = alpha3_pi_23 * std::exp(-a2r2);

  STOP_LOG("smoothed_force_exp()", "GGEMSystem");

  return g;
}

// ======================================================================
Real GGEMPoisson::green_tensor_unbounded(const Point       & x,
                                         const DeltaFunType& delta_type) const
{
  START_LOG("green_tensor_unbounded()", "GGEMPoisson");

  // Compute the values according to the delta_type
  Real GT;

  switch (delta_type)
  {
  case SINGULAR:
    GT = this->green_tensor_unbounded_singular(x);
    break;

  case SMOOTHED_EXP:
    GT = this->green_tensor_unbounded_smoothed(x, alpha);
    break;

  default:
    printf("*** error in GGEMPoisson::green_tensor: undefined DeltaFunType!\n");

    libmesh_error();
  }

  STOP_LOG("green_tensor_unbounded()", "GGEMPoisson");
  return GT;
}

// ======================================================================
Real GGEMPoisson::green_tensor_unbounded_singular(const Point& x
                                                  ) const
{
  START_LOG("green_tensor_unbounded_singular()", "GGEMPoisson");

  const Real r = x.norm();

  // output warning if x->0
  if (r < r_eps) {
    std::cout <<
      "Warning: warning in GGEMPoisson::green_tensor_exact, \
        r->0 which may lead to singular solution!" << std::endl;
  }

  Real G = 1 / r;

  STOP_LOG("green_tensor_unbounded_singular()", "GGEMPoisson");

  // done
  return G;
}

// ======================================================================
Real GGEMPoisson::green_tensor_unbounded_smoothed(const Point& x,
                                                  const Real & alpha_or_ksi) const
{
  START_LOG("green_tensor_unbounded_smoothed()", "GGEMPoisson");

  // DenseMatrix<Number> G(dim,dim);
  Real G;
  const Real r = x.norm();

  // check zero limit
  if (r < r_eps)
  {
    G = 2 * alpha_or_ksi / sqrt_pi;
  } // end if()
  else
  {
    G = std::erf(alpha_or_ksi * r) / r;
  } // end if-else
  STOP_LOG("green_tensor_unbounded_smoothed()", "GGEMPoisson");

  // done
  return G;
}

// ======================================================================
Real GGEMPoisson::green_tensor_local_singular(const Point& x) const
{
  START_LOG("green_tensor_local_singular()", "GGEMPoisson");

  Real G;

  const Real r = x.norm();

  // check if zero limit is true.
  if (r < r_eps)
  {
    std::cout <<
      "Warning: warning in GGEMPoisson::green_tensor_exact, \
        r->0 which may lead to singular solution!" << std::endl;
  }                               // end if()
  else
  {
    G = std::erfc(alpha * r) / r; // notice it's erfc function instead of erf
  }                               // end if-else

  STOP_LOG("green_tensor_local_singular()", "GGEMPoisson");

  // done
  return G;
}

// ======================================================================
Real GGEMPoisson::green_tensor_local_regularized(
  const Point& x) const
{
  START_LOG("green_tensor_local_regularized()", "GGEMPoisson");

  Real G;
  const Real r = x.norm();

  // check if zero limit is true.
  if (r < r_eps)
  {
    G = coeff_1; // G = 2 / sqrt_pi * (ksi - alpha)
  }              // end if()
  else
  {
    G = (std::erf(ksi * r) - std::erf(alpha * r)) / r;
  } // end if-else

  STOP_LOG("green_tensor_local_regularized()", "GGEMPoisson");

  // done
  return G;
}

// ======================================================================
Point GGEMPoisson::green_tensor_local_regularized_grad(
  const Point& x) const
{
  START_LOG("green_tensor_local_regularized_grad()", "GGEMPoisson");

  Point G_grad(0.);

  // in libmesh, x.norm() is evaluated by sqrt(r.norm_sq()), thus we do
  // this manually to avoid duplicate calculations
  const Real r2 = x.norm_sq();
  const Real r = std::sqrt(r2);

  // check if zero limit is true.
  if (r < r_eps)
  {
    // do nothing, since G_grad = (0., 0., 0.) if r->0
  }
  else
  {
    const Real r3 = r2 * r;
    // this tmp is the same for all directions
    const Real tmp =
      (two_alpha_sqrt_pi * std::exp(-ksi2 * r2) - two_alpha_sqrt_pi *std::exp(-alpha2 * r2)) / r
      + (std::erf(alpha * r) - std::erf(ksi * r)) / r3;

    // calculate gradient of green's function on all directions
    for (int i=0; i<3; i++)
      G_grad(i) = x(i) * tmp;
  } // end if-else

  STOP_LOG("green_tensor_local_regularized_grad()", "GGEMPoisson");

  // done
  return G_grad;
}



// ======================================================================
Real GGEMPoisson::local_solution_field(PointMesh<3>      *point_mesh,
                                  const Point      & ptx,
                                  const std::string& charge_type,
                                  const dof_id_type ptx_elem_id) const
{
  START_LOG("local_solution_field()", "GGEMPoisson");

  // create a reference to mesh
  const dof_id_type& n_elem = point_mesh->get_mesh().n_elem();

  // if elem_id (id of the element that contains ptx) of ptx isn't
  // given, we can build particle neighbor list using KD tree
  if ((ptx_elem_id<0) | (ptx_elem_id >= n_elem))
  {
    std::cout<<"Error: invalid element id for point ("<<ptx(0)<<","<<ptx(1)
      <<","<<ptx(2)<<"), element id = "<<ptx_elem_id<<std::endl;
    libmesh_error();
    //build the particle neighbor list around the given point \p ptx
    //    const bool is_sorted = false;
    //    std::vector<std::pair<std::size_t, Real> > IndicesDists;
    //    point_mesh->build_particle_neighbor_list(ptx, is_sorted, IndicesDists);
    //    for (std::size_t v=0; v<IndicesDists.size(); v++)
    //      point_nb_list.push_back(IndicesDists[v].first);
  }

  const std::vector<dof_id_type>& point_nb_list =
    point_mesh->get_elem_point_neighbor_list(ptx_elem_id);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     Loop over all the neighbor list beads, and compute the local potential:
     phi_l(i) = sum[v=1:Nl] G_v(x-x_v; i,j)*q_v(j)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  Real phi = 0.;
  Real GT;
  for (std::size_t v = 0; v < point_nb_list.size(); ++v)
  {
    // 0. particle id and position, vector x = ptx - pt0
    const dof_id_type& p_id = point_nb_list[v];
    const Point& pt0 = point_mesh->particles()[p_id]->point();
    const Point& x = point_mesh->pm_periodic_boundary()->point_vector(pt0, ptx);

    // 1. compute the Green function of particle-v
    if (charge_type == "regularized") {
      GT = this->green_tensor_local_regularized(x);
    }
    else {
      std::cout<<"Error: invalid charge_type. Exiting..."<<std::endl;
      libmesh_error();
//      libmesh_assert("GGEMPoisson::local_potential_field, wrong charge_type!");
    } // end if-else

    // 2. Get the charge of particle-v
    const Real q = point_mesh->particles()[p_id]->charge();

    // 3. compute phi
    phi += GT * q;
  } // end for v-loop

  STOP_LOG("local_solution_field()", "GGEMPoisson");

  return phi;
}

// ======================================================================
void GGEMPoisson::local_solution_field(PointMesh<3>*point_mesh,
                                      const Point      & ptx, /* a pt in space */
                                      const std::string& charge_type,
                                      const dof_id_type ptx_elem_id,
                                      const std::string& sol_option,
                                      std::map<Real, Point>& local_sol)
{

  // create a reference to mesh
  const dof_id_type& n_elem = point_mesh->get_mesh().n_elem();

  // if elem_id (id of the element that contains ptx) of ptx isn't
  // given, we can build particle neighbor list using KD tree
  if ((ptx_elem_id<0) | (ptx_elem_id >= n_elem))
  {
    std::cout<<"Error: invalid element id for point ("<<ptx(0)<<","<<ptx(1)
             <<","<<ptx(2)<<"), element id = "<<ptx_elem_id<<std::endl;
    libmesh_error();
    //build the particle neighbor list around the given point \p ptx
    //    const bool is_sorted = false;
    //    std::vector<std::pair<std::size_t, Real> > IndicesDists;
    //    point_mesh->build_particle_neighbor_list(ptx, is_sorted, IndicesDists);
    //    for (std::size_t v=0; v<IndicesDists.size(); v++)
    //      point_nb_list.push_back(IndicesDists[v].first);
  }

  const std::vector<dof_id_type>& point_nb_list =
    point_mesh->get_elem_point_neighbor_list(ptx_elem_id);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Loop over all the neighbor list beads, and
     and compute the local potential:
        phi_l(i) = sum[v=1:Nl] G_v(x-x_v; i,j)*q_v(j)
     and the gradient of the local potential
        phi_l_grad(i) = sum[v=1:Nl] grad(G_v(x-x_v; i,j)) * q_v(j)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  // clean local_sol before writing
  local_sol.clear();

  // Initialize variables
  Real phi = 0., GT = 0.;
  Point phi_grad, GT_grad;

  // loop over all neighbor points of this field point
  for (std::size_t v = 0; v < point_nb_list.size(); ++v)
  {
    // 0. particle id and position, vector x = ptx - pt0
    const dof_id_type& p_id = point_nb_list[v];
    const Point& pt0 = point_mesh->particles()[p_id]->point();
    const Point& x = point_mesh->pm_periodic_boundary()->point_vector(pt0, ptx);

    // 1. compute the Green function of particle-v
    if (charge_type == "regularized") {
      if (sol_option=="phi")
        GT = this->green_tensor_local_regularized(x);
      else if (sol_option=="grad")
        GT_grad = this->green_tensor_local_regularized_grad(x);
      else if (sol_option=="phi&grad"){
        GT = this->green_tensor_local_regularized(x);
        GT_grad = this->green_tensor_local_regularized_grad(x);
      }
      else{
        std::cout<<"Error: invalid solution option when calling "
                   "GGEMPoisson::local_solution_field(). sol_option="
                 <<sol_option<<std::endl;
        libmesh_error();
      }
    }
    else {
      std::cout<<"Error: invalid charge_type. Exiting..."<<std::endl;
      libmesh_error();
    } // end if-else

    // 2. Get the charge of particle-v
    const Real q = point_mesh->particles()[p_id]->charge();

    // 3. compute phi and phi_grad. Notice that if phi_grad is not needed,
    // then GT_grad will always be Point(0.), it's ok to sum it up here.
    phi += GT * q;
    phi_grad += (GT_grad * q);
  } // end for v-loop

  // insert phi and phi_grad to local_sol.
  local_sol.insert(std::make_pair(phi, phi_grad));

  STOP_LOG("local_solution_field()", "GGEMPoisson");
}

// ======================================================================
Real GGEMPoisson::local_solution_bead(PointMesh<3>      *point_mesh,
                                 const dof_id_type& bead_id,
                                 const std::string& charge_type) const
{
  START_LOG("local_solution_bead()", "GGEMPoisson");

  const std::vector<dof_id_type>& nb_list = point_mesh->particles()
    [bead_id]->neighbor_list();
  const std::vector<Point>& nb_vector = point_mesh->particles()
    [bead_id]->neighbor_vector();
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Loop over all the neighbor list beads, and compute the local potential:
     phi_l(i) = sum[v=1:Nl] G_v(x-x_v; i,j)*q_v(j)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
  Real phi = 0.;
  Real GT;
  for (std::size_t v = 0; v < nb_list.size(); ++v)
  {
    // 1. compute the Green function of particle-v
    if (charge_type == "regularized") {
      GT = this->green_tensor_local_regularized(nb_vector[v]);
    }
    else {
      std::cout<<"Error: invalid charge_type. Exiting..."<<std::endl;
      libmesh_error();
    } // end if-else

    // 2. Get the charge of particle-v
    const Real q = point_mesh->particles()[nb_list[v]]->charge();

    // 3. compute phi
    phi += GT * q;
  } // end for v-loop

  STOP_LOG("local_solution_bead()", "GGEMPoisson");

  return phi;
}

// ===================================================================
void GGEMPoisson::local_solution_bead(PointMesh<3>*point_mesh,
                         const dof_id_type& bead_id,
                         const std::string& charge_type,
                         const std::string& sol_option,
                         std::map<Real, Point>& local_sol)
{
  START_LOG("local_solution_bead()", "GGEMPoisson");

  const std::vector<dof_id_type>& nb_list = point_mesh->particles()
  [bead_id]->neighbor_list();
  const std::vector<Point>& nb_vector = point_mesh->particles()
  [bead_id]->neighbor_vector();
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Loop over all the neighbor list beads, and compute the local potential:
     phi_l(i) = sum[v=1:Nl] G_v(x-x_v; i,j)*q_v(j)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  // clean local_sol before writing
  local_sol.clear();

  Real phi = 0., GT = 0.;
  Point phi_grad, GT_grad;

  // loop over all neighbor points of this bead. Notice that this neighbor
  // list will not include the particle itself
  for (std::size_t v = 0; v < nb_list.size(); ++v)
  {
    // 1. compute the Green function of particle-v
    if (charge_type == "regularized") {
      if (sol_option=="phi")
        GT = this->green_tensor_local_regularized(nb_vector[v]);
      else if (sol_option=="grad")
        GT_grad = this->green_tensor_local_regularized_grad(nb_vector[v]);
      else if (sol_option=="phi&grad"){
        GT = this->green_tensor_local_regularized(nb_vector[v]);
        GT_grad = this->green_tensor_local_regularized_grad(nb_vector[v]);
      }
      else{
        std::cout<<"Error: invalid solution option when calling "
                   "GGEMPoisson::local_solution_bead(). sol_option="
                 <<sol_option<<std::endl;
        libmesh_error();
      }
    }
    else {
      std::cout<<"Error: invalid charge_type. Exiting..."<<std::endl;
      libmesh_error();
    } // end if-else

    // 2. Get the charge of particle-v
    const Real q = point_mesh->particles()[nb_list[v]]->charge();

    // 3. compute phi and phi_grad. Notice that if phi_grad is not needed,
    // then GT_grad will always be Point(0.), it's ok to sum it up here.
    phi += GT * q;
    phi_grad += GT_grad * q;
  } // end for v-loop

  // insert phi and phi_grad to local_sol.
  local_sol.insert(std::make_pair(phi, phi_grad));

  STOP_LOG("local_solution_bead()", "GGEMPoisson");
}

} // end of namespace
