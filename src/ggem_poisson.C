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
Real GGEMPoisson::green_tensor_local_regularized(const Point& x
                                                 ) const
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
    const Real r = x.norm();
    G = (std::erf(ksi * r) - std::erf(alpha * r)) / r;
  } // end if-else

  STOP_LOG("green_tensor_local_regularized()", "GGEMPoisson");

  // done
  return G;
}

// ======================================================================
Real GGEMPoisson::local_potential_field(PointMesh<3>      *point_mesh,
                                        const Point      & ptx,
                                        const std::string& charge_type) const
{
  START_LOG("local_potential_field()", "GGEMPoisson");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     build the particle neighbor list around the given point \p ptx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  const bool is_sorted = false;
  std::vector<std::pair<std::size_t, Real> > IndicesDists;
  point_mesh->build_particle_neighbor_list(ptx, is_sorted, IndicesDists);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     Loop over all the neighbor list beads, and compute the local potential:
     phi_l(i) = sum[v=1:Nl] G_v(x-x_v; i,j)*q_v(j)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  Real phi = 0;
  Real GT;

  for (std::size_t v = 0; v < IndicesDists.size(); ++v)
  {
    // 0. particle id and position, vector x = ptx - pt0
    const std::size_t p_id = IndicesDists[v].first;
    const Point pt0        = point_mesh->particles()[p_id]->point();
    const Point x          = point_mesh->pm_periodic_boundary()->point_vector(pt0,
                                                                              ptx);

    // 1. compute the Green function of particle-v
    if (charge_type == "regularized") {
      GT = this->green_tensor_local_regularized(x);
    }
    else {
      libmesh_assert("GGEMPoisson::local_potential_field, wrong charge_type!");
    } // end if-else

    // 2. Get the charge of particle-v
    const Real q = point_mesh->particles()[p_id]->charge();

    // 3. compute phi
    phi += GT * q;
  } // end for v-loop

  STOP_LOG("local_potential_field()", "GGEMPoisson");

  return phi;
}

// ======================================================================
Real GGEMPoisson::local_potential_field(PointMesh<3>      *point_mesh,
                                        const Elem        *elem,
                                        const Point      & ptx,
                                        const std::string& charge_type) const
{
  START_LOG("local_potential_field()", "GGEMPoisson");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     build the particle neighbor list around the given point \p ptx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  const bool is_sorted                    = false;
  std::vector<std::size_t> elem_neighbors = point_mesh->elem_neighbor_list(elem);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     Loop over all the neighbor list beads, and compute the local potential:
     phi_l(i) = sum[v=1:Nl] G_v(x-x_v; i,j)*q_v(j)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  Real phi = 0.;
  Real GT;

  for (std::size_t v = 0; v < elem_neighbors.size(); ++v)
  {
    // 0. particle id and position, vector x = ptx - pt0
    const std::size_t p_id = elem_neighbors[v];
    const Point pt0        = point_mesh->particles()[p_id]->point();
    const Point x          = point_mesh->pm_periodic_boundary()->point_vector(pt0,
                                                                              ptx);

    if (charge_type ==
        "regularized") GT = this->green_tensor_local_regularized(x);
    else libmesh_assert("GGEMPoisson::local_potential_field, wrong charge_type!");

    // 2. Get the charge of particle-v
    const Real q = point_mesh->particles()[p_id]->charge();

    // 3. compute phi due to this particle
    phi += GT * q;
  } // end for v-loop

  STOP_LOG("local_velocity_fluid()", "GGEMPoisson");

  return phi;
}

// ======================================================================
Real GGEMPoisson::local_potential_bead(PointMesh<3>      *point_mesh,
                                       const std::size_t& pid0,
                                       const std::string& charge_type) const
{
  START_LOG("local_velocity_bead()", "GGEMPoisson");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     Find the bead position ptx and its neighbor list, and its point type
     NOTE: for a given bead/tracking pt, its neighbor list does NOT include
       itself!
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  const Point& ptx            = point_mesh->particles()[pid0]->point();
  const PointType point_type0 = point_mesh->particles()[pid0]->point_type();
  std::vector<std::pair<std::size_t, Real> > IndicesDists;
  IndicesDists = point_mesh->particles()[pid0]->neighbor_list();
  const std::vector<Point>& neighbor_vector =
    point_mesh->particles()[pid0]->neighbor_vector();

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    -
     Loop over all the neighbor list beads, and compute the local potential:
     u_l(i) = sum[v=1:Nl] G_v(x-x_v; i,j)*q_v(j)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       - */
  Real phi = 0.;
  Real GT;

  // const Real ksi = this->regularization_parameter(hmin, ibm_beta, br0,
  // point_type0);
  for (std::size_t v = 0; v < IndicesDists.size(); ++v)
  {
    // 0. particle id and position, vector x = ptx - pt0
    const Point x           = neighbor_vector[v];
    const unsigned int p_id = IndicesDists[v].first;

    // 1. compute the Green function (Oseen Tensor) of particle-v
    if (charge_type ==
        "regularized") GT = this->green_tensor_local_regularized(x);
    else libmesh_assert("GGEMPoisson::local_potential_bead, wrong charge_type");

    // 2. Get the charge of particle-v
    const Real q = point_mesh->particles()[p_id]->charge();

    // 3. compute u due to this particle
    phi += GT * q;
  } // end for v-loop
  STOP_LOG("local_velocity_bead()", "GGEMPoisson");

  return phi;
}
} // end of namespace
