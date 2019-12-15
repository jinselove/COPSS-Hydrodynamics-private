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
#include "ggem_stokes.h"


namespace libMesh {
// ======================================================================
GGEMStokes::GGEMStokes()
  : GGEMSystem() {}


// ======================================================================
GGEMStokes::~GGEMStokes()
{
  // do nothing
}

// ======================================================================
void GGEMStokes::set_mu(const Real& _mu)
{
  START_LOG("set_mu()", "GGEMStokes");

  mu              = _mu;
  one_eight_pi_mu = 1. / (8. * PI * mu);

  STOP_LOG("set_mu()", "GGEMStokes");
}

// ======================================================================
void GGEMStokes::set_alpha(const Real& _alpha)
{
  START_LOG("set_alpha()", "GGEMStokes");

  alpha             = _alpha;
  alpha2            = alpha * alpha;
  alpha3            = alpha2 * alpha;
  alpha3_pi_23      = alpha3 / pi_23;
  two_alpha_sqrt_pi = 2. * alpha / sqrt_pi;
  coeff_1           = one_eight_pi_mu * 4.0 * alpha / sqrt_pi;

  STOP_LOG("set_alpha()", "GGEMStokes");
}

// ======================================================================
void GGEMStokes::set_ksi()
{
  START_LOG("set_ksi()", "GGEMStokes");

  if ((point_type == NOT_DEFINED_POINTTYPE) || (br0 == NAN) || (hmin == NAN))
  {
    std::cout <<
      "Error: error in GGEMStokes::set_ksi(), required parameters\
                      are not initialized yet." << std::endl;
    libmesh_error();
  }

  switch (point_type)
  {
  case POLYMER_BEAD: // Case I: bead type point
  {
    ksi = sqrt_pi / (3. * br0);
    break;
  }

  case POINT_PARTICLE: // Case III: point particle
  {
    ksi = sqrt_pi / (3. * br0);
    break;
  }

  case LAGRANGIAN_POINT:           // Case II: Lagrangian point of immersed body
  {
    ksi = 1.0 / (ibm_beta * hmin); // default value of ibm_beta is 0.75, which
                                   // needs to be optimized for different number
                                   // of nodes
    break;
  }

  default:
  {
    std::cout << "Error: error in GGEMStokes::set_ksi(): Undefined PointType!" <<
        std::endl;
    libmesh_error();
  }
  }

  // init related helper parameter
  ksi2            = ksi * ksi;
  ksi3            = ksi2 * ksi;
  ksi3_pi_23      = ksi3 / pi_23;
  two_ksi_sqrt_pi = 2. * ksi / sqrt_pi;
  coeff_2         = one_eight_pi_mu * 4.0 * (ksi - alpha) / sqrt_pi;
  STOP_LOG("set_ksi", "GGEMStokes");
}

// ======================================================================
Real GGEMStokes::smoothed_force_exp(const Real& r) const
{
  START_LOG("smoothed_force_exp()", "GGEMSystem");

  const Real r2   = r * r;
  const Real a2r2 = alpha2 * r2;
  Real g          = alpha3_pi_23 * std::exp(-a2r2) * (2.5 - a2r2);

  STOP_LOG("smoothed_force_exp()", "GGEMSystem");

  return g;
}

// ======================================================================
DenseMatrix<Number>GGEMStokes::green_tensor_unbounded(const Point       & x,
                                                      const DeltaFunType& delta_type)
const
{
  START_LOG("green_tensor_unbounded()", "GGEMStokes");

  // Compute the values according to the delta_type
  DenseMatrix<Number> GT(dim, dim);

  switch (delta_type)
  {
  case SINGULAR:
    GT = this->green_tensor_unbounded_singular(x);
    break;

  case SMOOTHED_EXP:
    GT = this->green_tensor_unbounded_smoothed(x, alpha);
    break;

  default:
    printf("*** error in GGEMStokes::green_tensor: undefined DeltaFunType!\n");

    libmesh_error();
  }

  STOP_LOG("green_tensor_unbounded()", "GGEMStokes");
  return GT;
}

// ======================================================================
DenseMatrix<Number>GGEMStokes::green_tensor_unbounded_singular(const Point& x
                                                               ) const
{
  START_LOG("green_tensor_unbounded_singular()", "GGEMStokes");

  const Real r = x.norm(), r2  = r * r;

  // output warning if x->0
  if (r < r_eps) {
    std::cout <<
      "Warning: warning in GGEMStokes::green_tensor_exact, \
        r->0 which may lead to singular solution!" << std::endl;
  }

  // init G tensor
  DenseMatrix<Number> G(dim, dim);

  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j) {
      // G(i,j) = c0*( _kronecker_delta[i][j] + x(i)*x(j)/r2 );
      G(i, j) = one_eight_pi_mu / r * (_kronecker_delta[i][j] + x(i) * x(j) / r2); //
                                                                                   // equation
                                                                                   // (2)
                                                                                   // in
                                                                                   // Juan.
                                                                                   // Hernandez-Ortiz,
                                                                                   // PRL
                                                                                   // 98.140602
                                                                                   // (2007)
    }
  }

  STOP_LOG("green_tensor_unbounded_singular()", "GGEMStokes");

  // done
  return G;
}

// ======================================================================
DenseMatrix<Number>GGEMStokes::green_tensor_unbounded_smoothed(const Point& x,
                                                               const Real & alpha_or_ksi
                                                               ) const
{
  START_LOG("green_tensor_unbounded_smoothed()", "GGEMStokes");

  DenseMatrix<Number> G(dim, dim);

  const Real r = x.norm();

  // check zero limit
  if (r < r_eps)
  {
    Real c0 = one_eight_pi_mu * 4.0 * alpha_or_ksi / sqrt_pi;

    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) G(i, j) = c0 * _kronecker_delta[i][j];
  } // end if()
  else
  {
    // coefficients
    const Real r2 = r * r, alpha_or_ksi2 = alpha_or_ksi * alpha_or_ksi;
    const Real c1 = std::erf(alpha_or_ksi * r) / r;
    const Real c2 = 2. * alpha_or_ksi / sqrt_pi * std::exp(-alpha_or_ksi2 * r2);

    // compute G tensor; since G is symmetric, we will built the upper half and
    // the diagonal element first, and do a mirror
    Real tmp1 = 0, tmp2 = 0, xixj_r2 = 0.;

    for (int i = 0; i < dim; ++i)
    {
      // upper half and diagonal
      for (int j = i; j < dim; ++j)
      {
        xixj_r2 = x(i) * x(j) / r2;
        tmp1    = (_kronecker_delta[i][j] + xixj_r2) * c1;
        tmp2    = (_kronecker_delta[i][j] - xixj_r2) * c2;
        G(i, j) = tmp1 + tmp2;
      } // end for j-loop

      for (int j = i - 1; j >= 0; --j) {
        G(i, j) = G(j, i);
      }
    } // end for i-loop
    G *= one_eight_pi_mu;
  }   // end if-else
  STOP_LOG("green_tensor_unbounded_smoothed()", "GGEMStokes");

  // done
  return G;
}

// ======================================================================
DenseMatrix<Number>GGEMStokes::green_tensor_local_singular(const Point& x) const
{
  START_LOG("green_tensor_local_singular()", "GGEMStokes");

  DenseMatrix<Number> G(dim, dim);

  const Real r = x.norm();

  // check if zero limit is true.
  if (r < r_eps)
  {
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) G(i, j) = coeff_1 * _kronecker_delta[i][j];
  } // end if()
  else
  {
    // coefficients
    const Real r2 = r * r;
    const Real c1 = std::erfc(alpha * r) / r;
    const Real c2 = two_alpha_sqrt_pi * std::exp(-alpha2 * r2);

    // compute G tensor; since G is symmetric, we will built the upper half and
    // the diagonal element first, and do a mirror
    Real tmp1 = 0, tmp2 = 0, xixj_r2 = 0.;

    for (int i = 0; i < dim; ++i)
    {
      // upper half and diagonal
      for (int j = i; j < dim; ++j)
      {
        xixj_r2 = x(i) * x(j) / r2;
        tmp1    = (_kronecker_delta[i][j] + xixj_r2) * c1;
        tmp2    = (_kronecker_delta[i][j] - xixj_r2) * c2;
        G(i, j) = tmp1 - tmp2;
      } // end for j-loop

      for (int j = i - 1; j >= 0; --j) {
        G(i, j) = G(j, i);
      }
    } // end for i-loop
    G *= one_eight_pi_mu;
  }   // end if-else

  STOP_LOG("green_tensor_local_singular()", "GGEMStokes");

  // done
  return G;
}

// ======================================================================
DenseMatrix<Number>GGEMStokes::green_tensor_local_regularized(const Point& x)
const
{
  START_LOG("green_tensor_local_regularized()", "GGEMStokes");

  DenseMatrix<Number> G(dim, dim);

  const Real r = x.norm();

  // check if zero limit is true.
  if (r < r_eps)
  {
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        G(i, j) = coeff_2 * _kronecker_delta[i][j];
      }
  } // end if()
  else
  {
    // coefficients
    const Real r2 = r * r;
    const Real kr = ksi * r,   kr2 = kr * kr;
    const Real ar = alpha * r, ar2 = ar * ar;
    const Real c1 = (std::erf(kr) - std::erf(ar)) / r;
    const Real c2 = two_ksi_sqrt_pi * std::exp(-kr2) - two_alpha_sqrt_pi *
                    std::exp(-ar2);

    // compute G tensor; since G is symmetric, we will built the upper half and
    // the diagonal element first, and do a mirror
    Real tmp1 = 0, tmp2 = 0, xixj_r2 = 0.;

    for (int i = 0; i < dim; ++i)
    {
      for (int j = i; j < dim; ++j)
      {
        xixj_r2 = x(i) * x(j) / r2;
        tmp1    = (_kronecker_delta[i][j] + xixj_r2) * c1;
        tmp2    = (_kronecker_delta[i][j] - xixj_r2) * c2;
        G(i, j) = tmp1 + tmp2;
      } // end for j-loop

      for (int j = i - 1; j >= 0; --j)
      {
        G(i, j) = G(j, i);
      }
    } // end for i-loop
    G *= one_eight_pi_mu;
  }   // end if-else

  STOP_LOG("green_tensor_local_regularized()", "GGEMStokes");

  // done
  return G;
}

// ======================================================================
Point GGEMStokes::local_solution_field(PointMesh<3>      *point_mesh,
                                       const Point      & ptx,
                                        const std::string&force_type,
                                const dof_id_type& ptx_elem_id) const
{
  START_LOG("local_solution_field()", "GGEMStokes");

  // if elem_id (id of the element that contains ptx) of ptx isn't
  // given, we can build particle neighbor list using KD tree
  if ((ptx_elem_id<0) | (ptx_elem_id >= point_mesh->get_mesh().n_elem()))
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
     Loop over all the neighbor list beads, and compute the local velocity:
     u_l(i) = sum[v=1:Nl] G_v(x-x_v; i,j)*f_v(j)
     zero_limit = false, if the given point is 'fluid' (not a bead/tracking pt.)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */

  // const Real br0    = 1.0;    // normalized bead radius.
  //  const bool  zero_limit  = false;  // for fluid, we let it be false, and
  // allow singularity
  Point u;
  DenseMatrix<Number> GT;
  for (std::size_t v = 0; v < point_nb_list.size(); ++v)
  {
    // 0. particle id and position, vector x = ptx - pt0
    const std::size_t& p_id = point_nb_list[v];
    const Point& pt0 = point_mesh->particles()[p_id]->point();
    const Point& x = point_mesh->pm_periodic_boundary()->point_vector(pt0, ptx);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // 1. compute the Green function (Oseen Tensor) of particle-v
    if (force_type == "regularized") {
      GT = this->green_tensor_local_regularized(x);
    }
    else {
      std::cout<<"Error: invalid force_type:"<<force_type<<std::endl;
      libmesh_error();
    } // end if-else

    // 2. compute the force vector of particle-v
    const Point fv = point_mesh->particles()[p_id]->particle_force();

    // 3. compute u due to this particle
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        u(i) += GT(i, j) * fv(j);
  } // end for v-loop

  STOP_LOG("local_solution_field()", "GGEMStokes");
  return u;
}

// =======================================================================
Point GGEMStokes::local_solution_bead(PointMesh<3>      *point_mesh,
                                      const std::size_t& bead_id,
                                      const std::string&force_type) const
{
  START_LOG("local_solution_bead()", "GGEMStokes");

  const std::vector<dof_id_type>& nb_list = point_mesh->particles()
    [bead_id]->neighbor_list();
  const std::vector<Point>& nb_vector = point_mesh->particles()
    [bead_id]->neighbor_vector();
  const PointType& point_type = point_mesh->particles()[bead_id]->point_type();
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Loop over all the neighbor list beads, and compute the local velocity:
     u_l(i) = sum[v=1:Nl] G_v(x-x_v; i,j)*f_v(j)
     zero_limit = false, if the given point is 'fluid' (not a bead/tracking pt.)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  // const Real br0    = 1.0;    // normalized bead radius.
  //  const bool  zero_limit  = false;  // for fluid, we let it be false, and
  // allow singularity
  Point u;
  DenseMatrix<Number> GT;
  for (std::size_t v = 0; v < nb_vector.size(); ++v)
  {
    // 0. particle id and position, vector x = ptx - pt0
    const dof_id_type& p_id = nb_list[v];
    const Point& x = nb_vector[v];

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // 1. compute the Green function (Oseen Tensor) of particle-v
    if (force_type == "regularized") {
      GT = this->green_tensor_local_regularized(x);
    }
    else {
      std::cout<<"Error: invalid force_type:"<<force_type<<std::endl;
      libmesh_error();
    } // end if-else

    // 2. compute the force vector of particle-v
    const Point& fv = point_mesh->particles()[p_id]->particle_force();

    // 3. compute u due to this particle
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        u(i) += GT(i, j) * fv(j);
  } // end for v-loop

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  -
   The above velocity does NOT include the contribution from the bead itself,
     since
   it is not in its own neighbor list. This is true for a polymer bead.
     However,
   if this point is a tracking point of an immersed body in Immersed Boundary
     Method,
   this contribution should be included by letting x-->0
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     - */
  if (point_type == LAGRANGIAN_POINT) // for Lagrangian tracking point type
    // only
  {
    // 1. compute the Green function (Oseen Tensor) when x-->0
    Point self_vector(0.);
    GT = this->green_tensor_local_regularized(self_vector);

    // 2. compute the force vector of this particle
    const Point& fv = point_mesh->particles()[bead_id]->particle_force();

    // 3. compute u due to this particle
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        u(i) += GT(i, j) * fv(j);
  }

  STOP_LOG("local_solution_bead()", "GGEMStokes");
  return u;
}

// ======================================================================
Point GGEMStokes::global_self_exclusion(PointMesh<3>      *point_mesh,
                                        const std::size_t& pid0
                                        ) const
{
  START_LOG("self_exclusion()", "GGEMStokes");

  // 1. compute the force vector of the particle
  const Point& fv = point_mesh->particles()[pid0]->particle_force();

  // 2. compute the self exclusion term at the position of the i-th bead
  Point self_v;

  for (int i = 0; i < dim; ++i) self_v(i) = coeff_1 * fv(i);

  STOP_LOG("self_exclusion()", "GGEMStokes");

  return self_v; // done and return
}

// // ======================================================================
// DenseMatrix<Number> GGEMStokes::rpy_tensor(const Point& x,  /* vector x = pt1
// - pt0 */
//                                            const Real& mu,  /* viscosity */
//                                            const Real& a,   /* bead radius */
//                                            const std::size_t& dim) const
// /*dim==3*/
// {
//   START_LOG ("rpy_tensor()", "GGEMStokes");
//
//   const Real  r = x.norm(), r2  = r*r, a2 = a*a;
//   const Real t1 = 1.0/(8.*PI*mu*r);
//   const Real t2 = 1.0/(6.*PI*mu*r);
//   Real c0 = 0.0, C1 = 0., C2 = 0.;
//
//   // init G tensor
//   DenseMatrix<Number> G(dim,dim);
//
//   // avoid singularity
//   if(r >= r_eps)
//   {
//     for (int i=0; i<dim; ++i)
//     {
//       for (int j=0; j<dim; ++j)
//       {
//         if( r>=2.*a )
//         {
//           c0 = t1;
//           C1 = 1.0 + (2.0*a2)/(3.0*r2);
//           C2 = 1.0 - (2.0*a2)/r2;
//         }
//         else
//         {
//           c0 = t2;
//           C1 = 1.0 - (9.0*r)/(32.0*a);
//           C2 = (3.0*r)/(32.0*a);
//         }
//         const Real dij = _kronecker_delta[i][j];
//         G(i,j) = c0*( C1*dij + C2*x(i)*x(j)/r2 );
//       } // end j-loop
//     } // end i-loop
//   } // end if
//
//   STOP_LOG ("rpy_tensor()", "GGEMStokes");
//   return G;
// }
//
//
//
// // ======================================================================
// DenseMatrix<Number> GGEMStokes::mobility_tensor(const Point& x,     /* vector
// x = pt1 - pt0 */
//                                                 const Real& mu,     /*
// kinematic viscosity */
//                                                 const Real& a,      /* bead
// radius */
//                                                 const std::size_t& dim) const
// /*dim==3*/
// {
//   START_LOG ("mobility_tensor()", "GGEMStokes");
//
//   const Real  r = x.norm();
//   const Real t1 = 1.0/(6.*PI*mu*a);
//
//   // init M & G tensor
//   DenseMatrix<Number> M(dim,dim), G(dim,dim);
//   G = this->rpy_tensor(x, mu, a, dim);
//
//   if( r<r_eps ) // rk = rl <=> r_kl = 0   => rpy_tensor = 0
//   {
//     for (int i=0; i<dim; ++i) {
//       for (int j=0; j<dim; ++j) {
//         const Real dij = _kronecker_delta[i][j];
//         M(i,j) = t1*dij;
//       }
//     }
//   }
//   else
//   {
//     for (int i=0; i<dim; ++i) {
//       for (int j=0; j<dim; ++j) {
//         M(i,j) = G(i,j);
//       }
//     }
//   }
//
//
//   STOP_LOG ("mobility_tensor()", "GGEMStokes");
//   return M;
// }


// // ======================================================================
// DenseMatrix<Number> GGEMStokes::green_tensor_diff(const Point& x,
//                                                   const Real& mu,
//                                                   const std::size_t& dim,
//                                                   const bool& zero_limit,
//                                                   const DeltaFunType&
// delta_type_a,
//                                                   const DeltaFunType&
// delta_type_k) const
// {
//   START_LOG ("green_tensor_diff()", "GGEMStokes");
//
//   // Initialization
//   DenseMatrix<Number> GTa(dim,dim),  GTk(dim,dim);
//
//   // compute GT_a and GT_k
//   GTa = this->green_tensor(x, alpha, mu, dim, zero_limit, delta_type_a);
//   GTk = this->green_tensor(x, reg_param,   mu, dim, zero_limit,
// delta_type_k);
//   GTk -= GTa;
//
//   STOP_LOG ("green_tensor_diff()", "GGEMStokes");
//
//   // done
//   return GTk;
// }


//
//// ======================================================================
// std::vector<Real> GGEMStokes::local_velocity(ParticleMesh<3>* point_mesh,
//                                             const Elem* elem,
//                                             const Point& ptx,
//                                             const Real& alpha,
//                                             const Real& mu,
//                                             const Real& ksi,
//                                             const std::size_t& dim,
//                                             const std::string& option)
// {
//  START_LOG ("local_velocity(Elem)", "GGEMStokes");
//
//  // build the particle neighbor list around the element \p elem
//  // In fact, it is not necessary to use 'sorted' neighbor list here!
//  //const bool is_sorted = point_mesh->is_sorted();
//  const bool is_sorted = false;
//  std::vector<std::size_t> n_list;
//  point_mesh->build_elem_neighbor_list(elem,is_sorted,n_list);
//
//  // the local velocity is  u_l(i) = sum[v=1:Nl] G_v(x-x_v; i,j)*f(j)
//  // Loop over all the neighbor list particles
//  std::vector<Real> u(dim,0.0);
//  for (std::size_t v=0; v<n_list.size(); ++v)
//  {
//    const std::size_t p_id = n_list[v];
//    const Point pt0 = point_mesh->particles()[p_id]->point();
//
//    // 0. Exclude the self term of local velocity
//    const Point x = point_mesh->pm_periodic_boundary()->point_vector(pt0,ptx);
//    if ( x.size()<r_eps ) continue;
//
//    // 1. compute the Green function (Oseen Tensor) of particle-v
//    DenseMatrix<Number> GT;
//    if (option=="regularized")
//      GT = this->green_tensor_regularized(x,alpha,mu,ksi,dim);
//    else if (option=="smooth")
//      GT = this->green_tensor_smooth(x,alpha,mu,dim);
//    else
//      libmesh_assert ("GGEMStokes::local_velocity, wrong option!");
//    // end if-else
//
//    // 2. compute the force vector of particle-v
//    const std::vector<Real> fv =
// point_mesh->particles()[p_id]->particle_force();
//
//    // 3. compute u due to this particle
//    for (int i=0; i<dim; ++i)
//      for (int j=0; j<dim; ++j)
//        u[i] += GT(i,j)*fv[j];
//  } // end for v-loop
//
//  // ----------------------- test output -------------------------
//  Real u_mag = 0.0;
//  for (int i=0; i<dim; ++i) u_mag += u[i]*u[i];
//  u_mag = std::sqrt( u_mag );
//  if ( n_list.size()>0 )// && u_mag>1e-10
//  {
//    printf("------------ GGEMStokes::local_velocity(elem) -------------\n");
//    printf("There are %lu neighbor particles around the element within r =
// %f\n",
//           n_list.size(), point_mesh->search_radius("e") );
//
//    for (std::size_t i=0; i<n_list.size(); ++i)
//      printf("%lu ",n_list[i]);
//    printf("\n");
//
//    printf("vel_bc = ");
//    for (int i=0; i<dim; ++i)
//      printf("%f, ", u[i] );
//    printf(" at node (%f, %f, %f)\n\n", ptx(0),ptx(1),ptx(2) );
//  }
//  // -------------------------------------------------------------
//
//  STOP_LOG ("local_velocity(Elem)", "GGEMStokes");
//  return u;
// }
} // end of namespace
