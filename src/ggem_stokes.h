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
#pragma once


#include <stdio.h>
#include <cmath>
#include <limits>

#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/point.h"
#include "point_mesh.h"
#include "ggem_system.h"

using namespace libMesh;


namespace libMesh {
/*
 * Remember that:
 * The PMSystemStokes is designed to solve
 * the Stokes equation with regularized Gaussian
 * point forces due to particles using mixed FEM, which
 * is usually known as 'global' part of the solution
 *
 * This class GGEMStokes will provide the main functionalities of
 * GGEM (General Geometry Ewald-like Method), and it will provide
 * the 'local' solution. Together with 'PMSystemStokes', this will
 * provide the complete solution for the Stokes system with
 * multiple point forces.
 *
 *
 * The current Green's functions are only for 3D cases!
 * Green's function for 2D situation is different.
 */

// class GGEMStokes :  public ReferenceCountedObject<GGEMStokes>
class GGEMStokes : public GGEMSystem {
public:

  // Constructor
  GGEMStokes();


  // Destructor
  ~GGEMStokes();

  /*
   * Modified smoothed Dirac delta function with exponent form
   * Reference: Eqn (27) in J Chem Phys. 136, 014901(2012), Yu Zhang, de Pablo
   *and Graham.
   */
  Real smoothed_force_exp(const Real& r) const;


  /*
   * Green tensors for different types of Dirac Delta function in unbounded
   *domain
   * if delta_type==SINGULAR: will return this ->
   *green_tensor_unbounded_singular
   * if delta_type==SMOOTHED_EXP: will return this ->
   *green_tensor_unbounded_smoothed
   */
  DenseMatrix<Number>green_tensor_unbounded(const Point       & x, /* vector x =
                                                                      pt1 - pt0
                                                                      */
                                            const DeltaFunType& delta_type) const;


  /*
   * Free-space Green's function of the Stokes equation, known as the stokeslet
   *if
   * force density is written as dirac function, i.e., rho(x) = sum[v=1:N]
   *Dirac(x, x_v)
   * This stokeslet is singular at x = x_v.
   * Reference: equation (2) in PRL 98,140602 (2007)
   */
  DenseMatrix<Number>green_tensor_unbounded_singular(const Point& x /* vector x
                                                                       = pt1 -
                                                                       pt0 */
                                                     ) const;


  /*
   * Green tensor for unbounded domain due to a smoothed force defined by a
   *Gaussian
   * force density: rho = sum[v=1:N] f_v*( g ).
   * If alpha_or_ksi == ksi, then this tensor gives the solution to the total
   *local_velocity
   * at the boundary of a periodic box.
   * Reference: The ksi part of Equation (31) in J Chem Phys. 136, 014901(2012),
   *Yu Zhang, de Pablo and Graham.
   * The reason why we only take the ksi part is because Equation (31) is the
   *local part of the green
   * tensor, i.e, the green tensor because of the local force, i.e., g_ski -
   *g_alpha. However, here we need
   * the total green tensor because of the total force, i.e., g_ksi
   */
  DenseMatrix<Number>green_tensor_unbounded_smoothed(const Point& x, /* vector x
                                                                        = pt1 -
                                                                        pt0 */
                                                     const Real & alpha_or_ksi)
  const;


  /*
   * Green function for the local force density, which is written as a singular
   *point
   * force (delta) - a 'local' smoothed force (g)
   * with exp form: rho = sum[v=1:N] f_v*( delta - g ).
   * The solution is u_loc(i) = sum[v=1:N] G_v(x-x_v; i,j)*f_v(j)
   * Eqn (3) in PRL 98, 140602(2007), JP Hernandez-Ortiz et al.
   * or Eqn (30) in J Chem Phys. 136, 014901(2012), Yu Zhang, de Pablo and
   *Graham.
   */
  DenseMatrix<Number>green_tensor_local_singular(const Point& x /* vector x =
                                                                   pt1 - pt0 */
                                                 ) const;


  /*
   * Green function for the local force density. The local force density is
   *regularized
   * in order to remove the singularity. For ksi^-1 = 3a/sqrt(PI), the max fluid
   *velocity
   * is equal to that of a particle with radius a and the pair mobility remains
   *positive definite.
   * See Eqn (4) in PRL 98, 140602(2007), JP Hernandez-Ortiz et al.
   * or Eqn (31)(32) in J Chem Phys. 136, 014901(2012), Yu Zhang, de Pablo and
   *Graham.
   */
  DenseMatrix<Number>green_tensor_local_regularized(const Point& x /* vector x =
                                                                      pt1 - pt0
                                                                      */
                                                    ) const;


  /* Compute the local vel-solution of the fluid at a given point ptx in the
   * field. This point cannot be a particle.
   * due to smoothed/regularized point forces.
   * force_type: regularized
   *
   * NOTE: due to the fast convergence of gaussian function, only a small group
   *of
   * particles within the neighbor list are considered. There are two ways
   * to construct this neighbor list:
   * (1) element independent: directly search particles near the given point
   *using KDTree;
   * (2) element dependent: directly use the neighbor list of the parent
   *element;
   * this function implement method (1), which contains a short neighbor list
   *
   * Eqn (33) in J Chem Phys. 136, 014901(2012), Yu Zhang, de Pablo and Graham.
   *
   * Notice that this solution only contains velocity
   */
  Point local_solution_field(PointMesh<3>      *point_mesh,
        const Point      & ptx, /* a pt in space*/
        const std::string& force_type,
        const dof_id_type& ptx_elem_id)const;

  /* Compute the local vel-solution of the fluid at the position of bead
   * due to smoothed/regularized point forces.
   * force_type: regularized
   * Notice that this solution only contains velocity
   */
  Point local_solution_bead(PointMesh<3>      *point_mesh,
                            const std::size_t& bead_id,
                             const std::string& force_type)const;


    /*
   * Self-exclusion term for the GLOBAL velocity at the i-th bead
   */
  Point global_self_exclusion(PointMesh<3>      *point_mesh,
                              const std::size_t& pid0
                              ) const;


  // ! set minimum mesh size
  void set_hmin(const Real& _hmin) {
    hmin = _hmin;
  }

  Real get_hmin() {
    return hmin;
  }

  // ! set Immersed boundary method regularization parameter, ibm_beta
  void set_ibm_beta(const Real& _ibm_beta) {
    ibm_beta = _ibm_beta;
  }

  Real get_ibm_beta() {
    return ibm_beta;
  }

  // ! set kinematic viscosity
  void set_mu(const Real& _mu);
  Real get_mu() {
    return mu;
  }

  // ! set GGEM alpha
  void set_alpha(const Real& _alpha);

  // ! set regularization parameter
  void set_ksi();


  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Difference between green_tensor_local(Gl) and green_tensor_regularized(Gr):
     G_reg(x) = G_exp(x,ksi) - G_exp(x,alpha)
     G_loc(x) = G_exact(x) - G_exp(x,alpha)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       -*/

  // DenseMatrix<Number> green_tensor_diff(const Point& x,     /* vector x = pt1
  // - pt0 */
  //                                       // const Real& alpha,  /* alpha
  // parameter      */
  //                                       const Real& mu,     /* kinematic
  // viscosity  */
  //                                       // const Real& ksi,    /*
  // regularization parameter */
  //                                       const std::size_t& dim,   /* dim==3
  // */
  //                                       const bool& zero_limit,   /* x-->0
  //  */
  //                                       const DeltaFunType& delta_type_a, /*
  // delta type of alpha */
  //                                       const DeltaFunType& delta_type_k)
  // const; /* delta type of ksi */
  //


  /* Compute the local v-solution at a point ptx in an unbounded domain
   * due to smoothed/regularized point forces.
   * NOTE: due to the fast convergence of gaussian function, only a small group
   *of
   * particles within the neighbor list are considered. There are two ways
   * to construct this neighbor list:
   * (1) element independent: directly search particles near the given point
   *using KDTree;
   * (2) element dependent: directly use the neighbor list of the parent
   *element;
   * this function implement method (2), which contains a longer neighbor list
   */

  //  std::vector<Real> local_velocity(ParticleMesh<3>* particle_mesh,
  //                                   const Elem* elem,   /* the parent element
  // */
  //                                   const Point& ptx,   /* a pt in space */
  //                                   const Real& alpha,  /* alpha parameter */
  //                                   const Real& mu,     /* kinematic
  // viscosity */
  //                                   const Real& ksi,    /* ksi */
  //                                   const std::size_t& dim, /*dim==3*/
  //                                   const std::string& option);


  /*
   * Rotne-Prager-Yamakawa(RPY) tensor
   * Reference: Jendrejack, Schwartz, de Pablo and Graham, J Chem Phys (2003)
   * Eqn (5) - (9)
   */

  // DenseMatrix<Number> rpy_tensor(const Point& x,     /* vector x = pt1 - pt0
  // */
  //                                const Real& mu,     /* kinematic viscosity
  // */
  //                                const Real& a,      /* bead radius */
  //                                const std::size_t& dim) const; /*dim==3*/


  /*
   * Mobility tensor in an infinite domain(no walls) using RPY tensor
   * Reference: Jendrejack, Schwartz, de Pablo and Graham, J Chem Phys (2003)
   * Eqn (3) & (5)
   * Note, diffusion tensor eqn(3) = this mobility tensor x kB*T.
   */

  // DenseMatrix<Number> mobility_tensor(const Point& x,     /* vector x = pt1 -
  // pt0 */
  //                                     const Real& mu,     /* kinematic
  // viscosity */
  //                                     const Real& a,      /* bead radius */
  //                                     const std::size_t& dim) const;
  // /*dim==3*/
  //

private:

  // ! kinematic viscosity
  Real mu              = NAN;
  Real one_eight_pi_mu = NAN; // one_eight_pi_mu = 1 / (8 * pi * mu)


  // ! alpha
  Real alpha2            = NAN; // alpha2 = alpha * alpha
  Real alpha3            = NAN; // alpha3 = alpha * alpha * alpha
  Real alpha3_pi_23      = NAN; // alpha3_pi_23 = alpha^3. / pi^(3./2.)
  Real two_alpha_sqrt_pi = NAN; // two_alpha_sqrt_pi = 2. * alpha / pi^(1./2.)
  Real coeff_1           = NAN; // coeff_1 = one_eight_pi_mu * 4.0 * alpha /
                                // sqrt_pi


  // ! regularization parameter
  Real ksi2            = NAN; // ksi2 = ksi * ksi
  Real ksi3            = NAN; // ksi3 = ksi^3.
  Real ksi3_pi_23      = NAN; // ksi3_pi_23 = ksi^3. / pi^(3./2.)
  Real two_ksi_sqrt_pi = NAN; // two_ksi_sqrt_pi = 2. * ksi / pi^(1./2.)
  Real coeff_2         = NAN; // coeff_2 = one_eight_pi_mu * 4.0 * (ksi - alpha)
                              // / sqrt_pi


  // ! minimum mesh size in the system
  Real hmin = NAN;

  // ! ibm_beta for Immersed boundary method
  Real ibm_beta = NAN;
};
} // end of namespace
