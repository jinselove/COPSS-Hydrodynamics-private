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
#include "libmesh/equation_systems.h"

#include "../particle_mesh.h"
#include "../point_mesh.h"
#include "../pm_periodic_boundary.h"
#include "fix_base.h"


using namespace libMesh;


// ======================================================================
FixBase::FixBase()
{
  // Do nothing
}

// ======================================================================
// FixBase::FixBase(PMLinearImplicitSystem& pm_sys)
// : _pm_system(pm_sys)
// {
//  // initialize the private memebers
// }


// ======================================================================
FixBase::~FixBase()
{
  // do nothing
}

// ======================================================================
Point FixBase::spring_force_wls(const Point& pt_ij,
                                const Real & c1,
                                const Real & Ls) const
{
  START_LOG("spring_force_wls(pt_ij)", "FixBase");

  // The spring force: pt_ij = ptj - pti(with periodicity), R_ij = |pt_ij|
  const Real R_ij  = pt_ij.norm();
  const Real ratio = R_ij / Ls;
  Real val0        = 1.0 - ratio;

  // avoid singularity
  if (std::abs(val0) <  tol) val0 = val0 / std::abs(val0) * tol;

  // compute force
  const Real tmp = 1. / (val0 * val0);
  Point F_ij     = c1 * (tmp - 1. + 4. * ratio) * pt_ij / R_ij;

  STOP_LOG("spring_force_wls(pt_ij)", "FixBase");
  return F_ij;
}

// ======================================================================
Point FixBase::spring_force_fene(const Point& pt_ij,
                                 const Real & c1,
                                 const Real & Ls) const
{
  START_LOG("spring_force_fene(pt_ij)", "FixBase");

  // The spring force: pt_ij = ptj - pti(with periodicity), R_ij = |pt_ij|
  const Real R_ij  = pt_ij.norm();
  const Real ratio = R_ij / Ls;
  Real val0        =  1.0 - ratio * ratio;

  // avoid singularity
  if (std::abs(val0) <  tol) val0 = val0 / std::abs(val0) * tol * tol;

  // compute force
  Point F_ij = c1 * ratio / val0 * pt_ij / R_ij;

  STOP_LOG("spring_force_fene(pt_ij)", "FixBase");
  return F_ij;
}

// ======================================================================
Point FixBase::spring_force_ud(const Point& pt_ij,
                               const Real & c1,
                               const Real & Ls) const
{
  START_LOG("spring_force_ud(pt_ij)", "FixBase");

  // The spring force: pt_ij = ptj - pti(with periodicity), R_ij = |pt_ij|
  const Real R_ij  = pt_ij.norm(); // vector length
  const Real ratio = R_ij / Ls;
  const Real r2    = ratio * ratio;
  const Real c2    = c1 * c1;

  const Real a1 = 1.0;
  const Real a2 = -7.0 * c1;
  const Real a3 = 3.0 / 32.0 - 0.75 * c1 - 6. * c2;
  const Real a4 = (13.0 / 32.0 + 0.8172 * c1 - 14.79 * c2) /
                  (1.0 - 4.225 * c1 + 4.87 * c2);


  const Real T4 = (1.0 - r2);
  const Real T2 = 1.0 / T4;
  const Real T1 = T2 / T4; // = 1/(T4*T4)
  Point F_ij    = (a1 * T1 + a2 * T2 + a3 + a4 * T4) * pt_ij / R_ij;

  STOP_LOG("spring_force_ud(pt_ij)", "FixBase");
  return F_ij;
}

// ======================================================================
Point FixBase::spring_force_lhs(const Point& pt_ij,
                                const Real & l0,
                                const Real & k0) const
{
  START_LOG("spring_force_lhs(pt_ij)", "FixBase");

  // f_ij = k0*( |R_ij| - l0 )  * R_ij/|R_ij|

  // The spring force: pt_ij = ptj - pti(with periodicity), R_ij = |pt_ij|
  const Real R_ij = pt_ij.norm();

  if (std::abs(R_ij) < 1E-6)
  {
    printf("--->test in ForceField::spring_force_lhs(): Rij = %f\n", R_ij);
  }
  const Real F0 = k0 * (R_ij - l0);  // force magnitude
  Point F_ij    = F0 * pt_ij / R_ij; // force direction

  STOP_LOG("spring_force_lhs(pt_ij)", "FixBase");
  return F_ij;
}

// ======================================================================
Point FixBase::gaussian_force(const Point& r_ij,
                              const Real & c1,
                              const Real & c2) const
{
  START_LOG("gaussian_force(pt_ij)", "FixBase");

  // f_ij = c1*c2* exp( -c2*|r_ij|^2 ) * r_ij
  // The spring force: pt_ij = ptj - pti(with periodicity), R_ij = |pt_ij|
  const Real r_ij_size = r_ij.norm();
  Point force = -c1 *c2 *std::exp(-c2 * r_ij_size * r_ij_size) * r_ij;

  STOP_LOG("gaussian_force(pt_ij)", "FixBase");
  return force;
}

// ======================================================================
Point FixBase::lj_force(const Point& r_ij,        // direction vector
                        const Real & epsilon,     // energy coefficient
                        const Real & sigma) const // distance coefficient
{
  START_LOG("lj_force(&r_ij, &epsilon, &sigma)", "FixBase");

  // f_ij = -24 * epsilon * (2*(sigma/|r_ij|)^12 - (sigma/|r_ij|)^6 ) * r_ij /
  // |r_ij|^2

  Real  r     = r_ij.norm();
  Real  r2    = r * r;
  Real  sr2   = sigma * sigma / r2;
  Real  sr6   = sr2 * sr2 * sr2;
  Real  sr12  = sr6 * sr6;
  Point force = -24. * epsilon * (2. * sr12 - sr6) / r2 * r_ij;

  STOP_LOG("lj_force(&r_ij), &epsilon, &sigma", "FixBase");
  return force;
}

// ======================================================================
Point FixBase::harmonic_force(const Point& r_ij,
                              const Real & k,
                              const Real & r0) const
{
  START_LOG("harmonic_force(&r_ij, &r0, &k)", "FixBase");

  Point force = k * (r_ij.norm() - r0) * r_ij.unit();

  STOP_LOG("harmonic_force(&r_ij, &r0, &k)", "FixBase");
  return force;
}

// ======================================================================
Point FixBase::polymer_wall_empirical_force(const Point& r_ij,     // vector
                                                                   // from
                                                                   // particle
                                                                   // to wall,
                                                                   // r_j-r_i
                                            const Real & c0,       // constant
                                                                   // 1:
                                            const Real & d0) const // constant
                                                                   // 2:
{
  START_LOG("polymer_wall_empirical_force(pt_ij)", "FixBase");

  /*
   * f_i = c0*( 1 - y0/d0 )^2, must be > 0
   * Theoretically, we require that y0 > 0, and d0 > 0;
   */
  Point fij(0.);

  if (std::abs(r_ij.norm()) < std::abs(d0)) // if distance to wall is less than
                                            // cutoff radius
  {
    fij = -c0 * (1 - r_ij.norm() / d0) * (1 - r_ij.norm() / d0) * r_ij.unit();
  }
  STOP_LOG("polymer_wall_empiricalforce(pt_ij)", "FixBase");
  return fij;
}

// ======================================================================
Point FixBase::friction_force(const Point            & bead_1,
                              const Point            & bead_2,
                              const std::vector<Real>& v1,
                              const std::vector<Real>& v2,
                              const std::vector<Real>& fxv_12,
                              const Real             & Hf,
                              const Real             & dmin) const
{
  START_LOG("friction_force()", "FixBase");

  Point fij(0.);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     First, check the distance of these two beads, if they are far away from
     each other, there is no frcition force
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Point v_12 = bead_1 - bead_2;

  if (v_12.norm() > dmin) return fij;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     If these two beads are close enough, we compute their friction force
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Real f_12_norm = 0.0;

  for (std::size_t i = 0; i < 3; ++i)
  {
    v_12(i)    = v2[i] - v1[i];
    f_12_norm += fxv_12[i] * fxv_12[i];
  }
  f_12_norm = std::sqrt(f_12_norm);
  const Real v_12_norm = v_12.norm();
  fij = Hf * f_12_norm * v_12 / v_12_norm;

  STOP_LOG("friction_force()", "FixBase");
  return fij;
}
