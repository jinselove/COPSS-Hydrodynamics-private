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

#include "libmesh/reference_counted_object.h"



namespace libMesh
{
  
  /*
   * This is a base class to define the force field of immersed structures.
   * Users can define their own force field via derived classes!
   *
   */
  
class FixBase
{
public:
  // Constructor
  FixBase();
  
  
  // Destructor
  ~FixBase();
  
  
  // 0. define parameters of the force field
  // *** This is illegal in C++ (see Practical C++ programming 2nd ed P215)
  const Real PI  = libMesh::pi;  //3.1415926;
  const Real tol = 1E-3;
  
  
  /* --------------------------------------------------------------
   * The following forces are computed according to the PRL paper:
   * DNA dynamics in a microchannel, R.M. Jendrejack et al. (2003)
   * See its Ref. [38], [39], [40] therein, and 
   * Doyle&Underhill (2005) HandBook of Material modeling
   -------------------------------------------------------------- */
  
  
  /*
   * Compute the spring force between the bead i and j, the direction is from i to j
   * WLS:  Worm-like springs model. ( vector R_ij = Rj - Ri )
   * f_ij = c1*[ (1-|R_ij|/Ls)^(-2) -1 + 4*|R_ij|/Ls ] * R_ij/|R_ij|
   *
   * Originally, c1 = kB*T/(2bk), and Ls = q0 = N_ks*bk
   * If force is normalized by bead radius a, we have 
   * c1 = a/(2bk), and c2 = q0/a = Nks*bk/a
   */
  std::vector<Real> spring_force_wls(const Point& R_ij,      // direction vector
                                             const Real&  c1,        // constant 1
                                             const Real&  Ls) const; // constant 2
  
  
  /*
   * Compute the force on the bead i, whose direction is from i to j
   * FENE: finitely extensible nonlinear elastic spring model. ( vector R_ij = Rj - Ri )
   *                 |R_ij|/Ls
   * f_ij = c1*[ ------------------- ] * R_ij/|R_ij|
   *             1 - ( |R_ij|/Ls )^2
   *
   * Originally, c1 = 3*kB*T/(bk), and Ls = q0
   * If force is normalized by bead radius a, we have
   * c1 = 3a/bk, and c2 = q0/a = Nks*bk/a
   */
  std::vector<Real> spring_force_fene(const Point& R_ij,      // direction vector
                                              const Real&  c1,        // constant 1:
                                              const Real&  Ls) const; // constant 2:
  
  
  
  /*
   * Compute the force on the bead i, whose direction is from i to j
   * Underhill-Doyle(UD): spring model. ( vector R_ij = Rj - Ri )
   *
   * f_ij = [ a1*(1-r2)^-2 + a2*(1-r2)^-1 + a3 + a4*(1-r2) ] * R_ij/|R_ij|
   * where
   *     r2 = (|R_ij|/Ls)^2
   *     a1 ~ a4 can be found in Kounovsky-Shafer et al. Macromolecules(2013) 46 p.8356
   *     c1 input is chi = 1/Nks, 
   *     Ls is q0 (max spring extension) or normalized q0 = q0/a = Nks*bk/a
   */
  std::vector<Real> spring_force_ud(const Point& R_ij,      // direction vector
                                            const Real&  c1,        // constant 1:
                                            const Real&  Ls) const; // constant 2:
  
  
  /*
   * Compute the force on the bead i, whose direction is from i to j
   * LHS:  Linear Hookean Spring
   *
   * f_ij = k0*( |R_ij| - l0 )  * R_ij/|R_ij|
   */
  std::vector<Real> spring_force_lhs(const Point& R_ij,      // direction vector
                                             const Real&  l0,        // equilibrium distance
                                             const Real&  k0) const; // spring constant
  
  
  
  /*
   * Compute the force on the bead i, whose direction is from i to j
   * Dimensionless form:
   * potential: u_gaussian  = 1/2 * c1 * exp( -c2*|r_ij|^2 )
   * where r_ij is dimensionless : r_ij = R_ij / a and r_ij = rj - ri 
   * where c1 > 0, c2 > 0
   * force    : f_ij = -du_gaussian/dr_ij = -c1*c2* exp( -c2*|r_ij|^2 ) * r_ij
   * there is a negative sign before the force because r_ij = r_j - r_i, which has
   * an opposite direction as f_ij
   */
  std::vector<Real> gaussian_force(const Point& r_ij,      // direction vector
                                           const Real&  c1,        // constant 1: coefficient
                                           const Real&  c2) const; // constant 2: exp coef
  /*
   * Compute the force acting on bead i
   * Dimensionless form:
   * potential: u_lj = 4 * epsilon *[(sigma / |r_ij|)^12 - (sigma / |r_ij|)^6]
     // f_ij = -24 * epsilon * (2*(sigma/|r_ij|)^12 - (sigma/|r_ij|)^6 ) * r_ij / |r_ij|^2
   * where r_ij (vector) is dimensionless : r_ij = R_ij / a and r_ij = rj - ri
   * sigma is dimensionless : sigma = SIGMA / a
   * epsilon is dimensionless : epsilon = EPSILON / kBT
   * @ http://www.physics.buffalo.edu/phy516/jan28.pdf
   */
  std::vector<Real> lj_force(const Point& r_ij, // direction vector
                                         const Real& epsilon, // energy coefficient
                                         const Real& sigma) const; // distance coefficient

  /*
   * Compute the force acting on bead i
   * Dimensionless form
   * potential uij = 1/2 * k * (|r_ij| - |r0|)^2
   * force fij = k (|r_ij| - r0) * r_ij / |r_ij|
   * where r_ij (vector) is dimensionless: r_ij = R_ij / a and r_ij = rj - ri
   * there is no minus sign before fij because r_ij = rj - ri
   * r0 is equilibrim distance
   * k is dimensionless: k = K /kBT
   * @ Edmond Chow and Jeffrey Skolnick (2015), www.pnas.org/cgi/doi/10.1073/pnas.1514757112
   */
  std::vector<Real> harmonic_force(const Point& r_ij, // direction vector
                                           const Real& k, // equilibrium distance
                                           const Real& r0) const; // energy coefficient
  
  
  /*
   * Compute the force on the bead i, whose direction is from i to j
   * Particle-wall repulsive force
   * potential: Uw   = -c0/3 *dwall * ( 1 - y/dwall ), where y is the distance to a wall
   * force    : f_i = -dUw/dy = -c0*( 1 - y/dwall )^2 * r_ij.unit(), where c0 > 0 and r_ij = rj - ri
   * Jendrejack, R. M., Schwartz, D. C., Graham, M. D., & de Pablo, J. J. (2003). Effect of confinement on DNA dynamics in microfluidic devices. The Journal of Chemical Physics, 119(2), 1165–10. http://doi.org/10.1063/1.1575200
   */
  std::vector<Real> polymer_wall_empirical_force(const Point& r_ij,    // vector from particle to wall (or particle i to particle j, r_j-r_i)
                                                const Real&  c0,        // constant 1:
                                                const Real&  d0) const; // constant 2:
  
  
  
  /*
   * Compute the friction force between two moving beads using a Coulombic law.
   * Ref:
   * Londono-Hurtado et al. J Reinforced Plastic&Composites 30(9) 781-790(2011)
   */
  std::vector<Real> friction_force(const Point& bead_1,       // position of bead 1
                                           const Point& bead_2,       // position of bead 2
                                           const std::vector<Real>& v1, // velocity of bead 1
                                           const std::vector<Real>& v2, // velocity of bead 2
                                           const std::vector<Real>& fxv_12, // excluded vol force
                                           const Real& Hf,            // friction coef
                                           const Real& dmin) const;   // minimum distance
  
  
  
  
  /*
   * This virtual function needs to be defined by the user in the derived class.
   */
//  static void reinit_force_field() = 0;

  
  

  
};  // end of class
  

} // end of namespace


