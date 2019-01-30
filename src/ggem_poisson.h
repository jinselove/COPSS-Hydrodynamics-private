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


namespace libMesh
{
  
/*
 * Remember that:
 * The PMSystemPoisson is designed to solve
 * the Poisson equation with regularized Gaussian
 * point charges due to particles using mixed FEM, which 
 * is usually known as 'global' part of the solution
 *
 * This class GGEMPoisson will provide the main functionalities of
 * GGEM (General Geometry Ewald-like Method), and it will provide
 * the 'local' solution. Together with 'PMSystemPoisson', this will
 * provide the complete solution for the Poisson system with
 * multiple point charges.
 *
 *
 * The current Green's functions are only for 3D cases!
 * Green's function for 2D situation is different.
 
 */
  
// class GGEMPoisson :  public ReferenceCountedObject<GGEMPoisson>
class GGEMPoisson : public GGEMSystem
{
public:

    // Constructor
    GGEMPoisson();


    // Destructor
    ~GGEMPoisson();
    

    /*
    * Modified smoothed Dirac delta function with exponent form
    * Reference: equation (13) in JCP, 143. 014108 (2015)
    */
    Real smoothed_charge_exp(const Real& r) const;


    /*
    * Green tensors for different types of Dirac Delta function
    * if delta_type==SINGULAR: will return this -> green_tensor_unbounded_singular
    * if delta_type==SMOOTHED_EXP: will return this -> green_tensor_unbounded_smoothed
    */
    Real green_tensor_unbounded(const Point& x,     /* vector x = pt1 - pt0 */
                                const DeltaFunType& delta_type) const;


    /*
    * Free-space Green's function of the Poisson equation if the charge density 
    * is written as dirac function, i.e., rho(x) = sum[v=1:N] Dirac(x, x_v) 
    * This tensor is singular at x = x_v.
    * G = 1/(4*pi*epsilon*epsilon_0) * 1 / r
    * The dimensionless unit of G is 1 / (4*pi*epsilon*epsilon_0)
    * In dimensionless form, G = 1 / r, the later r is the first r non-dimensionalized 
    * by bead radius a. The resulting potential is non-dimensionalized
    * by e / (4 * pi * epsilon * epsilon_0 * a)
    */
    Real green_tensor_unbounded_singular(const Point& x /* vector x = pt1 - pt0 */
                                                        ) const;  
                                        
                                        
    /*
     * Green tensor for unbounded domain due to a smoothed charge defined by a Gaussian
     * charge density: rho = sum[v=1:N] z_v *( g ).
     * If alpha_or_ksi == ksi, then this tensor gives the solution to the total local_velocity
     * at the boundary of a periodic box.
     * Reference: The ksi part of Equation (17) in Membrane paper, Jarol E Molina, J de Pablo, J.Hernandez-Ortiz
     * The reason why we only take the ksi part is because Equation (17) is the local part of the green
     * tensor, i.e, the green tensor because of the local charge, i.e., g_ski - g_alpha. However, here we need
     * the total green tensor because of the total charge, i.e., g_ksi 
     * The dimensionless unit of G is 1 / (4*pi*epsilon*epsilon_0)
     * In dimensionless form, G = 1 / r, the later r is the first r non-dimensionalized 
     * by bead radius a. The resulting potential is non-dimensionalized
     * by e / (4 * pi * epsilon * epsilon_0 * a)
     */
     Real green_tensor_unbounded_smoothed(const Point& x,     /* vector x = pt1 - pt0 */
                                         const Real& alpha_or_ksi
                                         ) const;   


    /*
    * Green tensor for the local charge density. The charge density is written as a singular
    * point charge (delta) - a "local" smoothed charge density (g)
    * with exp form: rho = sum[v=1:N] f_v*( delta - g ).
    * The solution is u_loc(i) = sum[v=1:N] G_v(x-x_v; i,j)*f_v(j)
    * Reference: Equation (14) in J. Chem. Theory Comput. 2018, 14, 4901-4913
    * The dimensionless unit of G is 1 / (4*pi*epsilon*epsilon_0)
    */
    Real green_tensor_local_singular(const Point& x     /* vector x = pt1 - pt0 */
                                    ) const;    


    /*
    * Regularized Green function in order to remove the singularity of v
    * For ksi^-1 = a/3. The charge is distributed through a gaussian function with a
    * valence relative to the particle/bead diameter. 
    * Reference: Equation (17) in Membrane paper, Jarol E Molina, J de Pablo, J.Hernandez-Ortiz
    * The dimensionless unit of G is 1 / (4*pi*epsilon*epsilon_0)
    */
    Real green_tensor_local_regularized(const Point& x     /* vector x = pt1 - pt0 */
                                       ) const; 


    /* Compute the local electrostatic potential field at a given point ptx
    * due to smoothed/regularized point charges.
    * charge_type: regularized
    *
    * NOTE: due to the fast convergence of gaussian function, only a small group of 
    * particles within the neighbor list are considered. There are two ways
    * to construct this neighbor list:
    * (1) element independent: directly search particles near the given point using KDTree;
    * (2) element dependent: directly use the neighbor list of the parent element;
    * this function implement method (1), which contains a short neighbor list
    *
    * Eqn (16) in Membrane paper, Jarol E Molina, J de Pablo, J.Hernandez-Ortiz
    */
    Real local_potential_field(PointMesh<3>*  point_mesh,
                               const Point&   ptx,      /* a pt in space */
                               const std::string& charge_type) const;

    /* Compute the local electrostatic potential field at a given point ptx
    * due to smoothed/regularized point charges.
    * charge_type: regularized
    *
    * NOTE: due to the fast convergence of gaussian function, only a small group of 
    * particles within the neighbor list are considered. There are two ways
    * to construct this neighbor list:
    * (1) element independent: directly search particles near the given point using KDTree;
    * (2) element dependent: directly use the neighbor list of the parent element;
    * this function implement method (2), which contains a short neighbor list
    *
    * Eqn (16) in Membrane paper, Jarol E Molina, J de Pablo, J.Hernandez-Ortiz
    */
    Real local_potential_field(PointMesh<3>*  point_mesh,
                             const Elem* elem,
                             const Point&   ptx,      /* a pt in space */
                             const std::string& charge_type) const;


    /*
    * Compute the local electrostatic potential at a point/bead with point_id = pid0.
    * charge_type: regularized
    *
    * Eqn (16) in Membrane paper, Jarol E Molina, J de Pablo, J.Hernandez-Ortiz
    * This potential should be include the potential field induced by the bead itself.
    */
    Real local_potential_bead(PointMesh<3>*  point_mesh,
                              const std::size_t& pid0,  /* point id */
                              const std::string& charge_type) const;


    //! set PointType 
    void set_point_type(const PointType& _point_type) {point_type = _point_type;};
    PointType get_point_type() {return point_type;};


    //! set GGEM alpha
    void set_alpha(const Real& _alpha);
    Real get_alpha() {return alpha;};


    //! set regularization parameter
    void set_ksi();
    Real get_ksi() {return ksi;};


    //! set bead radius (non-dimensional), it should be 1.
    void set_br0(const Real& _br0) {br0 = _br0;};
    Real get_br0() {return br0;};

    //! set phi0
    void set_phi0(const Real& _phi0) {phi0 = _phi0;};
    Real get_phi0() {return phi0;};

  
private:

  //! alpha
  Real alpha2 = NAN; // alpha2 = alpha * alpha
  Real alpha3 = NAN; // alpha3 = alpha * alpha * alpha
  Real alpha3_pi_23 = NAN; // alpha3_pi_23 = alpha^3. / pi^(3./2.)
  Real two_alpha_sqrt_pi = NAN; // two_alpha_sqrt_pi = 2. * alpha / pi^(1./2.)
  
  //! regularization parameter
  Real ksi2 = NAN; // ksi2 = ksi * ksi
  Real ksi3 = NAN; // ksi3 = ksi^3.
  Real ksi3_pi_23 = NAN; // ksi3_pi_23 = ksi^3. / pi^(3./2.)
  Real two_ksi_sqrt_pi = NAN; // two_ksi_sqrt_pi = 2. * ksi / pi^(1./2.)
  Real coeff_1 = NAN; // coeff_1 = 2 / pi^(1./2.) * (ksi - alpha)
    
  //! potential field dimensionless unit, i.e., the value of e / (4*pi*epsilon*epsilon_0) 
  Real phi0 = NAN;
};

  
} // end of namespace
