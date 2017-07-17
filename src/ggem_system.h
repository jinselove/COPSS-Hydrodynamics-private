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



#include <stdio.h>
#include <cmath>

#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/point.h"


#include "point_mesh.h"
//#include "particle_mesh.h"


using namespace libMesh;



namespace libMesh
{
  
  
  /*
   * Type of a Dirac Delta function(for a point force)
   */
  enum DeltaFunType {
    SINGULAR      = 0,    // Original Dirac delta function with singularity
    SMOOTHED_EXP  = 1,    // Smoothed force(Exponent/Gaussian form)
    SMOOTHED_POL1 = 2,    // Smoothed force(Polynomial type 1)
    SMOOTHED_POL2 = 3     // Smoothed force(Polynomial type 2)
  };
  
  
  
  
/*
 * Remember that:
 * The PMLinearImplicitSystem is designed to solve
 * the Stokes equation with regularized Gaussian
 * point forces due to particles using mixed FEM, which 
 * is usually known as 'global' part of the solution
 *
 * This class GGEMSystem will provide the main functionalities of
 * GGEM (General Geometry Ewald-like Method), and it will provide
 * the 'local' solution. Together with 'PMLI_system', this will
 * provide the complete solution for the Stokes system with
 * multiple point forces.
 *
 *
 * The current Green's functions are only for 3D cases!
 * Green's function for 2D situation is different.
 */
  
class GGEMSystem :  public ReferenceCountedObject<GGEMSystem>
{
public:

  // Constructor
  GGEMSystem();
  
  
  // Destructor
  ~GGEMSystem();
  
  
  // PI constant
  const Real  PI      = libMesh::pi;
  const Real  sqrt_pi = std::sqrt(PI);
  const Real pi_23 = std::pow(PI, 3./2.);
  const Real  r_eps   = 1E-6;         // a small value
  
  
  /* 
   * Kronecker delta function 
   */
  inline Real kronecker_delta(const std::size_t i,
                              const std::size_t j) const
  { return  i==j ? 1.0 : 0.0; };
  

  /*
   * Determine the regularization parameter according to the point type:
   * regularization parameter ksi=sqrt(PI)/(3*R0) for a polymer bead;
   * and we ksi ~ hc^-1 for tracking points in IBM
   * See Phys. Fluids 22, 123103(2010) Pranay et al., in which they took (ksi*hc)^-1=0.75
   */
  Real regularization_parameter(const Real& hc,
                                const Real& R0,
                                const PointType& point_type) const;
  
  
  
  /*
   * Modified smoothed Dirac delta function with exponent form
   * 'alpha' is a regularization parameter
   */
  Real smoothed_force_exp(const Real& r,
                          const Real& alpha) const;


  
  /*
   * Green tensors for different types of Dirac Delta function
   */
  DenseMatrix<Number> green_tensor(const Point& x,     /* vector x = pt1 - pt0 */
                                   const Real& alpha,  /* smoothing parameter  */
                                   const Real& mu,     /* kinematic viscosity  */
                                   const std::size_t& dim,           /* dim==3 */
                                   const bool& zero_limit,           /* x-->0  */
                                   const DeltaFunType& delta_type) const;
  
  
  
  /*
   * Original Green's function for a point force, which is singular at x = 0.
   */
  DenseMatrix<Number> green_tensor_exact(const Point& x, /* vector x = pt1 - pt0 */
                                         const Real& mu,  /* kinematic viscosity */
                                         const std::size_t& dim) const;  /*dim==3*/
  
  
  
  /* 
   * Green function for unbounded Stokes due to a 'local' smoothed force with exp form
   * rho = sum[v=1:N] f_v*( g ).
   * The solution is u_loc(i) = sum[v=1:N] G_v(x-x_v; i,j)*f_v(j)
   * Eqn (2) - (3) in PRL 98, 140602(2007), JP Hernandez-Ortiz et al.
   */
  DenseMatrix<Number> green_tensor_exp(const Point& x,     /* vector x = pt1 - pt0 */
                                       const Real& alpha,  /* alpha parameter      */
                                       const Real& mu,     /* kinematic viscosity  */
                                       const std::size_t& dim,           /* dim==3 */
                                       const bool& zero_limit) const;    /* x-->0  */
  
  
  /*
   * Green function for unbounded Stokes due to a point force + a 'local' smoothed force
   * with exp form: rho = sum[v=1:N] f_v*( delta - g ).
   * The solution is u_loc(i) = sum[v=1:N] G_v(x-x_v; i,j)*f_v(j)
   * Eqn (3) in PRL 98, 140602(2007), JP Hernandez-Ortiz et al.
   * or Eqn (30) in J Chem Phys. 136, 014901(2012), Yu Zhang, de Pablo and Graham.
   */
  DenseMatrix<Number> green_tensor_local(const Point& x,     /* vector x = pt1 - pt0 */
                                         const Real& alpha,  /* alpha parameter      */
                                         const Real& mu,     /* kinematic viscosity  */
                                         const std::size_t& dim,           /* dim==3 */
                                         const bool& zero_limit) const;    /* x-->0  */
  
  
  /*
   * Regularized Green function in order to remove the singularity of v
   * For ksi^-1 = 3a/sqrt(PI), the max fluid velocity is equal to that
   * of a particle with radius a and the pair mobility remains positive definite.
   *
   * See Eqn (4) in PRL 98, 140602(2007), JP Hernandez-Ortiz et al.
   * or Eqn (31)(32) in J Chem Phys. 136, 014901(2012), Yu Zhang, de Pablo and Graham.
   */
  DenseMatrix<Number> green_tensor_regularized(const Point& x,     /* vector x = pt1 - pt0 */
                                               const Real& alpha,  /* alpha parameter      */
                                               const Real& mu,     /* kinematic viscosity  */
                                               const Real& ksi,    /* regularization parameter */
                                               const std::size_t& dim,          /* dim==3 */
                                               const bool& zero_limit) const;   /* x-->0  */
  
  
  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Difference between green_tensor_local(Gl) and green_tensor_regularized(Gr):
   G_reg(x) = G_exp(x,ksi) - G_exp(x,alpha)
   G_loc(x) = G_exact(x) - G_exp(x,alpha)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  DenseMatrix<Number> green_tensor_diff(const Point& x,     /* vector x = pt1 - pt0 */
                                        const Real& alpha,  /* alpha parameter      */
                                        const Real& mu,     /* kinematic viscosity  */
                                        const Real& ksi,    /* regularization parameter */
                                        const std::size_t& dim,   /* dim==3 */
                                        const bool& zero_limit,   /* x-->0  */
                                        const DeltaFunType& delta_type_a, /* delta type of alpha */
                                        const DeltaFunType& delta_type_k) const; /* delta type of ksi */
  
  

  /* Compute the local vel-solution of the fluid at a given point ptx
   * due to smoothed/regularized point forces.
   * force_type: regularized, smooth
   *
   * NOTE: due to the fast convergence of gaussian function, only a small group of 
   * particles within the neighbor list are considered. There are two ways
   * to construct this neighbor list:
   * (1) element independent: directly search particles near the given point using KDTree;
   * (2) element dependent: directly use the neighbor list of the parent element;
   * this function implement method (1), which contains a short neighbor list
   *
   * Eqn (33) in J Chem Phys. 136, 014901(2012), Yu Zhang, de Pablo and Graham.
   */
  std::vector<Real> local_velocity_fluid(PointMesh<3>*  point_mesh,
                                         const Point&   ptx,      /* a pt in space */
                                         const Real&    alpha,    /* alpha parameter */
                                         const Real&    mu,       /* kinematic viscosity */
                                         const Real&    br0,      /* normalized bead radius */
                                         const Real&    hs,       /* mesh size of solids */
                                         const std::size_t& dim,  /*dim==3*/
                                         const std::string& force_type) const;
  
  
  /*
   * Compute the local velocity of a point/bead with point_id = pid0.
   * force_type: regularized, smooth
   *
   * Eqn (34)&(35) in J Chem Phys. 136, 014901(2012), Yu Zhang, de Pablo and Graham.
   */
  std::vector<Real> local_velocity_bead(PointMesh<3>*  point_mesh,
                                        const std::size_t& pid0,  /* point id */
                                        const Real&    alpha,     /* alpha parameter */
                                        const Real&    mu,        /* kinematic viscosity */
                                        const Real&    br0,       /* normalized bead radius */
                                        const Real&    hs,        /* mesh size of solids*/
                                        const std::size_t& dim,   /*dim==3*/
                                        const std::string& force_type) const;
  
  
  
  /* Compute the local v-solution at a point ptx in an unbounded domain
   * due to smoothed/regularized point forces.
   * NOTE: due to the fast convergence of gaussian function, only a small group of
   * particles within the neighbor list are considered. There are two ways
   * to construct this neighbor list:
   * (1) element independent: directly search particles near the given point using KDTree;
   * (2) element dependent: directly use the neighbor list of the parent element;
   * this function implement method (2), which contains a longer neighbor list
   */
//  std::vector<Real> local_velocity(ParticleMesh<3>* particle_mesh,
//                                   const Elem* elem,   /* the parent element */
//                                   const Point& ptx,   /* a pt in space */
//                                   const Real& alpha,  /* alpha parameter */
//                                   const Real& mu,     /* kinematic viscosity */
//                                   const Real& ksi,    /* ksi */
//                                   const std::size_t& dim, /*dim==3*/
//                                   const std::string& option);
  
  
  
  /*
   * Self-exclusion term for the GLOBAL velocity at the i-th bead
   */
  std::vector<Real> global_self_exclusion(PointMesh<3>* point_mesh,
                                          const std::size_t&  pid0,
                                          const Real& alpha,
                                          const Real& mu,
                                          const std::size_t& dim) const;
  
  
  
  /*
   * Rotne-Prager-Yamakawa(RPY) tensor
   * Reference: Jendrejack, Schwartz, de Pablo and Graham, J Chem Phys (2003)
   * Eqn (5) - (9)
   */
  DenseMatrix<Number> rpy_tensor(const Point& x,     /* vector x = pt1 - pt0 */
                                 const Real& mu,     /* kinematic viscosity */
                                 const Real& a,      /* bead radius */
                                 const std::size_t& dim) const; /*dim==3*/
  
  

  /*
   * Mobility tensor in an infinite domain(no walls) using RPY tensor
   * Reference: Jendrejack, Schwartz, de Pablo and Graham, J Chem Phys (2003)
   * Eqn (3) & (5)
   * Note, diffusion tensor eqn(3) = this mobility tensor x kB*T.
   */
  DenseMatrix<Number> mobility_tensor(const Point& x,     /* vector x = pt1 - pt0 */
                                      const Real& mu,     /* kinematic viscosity */
                                      const Real& a,      /* bead radius */
                                      const std::size_t& dim) const; /*dim==3*/
  
  
};

  
} // end of namespace
