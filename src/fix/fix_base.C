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
//FixBase::FixBase(PMLinearImplicitSystem& pm_sys)
//: _pm_system(pm_sys)
//{
//  // initialize the private memebers
//}



// ======================================================================
FixBase::~FixBase()
{
  // do nothing
}



// ======================================================================
std::vector<Real> FixBase::spring_force_wls(const Point& pt_ij,
                                                   const Real&  c1,
                                                   const Real&  Ls) const
{
  START_LOG ("spring_force_wls(pt_ij)", "FixBase");

  std::vector<Real> fij(3);
    
  // The spring force: pt_ij = ptj - pti(with periodicity), R_ij = |pt_ij|
  const Real  R_ij  = pt_ij.norm();
  const Real ratio  = R_ij/Ls;
  Real        val0  = 1.0 - ratio;
  
  // avoid singularity
  if( std::abs(val0) <  tol) val0 = val0/std::abs(val0)*tol;
  
  // compute force
  const Real  tmp   = 1./( val0*val0 );
  const Point F_ij  = c1*( tmp - 1. + 4.*ratio )*pt_ij/R_ij;
  for (std::size_t j=0; j<3; ++j){
    fij[j] = F_ij(j);
  }
  
  STOP_LOG ("spring_force_wls(pt_ij)", "FixBase");
  return fij;
}



// ======================================================================
std::vector<Real> FixBase::spring_force_fene(const Point& pt_ij,
                                                    const Real&  c1,
                                                    const Real&  Ls) const
{
  START_LOG ("spring_force_fene(pt_ij)", "FixBase");
  std::vector<Real> fij(3);  
  
  // The spring force: pt_ij = ptj - pti(with periodicity), R_ij = |pt_ij|
  const Real  R_ij  = pt_ij.norm();
  const Real ratio  = R_ij/Ls;
  Real        val0  =  1.0 - ratio*ratio;
  
  // avoid singularity
  if( std::abs(val0) <  tol) val0 = val0/std::abs(val0)*tol*tol;
  
  // compute force
  const Point F_ij  = c1*ratio/val0 *pt_ij/R_ij;
  for (std::size_t j=0; j<3; ++j){
    fij[j] = F_ij(j);
  }
  
  STOP_LOG ("spring_force_fene(pt_ij)", "FixBase");
  return fij;
}



// ======================================================================
std::vector<Real> FixBase::spring_force_ud(const Point& pt_ij,
                                                  const Real&  c1,
                                                  const Real&  Ls) const
{
  START_LOG ("spring_force_ud(pt_ij)", "FixBase");
  
  std::vector<Real> fij(3);
  // The spring force: pt_ij = ptj - pti(with periodicity), R_ij = |pt_ij|
  const Real  R_ij  = pt_ij.norm(); // vector length
  const Real  ratio = R_ij/Ls;
  const Real  r2    = ratio*ratio;
  const Real  c2    = c1*c1;
  
  const Real a1 = 1.0;
  const Real a2 = -7.0*c1;
  const Real a3 = 3.0/32.0 - 0.75*c1 - 6.*c2;
  const Real a4 = (13.0/32.0 + 0.8172*c1 - 14.79*c2)/(1.0 - 4.225*c1 + 4.87*c2);
  
  
  const Real T4 = (1.0 - r2);
  const Real T2 = 1.0/T4;
  const Real T1 = T2/T4;  // = 1/(T4*T4)
  const Point F_ij  = (a1*T1 + a2*T2 + a3 + a4*T4) * pt_ij/R_ij;
  for (std::size_t j=0; j<3; ++j){
    fij[j] = F_ij(j);
  }
  
  STOP_LOG ("spring_force_ud(pt_ij)", "FixBase");
  return fij;
}



// ======================================================================
std::vector<Real> FixBase::spring_force_lhs(const Point& pt_ij,
                                                   const Real&  l0,
                                                   const Real&  k0) const
{
  START_LOG ("spring_force_lhs(pt_ij)", "FixBase");
  
  //f_ij = k0*( |R_ij| - l0 )  * R_ij/|R_ij|
  std::vector<Real> fij(3);

  // The spring force: pt_ij = ptj - pti(with periodicity), R_ij = |pt_ij|
  const Real  R_ij  = pt_ij.norm();
  if(std::abs(R_ij)<1E-6 )
  {
    printf("--->test in ForceField::spring_force_lhs(): Rij = %f\n",R_ij);
  }
  
  const Real  F0   = k0*( R_ij - l0 );  // force magnitude
  const Point F_ij = F0*pt_ij/R_ij;     // force direction
  for (std::size_t j=0; j<3; ++j){
    fij[j] = F_ij(j);
  }
  
  STOP_LOG ("spring_force_lhs(pt_ij)", "FixBase");
  return fij;
}



// ======================================================================
std::vector<Real> FixBase::gaussian_force(const Point& r_ij,
                                                        const Real&  c1,
                                                        const Real&  c2) const
{
  START_LOG ("gaussian_force(pt_ij)", "FixBase");
 
  // f_ij = c1*c2* exp( -c2*|r_ij|^2 ) * r_ij
  std::vector<Real> fij(3);  
  // The spring force: pt_ij = ptj - pti(with periodicity), R_ij = |pt_ij|
  const Real  r_ij_size = r_ij.norm();
  const Point force = -c1*c2*std::exp(-c2*r_ij_size*r_ij_size)*r_ij; 
  for (std::size_t j=0; j<3; ++j){
    fij[j] = force(j);
  }
  
  STOP_LOG ("gaussian_force(pt_ij)", "FixBase");
  return fij;
}

// ======================================================================
std::vector<Real> FixBase::lj_force(const Point& r_ij, // direction vector
                                         const Real& epsilon, // energy coefficient
                                         const Real& sigma) const // distance coefficient
{
  START_LOG("lj_force(&r_ij, &epsilon, &sigma)","FixBase");
  // f_ij = -24 * epsilon * (2*(sigma/|r_ij|)^12 - (sigma/|r_ij|)^6 ) * r_ij / |r_ij|^2
  std::vector<Real> fij(3);
 
  Real r = r_ij.norm();
  Real r2 = r * r;
  Real sr2 = sigma*sigma / r2;
  Real sr6 = sr2 * sr2 * sr2;
  Real sr12 = sr6 * sr6;
  const Point force = - 24. * epsilon * (2.* sr12 - sr6) / r2 * r_ij ;
  for (std::size_t j=0; j<3; ++j){
    fij[j] = force(j);
  }
  STOP_LOG("lj_force(&r_ij), &epsilon, &sigma", "FixBase");
  return fij;
}

// ======================================================================
std::vector<Real> FixBase::harmonic_force(const Point& r_ij,
                                                 const Real& k,
                                                 const Real& r0) const
{
  START_LOG("harmonic_force(&r_ij, &r0, &k)","FixBase");

  std::vector<Real> fij(3);

  const Point force = k * (r_ij.norm() - r0) * r_ij.unit();
  for (std::size_t j=0; j<3; ++j){
    fij[j] = force(j);
  }
  STOP_LOG("harmonic_force(&r_ij, &r0, &k)","FixBase");
  return fij;
}

// ======================================================================
std::vector<Real> FixBase::polymer_wall_empirical_force(const Point&  r_ij,   //vector from particle to wall, r_j-r_i
                                                      const Real&  c0,        // constant 1:
                                                      const Real&  d0) const // constant 2:
{
  START_LOG ("polymer_wall_empirical_force(pt_ij)", "FixBase");
  
  /* 
   * f_i = c0*( 1 - y0/d0 )^2, must be > 0
   * Theoretically, we require that y0 > 0, and d0 > 0;
   */
  std::vector<Real> fij(3);
  if( std::abs(r_ij.norm()) < std::abs(d0) ) // if distance to wall is less than cutoff radius
  {
    const Point f_i = -c0*( 1 - r_ij.norm()/d0 )*( 1 - r_ij.norm()/d0 ) * r_ij.unit();
    for(std::size_t _dim = 0; _dim < 3; _dim++){
      fij[_dim] = f_i(_dim);
    }
  }  
  STOP_LOG ("polymer_wall_empiricalforce(pt_ij)", "FixBase");
  return fij;
}



// ======================================================================
std::vector<Real> FixBase::friction_force(const Point& bead_1,
                                                 const Point& bead_2,
                                                 const std::vector<Real>& v1,
                                                 const std::vector<Real>& v2,
                                                 const std::vector<Real>& fxv_12,
                                                 const Real& Hf,
                                                 const Real& dmin) const
{
  START_LOG ("friction_force()", "FixBase");

  std::vector<Real> fij(3);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   First, check the distance of these two beads, if they are far away from
   each other, there is no frcition force
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Point v_12 = bead_1 - bead_2;
  if( v_12.norm() > dmin ) return fij;
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   If these two beads are close enough, we compute their friction force
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Real f_12_norm = 0.0;
  for(std::size_t i=0; i<3; ++i)
  {
    v_12(i) = v2[i] - v1[i];
    f_12_norm += fxv_12[i]*fxv_12[i];
  }
  f_12_norm = std::sqrt(f_12_norm);
  const Real v_12_norm = v_12.norm();
  
  for(std::size_t i=0; i<3; ++i)
  {
    fij[i] = Hf*f_12_norm*v_12(i)/v_12_norm;
  }
  
  STOP_LOG ("friction_force()", "FixBase");
  return fij;
}







//// ======================================================================
//std::vector<Real> FixBase::particle_wall_force(const Point& pti,       // bead i
//                                                      const Point& dist,      // distance to the wall
//                                                      const Real&  c0,        // constant 1:
//                                                      const Real&  d0) const  // constant 2:
//{
//  START_LOG ("particle_wall_force(pti, dist)", "FixBase");
//  
//  //  retrieve the box boundary
//  const Point box_min = _pm_system.point_mesh()->pm_periodic_boundary()->box_min();
//  const Point box_max = _pm_system.point_mesh()->pm_periodic_boundary()->box_max();
//  
//  // compute the particle-wall interaction force
//  const Real bead_r0 = 1.0;
//  std::vector<Real> fij(_dim,0.0);
//  for (std::size_t j=0; j<_dim; ++j)  // loop over each direction
//  {
//    // compute the distance to the wall if NO periodic boundary.
//    bool periodic = _pm_system.point_mesh()->pm_periodic_boundary()->periodic_direction(j);
//    if ( periodic==false )
//    {
//      Real dist_wall_a = pti(j) - box_min(j);
//      Real dist_wall_b = pti(j) - box_max(j);
//      
//      // FIXME, if dist_wall==0, fix it! because it has no physical meaning.
//      //if(dist_wall_a<bead_r0)  dist_wall_a = +bead_r0;
//      //if(dist_wall_b>-bead_r0) dist_wall_b = -bead_r0;
//      const Real dist_wall = std::abs(dist_wall_a)<std::abs(dist_wall_b) ? dist_wall_a : dist_wall_b;
//      
//      // It is non-zero only near a thin layer near the wall.
//      if ( std::abs(dist_wall) < d0 )
//      {
//        const Real unit_v = dist_wall/std::abs(dist_wall);
//        const Real tmp = (1.0 - dist_wall/d0);
//        fij[j] = c0*tmp*tmp*unit_v;
//      } // end if
//      else
//      {
//        fij[j] = 0.0;
//      } // end if-else
//      
//      // if the particle is right on the centraline between the wall,
//      // the interaction force along this direction cancels each other and it is 0.
//      if ( std::abs(dist_wall_a + dist_wall_b) < 1E-8)
//        fij[j] = 0.0;
//      
//    } // end if
//    
//  } // end for j-loop
//  
//  STOP_LOG ("particle_wall_force(pti, dist)", "FixBase");
//  return fij;
//}

