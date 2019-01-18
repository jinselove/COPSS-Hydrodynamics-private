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



namespace libMesh
{


// ======================================================================
GGEMStokes::GGEMStokes()
{
  _kronecker_delta.resize(3, std::vector<Real>(3,0.));
  for(int i=0; i<3; i++) _kronecker_delta[i][i] = 1.;
  // do nothing
}


// ======================================================================
GGEMStokes::~GGEMStokes()
{
  // do nothing
}



// ======================================================================
Real GGEMStokes::regularization_parameter(const Real& hc,
                                          const Real& ibm_beta,
                                          const Real& R0,
                                          const PointType& point_type) const
{
  START_LOG ("regularization_parameter()", "GGEMStokes");
  
  Real val = 0.0;
  switch (point_type)
  {
    case POLYMER_BEAD:      // Case I: bead type point
    {
      val = sqrt_pi/(3.*R0);
      break;
    }
    case LAGRANGIAN_POINT:  // Case II: Lagrangian point of immersed body
    {
      val = 1.0/(ibm_beta*hc);  // default value of ibm_beta is 0.75, which needs to be optimized for different number of nodes
      break;
    }
    case POINT_PARTICLE:    // Case III: point particle
    {
      val = sqrt_pi/(3.*R0);
      break;
    }
    default:
    {
      printf("*** error in GGEMStokes::regularization_parameter(): Undefined PointType!\n");
      libmesh_error();
    }
  }
  
  STOP_LOG ("regularization_parameter()", "GGEMStokes");
  //printf("--->test in regularization_parameter(): ksi = %f\n", val);
  return val;
}



// ======================================================================
Real GGEMStokes::smoothed_force_exp(const Real& r,
                                         const Real& alpha) const
{
  START_LOG ("smoothed_force_exp()", "GGEMStokes");
  
  const Real a2   = alpha*alpha,  a3 = a2*alpha;
  const Real r2   = r*r;
  const Real a2r2 = a2*r2;
  Real  force = ( a3/(pi_23) )*(2.5 - a2r2)*std::exp( -a2r2 );
  
  STOP_LOG ("smoothed_force_exp()", "GGEMStokes");
  
  return force;
}


  
// ======================================================================
DenseMatrix<Number> GGEMStokes::green_tensor(const Point& x,
                                             const Real& alpha,
                                             const Real& mu,
                                             const std::size_t& dim,
                                             const bool& zero_limit,
                                             const DeltaFunType& delta_type) const
{
  START_LOG ("green_tensor()", "GGEMStokes");
  
  // Compute the values according to the delta_type
  DenseMatrix<Number> GT(dim,dim);
  switch (delta_type)
  {
    case SINGULAR :
      GT = this->green_tensor_exact(x, mu, dim);
      break;
    case SMOOTHED_EXP :
      GT = this->green_tensor_exp(x, alpha, mu, dim, zero_limit);
      break;
    default:
      printf("*** error in GGEMStokes::green_tensor: undefined DeltaFunType!\n");
      libmesh_error();
  }
  
  STOP_LOG ("green_tensor()", "GGEMStokes");
  return GT;
}
  
  
  
// ======================================================================
DenseMatrix<Number> GGEMStokes::green_tensor_exact(const Point& x,
                                                   const Real& mu,
                                                   const std::size_t& dim) const
{
  START_LOG ("green_tensor_exact()", "GGEMStokes");
  
  const Real  r = x.norm(), r2  = r*r;
  const Real c0 = 1.0/(8.*PI*mu*r);
  
  // output warning if x->0
  if(r < r_eps) {
    printf("*** warning in GGEMStokes::green_tensor_exact, \
           r->0 which may lead to singular solution!\n");
  }
  
  // init G tensor
  DenseMatrix<Number> G(dim,dim);
  for (std::size_t i=0; i<dim; ++i){
    for (std::size_t j=0; j<dim; ++j){
      G(i,j) = c0*( _kronecker_delta[i][j] + x(i)*x(j)/r2 );
    }
  }
  
  STOP_LOG ("green_tensor_exact()", "GGEMStokes");
  
  // done
  return G;
}
  
  
  

// ======================================================================
DenseMatrix<Number> GGEMStokes::green_tensor_exp(const Point& x,
                                                 const Real& alpha,
                                                 const Real& mu,
                                                 const std::size_t& dim,
                                                 const bool& zero_limit) const
{
  START_LOG ("green_tensor_exp()", "GGEMStokes");
  
  const Real c0 = 1.0/(8.*PI*mu);
  DenseMatrix<Number> G(dim,dim);
  
  // check if zero limit is true.
  if (zero_limit)
  {
    for (std::size_t i=0; i<dim; ++i)
      for (std::size_t j=0; j<dim; ++j)
        G(i,j) = c0*4.0*alpha/sqrt_pi*_kronecker_delta[i][j];
  } // end if()
  else
  {
    // coefficients
    const Real  r = x.norm(), r2  = r*r, a2 = alpha*alpha;
    const Real c1 = std::erf(alpha*r)/r;
    const Real c2 = 2.*alpha/sqrt_pi*std::exp(-a2*r2);
    
    // compute G tensor
    for (std::size_t i=0; i<dim; ++i)
    {
      for (std::size_t j=0; j<dim; ++j)
      {
        Real tmp1 = ( _kronecker_delta[i][j] + x(i)*x(j)/r2 )*c1;
        Real tmp2 = ( _kronecker_delta[i][j] - x(i)*x(j)/r2 )*c2;
        G(i,j) = c0*( tmp1 + tmp2 );
      } // end for j-loop
    } // end for i-loop
    
  } // end if-else
  
  STOP_LOG ("green_tensor_exp()", "GGEMStokes");
  
  // done
  return G;
}


  
// ======================================================================
DenseMatrix<Number> GGEMStokes::green_tensor_local(const Point& x,
                                                   const Real& alpha,
                                                   const Real& mu,
                                                   const std::size_t& dim,
                                                   const bool& zero_limit) const
{
  START_LOG ("green_tensor_exp()", "GGEMStokes");
  
  const Real c0 = 1.0/(8.*PI*mu);
  DenseMatrix<Number> G(dim,dim);
  
  // check if zero limit is true.
  if (zero_limit)
  {
    for (std::size_t i=0; i<dim; ++i)
    for (std::size_t j=0; j<dim; ++j)
    G(i,j) = c0*4.0*alpha/sqrt_pi*_kronecker_delta[i][j];
  } // end if()
  else
  {
    // coefficients
    const Real  r = x.norm(), r2  = r*r, a2 = alpha*alpha;
    const Real c1 = std::erfc(alpha*r)/r;
    const Real c2 = 2.*alpha/sqrt_pi*std::exp(-a2*r2);
    
    // compute G tensor
    for (std::size_t i=0; i<dim; ++i)
    {
      for (std::size_t j=0; j<dim; ++j)
      {
        Real tmp1 = ( _kronecker_delta[i][j] + x(i)*x(j)/r2 )*c1;
        Real tmp2 = ( _kronecker_delta[i][j] - x(i)*x(j)/r2 )*c2;
        G(i,j) = c0*( tmp1 - tmp2 );
      } // end for j-loop
    } // end for i-loop
    
  } // end if-else
  
  STOP_LOG ("green_tensor_exp()", "GGEMStokes");
  
  // done
  return G;
}

// ======================================================================
DenseMatrix<Number> GGEMStokes::green_tensor_regularized(const Point& x,
                                                         const Real& alpha,
                                                         const Real& mu,
                                                         const Real& ksi,
                                                         const std::size_t& dim,
                                                         const bool& zero_limit) const
{
  START_LOG ("green_tensor_regularized()", "GGEMStokes");
  
  const Real c0 = 1.0/(8.*PI*mu);
  DenseMatrix<Number> G(dim,dim);
  
  // check if zero limit is true.
  if (zero_limit)
  {
    for (std::size_t i=0; i<dim; ++i)
      for (std::size_t j=0; j<dim; ++j)
        G(i,j) = c0*4.0/sqrt_pi*(ksi - alpha)*_kronecker_delta[i][j];
  } // end if()
  else
  {
    // coefficients
    const Real  r   = x.norm(), r2  = r*r;
    const Real  kr  = ksi*r,   kr2 = kr*kr;
    const Real  ar  = alpha*r, ar2 = ar*ar;
    const Real  c1  = ( std::erf(kr) - std::erf(ar) )/r;
    const Real  c2  = 2./sqrt_pi*( ksi*std::exp(-kr2) - alpha*std::exp(-ar2) );
    
    // compute G tensor
    for (std::size_t i=0; i<dim; ++i)
    {
      for (std::size_t j=0; j<dim; ++j)
      {
        Real tmp1 = ( _kronecker_delta[i][j] + x(i)*x(j)/r2 )*c1;
        Real tmp2 = ( _kronecker_delta[i][j] - x(i)*x(j)/r2 )*c2;
        G(i,j) = c0*( tmp1 + tmp2 );
      } // end for j-loop
    } // end for i-loop
    
  } // end if-else
  
  STOP_LOG ("green_tensor_regularized()", "GGEMStokes");
  
  // done
  return G;
}
  
  
  
// ======================================================================
DenseMatrix<Number> GGEMStokes::green_tensor_diff(const Point& x,
                                                  const Real& alpha,
                                                  const Real& mu,
                                                  const Real& ksi,
                                                  const std::size_t& dim,
                                                  const bool& zero_limit,
                                                  const DeltaFunType& delta_type_a,
                                                  const DeltaFunType& delta_type_k) const
{
  START_LOG ("green_tensor_diff()", "GGEMStokes");
  
  // Initialization
  DenseMatrix<Number> GTa(dim,dim),  GTk(dim,dim);

  // compute GT_a and GT_k
  GTa = this->green_tensor(x, alpha, mu, dim, zero_limit, delta_type_a);
  GTk = this->green_tensor(x, ksi,   mu, dim, zero_limit, delta_type_k);
  GTk -= GTa;
  
  STOP_LOG ("green_tensor_diff()", "GGEMStokes");
  
  // done
  return GTk;
}

  

// ======================================================================
std::vector<Real> GGEMStokes::local_velocity_fluid(PointMesh<3>*  point_mesh,
                                                   const Point& ptx,
                                                   const Real& alpha,
                                                   const Real& ibm_beta,
                                                   const Real& mu,
                                                   const Real& br0,
                                                   const Real& hmin,
                                                   const std::size_t& dim,
                                                   const std::string& force_type) const
{
  START_LOG ("local_velocity_fluid()", "GGEMStokes");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   build the particle neighbor list around the given point \p ptx
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const bool is_sorted  = false;
  std::vector<std::pair<std::size_t,Real> > IndicesDists;
  point_mesh->build_particle_neighbor_list(ptx, is_sorted, IndicesDists);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over all the neighbor list beads, and compute the local velocity:
   u_l(i) = sum[v=1:Nl] G_v(x-x_v; i,j)*f_v(j)
   zero_limit = false, if the given point is 'fluid' (not a bead/tracking pt.)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //const Real br0    = 1.0;    // normalized bead radius.
//  const bool  zero_limit  = false;  // for fluid, we let it be false, and allow singularity
  std::vector<Real> u(dim,0.0);
  for (std::size_t v=0; v<IndicesDists.size(); ++v)
  {
    // 0. particle id and position, vector x = ptx - pt0
    const std::size_t p_id = IndicesDists[v].first;
    const Point pt0 = point_mesh->particles()[p_id]->point();
    const Point x   = point_mesh->pm_periodic_boundary()->point_vector(pt0,ptx);
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    bool  zero_limit  = false;
    Point dpt = pt0 - ptx;
    if(dpt.norm()<r_eps) zero_limit  = true;
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    // Determine the regularization parameter ksi according to the point type.
    const PointType point_type = point_mesh->particles()[p_id]->point_type();
    const Real ksi     = this->regularization_parameter(hmin,ibm_beta,br0,point_type);
    
    // 1. compute the Green function (Oseen Tensor) of particle-v
    DenseMatrix<Number> GT; // Green function Tensor has the size: dimxdim
    if(force_type == "regularized")
      GT = this->green_tensor_regularized(x,alpha,mu,ksi,dim,zero_limit);
    else
      libmesh_assert ("GGEMStokes::local_velocity_fluid, wrong force_type!");

    // end if-else
    
    // 2. compute the force vector of particle-v
    const Point fv = point_mesh->particles()[p_id]->particle_force();
    
    // 3. compute u due to this particle
    for (std::size_t i=0; i<dim; ++i)
      for (std::size_t j=0; j<dim; ++j)
        u[i] += GT(i,j)*fv(j);
    
    // ================ test output of GT, fv, and u ===========
//    printf("--------------------------- output matrix GT ---------------------------\n");
//    for(unsigned int i=0; i<GT.m(); ++i)
//    {
//      for(unsigned int j=0; j<GT.n(); ++j)    printf("%f  ", GT(i,j) );
//      printf(", %f  \n", fv[i]);
//    }
//    printf("--->test: local velocity vector u = (%f %f %f)\n", u[0],u[1],u[2]);
//    printf("--------------------------- end of matrix GT ---------------------------\n");
    // ==============================================
    
  } // end for v-loop

  
  STOP_LOG ("local_velocity_fluid()", "GGEMStokes");
  return u;
}

// ======================================================================
std::vector<Real> GGEMStokes::local_velocity_fluid(PointMesh<3>*  point_mesh,
                                                   const Elem* elem,
                                                   const Point& ptx,
                                                   const Real& alpha,
                                                   const Real& ibm_beta,
                                                   const Real& mu,
                                                   const Real& br0,
                                                   const Real& hmin,
                                                   const std::size_t& dim,
                                                   const std::string& force_type) const
{
  START_LOG ("local_velocity_fluid()", "GGEMStokes");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   build the particle neighbor list around the given point \p ptx
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const bool is_sorted  = false;
  std::vector<std::size_t> elem_neighbors = point_mesh->elem_neighbor_list(elem);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over all the neighbor list beads, and compute the local velocity:
   u_l(i) = sum[v=1:Nl] G_v(x-x_v; i,j)*f_v(j)
   zero_limit = false, if the given point is 'fluid' (not a bead/tracking pt.)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //const Real br0    = 1.0;    // normalized bead radius.
//  const bool  zero_limit  = false;  // for fluid, we let it be false, and allow singularity
  std::vector<Real> u(dim,0.0);
  for (std::size_t v=0; v<elem_neighbors.size(); ++v)
  {
    // 0. particle id and position, vector x = ptx - pt0
    const std::size_t p_id = elem_neighbors[v];
    const Point pt0 = point_mesh->particles()[p_id]->point();
    const Point x   = point_mesh->pm_periodic_boundary()->point_vector(pt0,ptx);
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    bool  zero_limit  = false;
    Point dpt = pt0 - ptx;
    if(dpt.norm()<r_eps) zero_limit  = true;
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    // Determine the regularization parameter ksi according to the point type.
    const PointType point_type = point_mesh->particles()[p_id]->point_type();
    const Real ksi     = this->regularization_parameter(hmin,ibm_beta,br0,point_type);
    
    // 1. compute the Green function (Oseen Tensor) of particle-v
    DenseMatrix<Number> GT; // Green function Tensor has the size: dimxdim
    if (force_type=="regularized")
      GT = this->green_tensor_regularized(x,alpha,mu,ksi,dim,zero_limit);
    else
      libmesh_assert ("GGEMStokes::local_velocity_fluid, wrong force_type!");
    // end if-else
    
    // 2. compute the force vector of particle-v
    const Point fv = point_mesh->particles()[p_id]->particle_force();
    
    // 3. compute u due to this particle
    for (std::size_t i=0; i<dim; ++i)
      for (std::size_t j=0; j<dim; ++j)
        u[i] += GT(i,j)*fv(j);
    
    // ================ test output of GT, fv, and u ===========
//    printf("--------------------------- output matrix GT ---------------------------\n");
//    for(unsigned int i=0; i<GT.m(); ++i)
//    {
//      for(unsigned int j=0; j<GT.n(); ++j)    printf("%f  ", GT(i,j) );
//      printf(", %f  \n", fv[i]);
//    }
//    printf("--->test: local velocity vector u = (%f %f %f)\n", u[0],u[1],u[2]);
//    printf("--------------------------- end of matrix GT ---------------------------\n");
    // ==============================================
    
  } // end for v-loop

  
  STOP_LOG ("local_velocity_fluid()", "GGEMStokes");
  return u;
}

  
// ======================================================================
Point GGEMStokes::local_velocity_bead(PointMesh<3>*  point_mesh,
                                      const std::size_t& pid0,
                                      const Real&   alpha,
                                      const Real&   ibm_beta,
                                      const Real&   mu,
                                      const Real&   br0,
                                      const Real&   hmin,
                                      const std::size_t& dim,
                                      const std::string& force_type) const
{
  START_LOG ("local_velocity_bead()", "GGEMStokes");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Find the bead position ptx and its neighbor list, and its point type
   NOTE: for a given bead/tracking pt, its neighbor list does NOT include itself!
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const Point& ptx      = point_mesh->particles()[pid0]->point();
  const PointType point_type0 = point_mesh->particles()[pid0]->point_type();
  std::vector<std::pair<std::size_t,Real> > IndicesDists;
  IndicesDists          = point_mesh->particles()[pid0]->neighbor_list();
  const std::vector<Point>& neighbor_vector= point_mesh->particles()[pid0]->neighbor_vector();    
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over all the neighbor list beads, and compute the local velocity:
   u_l(i) = sum[v=1:Nl] G_v(x-x_v; i,j)*f_v(j)
   zero_limit = false, because the neighbor list does NOT include itself.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  bool  zero_limit    = false;    // this will be changed to "true" for tracking points
  Point u;
  DenseMatrix<Number> GT;         // Green function Tesnosr has the size: dim x dim
  const Real ksi = this->regularization_parameter(hmin, ibm_beta, br0, point_type0);
  for (std::size_t v=0; v<IndicesDists.size(); ++v)
  {
    // 0. particle id and position, vector x = ptx - pt0
    const Point x = neighbor_vector[v];
    const unsigned int p_id = IndicesDists[v].first;
    // Determine the regularization parameter ksi according to the point type.
    // const Real ksi     = this->regularization_parameter(hmin,br0,point_type);    
    // 1. compute the Green function (Oseen Tensor) of particle-v
    if(force_type=="regularized")
      GT = this->green_tensor_regularized(x,alpha,mu,ksi,dim,zero_limit);
    else
      libmesh_assert("GGEMStokes::local_velocity_bead, wrong force_type");
    // end if-else
    
    // 2. compute the force vector of particle-v
    const Point fv = point_mesh->particles()[p_id]->particle_force();
    
    // 3. compute u due to this particle
    for (std::size_t i=0; i<dim; ++i)
      for (std::size_t j=0; j<dim; ++j)
        u(i) += GT(i,j)*fv(j);
  } // end for v-loop
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   The above velocity does NOT include the contribution from the bead itself, since
   it is not in its own neighbor list. This is true for a polymer bead. However,
   if this point is a tracking point of an immersed body in Immersed Boundary Method,
   this contribution should be included by letting x-->0
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if( point_type0==LAGRANGIAN_POINT )  // for Lagrangian tracking point type only
  {
    // 1. compute the Green function (Oseen Tensor) when x-->0
    zero_limit      = true;
    const Real ksi0 = this->regularization_parameter(hmin,ibm_beta,br0,point_type0);
    GT = this->green_tensor_regularized(ptx,alpha,mu,ksi0,dim,zero_limit);
    // end if-else
    
    // 2. compute the force vector of this particle
    const Point fv = point_mesh->particles()[pid0]->particle_force();
    
    // 3. compute u due to this particle
    for (std::size_t i=0; i<dim; ++i)
      for (std::size_t j=0; j<dim; ++j)
        u(i) += GT(i,j)*fv(j);
  }
  
  
  STOP_LOG ("local_velocity_bead()", "GGEMStokes");
  return u;
}

  

// ======================================================================
Point GGEMStokes::global_self_exclusion(PointMesh<3>* point_mesh,
                                                    const std::size_t&  pid0,
                                                    const Real& alpha,
                                                    const Real& mu,
                                                    const std::size_t& dim) const
{
  START_LOG ("self_exclusion()", "GGEMStokes");
  
  // 1. compute the force vector of the particle
  const Point fv = point_mesh->particles()[pid0]->particle_force();
  const Real   c0 = 1.0/(8.*PI*mu);
  
  // 2. compute the self exclusion term at the position of the i-th bead
  Point self_v;
  for (std::size_t i=0; i<dim; ++i)
    self_v(i) = c0*4.0*alpha/sqrt_pi*fv(i);

  STOP_LOG ("self_exclusion()", "GGEMStokes");
  
  return self_v;  // done and return
}


// ======================================================================
DenseMatrix<Number> GGEMStokes::rpy_tensor(const Point& x,  /* vector x = pt1 - pt0 */
                                           const Real& mu,  /* viscosity */
                                           const Real& a,   /* bead radius */
                                           const std::size_t& dim) const /*dim==3*/
{
  START_LOG ("rpy_tensor()", "GGEMStokes");
  
  const Real  r = x.norm(), r2  = r*r, a2 = a*a;
  const Real t1 = 1.0/(8.*PI*mu*r);
  const Real t2 = 1.0/(6.*PI*mu*r);
  Real c0 = 0.0, C1 = 0., C2 = 0.;
  
  // init G tensor
  DenseMatrix<Number> G(dim,dim);
  
  // avoid singularity
  if(r >= r_eps)
  {
    for (std::size_t i=0; i<dim; ++i)
    {
      for (std::size_t j=0; j<dim; ++j)
      {
        if( r>=2.*a )
        {
          c0 = t1;
          C1 = 1.0 + (2.0*a2)/(3.0*r2);
          C2 = 1.0 - (2.0*a2)/r2;
        }
        else
        {
          c0 = t2;
          C1 = 1.0 - (9.0*r)/(32.0*a);
          C2 = (3.0*r)/(32.0*a);
        }
        const Real dij = _kronecker_delta[i][j];
        G(i,j) = c0*( C1*dij + C2*x(i)*x(j)/r2 );
      } // end j-loop
    } // end i-loop
  } // end if
  
  STOP_LOG ("rpy_tensor()", "GGEMStokes");
  return G;
}


  
// ======================================================================
DenseMatrix<Number> GGEMStokes::mobility_tensor(const Point& x,     /* vector x = pt1 - pt0 */
                                                const Real& mu,     /* kinematic viscosity */
                                                const Real& a,      /* bead radius */
                                                const std::size_t& dim) const /*dim==3*/
{
  START_LOG ("mobility_tensor()", "GGEMStokes");
  
  const Real  r = x.norm();
  const Real t1 = 1.0/(6.*PI*mu*a);
  
  // init M & G tensor
  DenseMatrix<Number> M(dim,dim), G(dim,dim);
  G = this->rpy_tensor(x, mu, a, dim);
  
  if( r<r_eps ) // rk = rl <=> r_kl = 0   => rpy_tensor = 0
  {
    for (std::size_t i=0; i<dim; ++i) {
      for (std::size_t j=0; j<dim; ++j) {
        const Real dij = _kronecker_delta[i][j];
        M(i,j) = t1*dij;
      }
    }
  }
  else
  {
    for (std::size_t i=0; i<dim; ++i) {
      for (std::size_t j=0; j<dim; ++j) {
        M(i,j) = G(i,j);
      }
    }
  }
  
  
  STOP_LOG ("mobility_tensor()", "GGEMStokes");
  return M;
}
  
  

// ======================================================================


  
  
  
  
  
  //
  //// ======================================================================
  //std::vector<Real> GGEMStokes::local_velocity(ParticleMesh<3>* point_mesh,
  //                                             const Elem* elem,
  //                                             const Point& ptx,
  //                                             const Real& alpha,
  //                                             const Real& mu,
  //                                             const Real& ksi,
  //                                             const std::size_t& dim,
  //                                             const std::string& option)
  //{
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
  //    const std::vector<Real> fv = point_mesh->particles()[p_id]->particle_force();
  //
  //    // 3. compute u due to this particle
  //    for (std::size_t i=0; i<dim; ++i)
  //      for (std::size_t j=0; j<dim; ++j)
  //        u[i] += GT(i,j)*fv[j];
  //  } // end for v-loop
  //
  //  // ----------------------- test output -------------------------
  //  Real u_mag = 0.0;
  //  for (std::size_t i=0; i<dim; ++i) u_mag += u[i]*u[i];
  //  u_mag = std::sqrt( u_mag );
  //  if ( n_list.size()>0 )// && u_mag>1e-10
  //  {
  //    printf("------------ GGEMStokes::local_velocity(elem) -------------\n");
  //    printf("There are %lu neighbor particles around the element within r = %f\n",
  //           n_list.size(), point_mesh->search_radius("e") );
  //
  //    for (std::size_t i=0; i<n_list.size(); ++i)
  //      printf("%lu ",n_list[i]);
  //    printf("\n");
  //
  //    printf("vel_bc = ");
  //    for (std::size_t i=0; i<dim; ++i)
  //      printf("%f, ", u[i] );
  //    printf(" at node (%f, %f, %f)\n\n", ptx(0),ptx(1),ptx(2) );
  //  }
  //  // -------------------------------------------------------------
  //  
  //  STOP_LOG ("local_velocity(Elem)", "GGEMStokes");
  //  return u;
  //}
  

  
} // end of namespace
