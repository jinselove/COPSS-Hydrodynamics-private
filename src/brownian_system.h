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
#include <chrono>

#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_config.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/slepc_macro.h"
#include "libmesh/parallel_object.h"
#include "libmesh/reference_counted_object.h"



#include "random_generator.h"
#include "pm_linear_implicit_system.h"
#include "fix/fix.h"
// include SLEPc EPS solver, "libmesh/slepc_macro.h" must be included before this!
EXTERN_C_FOR_SLEPC_BEGIN
# include <slepceps.h>
EXTERN_C_FOR_SLEPC_END


using namespace libMesh;


/**
 * Internal function if shell matrix mode is used, this just
 * calls the shell matrix's matrix multiplication function.
 * In this class, it will call Stokes solver. u = M*f
 * NOTE that the pm_system should be reinit() before solving Stokes
 */
extern "C"
{
  PetscErrorCode _MatMult_Stokes(Mat M,Vec f,Vec u);
}



namespace libMesh
{

  
class EquationSystems;
  
  
  
  /*
   * The BrownianSystem is designed to deal with the stochastic
   * process of Brownian Dynamics(BD)
   */
  
  
class BrownianSystem :  public ReferenceCountedObject<BrownianSystem>,
                        public ParallelObject
{
public:
  // Constructor
  BrownianSystem(EquationSystems& es);


  // Destructor
  ~BrownianSystem();



  /*
   * Return a constant reference of equation system
   */
  const EquationSystems& get_equation_systems() const
  { return _equation_systems; };

  
  /*
   * Return a reference of equation system
   */
  EquationSystems& get_equation_systems()
  { return _equation_systems; };

 
  /*
   * Manage PETSc Vec's scatters and gathers to all process
   * NOTE: We can't put VecScatter inside the function,
   * VecScatterCreateToAll() should be outside!
   */
  PetscErrorCode petsc_vec_scatter_to_all(Vec f,    // vin
                                          Vec vf,   // vout
                                          VecScatter scatter,
                                          const std::string& mode);

  
  /*
   * Extract particle coordinates as a PETSc vector     : mode = "extract"
   * Assign the paticle coordinates from a PETSc vector : mode = "assign"
   * vec_type can be either "coordinate" or "force"
   */
  PetscErrorCode extract_particle_vector(Vec* x,
                                         const std::string& vec_type,
                                         const std::string& mode);
  
  
  /*
   * Transformation from a std::vector<>  to a PETSc Vec
   * mode == "forward":   std_vec   -> petsc_vec
   * mode == "backward":  petsc_vec -> std_vec
   * NOTE: the PETSc vector must be initialized with some parallel distributions
   */
  PetscErrorCode vector_transform(std::vector<Real>& std_vec,
                                  Vec* petsc_vec,
                                  const std::string& mode) const;
  
  
  // Random generator from C++ std library
  RandomGenerator& random_generator() { return _random_generator; }
  
  
  /*
   * set the random seed for standard random vectors
   */
  void set_std_random_seed(const std::size_t seed_val)
  {  _random_generator.set_seed(seed_val);  }
  
  
  /*
   * Generate a N-components random vector with gaussian/uniform distribution
   */
  PetscErrorCode std_random_vector(const Real& a,      // a (mean)
                                   const Real& b,      // b (std deviation)
                                   const std::string& random_type,
                                   Vec* u);
  
  
  /*
   * Generate a N-components random vector with uniform or gaussian(normal) distribution
   *
   * unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
   * may cause problems, because the seed (time) on each CPU may be different!
   */
  std::vector<Real> std_random_vector(const std::size_t N,// length
                                      const Real& a,      // a (mean)
                                      const Real& b,      // b (std deviation)
                                      const std::string& random_type);
  
  
  /*
   * Initialize PETSc random
   * The init and generate parts have to be separated, because in the real 
   * simulation, we only need to init once and generate rand number at each step.
   */
  PetscErrorCode init_petsc_random(PetscRandom* rand_ctx);
  
  
  /*
   * Generate PETSc random vector u with N components
   * This is only for the uniform distribution!
   */
  PetscErrorCode petsc_random_vector(const PetscInt N,
                                     PetscRandom* rand_ctx,
                                     Vec* u);
  
  
  /*
   * Return the smallest or largest eigenvalue of diffusion matrix D
   * option == "smallest" or "largest"
   */
  Real compute_eigenvalue(const std::string& option,
                          const Real tol);
  
  
  /*
   * Return the smallest or largest eigenvalue of diffusion matrix D
   * option == "smallest" or "largest"
   * The difference is that this allows EPS solver to use initial space
   */
  Real compute_eigenvalue(const std::string& option,
                          const Real tol,
                          const bool use_init_space,
                          Vec* v0);
  
  
  /*
   * Return both the smallest and largest eigenvalues
   */
  PetscErrorCode compute_eigenvalues(Real& eigen_min,
                                     Real& eigen_max,
                                     const Real& tol);
  
  /*
   * Compute the largest eigenvalue using the power_iteration
   * This is only for the purpose of test, and it is also supportted
   * by the SLEPc library.
   */
  Real power_iteration();
  
  
  
  /*
   * The transforming matrix from physical space to Chebyshev tranform space
   * using Gauss-Lobatto method
   */
  DenseMatrix<Number> chebyshev_transform_matrix(const std::size_t N) const;
  
  
  /*
   * Quadrature points in [-1,1] for the chebyshev transformation
   * (Gauss-Lobatto method)
   */
  DenseVector<Number> chebyshev_quadrature_point(const std::size_t N) const;
  
  
  /*
   * Chebyshev polynomial approximation:
   * y =  B*dw  or  y =  B^-1 * dw
   * where, B = sqrt(D),  D = kB*T*M, and dw is a random force vector.
   */
  bool chebyshev_polynomial_approximation(const std::size_t N,
                                          const Real eigen_min,
                                          const Real eigen_max,
                                          const Real tol_cheb,
                                          Vec* dw);
  
  /*
   * A wrap-up Stokes solver that simply call _MatMult_Stokes() to compute
   * the particle velocity vector provided that the force vector is given.
   */
  PetscErrorCode hi_ewald(Mat M, Vec f, Vec u)
  {   return _MatMult_Stokes(M, f, u);  }
  
  
  
  // this generate a diagonal matrix, which is used only for test purpose
  // u = M*f, where M(i,i) = i.
  PetscErrorCode _test_diagonal_mat(Vec f, Vec u);
  
  // This is another test function that calculates a vector's mean and variance
  PetscErrorCode _vector_mean_variance(Vec u,
                                       Real& mean,
                                       Real& variance) const;
  PetscErrorCode _vector_mean_variance(const std::vector<Real>& u) const;
  
  
  
  // create a shell matrix with dimension N x N
  PetscErrorCode _create_shell_mat(const std::size_t N, Mat* M);
  
  // create a vector with dimension N which has the same partition as M
  PetscErrorCode _create_petsc_vec(const std::size_t N, Vec* V);
  
  

private:
  // Reference of Equation systems
  EquationSystems & 	_equation_systems;
  
  // total number of points
  unsigned int _n_points;

  // point mesh
  PointMesh<3>* _point_mesh;

  // fixes
  std::vector<Fix*> _fixes;

  // dimension of mesh
  unsigned int _dim;  
  
  // Random generator from C++ std library
  RandomGenerator _random_generator;
  

  // friend class/function, so that it can reach the members in the current class.
  friend PetscErrorCode _MatMult_Stokes(Mat M, Vec f, Vec u);

  
};  // end of namespace
  
  
  
  
} // end of namespace
