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
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <random>

#include "libmesh/libmesh_logging.h"
#include "libmesh/equation_systems.h"


#include "rigid_particle.h"
#include "particle_mesh.h"
#include "point_particle.h"
#include "point_mesh.h"
#include "brownian_system.h"
#include "chebyshev.h"



// ========================================================================================
// NOTE, the input matrix and vectors are required to be initialized parallelly before use.
// ========================================================================================
PetscErrorCode _MatMult_Stokes(Mat M,Vec f,Vec u)
{
  void              *shell_ctx;
  PetscErrorCode    ierr;
  PetscFunctionBeginUser;
  START_LOG ("_MatMult_Stokes()", "BrownianSystem");
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Extract the system from ctx and get the particle mesh system
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = MatShellGetContext(M,&shell_ctx);    CHKERRQ(ierr);
  BrownianSystem& brownian_sys      = *(BrownianSystem*)shell_ctx;
  PMLinearImplicitSystem & pm_system =
  brownian_sys.get_equation_systems().get_system<PMLinearImplicitSystem> ("Stokes");
  std::vector<PointParticle*> particles  = pm_system.point_mesh()->particles();
  unsigned int  _n_points     = pm_system.point_mesh()->num_particles();
  unsigned int _dim = pm_system.get_mesh().mesh_dimension();
  //std::vector<RigidParticle*> particles  = pm_system.particle_mesh()->particles();
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Distribute the force vector f to each particle.
   Here, f = [ f1x, f1y, f1z; f2x,f2y,f2z; ...; fNx, fNy, fNz]
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  brownian_sys.extract_particle_vector(&f,"force","assign");
  MPI_Barrier( PETSC_COMM_WORLD );
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   TEST: print out the particle force on each particle
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //if(pm_system.comm().rank()==0) printf("--->test: MatMult_Stokes() - 2\n");
  //ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  //for (std::size_t i=0; i<_n_points; ++i)
  //{
  //  std::vector<Real> pforce = particles[i]->particle_force();
  //  printf("the %lu-th particle force = (%f, %f, %f) on the rank = %d\n",
  //         i, pforce[0], pforce[1], pforce[2],rank );
  //}
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Solve the Stokes equation to obtain the particle velocity vector pv
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const bool reinit_stokes = !( pm_system.stokes_solver().is_ksp_initialized() );
  pm_system.solve_stokes("disturbed", reinit_stokes);
  std::vector<Real> pvelocity(_dim*_n_points,0);
  pm_system.compute_point_velocity("disturbed", pvelocity);
  MPI_Barrier( PETSC_COMM_WORLD );
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Transform the velocity vector u from the std::vector form to the PETSc Vec form
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  brownian_sys.vector_transform(pvelocity,&u,"forward");
  MPI_Barrier( PETSC_COMM_WORLD );
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   TEST: view Vec u
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//  int vsize; ierr = VecGetLocalSize(u,&vsize); CHKERRQ(ierr);
//  printf("--->test: size(u) = %d on the rank %d\n",vsize, rank);CHKERRQ(ierr);
//  MPI_Barrier( PETSC_COMM_WORLD );
//  if(pm_system.comm().rank()==0) printf("--->test: MatMult_Stokes() - 4\n");
//  ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); // View the vector */
//  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n\n\n");CHKERRQ(ierr);

  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Return
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  STOP_LOG ("_MatMult_Stokes()", "BrownianSystem");
  PetscFunctionReturn(0);
}
// ======================================================================================


//namespace libMesh
//{

// ======================================================================================
BrownianSystem::BrownianSystem(EquationSystems& es)
: ParallelObject (es),
  _equation_systems(es)
{
  PMLinearImplicitSystem & _pm_system = es.get_system<PMLinearImplicitSystem>("Stokes");
  _point_mesh = _pm_system.point_mesh();
  _n_points     = _point_mesh->num_particles();
  _fixes = _pm_system.fixes();
  _dim = _pm_system.get_mesh().mesh_dimension();
}


// ======================================================================================
BrownianSystem::~BrownianSystem()
{
  // do nothing
}



// ======================================================================================
PetscErrorCode BrownianSystem::petsc_vec_scatter_to_all(Vec f,    // vin
                                                        Vec vf,   // vout
                                                        VecScatter scatter,
                                                        const std::string& mode)
{
#ifdef LIBMESH_HAVE_PETSC
  PetscErrorCode    ierr;
  PetscFunctionBeginUser;
  START_LOG ("petsc_vec_scatter_to_all()", "BrownianSystem");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     A generalized scatter from one vector to another.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if(mode == "forward") { // scatter: f -> vf
    ierr = VecScatterBegin(scatter,f,vf,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd  (scatter,f,vf,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  }
  else if(mode == "reverse") {  // gather: vf -> f
    ierr = VecScatterBegin(scatter,vf,f,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
    ierr = VecScatterEnd  (scatter,vf,f,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
  }
  else {
    libmesh_assert("*** error in BrownianSystem::petsc_vec_scatter(): WRONG ScatterMode!");
    libmesh_error();
  }
  
  STOP_LOG ("petsc_vec_scatter_to_all()", "BrownianSystem");
  PetscFunctionReturn(0);
#endif
}



// ======================================================================================
PetscErrorCode BrownianSystem::extract_particle_vector(Vec* x,
                                                       const std::string& vec_type,
                                                       const std::string& mode)
{
#ifdef LIBMESH_HAVE_PETSC
  Vec             vx;
  VecScatter      scatter;
  PetscErrorCode  ierr;
  PetscScalar     *px;
  PetscFunctionBeginUser;
  START_LOG ("extract_particle_vector()", "BrownianSystem");
  

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Get the particle mesh system and its parameters
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  std::vector<PointParticle*> point_particles  = _point_mesh->particles();
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   form the vec with the desired distribution if(mode=="extract").
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if(mode=="extract")
  {
    const std::size_t n_vec = _dim*_n_points;
    ierr = this->_create_petsc_vec(n_vec, x);CHKERRQ(ierr);
  }
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Creates a scatter context that copies all values of x to vx on each processor
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = VecScatterCreateToAll(*x,&scatter,&vx);   CHKERRQ(ierr);
  ierr = this->petsc_vec_scatter_to_all(*x,vx,scatter,"forward"); // x->vx
  ierr = VecGetArray(vx,&px);                      CHKERRQ(ierr);
  
  // ---extract/assign the coordinates to vx = [x1,y1,z1; x2,y2,z2, ..., xN,yN,zN]
  if(vec_type=="coordinate")
  {
    for(std::size_t i=0; i<_n_points; ++i)
    {
      Point& pt = point_particles[i]->point();
      if(mode=="extract")
        for(std::size_t j=0; j<_dim; ++j) px[i*_dim + j] = pt(j);
      else if(mode=="assign"){
        for(std::size_t j=0; j<_dim; ++j) pt(j) = px[i*_dim + j];
        // apply periodic boundary condition
        _fixes[0]->check_walls();
      }
      else
        libmesh_error();
      // end if-else
    } // end for i-loop
  }
  // ---extract/assign the force to vf = [f1x,f1y,f1z; f2x,f2y,f2z; ...; fNx,fNy,fNz]
  else if(vec_type=="force")
  {
    for (std::size_t i=0; i<_n_points; ++i)
    {
      std::vector<Real> pforce(_dim);
      if(mode=="extract")
      {
        pforce = point_particles[i]->particle_force();
        for(std::size_t j=0; j<_dim; ++j) px[i*_dim+j] = pforce[j];
      }
      else if(mode=="assign")
      {
        for(std::size_t j=0; j<_dim; ++j) pforce[j] = px[i*_dim+j];
        point_particles[i]->set_particle_force(pforce);
      }
      else
        libmesh_error();
      // end if-else
    } // end for i-loop
  }
  else
  {
    libmesh_error();
  } // end if- else
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Scatter reversely from vx to x only when we want to extract values to x
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if(mode=="extract")
  { ierr = this->petsc_vec_scatter_to_all(*x,vx,scatter,"reverse");CHKERRQ(ierr); }
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Destroy objects and return
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = VecRestoreArray(vx,&px);     CHKERRQ(ierr);
  ierr = VecDestroy(&vx);             CHKERRQ(ierr);
  ierr = VecScatterDestroy(&scatter); CHKERRQ(ierr);
  
  STOP_LOG ("extract_particle_vector()", "BrownianSystem");
  PetscFunctionReturn(0);
#endif
}



// ======================================================================================
PetscErrorCode BrownianSystem::vector_transform(std::vector<Real>& std_vec,
                                                Vec* petsc_vec,
                                                const std::string& mode) const
{
#ifdef LIBMESH_HAVE_PETSC
  PetscInt          low, high, nlocal, vsize;
  PetscScalar       *pv;
  PetscErrorCode    ierr;
  PetscFunctionBeginUser;
  START_LOG ("vector_transform()", "BrownianSystem");

  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Get the Ownership range and local components, then copy!
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = VecGetOwnershipRange(*petsc_vec,&low,&high); CHKERRQ(ierr);
  ierr = VecGetLocalSize(*petsc_vec,&nlocal);         CHKERRQ(ierr);
  ierr = VecGetArray(*petsc_vec,&pv);                 CHKERRQ(ierr);
  
  
  if(mode=="forward")
  {
    for(int i=0; i<nlocal; ++i)
      pv[i] = std_vec[low+i];
  }
  else if(mode=="backward")
  {
    std_vec.clear();
    std_vec.resize(std::size_t(nlocal));
    for(int i=0; i<nlocal; ++i)
      std_vec[i] = pv[i];
    
    this->comm().allgather(std_vec);
  }
  else
  {
    libMesh::err << "*** illegal vector transform mode (forward or backward only)!\n";
    libmesh_here();
    libmesh_error();
  }

  // restore and return
  ierr = VecRestoreArray(*petsc_vec,&pv);             CHKERRQ(ierr);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Make sure the size of these two vectors are the same!
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = VecGetSize(*petsc_vec, &vsize);CHKERRQ(ierr);
  const std::size_t pv_size = (std::size_t)vsize;
  if ( pv_size != std_vec.size() )
  {
    libMesh::err << "*** vector_transform(): The size of vectors are not equal! \n";
    libmesh_error();
  }
  
  STOP_LOG ("vector_transform()", "BrownianSystem");
  PetscFunctionReturn(0);
#endif
}



// ======================================================================================
PetscErrorCode BrownianSystem::std_random_vector(const Real& a,      // a (mean)
                                                 const Real& b,      // b (std deviation)
                                                 const std::string& random_type,
                                                 Vec* u)
{
#ifdef LIBMESH_HAVE_PETSC
  PetscFunctionBeginUser;
  START_LOG ("std_random_vector()", "BrownianSystem");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Get the vector size and problem dimension.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const std::size_t       n_vec = _dim*_n_points;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Generate random vector and assign it to vu on each process
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  std::vector<Real> randv = this->std_random_vector(n_vec,a,b,random_type);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Transform to PETSc Vec (randv-->u) and destroy objects
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  this->_create_petsc_vec(n_vec,u);
  this->vector_transform(randv,u,"forward");
  
  STOP_LOG ("std_random_vector()", "BrownianSystem");
  PetscFunctionReturn(0);
#endif
}



// ======================================================================================
std::vector<Real> BrownianSystem::std_random_vector(const std::size_t N,// length
                                                    const Real& a,      // a (mean)
                                                    const Real& b,      // b (std deviation)
                                                    const std::string& random_type)
{
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Generate random numbers with Gaussian(normal) distribution.
   mean = 0, standard deviation = 1. (variance = deviation^2 = 1)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  std::vector<Real> rand_v;
  if (random_type=="normal" || random_type=="gaussian")
  {
    rand_v = _random_generator.random_vector_normal(N,a,b);
  }
  else if (random_type=="uniform")
  {
    rand_v = _random_generator.random_vector_uniform(N,a,b);
  }
  else
  {
    libMesh::err << "*** BrownianSystem:: illegal random type (uniform gaussian or normal)!\n";
    libmesh_error();
  }
  this->comm().barrier();
  
  /* - - - - - test code: check the output mean/variance/deviation - - - - - - - - */
  this->_vector_mean_variance(rand_v);
  
  //return
  return rand_v;
}


// ======================================================================================
PetscErrorCode BrownianSystem::init_petsc_random(PetscRandom* rand_ctx)
{
  PetscFunctionBeginUser;
  START_LOG ("init_petsc_random()", "BrownianSystem");
  
  // Create Random
  PetscRandomCreate(PETSC_COMM_WORLD, rand_ctx);
#if defined(PETSC_HAVE_DRAND48)
  PetscRandomSetType(*rand_ctx,PETSCRAND48);
#elif defined(PETSC_HAVE_RAND)
  PetscRandomSetType(*rand_ctx,PETSCRAND);
#endif
  
  // FIXME: should we put it here or outside the function?
  PetscRandomSetFromOptions(*rand_ctx);
  
  STOP_LOG ("init_petsc_random()", "BrownianSystem");
  PetscFunctionReturn(0);
}


// ======================================================================================
PetscErrorCode BrownianSystem::petsc_random_vector(const PetscInt N,
                                                   PetscRandom* rand_ctx,
                                                   Vec* u)
{
#ifdef LIBMESH_HAVE_PETSC
  PetscErrorCode  ierr;
  PetscFunctionBeginUser;
  START_LOG ("petsc_random_vector()", "BrownianSystem");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create a uniform-distributed random vector with desired parallel distribution
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = this->_create_petsc_vec(N, u);  CHKERRQ(ierr);
  ierr = VecSetRandom(*u,*rand_ctx);     CHKERRQ(ierr);
  
  /* - - -  Test: check the mean and variance of the random vector - - - - - - - - */
//  Real mean = 0.0, variance = 0.0;
//  this->_vector_mean_variance(*u, mean, variance);
  
  // return
  STOP_LOG ("petsc_random_vector()", "BrownianSystem");
  PetscFunctionReturn(0);
#endif
}



// ======================================================================================
Real BrownianSystem::compute_eigenvalue(const std::string& option,
                                        const Real tol)
{
#ifdef LIBMESH_HAVE_SLEPC
  Mat            M;               /* operator matrix, which is diffusion matrix M */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  EPSConvergedReason reason;
  PetscInt       N = _n_points*3, its, maxits = 20;
  PetscScalar    value;
  PetscErrorCode ierr;
  START_LOG ("compute_eigenvalue()", "BrownianSystem");
  
  //PetscPrintf(PETSC_COMM_WORLD,"--->test: BrownianSystem::compute_eigenvalue()\n");
  //if(option=="smallest")
  //  PetscPrintf(PETSC_COMM_WORLD,"--->    (1) compute minimum eigenvalues.\n");
  //else if(option=="largest")
  //  PetscPrintf(PETSC_COMM_WORLD,"--->    (2) compute maximum eigenvalues.\n");
  //else
  //  PetscPrintf(PETSC_COMM_WORLD,"*** warning: input option is neither largest nor smallest!\n");
  // end if-else
    
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Compute the operator matrix that defines the eigensystem, Ax=kx
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = this->_create_shell_mat(N,&M);       CHKERRQ(ierr);
  ierr = MatSetFromOptions(M);                CHKERRQ(ierr);
  ierr = MatShellSetOperation(M,MATOP_MULT,(void(*)())_MatMult_Stokes);CHKERRQ(ierr);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create the eigensolver and set various options
   Start from Non-Hermitian type: EPS_NHEP or EPS_HEP(symmetric)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);      CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,M,NULL);           CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);       CHKERRQ(ierr);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set the tolerance and maximum iteration for EPS solver and EPS type:
   EPSKRYLOVSCHUR(Default)/EPSARNOLDI/EPSPOWER/EPSLAPACK/EPSGD/EPSJD
   EPSBLOPEX/EPSRQCG/EPSLANCZOS/EPSSUBSPACE/
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSSetType(eps,EPSKRYLOVSCHUR);CHKERRQ(ierr);
  //if(option=="largest"){
  //  ierr = EPSSetType(eps,EPSPOWER);CHKERRQ(ierr);
  //}
  ierr = EPSSetTolerances(eps,  tol,  maxits);  CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"EPS tol = %f, maxits = %d\n",tol,maxits);CHKERRQ(ierr);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Select portion of spectrum.
   Note: EPS_LARGEST_MAGNITUDE should be used with EPSPOWER
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if(option=="smallest")  // LOBPCG for smallest eigenvalue problem!
  { ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);CHKERRQ(ierr);  }
  else if(option=="largest")
  { ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE); CHKERRQ(ierr);  }
  else
  { ierr = EPSSetFromOptions(eps);  CHKERRQ(ierr);  }
  // end if-else
  ierr = EPSSetFromOptions(eps);  CHKERRQ(ierr);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Solve the eigensystem and get the solution
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"EPS solve starts ...\n");
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Print the results on the screen
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"\n--->test: EPS solve results:\n");
  ierr = EPSGetEigenvalue(eps,0,&value,NULL); CHKERRQ(ierr);
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"          Solution method: %s\n",type);
  ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"          Number of iterations of the method: %D\n",its);
  ierr = EPSGetConvergedReason(eps,&reason);
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"          EPS converged reason : %D\n",reason);
  //if(option=="smallest") {
  //  ierr = PetscPrintf(PETSC_COMM_WORLD,"          The computed minimum eigenvalue: %f\n",value);
  //}
  //else if(option=="largest") {
  //  ierr = PetscPrintf(PETSC_COMM_WORLD,"          The computed maximum eigenvalue: %f\n",value);
  //}
  //else {
  //  ierr = PetscPrintf(PETSC_COMM_WORLD,"          The computed eigenvalue: %f\n",value);
  //}
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"          the computed eigenvector:\n");
//  ierr = VecView(*v0,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); // View the vector */
 // ierr = PetscPrintf(PETSC_COMM_WORLD,"          \n\n");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Distroy solution and clean up
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&M);  CHKERRQ(ierr);
  
  
  STOP_LOG ("compute_eigenvalue()", "BrownianSystem");
  return value;
#endif
}



// ======================================================================================
Real BrownianSystem::compute_eigenvalue(const std::string& option,
                                        const Real tol,
                                        const bool use_init_space,
                                        Vec* v0)
{
#ifdef LIBMESH_HAVE_SLEPC
  Mat            M;               /* operator matrix, which is diffusion matrix M */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  EPSConvergedReason reason;
  PetscInt       N = _n_points*3, its, maxits = 20;
  PetscScalar    value;
  PetscErrorCode ierr;
  START_LOG ("compute_eigenvalue()", "BrownianSystem");
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Output on the screen
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscPrintf(PETSC_COMM_WORLD,"--->test: BrownianSystem::compute_eigenvalue()\n");
  if(option=="smallest")
    PetscPrintf(PETSC_COMM_WORLD,"--->    (1) compute minimum eigenvalues.\n");
  else if(option=="largest")
    PetscPrintf(PETSC_COMM_WORLD,"--->    (2) compute maximum eigenvalues.\n");
  else
    PetscPrintf(PETSC_COMM_WORLD,"*** warning: input option is neither largest nor smallest!\n");
  // end if-else
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Compute the operator matrix that defines the eigensystem, Ax=kx
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = this->_create_shell_mat(N,&M);       CHKERRQ(ierr);
  ierr = MatSetFromOptions(M);                CHKERRQ(ierr);
  ierr = MatShellSetOperation(M,MATOP_MULT,(void(*)())_MatMult_Stokes);CHKERRQ(ierr);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create the eigensolver and set various options
   Start from Non-Hermitian type: EPS_NHEP or EPS_HEP(symmetric)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);      CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,M,NULL);           CHKERRQ(ierr);
  if(use_init_space) ierr = EPSSetInitialSpace(eps,1,v0);
  ierr = EPSSetProblemType(eps,EPS_HEP);       CHKERRQ(ierr);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set the tolerance and maximum iteration for EPS solver and EPS type:
   EPSKRYLOVSCHUR(Default)/EPSARNOLDI/EPSPOWER/EPSLAPACK/EPSGD/EPSJD
   EPSBLOPEX/EPSRQCG/EPSLANCZOS/EPSSUBSPACE/
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //if(option=="largest") ierr = EPSSetType(eps,EPSPOWER);
  ierr = EPSSetType(eps,EPSKRYLOVSCHUR);
  ierr = EPSSetTolerances(eps,  tol,  maxits);  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"EPS tol = %f, maxits = %d\n",tol,maxits);CHKERRQ(ierr);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Select portion of spectrum.
   Note: EPS_LARGEST_MAGNITUDE should be used with EPSPOWER
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if(option=="smallest")  // LOBPCG for smallest eigenvalue problem!
  { ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);CHKERRQ(ierr);  }
  else if(option=="largest")
  { ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE); CHKERRQ(ierr);  }
  else
  { ierr = EPSSetFromOptions(eps);  CHKERRQ(ierr);  }
  // end if-else
  ierr = EPSSetFromOptions(eps);  CHKERRQ(ierr);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Solve the eigensystem and get the solution
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"EPS solve starts ...\n");
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n--->test: EPS solve results:\n");
  ierr = EPSGetEigenvalue(eps,0,&value,NULL); CHKERRQ(ierr);
  ierr = EPSGetEigenvector(eps,0,*v0,NULL); CHKERRQ(ierr);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Output the results on the screen
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"          Solution method: %s\n",type);
  ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"          Number of iterations of the method: %D\n",its);
  ierr = EPSGetConvergedReason(eps,&reason);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"          EPS converged reason : %D\n",reason);
  if(option=="smallest") {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"          The computed minimum eigenvalue: %f\n",value);
  }
  else if(option=="largest") {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"          The computed maximum eigenvalue: %f\n",value);
  }
  else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"          The computed eigenvalue: %f\n",value);
  }
//  ierr = PetscPrintf(PETSC_COMM_WORLD,"          the computed eigenvector:\n");
//  ierr = VecView(*v0,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); // View the vector */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"          \n\n");
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Distroy solution and clean up
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&M);  CHKERRQ(ierr);
  
  
  STOP_LOG ("compute_eigenvalue()", "BrownianSystem");
  return value;
#endif
}





// ======================================================================================
PetscErrorCode BrownianSystem::compute_eigenvalues(Real& eig_min,
                                                   Real& eig_max,
                                                   const Real& tol)
{
  PetscFunctionBeginUser;
  
  //PetscPrintf(PETSC_COMM_WORLD,"--->test in BrownianSystem: compute min/max eigenvalues.\n");
  eig_min = this->compute_eigenvalue("smallest",tol);
  eig_max = this->compute_eigenvalue("largest", tol);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Using initial space that the solver starts to iterate for small eigvalue.
   It seems using initial space doesn't improve the speed obviously.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//  Vec    v0;
//  PetscInt       N = _n_points*3;
//  this->_create_petsc_vec(N,&v0);
//  
//  bool use_init_space = false;
//  eig_max = this->compute_eigenvalue("largest", tol, use_init_space, &v0);
//  use_init_space = true;
//  VecScale(v0,0.00);
//  eig_min = this->compute_eigenvalue("smallest",tol, use_init_space, &v0);
//  VecDestroy(&v0);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Check the positiveness of the matrix
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if( eig_min<=0. || eig_max<=0. )
  {
    PetscPrintf(PETSC_COMM_WORLD,
                "*** warning in BrownianSystem::compute_eigenvalues():eigenvalue <= 0!");
    libmesh_error();
  }
  
  PetscFunctionReturn(0);
}


// ======================================================================================
Real BrownianSystem::power_iteration()
{
  Vec            v0, v1, w0;
  Mat            M;
  PetscInt       N = _n_points*3;
  PetscReal      tol = 1E-3, norm, eig_value;
  PetscErrorCode ierr;
  PetscScalar    one = 1.0;
  START_LOG ("power_iteration()", "BrownianSystem");
  PetscPrintf(PETSC_COMM_WORLD,"--->test in BrownianSystem: power_interation.\n");
  PetscPrintf(PETSC_COMM_WORLD,"--->                        compute max eigenvalue.\n");

  // Create ShellMat and Vec, and Initialize the v0 such that ||v0|| = 1
  this->_create_shell_mat(N, &M);
  this->_create_petsc_vec(N, &v0);
  ierr = VecSetFromOptions(v0); CHKERRQ(ierr);
  ierr = VecDuplicate(v0,&v1);  CHKERRQ(ierr); // Duplicate vectors
  ierr = VecDuplicate(v0,&w0);  CHKERRQ(ierr);
  
  // normalize v0 as the inital value
  ierr = VecSet(v0, one);            CHKERRQ(ierr);
  ierr = VecNorm(v0, NORM_2,&norm);  CHKERRQ(ierr);
  ierr = VecScale(v0, 1.0/norm);     CHKERRQ(ierr);
  
  // start the iteration
  std::size_t  maxits = 50;
  for(std::size_t i=0; i<maxits; ++i)
  {
    // apply A on v0 = v(k-1):  w0 = M*v0
    _MatMult_Stokes( M,v0,w0 );   // w0 -> w
    
    // normalize w
    ierr = VecNorm(w0, NORM_2,&norm);  CHKERRQ(ierr);
    ierr = VecScale(w0, 1.0/norm);     CHKERRQ(ierr);   // now w0 -> vk
    ierr = VecCopy(w0, v1); CHKERRQ(ierr); // keep a copy of w0, v1 = w0 (= vk)
    
    // compare v0 and w0 or in other word: v(k-1) and v(k)
    ierr = VecAXPY(w0, -1.0, v0); CHKERRQ(ierr); // w0 = w0 + a*v0, w0 is changed here!
    ierr = VecNorm(w0, NORM_2,&norm);  CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD," |v1 - v0| at iteration step %lu is %f\n",i,norm);
    ierr = VecCopy(v1, v0);       CHKERRQ(ierr); // let v0 = v1, and go to next iteration!
    
    // check if the iteration converges
    if(norm<tol)
    {
      _MatMult_Stokes( M,v0,w0 );   // w0 -> A*vk = M*v0
      ierr = VecDot(v0, w0, &eig_value); CHKERRQ(ierr);   // norm -> lambda_k
      PetscPrintf(PETSC_COMM_WORLD," Eigenvalue of iteration step %lu is %f\n",i,eig_value);
      break;
    } // end if
  } // end for i-loop
  
  // Destroy
  ierr = VecDestroy(&v0);  CHKERRQ(ierr);
  ierr = VecDestroy(&v1);  CHKERRQ(ierr);
  ierr = VecDestroy(&w0);  CHKERRQ(ierr);
  ierr = MatDestroy(&M );  CHKERRQ(ierr);
  
  // return
  STOP_LOG ("power_iteration()", "BrownianSystem");
  return eig_value;
}



// ======================================================================================
DenseMatrix<Number> BrownianSystem::chebyshev_transform_matrix(const std::size_t N) const
{
  const Real RN = Real(N);
  DenseMatrix<Number> Tmat(N+1,N+1);
  for (std::size_t i=0; i<=N; ++i)
  {
    Real ci = 1.;
    if(i==0 || i==N ) ci = 2.;
    for (std::size_t j=0; j<=N; ++j)
    {
      Real cj = 1.;
      if(j==0 || j==N ) cj = 2.;
      Tmat(i,j) = 2./( ci*cj*RN )*std::cos( libMesh::pi*Real(i)*Real(j)/RN );
      if( std::abs(Tmat(i,j))<1E-10 ) Tmat(i,j) = 0.0;    // filter
    }
  }
  
  return Tmat;
}



// ======================================================================================
DenseVector<Number> BrownianSystem::chebyshev_quadrature_point(const std::size_t N) const
{
  DenseVector<Number> Vx(N+1);
  for (std::size_t i=0; i<=N; ++i)
  {
    Vx(i) = std::cos( libMesh::pi*Real(i)/Real(N) );
    if( std::abs(Vx(i))<1E-10 ) Vx(i) = 0.0;    // filter
  }
  
  return Vx;
}



// ======================================================================================
bool BrownianSystem::chebyshev_polynomial_approximation(const std::size_t N,
                                                        const Real eigen_min,
                                                        const Real eigen_max,
                                                        const Real tol_cheb,
                                                        Vec* dw)
{
  Vec            x0, x1, xmid, dw_mid;
  Mat            M;
  PetscScalar    aux1, aux2, aux3, aux4; // auxiliary vars for error estimators.
  PetscScalar    error_cheb1=0., error_cheb2=0.;
  PetscInt       vsize;
  PetscErrorCode ierr;
  bool convergence = false;
  PetscFunctionBeginUser;
  START_LOG ("chebyshev_polynomial_approximation()", "BrownianSystem");
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create the shell matrix with dimension vsize x vsize and vectors that 
   are distributed in the same manner as the input dw
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = VecGetSize(*dw, &vsize);CHKERRQ(ierr);
  ierr = VecDuplicate(*dw,&x0);  CHKERRQ(ierr); // Duplicate vectors x0
  ierr = VecDuplicate(*dw,&x1);  CHKERRQ(ierr); // Duplicate vectors x1
  ierr = VecDuplicate(*dw,&xmid);CHKERRQ(ierr); // Duplicate vectors xmid
  ierr = VecDuplicate(*dw,&dw_mid);CHKERRQ(ierr); // Duplicate vectors dw_mid
  this->_create_shell_mat(vsize, &M);

  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   compute the coefficients da and db according to max/min eigenvalues
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const PetscScalar d_eig   = eigen_max - eigen_min;
  const PetscScalar da_cheb = 2.0/d_eig;
  const PetscScalar db_cheb = -(eigen_max + eigen_min)/d_eig;
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Compute coefficients of the Chebyshev series expansion X_j
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//  const DenseMatrix<Real> C_kj = this->chebyshev_transform_matrix(N);
//  DenseVector<Real> X_j  = this->chebyshev_quadrature_point(N);
//  DenseVector<Real> f_j(N+1);
//  for (std::size_t j=0; j<=N; ++j) {
//    f_j(j) = std::sqrt( (X_j(j) - db_cheb)/da_cheb );
//  }
//  C_kj.vector_mult(X_j,f_j);  // X_j = C_kj*f_j
//  //X_j.print(libMesh::out);
  
  Chebyshev chebyshev; // Gauss_Lobatto/Gauss/Gauss_Radau
  DenseVector<Real> X_j  = chebyshev.chebyshev_coefficients(N,da_cheb,db_cheb,"Gauss");
  //PetscPrintf(PETSC_COMM_WORLD,"\nTotal number of Chebyshev coefficients is %lu \n",X_j.size());
  //PetscPrintf(PETSC_COMM_WORLD,"--->test: Starting iterations for Chebyshev polynomial ...\n");
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Chebyshev polynomial approximation i = 0,1: x0 = dw, x1 = da*D*dw + db*dw
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = VecCopy(*dw,x0);CHKERRQ(ierr);       CHKERRQ(ierr);  // x0 = dw
  ierr = _MatMult_Stokes(M, x0, x1);          CHKERRQ(ierr);  // x1 = M*x0
//  ierr = _test_diagonal_mat(M, x0, x1);
  ierr = VecDot(x1, x0, &aux2);               CHKERRQ(ierr);  // dw*D*dw remains const.
  ierr = VecAXPBY(x1, db_cheb, da_cheb, x0);  CHKERRQ(ierr);  // x1 = da*x1 + db*x0
  //this->comm().barrier();
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   pre-calculate dw = c0*x0 + c1*x1, and the error estimator
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = VecSet(*dw, 0.0);          CHKERRQ(ierr); // set y = dw = 0
  ierr = VecAXPY(*dw,X_j(0),x0);    CHKERRQ(ierr); // y = dw += c0*x0
  ierr = VecNorm(*dw,NORM_2,&aux4); CHKERRQ(ierr); // aux4 = ||dw0|| = ||y_k||
  ierr = VecAXPY(*dw,X_j(1),x1);    CHKERRQ(ierr); // dw += c1*x1
  ierr = VecNorm(x1,NORM_2,&aux3);  CHKERRQ(ierr); // aux3 = ||x1||
  aux3 *= std::abs( X_j(1) );                      // aux3 = ||y_k - y_k-1|| = ||ck*xk||
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Error estimator 1 and 2. These operations don't change dw!
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = _MatMult_Stokes(M, *dw, dw_mid); CHKERRQ(ierr); // dw_mid = M*dw = M*B^-1*dw0 = B*dw0
  ierr = VecDot(dw_mid,dw_mid,&aux1);     CHKERRQ(ierr); // estimator for B^-1*dw
  //ierr = VecDot(*dw,*dw,&aux1);         CHKERRQ(ierr); // estimator for B*dw
  if (aux2 != 0.0) error_cheb1 = std::sqrt( std::abs(aux1 - aux2) / aux2 );
  if (aux4 != 0.0) error_cheb2 = aux3 / aux4 ;
  //PetscPrintf(PETSC_COMM_WORLD,"--->test: error_ch1 = %E, error_ch2 = %E at iteration step %lu\n",
  //            error_cheb1,error_cheb2,1);
  //PetscPrintf(PETSC_COMM_WORLD,"--->test: A1 = %E, A2 = %E, error_ch1 = sqrt(|A1 - A2|/A2) = %E\n",
  //            aux1,aux2,error_cheb1);
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   loop starts ...
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //std::size_t N_iter = 20;
  //N_iter = std::min(N, N_iter);
  for(std::size_t i=1; i<N; ++i)
  {
    ierr = VecNorm(*dw,NORM_2,&aux4);      CHKERRQ(ierr);  // aux4 = ||dw|| = ||y_k||
    ierr = _MatMult_Stokes(M, x1, xmid);   CHKERRQ(ierr);  // xmid = M*x1
//    ierr = _test_diagonal_mat(M, x1, xmid);CHKERRQ(ierr);
    
    // x2 = 2*da*xmid + 2*db*x1 - x0;
    ierr = VecAXPBY(xmid,2.*db_cheb,2.*da_cheb, x1);CHKERRQ(ierr);// xmid = 2*da*xmid + 2*db*x1
    ierr = VecAXPY(xmid, -1.0, x0);        CHKERRQ(ierr); // xmid = xmid - x0
    ierr = VecAXPY(*dw,X_j(i+1),xmid);     CHKERRQ(ierr); // dw += ck*xk
    ierr = VecNorm(xmid,NORM_2,&aux3);     CHKERRQ(ierr); // aux3 = ||xk||
    aux3 *= std::abs( X_j(i+1) );          // aux3 = ||y_k - y_k-1|| = ||ck*xk||
    
    /* Note: the following operation will not change dw, but only for error estimation */
    // check error: error estimator 1 - theoretical error
    ierr = _MatMult_Stokes(M, *dw, dw_mid); CHKERRQ(ierr); // dw_mid = M*dw = M*B^-1*dw0 = B*dw0
    ierr = VecDot(dw_mid,dw_mid,&aux1);     CHKERRQ(ierr); // estimator for B^-1*dw
    //ierr = VecDot(*dw,*dw,&aux1);         CHKERRQ(ierr); // estimator for B*dw
    if (aux2 != 0.0) error_cheb1 = std::sqrt( std::abs(aux1 - aux2) / aux2 );
    
    // check error: error estimator 2 - numerical relative error
    if (aux4 != 0.0) error_cheb2 = aux3 / aux4 ;
    //PetscPrintf(PETSC_COMM_WORLD,"--->test: error_ch1 = %E, error_ch2 = %E at iteration step %lu\n",
    //            error_cheb1,error_cheb2,i+1);
    //PetscPrintf(PETSC_COMM_WORLD,"--->test: A1 = %E, A2 = %E, error_ch1 = sqrt(|A1 - A2|/A2) = %E\n",
    //            aux1,aux2,error_cheb1);
    
    // Break when convergence occurs!
    if( error_cheb1<=tol_cheb )
    {
      convergence = true;
    //  PetscPrintf(PETSC_COMM_WORLD,"--->test: Chebyshev polynomial converges at step %lu:\n",i+1);
    //  PetscPrintf(PETSC_COMM_WORLD,"          error1 = %E, error2 = %E!\n",error_cheb1,error_cheb2);
      break;
    }
    
    // if reach the relative tol, break anyway! But theoretical error is large.
    if(error_cheb2<=1E-9) break;
    
    ierr = VecCopy(x1,x0);  CHKERRQ(ierr);  // x0 = x1 = x_k-1
    ierr = VecCopy(xmid,x1);CHKERRQ(ierr);  // x1 = xmid = x_k, which will be used to calculate x_k+1
  } // end for i-loop
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Print out the warning message if not converged!
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if(!convergence)
    PetscPrintf(PETSC_COMM_WORLD,"--->test: Chebyshev polynomial fails to converge in this case!\n");
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   destroy and return
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = VecDestroy(&x0);     CHKERRQ(ierr);
  ierr = VecDestroy(&x1);     CHKERRQ(ierr);
  ierr = VecDestroy(&xmid);   CHKERRQ(ierr);
  ierr = VecDestroy(&dw_mid); CHKERRQ(ierr);
  ierr = MatDestroy(&M);      CHKERRQ(ierr);
  
  STOP_LOG ("chebyshev_polynomial_approximation()", "BrownianSystem");
  PetscFunctionReturn(convergence);
}


// ======================================================================================
PetscErrorCode BrownianSystem::_test_diagonal_mat(Vec f, Vec u)
{
  PetscScalar    *pu;
  PetscInt       low, high, i, j=0;
  PetscMPIInt    size, rank ;
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  
  // copy u
  ierr = VecCopy(f,u);        CHKERRQ(ierr);
  ierr = VecGetArray(u,&pu);  CHKERRQ(ierr);
  
  // rank and size, local ownership
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(u,&low,&high);   CHKERRQ(ierr);
  //printf("low = %d, high = %d on the rank %d \n", low,high,rank);
  
  for(i=low; i<high; ++i) {
    pu[j] *= PetscScalar(i+1); j++;
    //printf("--->test in _test_diagonal_mat: i = %d on rank %d\n",i+1,rank);
  }
  ierr = MPI_Barrier(PETSC_COMM_WORLD);
  
  // restore pu
  ierr = VecRestoreArray(u,&pu); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}



// ======================================================================================
PetscErrorCode BrownianSystem::_vector_mean_variance(Vec u,
                                                     Real& mean_u,
                                                     Real& variance_u) const
{
  PetscErrorCode  ierr;
  PetscInt        N;
  Vec             v;
  PetscScalar     mean = 0., variance = 0.;
  PetscFunctionBeginUser;
  
  // Compute the mean and variance value of Vec u.
  ierr = VecGetSize(u,&N);   CHKERRQ(ierr);
  ierr = VecDuplicate(u,&v); CHKERRQ(ierr);  // init Vec v.
  ierr = VecSum(u,&mean);    CHKERRQ(ierr);  mean /= PetscScalar(N);
  ierr = VecSet(v,mean);     CHKERRQ(ierr);
  ierr = VecAYPX(v,-1.0,u);  CHKERRQ(ierr);
  ierr = VecNorm(v,NORM_2,&variance);CHKERRQ(ierr);
  ierr = VecDestroy(&v);      CHKERRQ(ierr);
  variance = variance*variance/PetscScalar(N);
  
  // Print out the values
//  PetscPrintf(PETSC_COMM_WORLD,"--->test in petsc_random_vector:       mean = %f, variance = %f\n",
//              mean, variance);
//  PetscPrintf(PETSC_COMM_WORLD,"Exact values for uniform distribution: mean = %f, variance = %f\n",
//              0.5, 1./12.);
  
  mean_u      = Real(mean);
  variance_u  = Real(variance);
  
  PetscFunctionReturn(0);
}



// ======================================================================================
PetscErrorCode BrownianSystem:: _vector_mean_variance(const std::vector<Real>& randn) const
{
  PetscFunctionBeginUser;
  Real mean = 0., variance = 0., deviation = 0.;
  const std::size_t N = randn.size();
  
  // mean value
  for (std::size_t i=0; i<N; ++i)
  {
    mean += randn[i];
  }
  mean /= Real(N);
  
  // variance & standard deviation
  for (std::size_t i=0; i<N; ++i) variance += (randn[i] - mean)*(randn[i] - mean);
  variance /= Real(N);
  deviation = std::sqrt(variance);
    
  // print out info
  //if(this->comm().rank()==0)
  //{
  //  printf("mean = %f, variance = %f, and deviation = %f on the rank %d\n",
  //         mean, variance, deviation, this->comm().rank());
  //}
  
  PetscFunctionReturn(0);
}



// ======================================================================================
PetscErrorCode BrownianSystem::_create_shell_mat(const std::size_t N,
                                                 Mat* M)
{
  // declare the variables
  PetscInt       n, dn, Nmat;
  PetscMPIInt    size, rank;
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create the shell matrix
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Nmat = PetscInt(N); // tranform from std::size_t to PetscInt
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  n  = Nmat/size + 1;
  dn = n*size - Nmat;
  if ( dn>0 && rank<dn ) n -= 1;
  ierr = MatCreateShell(PETSC_COMM_WORLD, n, n, Nmat, Nmat, this, M);  CHKERRQ(ierr);
  
//  ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);
//  printf("--->test in _create_shell_mat(): n = %d, N = %d, rank = %d\n",(int)n, (int)N, (int)rank);
  
  PetscFunctionReturn(0);
}



// ======================================================================================
PetscErrorCode BrownianSystem:: _create_petsc_vec(const std::size_t N, Vec* V)
{
  // declare the variables
  PetscInt       n, dn, Nmat;
  PetscMPIInt    size, rank;
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create the petsc vector
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Nmat = PetscInt(N); // tranform from std::size_t to PetscInt
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  n  = Nmat/size + 1;
  dn = n*size - Nmat;
  if ( dn>0 && rank<dn ) n -= 1;
  ierr = VecCreate(PETSC_COMM_WORLD, V);  CHKERRQ(ierr);
  ierr = VecSetSizes(*V, n, N); CHKERRQ(ierr);
  ierr = VecSetFromOptions(*V);          CHKERRQ(ierr);
//  printf("--->test in _create_petsc_vec(): n = %d, N = %d, rank = %d\n",(int)n, (int)N, (int)rank);
  
  PetscFunctionReturn(0);
}




//} //end of namespace
