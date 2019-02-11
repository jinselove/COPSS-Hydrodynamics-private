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


// std C++
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <time.h>

// libmesh headers
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmv_io.h"
#include "libmesh/vtk_io.h"

// particle-mesh header files
#include "pm_toolbox.h"
#include "brownian_system.h"
#include "pm_linear_implicit_system.h"


namespace libMesh
{

// ======================================================================================
PMLinearImplicitSystem::PMLinearImplicitSystem(EquationSystems& es,
                                               const std::string& name,
                                               const unsigned int number)
: LinearImplicitSystem (es, name, number),
  _point_mesh(NULL),
  _particle_mesh(NULL)
{
}



// ==================================================================================
PMLinearImplicitSystem::~PMLinearImplicitSystem()
{
  // Clear the parent data
  LinearImplicitSystem::clear();
}



// ==================================================================================
void PMLinearImplicitSystem::write_out_single_particle(const Point& coords,
                                                       const std::vector<Real>& vel,
                                                       const int i_step,
                                                       const Real time) const
{
  START_LOG("write_out_single_particle()", "PMLinearImplicitSystem");

  this->comm().barrier();
  std::size_t dim = this->get_mesh().mesh_dimension();
  std::string filename = "single_particle_history.txt";
  std::ofstream outfile;
  int o_width = 15, o_precision = 9;

  // the first step, ios_base::out, File open for writing
  if(i_step<=0)
  {
    outfile.open(filename,std::ios_base::out);
    outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
    outfile.precision(o_precision);   outfile.width(o_width/3);
    if( this->comm().rank()==0 ) outfile << i_step << "  " << time;

    for ( std::size_t i=0; i<dim; ++i )
    {
      outfile.setf(std::ios::right);  outfile.setf(std::ios::fixed);
      outfile.precision(o_precision); outfile.width(o_width);
      if( this->comm().rank()==0 )    outfile << coords(i);
    }
    for ( std::size_t i=0; i<dim; ++i )
    {
      outfile.setf(std::ios::right);  outfile.setf(std::ios::fixed);
      outfile.precision(o_precision); outfile.width(o_width);
      if( this->comm().rank()==0 )    outfile << vel[i];
    }
    if( this->comm().rank()==0 ) outfile << "\n";
  }

  // the following steps, ios_base::app, output operations happen at the end of the file,
  // appending to its existing contents.
  if(i_step>0)
  {
    outfile.open(filename,std::ios_base::app);
    outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
    outfile.precision(o_precision);   outfile.width(o_width);
    if( this->comm().rank()==0 ) outfile << i_step << "  " << time;

    for ( std::size_t i=0; i<dim; ++i )
    {
      outfile.setf(std::ios::right);  outfile.setf(std::ios::fixed);
      outfile.precision(o_precision); outfile.width(o_width);
      if( this->comm().rank()==0 )    outfile << coords(i);
    }
    for ( std::size_t i=0; i<dim; ++i )
    {
      outfile.setf(std::ios::right);  outfile.setf(std::ios::fixed);
      outfile.precision(o_precision); outfile.width(o_width);
      if( this->comm().rank()==0 )    outfile << vel[i];
    }
    if( this->comm().rank()==0 ) outfile << "\n";
  }

  // close the file
  outfile.close();
  this->comm().barrier();

  STOP_LOG("write_out_single_particle()", "PMLinearImplicitSystem");
}



// ==================================================================================
void PMLinearImplicitSystem::write_out_point_coordinate(Vec* vin,
                                                        const std::size_t istep,
                                                        const Real& time,
                                                        const std::string& filename,
                                                        const std::string& openmode) const
{
  START_LOG("write_out_particle_coordinate()", "PMLinearImplicitSystem");

  VecScatter ctx;
  Vec vout;
  const PetscScalar     *px;

  //  PetscPrintf(PETSC_COMM_WORLD,"--->test in write_out_particle_coordinate vin = \n");
  //  VecView(*vin,PETSC_VIEWER_STDOUT_WORLD);  // View the random vector

  VecScatterCreateToZero(*vin,&ctx,&vout);
  VecScatterBegin(ctx,*vin,vout,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,*vin,vout,INSERT_VALUES,SCATTER_FORWARD);
  VecGetArrayRead(vout,&px);

  //  PetscPrintf(PETSC_COMM_WORLD,"--->test in write_out_particle_coordinate vout = \n");
  //  VecView(vout,PETSC_VIEWER_STDOUT_WORLD);  // View the random vector

  //
  const std::size_t NP = _point_mesh->num_particles();
  const std::size_t dim = this->get_mesh().mesh_dimension();

  std::ofstream outfile;
  const int o_width = 5, o_precision = 9;
  if(this->comm().rank()==0)
  {
    if(openmode=="out"){
      outfile.open(filename,std::ios_base::out);
    }
    else if(openmode=="app"){
      outfile.open(filename,std::ios_base::app);
    }
    else{
      outfile.open(filename,std::ios_base::app);
    }
    // end of if-else

    if(istep==0) outfile << NP << "\n";
    outfile << "Step " << istep << " time = " <<time<< "\n";
    for (std::size_t i=0; i<NP; ++i)
    {
      outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
      outfile.precision(o_precision);   outfile.width(o_width);
      outfile << i << " ";
      for(std::size_t j=0; j<dim; ++j)
      {
        Real xyz = Real( px[i*dim+j] );
        outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
        outfile.precision(o_precision);   outfile.width(o_width);
        outfile << xyz << "  ";
      }
      outfile << "\n";
    }

    outfile << "End Step " << istep << "\n\n";
    outfile.close();
  }

  // restore vector and destroy ctx
  VecRestoreArrayRead(vout,&px);
  VecScatterDestroy(&ctx);
  VecDestroy(&vout);

  STOP_LOG("write_out_particle_coordinate()", "PMLinearImplicitSystem");
}



// ==================================================================================
void PMLinearImplicitSystem::write_point_csv(const std::string& filename,
                                             const std::vector<Real>& pv,
                                             const bool write_velocity) const
{
  START_LOG("write_point_csv()", "PMLinearImplicitSystem");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Write out the CSV file on the 0-th processor
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const std::size_t NP  = _point_mesh->num_particles();
  const std::size_t dim = this->get_mesh().mesh_dimension();
  std::ofstream outfile;
  const int o_width = 5, o_precision = 9;
  if(this->comm().rank()==0)
  {
    outfile.open(filename,std::ios_base::out);
    outfile << "Point #, "<< "X Coord, "<< "Y Coord, "<< "Z Coord";
    if(write_velocity)
      outfile << ",  Vx,      "<< "Vy,       "<< "Vz,      "<< "Vmag ";
    outfile << "\n";

    for (std::size_t i=0; i<NP; ++i)
    {
      // write the point id and coordinates
      outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
      outfile.precision(o_precision);   outfile.width(o_width);
      outfile << i << " ";
      for(std::size_t j=0; j<dim; ++j)
      {
        const Real xyz = _point_mesh->particles()[i]->point()(j);
        outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
        outfile.precision(o_precision);   outfile.width(o_width);
        outfile << ",  "<< xyz;
      } // end for j-loop

      // write the point velocity if needed
      if(write_velocity)
      {
        Real Vmag = 0.0;
        for(std::size_t j=0; j<dim; ++j)
        {
          const Real vxyz = pv[i*dim + j];
          outfile.setf(std::ios::right);    outfile.setf(std::ios::fixed);
          outfile.precision(o_precision);   outfile.width(o_width);
          outfile << ",  "<< vxyz;
          Vmag += vxyz*vxyz;
        }
        Vmag = std::sqrt(Vmag);
        outfile << ",  "<< Vmag;
      }

      outfile << "\n";
    }
    outfile.close();
  }

  STOP_LOG("write_point_csv()", "PMLinearImplicitSystem");
}



// ==================================================================================
PetscErrorCode PMLinearImplicitSystem::write_point_csv(const std::string& filename,
                                                       Vec * petsc_vec,
                                                       const bool write_velocity) const
{
  PetscInt          low, high, nlocal, vsize;
  PetscScalar       *pv;
  PetscErrorCode    ierr;
  PetscFunctionBeginUser;
  START_LOG ("write_point_csv()", "PMLinearImplicitSystem");

  std::vector<Real> std_vec;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Get the Ownership range and local components, then copy to a std::vector!
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (petsc_vec && write_velocity)
  {
    ierr = VecGetSize(*petsc_vec, &vsize);CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(*petsc_vec,&low,&high); CHKERRQ(ierr);
    ierr = VecGetLocalSize(*petsc_vec,&nlocal);         CHKERRQ(ierr);
    ierr = VecGetArray(*petsc_vec,&pv);                 CHKERRQ(ierr);

    std_vec.resize( (std::size_t)nlocal );
    for(int i=0; i<nlocal; ++i) std_vec[i] = pv[i];
    this->comm().allgather(std_vec);
    ierr = VecRestoreArray(*petsc_vec,&pv);             CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Call the member function to write out point CSV file
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  this->write_point_csv(filename, std_vec, write_velocity);

  STOP_LOG("write_point_csv()", "PMLinearImplicitSystem");
  PetscFunctionReturn(0);
}

} // end namespace libMesh
