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

// libmesh head files
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmv_io.h"
#include "libmesh/vtk_io.h"

// local header files
#include "pm_toolbox.h"
#include "ggem_poisson.h"
#include "brownian_system.h"
#include "analytical_solution.h"
#include "pm_system_poisson.h"


namespace libMesh
{

// ======================================================================================
PMSystemPoisson::PMSystemPoisson(EquationSystems& es,
                               const std::string& name,
                               const unsigned int number)
: PMLinearImplicitSystem (es, name, number),
  _poisson_solver(es)
{
  // Stokes equation assembly
  _assemble_poisson = ( new AssemblePoisson(es) );
}



// ==================================================================================
PMSystemPoisson::~PMSystemPoisson()
{
  // Clear data
  this->clear();
}



// ==================================================================================
void PMSystemPoisson::clear ()
{
  // delete the pointer
  if(_assemble_poisson) {
    delete _assemble_poisson;
  }
}



// ===========================================================
void PMSystemPoisson::assemble_matrix(const std::string& system_name,
             			     const std::string& option)
{
  libmesh_assert (this->matrix);
  libmesh_assert (this->matrix->initialized());

  START_LOG("assemble_matrix()", "PMSystemPoisson");

  STOP_LOG("assemble_matrix()", "PMSystemPoisson");
}



// ==================================================================================
void PMSystemPoisson::assemble_rhs(const std::string& system_name,
                                  const std::string& option)
{
  libmesh_assert (this->rhs);
  libmesh_assert (this->rhs->initialized());

  START_LOG("assemble_rhs()", "PMSystemPoisson");

  STOP_LOG("assemble_rhs()", "PMSystemPoisson");
}



// ==================================================================================
void PMSystemPoisson::solve(const std::string& option,
                            const bool& re_init)
{
  START_LOG("solve()", "PMSystemPoisson");

  STOP_LOG("solve()", "PMSystemPoisson");
}

} // end of namespace
