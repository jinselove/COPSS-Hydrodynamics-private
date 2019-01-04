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
#include "ggem_system.h"
#include "brownian_system.h"
#include "analytical_solution.h"
#include "pm_system_np.h"


namespace libMesh
{

// ======================================================================================
PMSystemNP::PMSystemNP(EquationSystems& es,
                               const std::string& name,
                               const unsigned int number)
: PMLinearImplicitSystem (es, name, number),
  _np_solver(es)
{
  // Stokes equation assembly
  _assemble_np = ( new AssembleNP(es) );
}



// ==================================================================================
PMSystemNP::~PMSystemNP()
{
  // Clear data
  this->clear();
}



// ==================================================================================
void PMSystemNP::clear ()
{
  // delete the pointer
  if(_assemble_np) {
    delete _assemble_np;
  }
}



// ===========================================================
void PMSystemNP::assemble_matrix(const std::string& system_name,
             			     const std::string& option)
{
  libmesh_assert (this->matrix);
  libmesh_assert (this->matrix->initialized());

  START_LOG("assemble_matrix()", "PMSystemNP");

  STOP_LOG("assemble_matrix()", "PMSystemNP");
}



// ==================================================================================
void PMSystemNP::assemble_rhs(const std::string& system_name,
                                  const std::string& option)
{
  libmesh_assert (this->rhs);
  libmesh_assert (this->rhs->initialized());

  START_LOG("assemble_rhs()", "PMSystemNP");

  STOP_LOG("assemble_rhs()", "PMSystemNP");
}



// ==================================================================================
void PMSystemNP::solve(const std::string& option,
                            const bool& re_init)
{
  START_LOG("solve()", "PMSystemNP");

  STOP_LOG("solve()", "PMSystemNP");
}

} // end of namespace
