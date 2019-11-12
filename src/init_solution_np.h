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

// C++ Includes
#include <math.h>

// Mesh library includes
#include "libmesh/libmesh_common.h"

// Bring in everything from the libMesh namespace
namespace libMesh
{
    /**
     * This function gives the initial solution, i.e., t=0., at point p
     *
     */
Real init_solution(const Point& p,
                   const Parameters& parameters)
{
  // where the instant source is at
  Point p0(0.);

  Real num = 0.;
  for (int dim_i=0; dim_i<3; dim_i++)
  {
    num += pow(p(dim_i) - p0(dim_i), 2.);
  }

  // den = D * (4. * t + 1.)
  const Real den = parameters.get<Real>("ion_diffusivity");

  // return exp(-num/den) / pow(4.*t + 1., 1.5)
  return exp(-num / den);
}
}

