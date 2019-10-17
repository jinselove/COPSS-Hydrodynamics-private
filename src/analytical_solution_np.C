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

// c++ headers
#include <cmath>

// libmesh headers
#include "libmesh/dense_matrix.h"
#include "libmesh/point.h"

// user defined headers
#include "pm_toolbox.h"
#include "analytical_solution_np.h"
#include "point_mesh.h"

// ======================================================================
AnalyticalSolutionNP::AnalyticalSolutionNP(const std::string& name)
{
  if (name != "NP") libmesh_error();
}

// ======================================================================
AnalyticalSolutionNP::~AnalyticalSolutionNP()
{
  // do nothing
}

// ======================================================================
void AnalyticalSolutionNP::attach_point_mesh(PointMesh<3> *point_mesh)
{
  START_LOG("attach_point_mesh()", "AnalyticalSolutionNP");

  _point_mesh = point_mesh;

  STOP_LOG("attach_point_mesh()", "AnalyticalSolutionNP")
}

// ======================================================================
PointMesh<3> * AnalyticalSolutionNP::get_point_mesh()
{
  START_LOG("get_point_mesh()", "AnalyticalSolutionNP");

  return _point_mesh;

  STOP_LOG("get_point_mesh()", "AnalyticalSolutionNP");
}

// ======================================================================
Real AnalyticalSolutionNP::exact_solution_infinite_domain(const Point& pt0,
  const Real& t)
const
{
  START_LOG("exact_solution_infinite_domain()", "AnalyticalSolutionNP");

  STOP_LOG("exact_solution_infinite_domain()", "AnalyticalSolutionNP");
  return 0;
}
