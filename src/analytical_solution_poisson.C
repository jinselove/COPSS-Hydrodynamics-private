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
#include "analytical_solution_poisson.h"
#include "point_mesh.h"

// ======================================================================
AnalyticalSolutionPoisson::AnalyticalSolutionPoisson(const std::string& name)
{
  if (name != "Poisson") libmesh_error();
}

// ======================================================================
AnalyticalSolutionPoisson::~AnalyticalSolutionPoisson()
{
  // do nothing
}

// ======================================================================
void AnalyticalSolutionPoisson::attach_point_mesh(PointMesh<3> *point_mesh)
{
  START_LOG("attach_point_mesh()", "AnalyticalSolutionPoisson");

  _point_mesh = point_mesh;

  STOP_LOG("attach_point_mesh()", "AnalyticalSolutionPoisson")
}

// ======================================================================
PointMesh<3> * AnalyticalSolutionPoisson::get_point_mesh()
{
  START_LOG("get_point_mesh()", "AnalyticalSolutionPoisson");

  return _point_mesh;

  STOP_LOG("get_point_mesh()", "AnalyticalSolutionPoisson");
}

// ======================================================================
Real AnalyticalSolutionPoisson::exact_solution_infinite_domain(
  GGEMPoisson& ggem_poisson,
  const Point& pt0)
const
{
  START_LOG("exact_solution_infinite_domain()", "AnalyticalSolutionPoisson");

  // std::vector<Real> UA(dim, 0.;
  Real phi = 0.;

  Real GT;

  // GGEM object and number of points in the system
  const std::size_t n_points = _point_mesh->num_particles();

  // loop over each point
  for (std::size_t i = 0; i < n_points; ++i)
  {
    const Point pti = _point_mesh->particles()[i]->point();
    const Point x   = pt0 - pti;

    // use ksi instead of alpha
    GT = ggem_poisson.green_tensor_unbounded_smoothed(x, ggem_poisson.get_ksi());
    Real q = _point_mesh->particles()[i]->charge();

    // 3. compute potential field due to this point charge
    phi += GT * q;
  } // end for i

  STOP_LOG("exact_solution_infinite_domain()", "AnalyticalSolutionPoisson");
  return phi;
}

// ======================================================================
Point AnalyticalSolutionPoisson::exact_solution_infinite_domain_grad(
  GGEMPoisson& ggem_poisson, const Point& pt0) const
{
  START_LOG("exact_solution_infinite_domain_grad()",
    "AnalyticalSolutionPoisson");

  // std::vector<Real> UA(dim, 0.;
  Point phi_grad, GT_grad;

  // GGEM object and number of points in the system
  const std::size_t n_points = _point_mesh->num_particles();

  // loop over each point
  for (std::size_t i = 0; i < n_points; ++i)
  {
    const Point pti = _point_mesh->particles()[i]->point();
    const Point x   = pti - pt0;

    // use ksi instead of alpha
    GT_grad = ggem_poisson.green_tensor_unbounded_smoothed_grad(x, ggem_poisson
    .get_ksi());
    Real q = _point_mesh->particles()[i]->charge();

    // 3. compute potential field due to this point charge
    phi_grad += GT_grad * q;
  } // end for i

  STOP_LOG("exact_solution_infinite_domain_grad()",
    "AnalyticalSolutionPoisson");
  return phi_grad;
}
