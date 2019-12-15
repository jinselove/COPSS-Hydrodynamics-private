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

#include "libmesh/reference_counted_object.h"
#include "libmesh/linear_implicit_system.h"

#include "point_mesh.h"
#include "ggem_poisson.h"

namespace libMesh {
/*! \brief Analytic solution to unbounded potential field

 * Define analytical solutions that are
 * available for some special cases.
 * (This is only used for validation and test purpose.)
 */

class AnalyticalSolutionPoisson : public ReferenceCountedObject<AnalyticalSolutionPoisson>

// public ParallelObject
{
public:

  /*! \brief Constructor

   */
  AnalyticalSolutionPoisson(const std::string& name);


  /*! \brief Destructor

   */
  ~AnalyticalSolutionPoisson();


  /*! \brief Exact solution for point forces in an unbounded domain

   */
  Real exact_solution_infinite_domain(GGEMPoisson& ggem_poisson,
                                      const Point& pt0) const;


  /**
   * calculate the gradient of the exact solution at pt0
   */
  Point exact_solution_infinite_domain_grad(GGEMPoisson& ggem_poisson,
                                            const Point& pt0) const;


  /*! \brief Attach point mesh to the class
   *
   */
  void attach_point_mesh(PointMesh<3> *point_mesh);


  /*! \brief Get the pointer to the class member, point_mesh
   *
   */
  PointMesh<3>* get_point_mesh();

private:

  // Initialization a null _point_mesh pointer
  PointMesh<3> *_point_mesh = nullptr;


  // System dimension
  const int dim = 3;
}; // end of class defination
}  // end of namespace
