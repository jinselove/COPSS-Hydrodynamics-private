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
#include "ggem_stokes.h"

namespace libMesh
{

/*! \brief Analytic solution to unbounded flow field

 * Define analytical solutions that are 
 * available for some special cases.
 * (This is only used for validation and test purpose.)
 */

class AnalyticalSolutionStokes : public ReferenceCountedObject<AnalyticalSolutionStokes>
//public ParallelObject
{
public:

  /*! \brief Constructor

  */
  AnalyticalSolutionStokes(const std::string& name);


  /*! \brief Destructor

  */
  ~AnalyticalSolutionStokes();
  
  
  /*! \brief Exact solution for point forces in an unbounded domain

  */
  std::vector<Real> exact_solution_infinite_domain(GGEMStokes& ggem_stokes,
                                                   const Point& pt0) const;
  
  
  /*! \brief correction factor for a particle in a cylinder: Bohlin approximation
   *
   */
  Real correction_factor_bohlin(const Real r_ratio) const;
  
  
  /*! \brief correction factor for a particle in a cylinder: Haberman approximation
  *
  */
  Real correction_factor_haberman(const Real r_ratio) const;
  
  
  /*! \brief Attach point mesh to the class
  *
  */
  void attach_point_mesh(PointMesh<3>* point_mesh);
  
  
  /*! \brief Get the pointer to the class member, point_mesh
  *
  */
  PointMesh<3>* get_point_mesh();
  
  
private:
    
    // Initialization a null _point_mesh pointer
    PointMesh<3>* _point_mesh = NULL;
    
    
    // System dimension
    const int dim = 3;

}; // end of class defination



}  // end of namespace
