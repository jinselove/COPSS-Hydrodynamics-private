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

namespace libMesh {
/*! \brief Analytic solution to a simple convection-diffusion system where an
 * instant ion concentration source happened at time 0 in an unconfined
 * domain of size 2*2*2, and there is a constant velocity that induces
 * convection.
 *
 * governing equation:
 *     \frac{\partial c}{\partial t} = D * \nabla^2c-\mathbf{u}\cdot\nabla c
 *
 * Initial condition:
 *     c(\mathbf{x},t=0)=c_{\text{exact}}(\mathbf{x},t=0)
 *
 * Boundary condition:
 *      c(\mathbf{x},t)|_{\mathbf{x}\in\partial \Omega}=c_{\text{exact}}(\mathbf{x},t)

 * c_{\text{exact}}(\mathbf{x}, t) is the analytical solution of the system:
 *       c_{\text{exact}}(\mathbf{x},t) = \frac{1}{(4t+1)^{3/2}}e^{\frac{-(\mathbf{x}-\mathbf{x_0}-\mathbf{u}t)^2}{D(4t+1)}}\text{, where }\mathbf{x_0}=(1,1,1)

 * where c(\mathbf{x}, t) is the ion concentration at location \mathbf{x} at time t;
 * \mathbf{u} is a constant fluid velocity;
 * D is the diffusion coefficient of this ion species;

 * Define analytical solutions that are
 * available for some special cases.
 * (This is only used for validation and test purpose.)
 */

class AnalyticalSolutionNP : public ReferenceCountedObject<AnalyticalSolutionNP>

// public ParallelObject
{
public:

    /*! \brief Constructor

     */
    AnalyticalSolutionNP(const std::string& name);


    /*! \brief Destructor

     */
    ~AnalyticalSolutionNP();


    /*! \brief Exact solution at a point in the domain at time t
     *
     */
    Real exact_solution_infinite_domain(const Point& pt,
                                        const Real& t,
                                        const Real& D_ion) const;


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
