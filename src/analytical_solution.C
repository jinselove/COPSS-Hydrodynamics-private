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
#include "analytical_solution.h"
#include "point_mesh.h"

// ======================================================================
AnalyticalSolution::AnalyticalSolution(const std::string& name)
{
    if (name != "Stokes") libmesh_error();
}


// ======================================================================
AnalyticalSolution::~AnalyticalSolution()
{
    //do nothing
}


// ======================================================================
void AnalyticalSolution::attach_point_mesh(PointMesh<3>* point_mesh) 
{
    START_LOG("attach_point_mesh()", "AnalyticalSolution");  
    
    _point_mesh = point_mesh;

    STOP_LOG("attach_point_mesh()", "AnalyticalSolution")  
}


// ======================================================================
PointMesh<3>* AnalyticalSolution::get_point_mesh() 
{
    START_LOG("get_point_mesh()", "AnalyticalSolution");  
    
    return _point_mesh;

    STOP_LOG("get_point_mesh()", "AnalyticalSolution");    
}


// ======================================================================
std::vector<Real> AnalyticalSolution::exact_solution_infinite_domain(GGEMStokes& ggem_stokes,
                                                                     const Point& pt0) const
{
  
  START_LOG("exact_solution_infinite_domain()", "AnalyticalSolution");
  
  std::vector<Real> UA(dim, 0.);
  
  DenseMatrix<Number> GT;
  
  // GGEM object and number of points in the system
  // GGEMStokes ggem_stokes;
  
  const std::size_t n_points = _point_mesh->num_particles();
  
  // loop over each point
  for(std::size_t i=0; i<n_points; ++i)
  {
    const Point pti = _point_mesh->particles()[i]->point();
    const Point x   = pt0 - pti;
    
    bool  zero_limit  = false;
    if(x.norm()<1E-6) zero_limit  = true;
    
    // use ksi instead of alpha
    GT = ggem_stokes.green_tensor_unbounded_smoothed(x, ggem_stokes.get_ksi(), zero_limit);
    const Point fv = _point_mesh->particles()[i]->particle_force();
    
    // 3. compute u due to this particle
    for (std::size_t k=0; k<dim; ++k){
      for (std::size_t l=0; l<dim; ++l){
        UA[k] += GT(k,l)*fv(l);
      } // end for l
    } // end for k
  } // end for i
  
  STOP_LOG("exact_solution_infinite_domain()", "AnalyticalSolution");
  return UA;
}


// ======================================================================
Real AnalyticalSolution::correction_factor_bohlin(const Real r_r0) const
{
  START_LOG("correction_factor_bohlin()", "AnalyticalSolution");
  
  Real factor = 1.0 - 2.10443*r_r0 + 2.08877*std::pow(r_r0,3)
              - 0.94813*std::pow(r_r0,5) - 1.372*std::pow(r_r0,6)
              + 3.87*std::pow(r_r0,8) - 4.19*std::pow(r_r0,10);
  
  STOP_LOG("correction_factor_bohlin()", "AnalyticalSolution");
  return 1.0/factor;
}


// ======================================================================
Real AnalyticalSolution::correction_factor_haberman(const Real r_r0) const
{
  START_LOG("correction_factor_haberman()", "AnalyticalSolution");
  
  Real f1 = 1.0 - 0.75857*std::pow(r_r0,5);
  Real f2 = 1.0 - 2.1050*r_r0 + 2.0865*std::pow(r_r0,3)
          - 1.7068*std::pow(r_r0,5) + 0.72603*std::pow(r_r0,6);
  
  STOP_LOG("correction_factor_haberman()", "AnalyticalSolution");
  return f1/f2;
}
