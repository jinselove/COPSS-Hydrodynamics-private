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
#include <vector>

#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/elem.h"


using libMesh::MeshTools::BoundingBox;
using libMesh::Point;

namespace libMesh
{

  
  /*
   * The class is designed for dealing with the periodic boundary
   * of the ParticleMesh system
   */
  
  
class PMPeriodicBoundary  :  public ReferenceCountedObject<PMPeriodicBoundary>
{
public:
  // Constructor
  PMPeriodicBoundary();
  
  // Destructor
  ~PMPeriodicBoundary();
  
  // Constructor
  PMPeriodicBoundary(const std::pair<Point, Point> &bbox_pts,
                     const std::vector<bool>& periodic_directions);
  
  // Constructor
  PMPeriodicBoundary(const Point &bbox_pmin,
                     const Point &bbox_pmax,
                     const std::vector<bool>& periodic_directions);

  // Constructor
  PMPeriodicBoundary(const Point &bbox_pmin,
		     const Point &bbox_pmax,
		     const std::vector<bool>& periodic_directions,
		     const std::vector<bool>& inlet_directions,
                     const std::vector<Real>& inlet_pressure);
  
  // Constructor
  PMPeriodicBoundary(const BoundingBox& bbox,
                     const std::vector<bool>& periodic_directions);
  
  // Copy Constructor
  PMPeriodicBoundary(const PMPeriodicBoundary & pb);
  
  
  // periodic direction
  bool periodic_direction(const std::size_t i) const;
  const std::vector<bool>& periodic_direction() const;
 
  // inlet direction
  bool inlet_direction (const std::size_t i) const;
  const std::vector <bool>& inlet_direction() const;
  const Real inlet_pressure(const std::size_t i) const;
  const std::vector<Real>& inlet_pressure() const; 
  
  // size of the Bounding box
  Real  box_length (const std::size_t i) const;
  Point box_length () const;
  
  
  // Box min and max boundaries and the center
  Point box_min() const {  return _bbox.min(); };
  Point box_max() const {  return _bbox.max(); };
  Point box_mid() const {  return 0.5*( _bbox.max() + _bbox.min() ); };
  
  
  // Return the bounding box of the domain
  const BoundingBox & bounding_box() const { return _bbox; }
  
  
  
  /* 
   * Compute the closest image point in x(y,z) - directions
   * If there is no image in the specified direction, return "false"
   * This is a simple version that is not used in this code.
   */
  bool get_image_point(const Point& pt0,
                       const std::size_t i, // i=0,1,2 for x,y,z direction
                       Point& im_pt) const;
  
  
  /* 
   * Compute the closest image point in xyz-directions
   * If there is no image in the specified direction, return "false"
   * Note that if distance of this point to the periodic boundary is
   * smaller than the search_radius, it is not necessary to get the image!
   *
   // 2D: i = 0, 1, 2 -> x, y, xy
   // 3D: i = 0, 1, 2 -> x, y, z;
   //     i = 3, 4, 5 -> xy, xz, yz
   //     i = 6       -> xyz
   *
   * this is used in particle_mesh.C to build particle/elem neighbor list.
   */
  bool get_image_point(const Point& pt0,
                       const Real& search_radius,
                       const std::size_t i,
                       Point& im_pt) const;
  
  /*
   * Correct the particle position 
   * (1) when a particle passes through the periodic boundary,
   *     it must re-enter the box through the opposite face.
   * (2) when a particle passes an inpenetrable wall by numerical accidents,
   *     we will have to pull it back in order to avoid losing particles.
   */
  void correct_position(Point& pt0) const;
  
  
  /*
   * Correct the particle-particle distance when computing the interaction 
   * force between a particle near the periodic boundary and its ghost images
   */
  Real point_distance(const Point& pt0,
                      const Point& pt1) const;
  
  
  /*
   * Correct the particle-particle distance vector when a particle is near 
   * the periodic boundary.
   * Return vector x =  pt1 - pt0
   */
  Point point_vector(const Point& pt0,
                     const Point& pt1) const;
  

  /*
   * space dimension of the periodic box
   */
  std::size_t dimension() const {  return _periodic_directions.size();   }
  
  
  
  /*
   * If an elem is "split" by the i-th periodic boundaries, it is an "image" element
   */
  std::vector<bool> image_elem(const Elem* elem) const;
  
  
  /*
   * If an elem is "split" by the periodic boundary, we need to construct an "image"
   * elem by moving some of its nodes to the other side of boundary
   */
  void build_image_elem(Elem* elem) const;
  
  
  /*
   * After building the image elem, the nodes are modified, so we
   * need to restore the original elem.
   *
   * To avoid elem change, this function has to be used together with the above.
   */
  void restore_image_elem(Elem* elem) const;
  
  
  
private:
  
  // the bounding box.
  BoundingBox _bbox;
  
  // the periodic direction: X, Y or Z - direction
  std::vector<bool> _periodic_directions;

  // the inlet direction: X, Y or Z - direction
  std::vector<bool> _inlet_directions;
  
  // inlet pressure
  std::vector<Real> _inlet_pressure;

  std::size_t _dim;
  
  Point _box_length;

}; // end of class
  

}  // end of namespace

