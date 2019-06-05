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


#include "libmesh/libmesh_logging.h"


#include "pm_periodic_boundary.h"


using namespace libMesh;


// ======================================================================
PMPeriodicBoundary::PMPeriodicBoundary(const std::pair<Point, Point>& bbox_pts,
                                       const std::vector<bool>& periodic_directions)
  : _bbox(bbox_pts),
  _periodic_directions(periodic_directions)
{
  _box_length = _bbox.max() - _bbox.min();
  _dim        = _periodic_directions.size();
}

// ======================================================================
PMPeriodicBoundary::PMPeriodicBoundary(const BoundingBox      & bbox,
                                       const std::vector<bool>& periodic_directions)
  : _bbox(bbox),
  _periodic_directions(periodic_directions)
{
  _box_length = _bbox.max() - _bbox.min();
  _dim        = _periodic_directions.size();
}

// ======================================================================
PMPeriodicBoundary::PMPeriodicBoundary(const Point            & bbox_pmin,
                                       const Point            & bbox_pmax,
                                       const std::vector<bool>& periodic_directions)
  : _bbox(bbox_pmin, bbox_pmax),
  _periodic_directions(periodic_directions)
{
  _box_length = _bbox.max() - _bbox.min();
  _dim        = _periodic_directions.size();
}

// ======================================================================
PMPeriodicBoundary::PMPeriodicBoundary(const Point            & bbox_pmin,
                                       const Point            & bbox_pmax,
                                       const std::vector<bool>& periodic_directions,
                                       const std::vector<bool>& inlet_directions,
                                       const std::vector<Real>& inlet_pressure)
  :  _bbox(bbox_pmin, bbox_pmax),
  _periodic_directions(periodic_directions),
  _inlet_directions(inlet_directions),
  _inlet_pressure(inlet_pressure)
{
  _box_length = _bbox.max() - _bbox.min();
  _dim        = _periodic_directions.size();
}

// ======================================================================
PMPeriodicBoundary::PMPeriodicBoundary(const PMPeriodicBoundary& pb)
{
  _bbox                = pb._bbox;
  _periodic_directions = pb._periodic_directions;

  _box_length = _bbox.max() - _bbox.min();
  _dim        = _periodic_directions.size();
}

// ======================================================================
PMPeriodicBoundary::~PMPeriodicBoundary()
{
  // do nothing
}

// ======================================================================
bool PMPeriodicBoundary::periodic_direction(const std::size_t i) const
{
  return _periodic_directions[i];
}

// ======================================================================
const std::vector<bool>& PMPeriodicBoundary::periodic_direction() const
{
  return _periodic_directions;
}

// ======================================================================
bool PMPeriodicBoundary::inlet_direction(const std::size_t i) const
{
  return _inlet_directions[i];
}

// ======================================================================
const std::vector<bool>& PMPeriodicBoundary::inlet_direction() const
{
  return _inlet_directions;
}

// ======================================================================
const Real PMPeriodicBoundary::inlet_pressure(std::size_t i) const
{
  return _inlet_pressure[i];
}

// ======================================================================
const std::vector<Real>& PMPeriodicBoundary::inlet_pressure() const
{
  return _inlet_pressure;
}

// ======================================================================
Real PMPeriodicBoundary::box_length(const std::size_t i) const
{
  return _box_length(i);
}

// ======================================================================
Point PMPeriodicBoundary::box_length() const
{
  return _box_length;
}

// ======================================================================
bool PMPeriodicBoundary::get_image_point(const Point     & pt0,
                                         const std::size_t i, // i=0,1,2 for
                                                              // x,y,z dir
                                         Point           & im_pt) const
{
  bool has_image = false;

  im_pt = pt0;

  // first need to check if the point is in the bounding domain
  const bool in_box = _bbox.contains_point(pt0);

  // If the pt is in the box, and there is also a periodic boundary
  // in the i-th direction, then the image will be computed!
  if (in_box && this->periodic_direction(i))
  {
    const Point box_min = _bbox.min();
    const Point box_max = _bbox.max();
    const Real  xsize   =  this->box_length(i);

    // the image in the i-th direction
    if (pt0(i) - box_min(i) < xsize / 2.0)
    { im_pt(i) += xsize;  has_image = true; }

    if (box_max(i) - pt0(i) <= xsize / 2.0)
    { im_pt(i) -= xsize;  has_image = true; }
  } // end if (in_box)

  return has_image;
}

// ======================================================================
bool PMPeriodicBoundary::get_image_point(const Point     & pt0,
                                         const Real      & search_radius,
                                         const std::size_t i,
                                         Point           & im_pt) const
{
  START_LOG("get_image_point()", "PMPeriodicBoundary");

  bool has_image = false;
  libmesh_example_requires(2 <= _dim, "2D/3D support for PMPeriodicBoundary");

  // first need to check if the point is in the bounding domain
  // If this is true, and there is also periodic boundary in the i-th direction,
  // then the image will be computed!
  // ---------------------------------------------------------------
  // 2D: i = 0, 1, 2 -> x, y, xy
  // 3D: i = 0, 1, 2 -> x, y, z;
  //     i = 3, 4, 5 -> xy, xz, yz
  //     i = 6       -> xyz
  // ---------------------------------------------------------------
  // Note the condition in if and else if can NOT be satisfied simultaneously.
  // Otherwise, the domain size used in this periodic direction is too small.
  const bool in_box = _bbox.contains_point(pt0);

  if (in_box)
  {
    // get the box geometric information
    im_pt = pt0;
    const Point& box_min  = this->box_min();
    const Point& box_max  = this->box_max();
    const Point& box_size = this->box_length();

    // first, check the image along the i-th direction.
    if ((i < _dim) && _periodic_directions[i])
    {
      if (pt0(i) - box_min(i) <= search_radius)
      { im_pt(i) += box_size(i);  has_image = true; }
      else if (box_max(i) - pt0(i) < search_radius)
      { im_pt(i) -= box_size(i);  has_image = true; }
    }

    // // check the image along the i-j corner direction
    // else if( i==_dim && _periodic_directions[0] && _periodic_directions[1] )
    // // xy direction
    // {
    //   // check x-dir
    //   bool xclose = false, yclose = false;
    //   if (pt0(0)-box_min(0) <= search_radius)
    //   { im_pt(0) += box_size(0);  xclose = true; }
    //   else if (box_max(0)-pt0(0) < search_radius)
    //   { im_pt(0) -= box_size(0);  xclose = true; }

    //   // check y-dir
    //   if (pt0(1)-box_min(1) <= search_radius)
    //   { im_pt(1) += box_size(1);  yclose = true; }
    //   else if (box_max(1)-pt0(1) < search_radius)
    //   { im_pt(1) -= box_size(1);  yclose = true; }

    //   if ( xclose && yclose ) has_image = true;
    // }
    // else if( i==_dim+1 && _periodic_directions[0] && _periodic_directions[2]
    // ) // xz direction
    // {
    //   // check x-dir
    //   bool xclose = false, zclose = false;
    //   if (pt0(0)-box_min(0) <= search_radius)
    //   { im_pt(0) += box_size(0);  xclose = true; }
    //   else if (box_max(0)-pt0(0) < search_radius)
    //   { im_pt(0) -= box_size(0);  xclose = true; }

    //   // check z-dir
    //   if (pt0(2)-box_min(2) <= search_radius)
    //   { im_pt(2) += box_size(2);  zclose = true; }
    //   else if (box_max(2)-pt0(2) < search_radius)
    //   { im_pt(2) -= box_size(2);  zclose = true; }

    //   if ( xclose && zclose ) has_image = true;
    // }
    // else if( i==_dim+2 && _periodic_directions[2] && _periodic_directions[1]
    // ) // yz direction
    // {
    //   // check z-dir
    //   bool zclose = false, yclose = false;
    //   if (pt0(2)-box_min(2) <= search_radius)
    //   { im_pt(2) += box_size(2);  zclose = true; }
    //   else if (box_max(2)-pt0(2) < search_radius)
    //   { im_pt(2) -= box_size(2);  zclose = true; }

    //   // check y-dir
    //   if (pt0(1)-box_min(1) <= search_radius)
    //   { im_pt(1) += box_size(1);  yclose = true; }
    //   else if (box_max(1)-pt0(1) < search_radius)
    //   { im_pt(1) -= box_size(1);  yclose = true; }

    //   if ( zclose && yclose ) has_image = true;
    // }
    // else if( i==_dim+3 &&                                        // xyz
    // direction
    //         _periodic_directions[0] && _periodic_directions[1] &&
    // _periodic_directions[2] )
    // {
    //   // check x-dir
    //   bool xclose = false, yclose = false, zclose = false;
    //   if (pt0(0)-box_min(0) <= search_radius)
    //   { im_pt(0) += box_size(0);  xclose = true; }
    //   else if (box_max(0)-pt0(0) < search_radius)
    //   { im_pt(0) -= box_size(0);  xclose = true; }

    //   // check y-dir
    //   if (pt0(1)-box_min(1) <= search_radius)
    //   { im_pt(1) += box_size(1);  yclose = true; }
    //   else if (box_max(1)-pt0(1) < search_radius)
    //   { im_pt(1) -= box_size(1);  yclose = true; }

    //   // check z-dir
    //   if (pt0(2)-box_min(2) <= search_radius)
    //   { im_pt(2) += box_size(2);  zclose = true; }
    //   else if (box_max(2)-pt0(2) < search_radius)
    //   { im_pt(2) -= box_size(2);  zclose = true; }

    //   if ( xclose && yclose && zclose) has_image = true;
    // }
    // // end if-else
  } // end if (inbox)


  STOP_LOG("get_image_point()", "PMPeriodicBoundary");

  return has_image;
}

// ======================================================================
void PMPeriodicBoundary::correct_position(Point& pt0) const
{
  START_LOG("correct_position()", "PMPeriodicBoundary");

  for (std::size_t i = 0; i < _periodic_directions.size(); ++i)
  {
    if (_periodic_directions[i]) // periodic boundary
    {
      if (pt0(i) <  _bbox.min()(i)) pt0(i) += this->box_length()(i);

      if (pt0(i) >= _bbox.max()(i)) pt0(i) -= this->box_length()(i);
    }

    // Only correct the position due to the periodicity, not the wall,
    // because the wall may have complex geometry different from the bbox.
    //    else        // non-periodic boundary (inpenetrable wall)
    //    {
    //      if( pt0(i) <  _bbox.min()(i) ) pt0(i) = 2.*_bbox.min()(i) - pt0(i);
    //      if( pt0(i) >  _bbox.max()(i) ) pt0(i) = 2.*_bbox.max()(i) - pt0(i);
    //    }
  } // end for i-loop

  STOP_LOG("correct_position()", "PMPeriodicBoundary");
}

// ======================================================================
Real PMPeriodicBoundary::point_distance(const Point& pt0,
                                        const Point& pt1) const
{
  START_LOG("point_distance()", "PMPeriodicBoundary");

  // the distance vector
  const Point x = this->point_vector(pt0, pt1);

  //  Real dist = 0.0;
  //  for(std::size_t i=0; i<3; ++i) {
  //    dist += x(i)*x(i);
  //  }

  STOP_LOG("point_distance()", "PMPeriodicBoundary");

  return x.norm();
}

// ======================================================================
Point PMPeriodicBoundary::point_vector(const Point& pt0,
                                       const Point& pt1) const
{
  START_LOG("point_vector()", "PMPeriodicBoundary");

  // the original distance vector
  Point dpt = pt1 - pt0;

  // Modify the value due to periodic boundaries.
  for (std::size_t i = 0; i < _dim; ++i)
  {
    if (_periodic_directions[i] or _inlet_directions[i])
    {
      // modify the value when |dpt(i)| > xsize
      const Real xsize = this->box_length()(i);

      if (dpt(i) <= -xsize / 2.0) dpt(i) += xsize;

      if (dpt(i) > xsize / 2.0) dpt(i) -= xsize;

      // printf("******************** I am here!
      // point_vector()*****************\n");
    }
  }

  STOP_LOG("point_vector()", "PMPeriodicBoundary");

  return dpt;
}

// ======================================================================
std::vector<bool>PMPeriodicBoundary::image_elem(const Elem *elem) const
{
  START_LOG("image_elem()", "PMPeriodicBoundary");

  // Get elem maximum size and the domain box size
  const Real   hmax     = elem->hmax();
  const Point& box_size = this->box_length();


  // Loop over each direction
  std::vector<bool> direction_flag(_dim, false);

  for (std::size_t i = 0; i < _dim; ++i)
  {
    if (_periodic_directions[i] && (hmax > box_size(i) / 2.)) {
      direction_flag[i] = true;
    }
  }


  STOP_LOG("image_elem()", "PMPeriodicBoundary");
  return direction_flag;
}

// ======================================================================
void PMPeriodicBoundary::build_image_elem(Elem *elem) const
{
  START_LOG("build_image_elem()", "PMPeriodicBoundary");

  // First, we need to determine if this elem is split in the p-direction
  const std::vector<bool> direction_flag = this->image_elem(elem);
  const Point& box_size                  = this->box_length();
  const Point& box_min                   = this->box_min();
  const Point& box_max                   = this->box_max();

  // Loop over each direction
  for (std::size_t i = 0; i < direction_flag.size(); ++i)
  {
    if (direction_flag[i])
    {
      // Loop over each node
      for (std::size_t j = 0; j < elem->n_nodes(); ++j)
      {
        const Point& pti = elem->point(j); // point ref of the i-th node

        // if this pt is close to side-A, move it to side-B
        if (box_max(j) - pti(j) > box_size(i) / 2.) {
          elem->point(j)(i) += box_size(i);
        }
      } // end for j-loop
    }   // end if
  }     // end for i-loop


  STOP_LOG("build_image_elem()", "PMPeriodicBoundary");
}

// ======================================================================
void PMPeriodicBoundary::restore_image_elem(Elem *elem) const
{
  START_LOG("restore_image_elem()", "PMPeriodicBoundary");

  // Loop over each node
  for (std::size_t j = 0; j < elem->n_nodes(); ++j)
  {
    Point& pti = elem->point(j); // point ref of the i-th node
    this->correct_position(pti);
  }

  STOP_LOG("restore_image_elem()", "PMPeriodicBoundary");
}
