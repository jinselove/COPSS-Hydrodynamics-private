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
#include <cmath>
#include <limits>

#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/point.h"
#include "point_mesh.h"

using namespace libMesh;


namespace libMesh {
/*
 * Type of a Dirac Delta function
 */
enum DeltaFunType {
  NOT_DEFINED_DELTAFUNTYPE = -1, // NOT DEFINED delta function
  SINGULAR                 = 0,  // Original Dirac delta function with
                                 // singularity
  SMOOTHED_EXP             = 1,  // Smoothed (Exponent/Gaussian form)
  SMOOTHED_POL1            = 2,  // (Polynomial type 1)
  SMOOTHED_POL2            = 3   // (Polynomial type 2)
};

class GGEMSystem : public ReferenceCountedObject<GGEMSystem> {
public:

  // Constructor
  GGEMSystem() {
    _kronecker_delta.resize(dim, std::vector<Real>(dim, 0.));

    for (int i = 0; i < 3; i++) _kronecker_delta[i][i] = 1.;
  }

  // Destructor
  ~GGEMSystem() {}


  /*
   * Kronecker delta function
   */
  inline Real kronecker_delta(const std::size_t i,
                              const std::size_t j) const
  {
    return i == j ? 1.0 : 0.0;
  }

  // ! set PointType
  void set_point_type(const PointType& _point_type) {
    point_type = _point_type;
  }

  PointType get_point_type() {
    return point_type;
  }

  // ! set GGEM alpha
  virtual void set_alpha(const Real& _alpha) = 0;
  Real         get_alpha() {
    return alpha;
  }

  // ! set regularization parameter
  virtual void set_ksi() = 0;
  Real         get_ksi() {
    return ksi;
  }

  // ! set bead radius (non-dimensional), it should be 1.
  void set_br0(const Real& _br0) {
    br0 = _br0;
  }

  Real get_br0() {
    return br0;
  }

protected:

  std::vector<std::vector<Real> >_kronecker_delta;

  // ! Pi = 3.1415926...
  const Real PI = libMesh::pi;

  // ! pi^(1/2)
  const Real sqrt_pi = std::sqrt(PI);

  // ! pi^(3/2)
  const Real pi_23 = std::pow(PI, 3. / 2.);

  // tolerance
  const Real r_eps = 1E-6;

  // dimensions
  const int dim = 3;

  // ! alpha
  Real alpha = NAN;

  // ! regularization parameter
  Real ksi = NAN;

  // ! point type: POLYMER_BEAD or LAGRANGIAN_POINT or POINT_PARTICLE
  PointType point_type = NOT_DEFINED_POINTTYPE;

  // ! bead radius
  Real br0 = NAN;
};
} // end of namespace
