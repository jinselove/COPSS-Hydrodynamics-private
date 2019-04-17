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
#include <vector>


#include "libmesh/libmesh.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"


namespace libMesh {
/*
 * The Chebyshev is designed for implementing
 * Chebyshev polynomial series
 */

// using libMesh::Real;


class Chebyshev {
public:

  // Constructor
  Chebyshev();


  // Destructor
  ~Chebyshev();


  /*
   * The function of the Chebyshev expansion
   */
  Real chebyshev_function(const Real x,
                          const Real da_cheb,
                          const Real db_cheb) const;


  /*
   * Chebyshev expansion coefficients
   */
  DenseVector<Number>chebyshev_coefficients(const std::size_t N,
                                            const Real        da_cheb,
                                            const Real        db_cheb,
                                            const std::string method) const;


  /*
   * The transforming matrix from physical space to Chebyshev tranform space
   * using Gauss-Lotatto method
   */
  DenseMatrix<Number>transform_matrix(const std::size_t N) const;


  /*
   * Quadrature points and weights to evaluate the Chebyshev coeffecients
   * (1) Chebyshev-Gauss: n = 0, 1, ... , N
   */
  void chebyshev_gauss(const std::size_t  N,        // # of expansion terms
                       std::vector<Real>& x,        // quadrature points
                       std::vector<Real>& w) const; // weights


  /*
   * Quadrature points and weights to evaluate the Chebyshev coeffecients
   * (2) Chebyshev-Gauss-Radau: n = 0, 1, ... , N
   */
  void chebyshev_gauss_radau(const std::size_t  N,        // # of expansion
                                                          // terms
                             std::vector<Real>& x,        // quadrature points
                             std::vector<Real>& w) const; // weights


  /*
   * Quadrature points and weights to evaluate the Chebyshev coeffecients
   * (3) Chebyshev-Gauss-Lobatto: n = 0, 1, ... , N
   */
  void chebyshev_gauss_lobatto(const std::size_t  N,        // # of expansion
                                                            // terms
                               std::vector<Real>& x,        // quadrature points
                               std::vector<Real>& w) const; // weights
};                                                          // end class
                                                            // Chebyshev
} // end namespace
