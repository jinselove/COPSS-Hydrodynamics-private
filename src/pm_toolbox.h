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

// C++ Includes
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cstring>
#include <math.h>


// LibMesh library includes
// #include "libmesh/petsc_macro.h"
// #include "libmesh/libmesh_common.h"
#include "libmesh/libmesh.h" // define PI
#include "libmesh/point.h"
#include "libmesh/mesh.h"
#include "libmesh/serial_mesh.h"

// dense matrix help to debug
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/parallel_object.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;


/*
   this class defines basic tools used in our codes
 */
class PMToolBox {
public:

  // quadratic function used for applying BC to avoid singularities at corners
  static Real quadratic_function_2d(const Real& y,
                                    const Real& YA,
                                    const Real& YB);

  static Real quadratic_function_3d(const Real& y,
                                    const Real& z,
                                    const Real& YA,
                                    const Real& YB,
                                    const Real& ZA,
                                    const Real& ZB);


  static void output_dense_matrix(const DenseMatrix<Number>& Ke);

  static void output_dense_matrix(const DenseMatrix<Number>& Ke,
                                  const unsigned int         m,
                                  const unsigned int         n);

  static void output_dense_vector(const DenseVector<Number>& Fe);

  static void output_dense_vector(const DenseVector<Number>& Fe,
                                  const unsigned int         n);

  static void output_subdense_matrix(const DenseSubMatrix<Number>& Ke,
                                     const unsigned int            m,
                                     const unsigned int            n);

  static void output_subdense_vector(const DenseSubVector<Number>& Fe,
                                     const unsigned int            n);

  template<typename T>
  static void output_std_vector(const std::vector<T>& std_v);


  static void zero_filter_dense_matrix(DenseMatrix<Number>& Ae,
                                       const Real           tol);
  static void zero_filter_dense_vector(DenseVector<Number>& Ve,
                                       const Real           tol);


  /*
   * Check if a file exists or not.
   */
  static bool file_exist(const std::string& filename);


  /*
   * Rotation matrix for rotating angle = [ alpha, beta, theta ]
   * for x, y and z directions
   */
  static void coordinate_rotation(Point      & pt,
                                  const Point& angles);

  /*
   * Output a message on the screen
   */
  static void output_message(const std::string           & msg,
                             const Parallel::Communicator& comm_in);

  /*
   * Compute the min/max element size of a mesh.
   */
  static std::vector<Real>mesh_size(const MeshBase& _mesh);


  /*
   * Magnify mesh in the x/y/z directions according to the mag_factor
   */
  static void magnify_serial_mesh(SerialMesh & mesh,
                                  const Point& mag_factor);


  /*
   * Rotate mesh in the x/y/z axis according to the angles
   */
  static void rotate_serial_mesh(SerialMesh & mesh,
                                 const Point& angles);


  /*
   * shift mesh in the x/y/z dirction according to a given distance
   */
  static void shift_serial_mesh(SerialMesh             & mesh,
                                const std::vector<Real>& dist);
};
