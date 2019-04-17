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


#include <sys/stat.h>
#include <unistd.h>


#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "pm_toolbox.h"


// =============================================================================================
Real PMToolBox::quadratic_function_2d(const Real& y,
                                      const Real& YA,
                                      const Real& YB)
{
  Real Y0 = (YA + YB) / 2.0;
  Real DY = YB - Y0; // note, this also equals to -(YA - Y0)
  Real V0 = DY / 2.; // maximum velocity magnitude
  Real A0 = V0 / (DY * DY);

  return V0 - A0 * (y - Y0) * (y - Y0);

  // a controls the magnitude of velocity, b is the mag of geometry(height of
  // channal)
  // to make sure u=0 at boundary, it requires b^2 = a*y^2
  //  Real a = 1.0, b = 0.5, y0 = 0.0;
  //  return b*b - a*(yb - y0)*(yb - y0);
}

// =============================================================================================
Real PMToolBox::quadratic_function_3d(const Real& y,
                                      const Real& z,
                                      const Real& YA,
                                      const Real& YB,
                                      const Real& ZA,
                                      const Real& ZB)
{
  Real Y0 = (YA + YB) / 2.0,  Z0 = (ZA + ZB) / 2.0;
  Real DY = YB - Y0,        DZ = ZB - Z0; // note, this also equals to -(YA -
                                          // Y0)
  Real VY = DY / 2.,          VZ = VY;    // maximum velocity magnitudes
  Real AY = VY / (DY * DY),     AZ = VZ / (DZ * DZ);

  Real value = (VY - AY * (y - Y0) * (y - Y0)) * (VZ - AZ * (z - Z0) * (z - Z0));

  // std::cout<<"quadratic_function_3d test: value = "<<value<<std::endl;
  return value;


  // a controls the magnitude of velocity, b is the mag of geometry(height of
  // channal)
  // make sure that ( b^2-a*y^2 )*( b^2 - z^2 )
  //  Real a = 1.0, b = 0.5;
  //  Real y0 = 0.0,z0 = 0.0;
  //  return ( b*b - a*(yb - y0)*(yb - y0) )*( b*b - a*(zb - z0)*(zb - z0) );
}

// =============================================================================================
void PMToolBox::output_dense_matrix(const DenseMatrix<Number>& Ke)
{
  output_dense_matrix(Ke, Ke.m(), Ke.n());
}

// =============================================================================================
void PMToolBox::output_dense_matrix(const DenseMatrix<Number>& Ke,
                                    const unsigned int         m,
                                    const unsigned int         n)
{
  std::cout << "--------------------------- output matrix " << m << " x " << n
            << " ---------------------------" << std::endl;

  for (unsigned int i = 0; i < m; ++i)
  {
    for (unsigned int j = 0; j < n; ++j)
    {
      std::cout << Ke(i, j) << "  " << std::setw(5);
    }
    std::cout << std::endl;
  }
  std::cout <<
  "--------------------------- end of matrix ---------------------------" <<
  std::endl;
}

// =============================================================================================
void PMToolBox::output_dense_vector(const DenseVector<Number>& Fe)
{
  std::cout <<
  "--------------------------- output vector ---------------------------" <<
  std::endl;

  for (unsigned int j = 0; j < Fe.size();
       ++j) std::cout << Fe(j) << "  " << std::setw(5);
  std::cout << std::endl;
  std::cout <<
  "--------------------------- end of vector ---------------------------" <<
  std::endl;
}

// =============================================================================================
void PMToolBox::output_dense_vector(const DenseVector<Number>& Fe,
                                    const unsigned int         n)
{
  std::cout <<
  "--------------------------- output vector ---------------------------" <<
  std::endl;

  for (unsigned int j = 0; j < n; ++j) std::cout << Fe(j) << "  " << std::setw(5);
  std::cout << std::endl;
  std::cout <<
  "--------------------------- end of vector ---------------------------" <<
  std::endl;
}

// =============================================================================================
void PMToolBox::output_subdense_matrix(const DenseSubMatrix<Number>& Ke,
                                       const unsigned int            m,
                                       const unsigned int            n)
{
  std::cout << "--------------------------- output matrix " << m << " x " << n
            << " ---------------------------" << std::endl;

  for (unsigned int i = 0; i < m; ++i)
  {
    for (unsigned int j = 0; j < n; ++j)
    {
      std::cout << Ke(i, j) << "  " << std::setw(5);
    }
    std::cout << std::endl;
  }
  std::cout <<
  "--------------------------- end of matrix ---------------------------" <<
  std::endl;
}

// =============================================================================================
void PMToolBox::output_subdense_vector(const DenseSubVector<Number>& Fe,
                                       const unsigned int            n)
{
  std::cout <<
  "--------------------------- output vector ---------------------------" <<
  std::endl;

  for (unsigned int j = 0; j < n; ++j) std::cout << Fe(j) << "  " << std::setw(5);

  std::cout <<
  "--------------------------- end of vector ---------------------------" <<
  std::endl;
}

// =============================================================================================
template<typename T>
void PMToolBox::output_std_vector(const std::vector<T>& std_v)
{
  std::cout <<
  "--------------------------- output vector ---------------------------" <<
  std::endl;

  for (unsigned int j = 0; j < std_v.size(); ++j)
  {
    std::cout << std::setw(5) << std_v[j] << "  ";
  }
  std::cout << std::endl;
  std::cout <<
  "--------------------------- end of vector ---------------------------" <<
  std::endl;
  std::cout << std::endl;
}

// =============================================================================================
void PMToolBox::zero_filter_dense_matrix(DenseMatrix<Number>& Ke, const Real tol)
{
  for (unsigned int i = 0; i < Ke.m(); ++i)
  {
    for (unsigned int j = 0; j < Ke.n(); ++j)
    {
      if (std::abs(Ke(i, j)) <= tol) Ke(i, j) = 0.0;
    } // end for j-loop
  }   // end for i-loop
}

// =============================================================================================
void PMToolBox::zero_filter_dense_vector(DenseVector<Number>& Ve, const Real tol)
{
  for (unsigned int j = 0; j < Ve.size(); ++j)
  {
    if (std::abs(Ve(j)) <= tol) Ve(j) = 0.0;
  } // end for j-loop
}

// =============================================================================================
bool PMToolBox::file_exist(const std::string& filename)
{
  struct stat buffer;

  return stat(filename.c_str(), &buffer) == 0;
}

// =============================================================================================
void PMToolBox::coordinate_rotation(Point      & pt,
                                    const Point& angles)
{
  START_LOG("coordinate_rotation()", "PMToolBox");

  // Init the rotation matrix
  DenseMatrix<Number> Rx(3, 3), Ry(3, 3), Rz(3, 3);
  DenseVector<Number> Vin(3), Vout(3);
  Real s_a = 0., c_a = 0.; // store sin & cos values
  std::vector<Real> theta(3);

  for (unsigned int i = 0; i < 3; ++i) {
    theta[i] = angles(i) / 180. * libMesh::pi;
  }

  // Rx
  s_a      = std::sin(theta[0]); c_a = std::cos(theta[0]); // sin/cos
  Rx(0, 0) = 1.0;    Rx(0, 1) = 0.0;    Rx(0, 2) = 0.0;
  Rx(1, 0) = 0.0;    Rx(1, 1) = c_a;    Rx(1, 2) = s_a;
  Rx(2, 0) = 0.0;    Rx(2, 1) = -s_a;    Rx(2, 2) = c_a;

  // Ry
  s_a      = std::sin(theta[1]); c_a = std::cos(theta[1]); // sin/cos
  Ry(0, 0) = c_a;    Ry(0, 1) = 0.0;    Ry(0, 2) = -s_a;
  Ry(1, 0) = 0.0;    Ry(1, 1) = 1.0;    Ry(1, 2) = 0.0;
  Ry(2, 0) = s_a;    Ry(2, 1) = 0.0;    Ry(2, 2) = c_a;

  // Rz
  s_a      = std::sin(theta[2]); c_a = std::cos(theta[2]); // sin/cos
  Rz(0, 0) = c_a;    Rz(0, 1) = s_a;    Rz(0, 2) = 0.0;
  Rz(1, 0) = -s_a;    Rz(1, 1) = c_a;    Rz(1, 2) = 0.0;
  Rz(2, 0) = 0.0;    Rz(2, 1) = 0.0;    Rz(2, 2) = 1.0;

  // Rx = Rx*Ry*Rz
  Rx.right_multiply(Ry);
  Rx.right_multiply(Rz);

  // Rx.print();

  // Vx = pt
  for (unsigned int i = 0; i < 3; ++i) {
    Vin(i) = pt(i);
  }

  // Vout = Rx*Vin
  Rx.vector_mult(Vout, Vin);

  // pt = Vin
  for (unsigned int i = 0; i < 3; ++i) {
    pt(i) = Vout(i);
  }

  STOP_LOG("coordinate_rotation()", "PMToolBox");
}

// =============================================================================================
void PMToolBox::output_message(const std::string           & msg,
                               const Parallel::Communicator& comm_in)
{
  // comm_in.barrier();
  if (comm_in.rank() == 0)
  {
    printf("\n");
    printf(
      "---------------------------------------------------------------------\n");
    printf("%s \n", msg.c_str());
    printf(
      "---------------------------------------------------------------------\n");
    printf("\n");
  }
}

// =============================================================================================
std::vector<Real>PMToolBox::mesh_size(const MeshBase& _mesh)
{
  START_LOG("mesh_size()", "PMToolBox");

  Real hmin = 0.0, hmax = 0.0;
  std::size_t count = 0;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Loop over all the elements in the mesh that live
     on the local processor, and compute the element size.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  MeshBase::const_element_iterator el           = _mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_elements_end();

  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently working on.
    const Elem *elem = *el;
    Real hmin_t      = elem->hmin();
    Real hmax_t      = elem->hmax();

    if (count == 0)
    {
      hmin = hmin_t;
      hmax = hmax_t;
    }
    else
    {
      hmin = std::min(hmin, hmin_t);
      hmax = std::min(hmax, hmax_t);
    }
    count++;
  } // end for el-loop

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the min/max values on all the processors.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  _mesh.comm().min(hmin);
  _mesh.comm().max(hmax);
  std::vector<Real> min_max(2);
  min_max[0] = hmin;
  min_max[1] = hmax;

  STOP_LOG("mesh_size()", "PMToolBox");
  return min_max;
}

// =============================================================================================
void PMToolBox::magnify_serial_mesh(SerialMesh & mesh,
                                    const Point& mag_factor)
{
  START_LOG("magnify_serial_mesh()", "PMToolBox");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     1. magnify
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  const unsigned int dim               = mesh.spatial_dimension();
  MeshBase::node_iterator nd           = mesh.active_nodes_begin();
  const MeshBase::node_iterator end_nd = mesh.active_nodes_end();

  for (; nd != end_nd; ++nd)
  {
    // Store a pointer to the element we are currently working on.
    Node *node = *nd;

    // magnify x&y coordinates by R0, and z by H0
    for (unsigned int i = 0; i < dim; ++i) {
      (*node)(i) =  (*node)(i) * mag_factor(i);
    }
  } // end for nd-loop

  STOP_LOG("magnify_serial_mesh()", "PMToolBox");
}

// =============================================================================================
void PMToolBox::rotate_serial_mesh(SerialMesh & mesh,
                                   const Point& angles)
{
  START_LOG("rotate_serial_mesh()", "PMToolBox");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Rotate according to the angles
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  MeshBase::node_iterator nd           = mesh.active_nodes_begin();
  const MeshBase::node_iterator end_nd = mesh.active_nodes_end();

  for (; nd != end_nd; ++nd)
  {
    // Store a pointer to the element we are currently working on.
    Node *node = *nd;
    coordinate_rotation(*node, angles);
  }

  STOP_LOG("rotate_serial_mesh()", "PMToolBox");
}

// =============================================================================================
void PMToolBox::shift_serial_mesh(SerialMesh             & mesh,
                                  const std::vector<Real>& dist)
{
  START_LOG("shift_serial_mesh()", "PMToolBox");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     shift the mesh by a distance according to the angles
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  const unsigned int dim               = mesh.spatial_dimension();
  MeshBase::node_iterator nd           = mesh.active_nodes_begin();
  const MeshBase::node_iterator end_nd = mesh.active_nodes_end();

  for (; nd != end_nd; ++nd)
  {
    // Store a pointer to the element we are currently working on.
    Node *node = *nd;

    // magnify x&y coordinates by R0, and z by H0
    for (unsigned int i = 0; i < dim; ++i) {
      (*node)(i) +=  dist[i];
    }
  }

  STOP_LOG("shift_serial_mesh()", "PMToolBox");
}
