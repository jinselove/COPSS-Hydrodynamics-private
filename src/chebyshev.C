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


#include <math.h>

#include "libmesh/libmesh_config.h"
#include "chebyshev.h"


using namespace libMesh;


// ======================================================================================
Chebyshev::Chebyshev()
{
  // do nothing
}

// ======================================================================================
Chebyshev::~Chebyshev()
{
  // do nothing
}

// ======================================================================================
Real Chebyshev::chebyshev_function(const Real x,
                                   const Real da_cheb,
                                   const Real db_cheb) const
{
  START_LOG("chebyshev_function()", "Chebyshev");

  Real f = std::sqrt((x - db_cheb) / da_cheb);

  STOP_LOG("chebyshev_function()", "Chebyshev");
  return 1.0 / f;
}

// ======================================================================================
DenseVector<Number>Chebyshev::chebyshev_coefficients(const std::size_t N,
                                                     const Real        da_cheb,
                                                     const Real        db_cheb,
                                                     const std::string method)
const
{
  START_LOG("chebyshev_coefficients()", "Chebyshev");

  // Chebyshev nodes and weights
  std::vector<Real> x, w;

  if (method == "Gauss_Lobatto") this->chebyshev_gauss_lobatto(N, x, w);
  else if (method == "Gauss_Radau") this->chebyshev_gauss_radau(N, x, w);
  else this->chebyshev_gauss(N, x, w);

  // end if


  DenseVector<Number> coef(N + 1);
  const Real RN =  Real(N);
  const Real PI = libMesh::pi;

  for (std::size_t k = 0; k <= N; ++k)
  {
    coef(k) = 0.0;

    for (std::size_t j = 0; j <= N; ++j)
    {
      Real tmpt = 0.0;

      if (method == "Gauss_Lobatto") tmpt = PI * Real(k * j) / RN;
      else if (method ==
               "Gauss_Radau") tmpt = 2. * PI * Real(k * j) / (2. * RN + 1.);
      else tmpt = PI * Real(2. * j + 1.) * Real(k) / (2. * RN + 2.);

      // end if

      const Real fj = this->chebyshev_function(x[j], da_cheb, db_cheb);
      coef(k) += w[j] * std::cos(tmpt) * fj;
    } // end for j-loop

    // multiply by a factor
    Real ck = 1.0;

    if (k == 0) ck = 2.0;
    coef(k) *= 2. / (PI * ck);
  } // end for k-loop

  STOP_LOG("chebyshev_coefficients()", "Chebyshev");
  return coef;
}

// ======================================================================================
DenseMatrix<Number>Chebyshev::transform_matrix(const std::size_t N) const
{
  const Real RN = Real(N);
  DenseMatrix<Number> Tmat(N + 1, N + 1);

  for (std::size_t i = 0; i <= N; ++i)
  {
    Real ci = 1.;

    if ((i == 0) || (i == N)) ci = 2.;

    for (std::size_t j = 0; j <= N; ++j)
    {
      Real cj = 1.;

      if ((j == 0) || (j == N)) cj = 2.;
      Tmat(i, j) = 2. / (ci * cj * RN) * std::cos(libMesh::pi * Real(i) * Real(
                                                    j) / RN);

      if (std::abs(Tmat(i, j)) < 1E-10) Tmat(i, j) = 0.0;  // filter
    }
  }
  return Tmat;
}

// ======================================================================================
void Chebyshev::chebyshev_gauss(const std::size_t  N,
                                std::vector<Real>& x,
                                std::vector<Real>& w) const
{
  // init the vector x and w
  x.resize(N + 1);  w.resize(N + 1);
  const Real PI = libMesh::pi;

  // compute the quadrature points and weights
  for (std::size_t i = 0; i <= N; ++i)
  {
    const Real j = (Real)i;
    x[i] = std::cos((2. * j + 1.) * PI / (2. * Real(N) + 2.));
    w[i] = PI / (Real(N) + 1.);
  }
}

// ======================================================================================
void Chebyshev::chebyshev_gauss_radau(const std::size_t  N,
                                      std::vector<Real>& x,
                                      std::vector<Real>& w) const
{
  // init the vector x and w
  x.resize(N + 1);  w.resize(N + 1);
  const Real PI = libMesh::pi;

  // compute the quadrature points and weights
  for (std::size_t i = 0; i <= N; ++i)
  {
    const Real j = (Real)i;
    x[i] = std::cos(2. * PI * j / (2. * Real(N) + 1.));
    w[i] = 2. * PI / (2. * Real(N) + 1.);
  }
  w[0] = PI / (2. * Real(N) + 1.);
}

// ======================================================================================
void Chebyshev::chebyshev_gauss_lobatto(const std::size_t  N,
                                        std::vector<Real>& x,
                                        std::vector<Real>& w) const
{
  // init the vector x and w
  x.resize(N + 1);  w.resize(N + 1);
  const Real PI = libMesh::pi;

  // compute the quadrature points and weights
  for (std::size_t i = 0; i <= N; ++i)
  {
    const Real j = (Real)i;
    x[i] = std::cos(PI * j / Real(N));
    w[i] = PI / Real(N);
  }
  w[0] = PI / (2. * Real(N));
  w[N] = PI / (2. * Real(N));
}

// ======================================================================================
// ======================================================================================
