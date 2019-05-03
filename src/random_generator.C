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


// #include <random>

#include "libmesh/libmesh_common.h"

#include "pm_toolbox.h"
#include "random_generator.h"


// ===========================================================================
RandomGenerator::RandomGenerator()
{
  // Do nothing
}

// ===========================================================================
RandomGenerator::~RandomGenerator()
{
  // Do nothing
}

// ===========================================================================
std::vector<Real>RandomGenerator::random_vector_normal(const std::size_t n,
                                                       const Real      & mean_val,
                                                       const Real      & dev_val)
{
  START_LOG("random_vector_normal()", "RandomGenerator");

  std::vector<Real> rv(n);

  // random engine
  // std::default_random_engine generator;
  std::normal_distribution<Real> distribution(mean_val, dev_val);
  distribution.reset();

  // std::cout << "--->test in RandomGenerator: Current random seed is "
  //          << this->generator() << std::endl;

  for (std::size_t i = 0; i < n; ++i)
  {
    const Real number = distribution(_generator);
    rv[i] = number;
  }

  STOP_LOG("random_vector_normal()", "RandomGenerator");
  return rv;
}

// ===========================================================================
std::vector<Real>RandomGenerator::random_vector_uniform(const std::size_t n,
                                                        const Real      & a,
                                                        const Real      & b)
{
  START_LOG("random_vector_uniform()", "RandomGenerator");

  std::vector<Real> rv(n);

  // random engine
  std::uniform_real_distribution<Real> distribution(a, b);
  distribution.reset();

  // std::cout << "--->test in RandomGenerator: Current random seed is "
  //          << this->generator() << std::endl;

  for (std::size_t i = 0; i < n; ++i)
  {
    const Real number = distribution(_generator);
    rv[i] = number;
  }

  STOP_LOG("random_vector_uniform()", "RandomGenerator");
  return rv;
}

// ===========================================================================
void RandomGenerator::set_seed(std::size_t seed_val)
{
  _generator.seed(seed_val);
}
