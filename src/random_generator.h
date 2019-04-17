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
#include <random>

#include "libmesh/reference_counted_object.h"
#include "libmesh/parallel_object.h"


/*
 * This class provide a random number generator
 * for stochastic modeling
 */

// using namespace libMesh;

using libMesh::Real;
using libMesh::ParallelObject;
using libMesh::ReferenceCountedObject;


class RandomGenerator : public ReferenceCountedObject<RandomGenerator> // ,
                                                                       // public
                                                                       // ParallelObject
{
public:

  // Constructor
  RandomGenerator();


  // Destructor
  ~RandomGenerator();


  /*
   * Generate a random vector with Gaussian(normal) distribution
   * mean (μ) with a specific standard deviation (σ)
   */
  std::vector<Real>random_vector_normal(const std::size_t n,
                                        const Real      & mean_val,
                                        const Real      & dev_val);


  /*
   * Generate a random vector with uniform distribution
   * in a range [a b]
   */
  std::vector<Real>random_vector_uniform(const std::size_t n,
                                         const Real      & a,
                                         const Real      & b);


  /*
   * reset the seed of the random generator.
   * Thus it will generate a new group of random numbers.
   */
  void set_seed(std::size_t seed_val);


  /*
   * random generator
   */
  std::default_random_engine& generator() {
    return _generator;
  }

private:

  /*
   * random engine should be put outside any function to
   * avoid generating the same group of random values every time.
   * std::default_random_engine generator(seed);
   */
  std::default_random_engine _generator;
}; // end of class
