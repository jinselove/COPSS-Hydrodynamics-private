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

// Local Includes
#include "assemble_poisson.h"
#include "solver_poisson.h"
#include "pm_linear_implicit_system.h"
#include "ggem_poisson.h"

namespace libMesh
{

/*
 * The PMSystemPoisson is designed to solve the Poisson
 * equation with point charge sources on a finite element mesh.
 */

class PMSystemPoisson : public PMLinearImplicitSystem
{
public:

  /**
   * Constructor.
   */
  PMSystemPoisson (EquationSystems& es,
                   const std::string& name,
                   const unsigned int number); // number of systems


  /**
   * Destructor.
   */
  virtual ~PMSystemPoisson ();


  /**
   * The type of system.
   */
  typedef PMSystemPoisson sys_type;


  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }


  /**
   * Clear all the data structures associated with the system.
   */
  void clear ();


  /**
   * Assemble the system matrix.
   */
  void assemble_matrix (const std::string& system_name);


  /**
   * Assemble the system rhs.
   */
  void assemble_rhs (const std::string& system_name);


  /*
   * Solve the system.
   * re_init = true => re-assemble the matrix and reinit the KSP solver.
   */
  void solve (const bool& re_init);


  /*
   * Add the local solution to the global solution for electrical potential
   */
  void add_local_solution();


  /*
   * Compute electrical potential at all beads' locations
   */
  void compute_point_potential(std::vector<Real>& pv);


  /*
   * Compute electric field (Ex, Ey, Ez) at all beads' locations
   */
  void compute_point_efield(std::vector<Real>& pv);


  /*
   * Add electrostatic forces for all beads
   */
  void add_electrostatic_forces();


  /*
   * Return the SolverPoisson
   */
  SolverPoisson& solver_poisson() { return _solver_poisson; }


  /**
   * Local electrical potential of a point in an unbounded space,
   * which is computed from Green's function. This function directly
   * search the neighbor list around Point &p.
   * charge_type: "regularized"
   */
  Real local_potential_field(const Point &p,
                             const std::string& charge_type) const;


  /**
   * Local electrical potential of a point in an unbounded space,
   * which is computed from Green's function. This function uses
   * already-built neighbor list of Elem* elem.
   * charge_type: "regularized"
   */
  Real local_potential_field(const Elem* elem,
                             const Point &p,
                             const std::string& charge_type) const;


  /**
   * Local electrical potential at a bead in an unbounded space,
   * which is computed from Green's function.
   * charge_type: "regularized"
   */
  Real local_potential_bead(const std::size_t& bead_id,
                            const std::string& charge_type) const;


private:

  // Poisson solver
  SolverPoisson _solver_poisson;

  // Assemble Poisson system
  AssemblePoisson* _assemble_poisson;

  // Get a pointer to GGEMPoisson
  GGEMPoisson* ggem_poisson;

};

} // end namespace libMesh
