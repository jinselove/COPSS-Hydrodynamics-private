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

#include "fix_rigid.h"
namespace libMesh
{

class FixRigidAnalyticPolarization : public FixRigid
{
public:
	FixRigidAnalyticPolarization(PMLinearImplicitSystem& pm_sys);

	~FixRigidAnalyticPolarization(){};

	void print_fix();

	void initParams();

	void compute();
	
	void compute_coulombic_force(std::vector<Point>& body_force);
	
	void compute_first_order_polarization(std::vector<Point>& body_force);
	
	void force_pol(std::vector<Point>& body_force, int i, int j, int k);
	
private:
	
	int order;
	
	double xlo= 0.0;
	
	double xhi = 1.0;
	
	const int ngauss = 5;
	
	const Real _xg0[5] = {-0.9061798459386640,-0.5384693101056831,0.00000000000000000,0.5384693101056831,0.9061798459386640};
	
	const Real _wg0[5] = {0.2369268850561891,0.4786286704993665,0.5688888888888889,0.4786286704993665,0.2369268850561891};

	Real eouter = 1.;

	Real einner = 1.;

	Real _e, _g;
	
	// helper variables for image method kernel function 
	Real Rxkj, Rykj, Rzkj, Rkj2, rkj;
	Real Rxij, Ryij, Rzij, Rij2, rij;
	Real ukj, vkj, wkj;
	Real aux1, aux2;
    Real auxv_x_integ, auxv_y_integ, auxv_z_integ, aux3_integ, aux3Sqrt_integ;
    Real auxv_x_delta, auxv_y_delta, auxv_z_delta, aux3_delta, aux3Sqrt_delta;
    Real radius;
};

}
