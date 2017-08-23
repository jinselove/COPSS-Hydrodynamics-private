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

#include "../copss.h"
#include "../polymer_chain.h"
#include "../include_fixes.h"

using std::cout;
using std::endl;
using std::string;

namespace libMesh{

class CopssPointParticleSystem : public Copss
{
public:
	
	CopssPointParticleSystem (CopssInit& init);

	~CopssPointParticleSystem();

	// integrator
	void run(EquationSystems& equation_systems) override;

protected:
	std::string point_particle_model;
	// ===========for beads and polymer chains
	unsigned int Nb; // total # of beads
	unsigned int Ns; // total # of springs per Chain
	unsigned int nBonds; // total # of springs/bonds
	// ===========for polymer chains
	unsigned int nChains; // total # of chains
	Real bk; // Kuhn length (um)
	Real Nks; // # of Kuhn length per Chain
	Real Ss2; // (um^2)
	Real q0; //maximum spring length (um)
	Real chain_length; // contour length of the spring (um)
	Real Dc; // Diffusivity of the chain (um^2/s)

	// extra parameters for dynamic process
	Real max_spring_len;
	bool chain_broken;

	//void read_data(std::string control_file){};
	// override read_particle_parameters() function in Copss class
	void read_particle_info () override;

	// create objects, polymer chains
	void create_object() override;

	// create object mesh
	void create_object_mesh() override;

	// attach object mesh to system
	void attach_object_mesh(PMLinearImplicitSystem& system) override;

	// set parameters for equations systems
	void set_parameters(EquationSystems& equation_systems) override;

	// update object due to PBC after check_wall()
	void update_object(std::string stage) override;




private:

	PolymerChain* polymer_chain;
	//std::unique_ptr<PolymerChain> polymer_chain;

	//std::unique_ptr<PointMesh<3> > point_mesh; 



};


}
