#pragma once

#include "fix.h"
#include "fix_point_LJCut.h"
#include "fix_point_Gaussian.h"
#include "fix_point_Gaussian_PolymerChain.h"
#include "fix_point_WormLikeSpring.h"
#include "../pm_linear_implicit_system.h"

namespace libMesh
{
class FixFactory
{
public:
	// constructor
	FixFactory() {};
	// destructor
	~ FixFactory() {}


// add else branch if you need to implment new force field 
	Fix* GetFix(std::string fix_name, PMLinearImplicitSystem& system)
	{
		if(fix_name == "point_LJCut"){
			return new FixPointLJCut(system);
		}
		else if (fix_name == "point_Gaussian"){
			return new FixPointGaussian(system);
		}
		else if (fix_name == "point_Gaussian_PolymerChain"){
			return new FixPointGaussianPolymerChain(system);
		}
		else if (fix_name == "point_WormLikeSpring"){
			return new FixPointWormLikeSpring(system);
		}
		else{
			std::cout <<"Error: undefined force type: " << fix_name	<< std::endl;
			libmesh_error();
		}
	}
};

}
