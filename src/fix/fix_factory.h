#pragma once

#include "fix.h"
#include "fix_point_lj_cut.h"
#include "fix_point_gaussian.h"
#include "fix_point_gaussian_dna.h"
#include "fix_point_wls.h"
#include "fix_point_slitWall_lj_cut.h"
#include "fix_point_sphereWall_lj_cut.h"
#include "fix_point_discretizedWall_lj_cut.h"
#include "fix_point_slitWall_empirical_dna.h"
#include "fix_point_sphereWall_empirical_dna.h"
#include "fix_point_constant.h"
#include "fix_rigid_surface_constraint.h"
#include "fix_rigid_sedimentation.h"
#include "fix_rigid_analytic_polarization.h"
#include "fix_rigid_lj_cut.h"
#include "fix_rigid_slitWall_lj_cut.h"
#include "fix_rigid_sphereWall_lj_cut.h"
#include "../pm_linear_implicit_system.h"

namespace libMesh
{
class FixFactory
{
public:
	// constructor
	FixFactory(){};
	// destructor
	virtual ~FixFactory(){};


// add else branch if you need to implment new force field 
	Fix* buildFix(std::string& fix_name, PMLinearImplicitSystem& pm_system);

private:
	std::string particle_type;
	std::string wall_type;
};

}
