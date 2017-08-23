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
#include "../rigid_particle.h"
#include "../particle_mesh.h"

using std::cout;
using std::endl;
using std::string;

namespace libMesh{

class CopssRigidParticleSystem : public Copss
{
public:
	
	CopssRigidParticleSystem (CopssInit& init);

	~CopssRigidParticleSystem();

	// integrator
	void run(EquationSystems& equation_systems) override;

protected:
	// override read_particle_info() function in Copss class
	void read_particle_info () override;

	// create objects, polymer chains
	void create_object() override;

	// attach mesh spring network to particle mesh
	void attach_mesh_spring_network();

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Build MeshSpringNetwork according to the particle's mesh, which will be used to
   apply the rigid-body constraint force.

   Note: if the particles use different meshes, or have different sizes,
   we need to build them for each of particles!
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	void create_object_mesh() override;

	// attach object mesh to system
	void attach_object_mesh(PMLinearImplicitSystem& system) override;

	// set parameters for equations systems
	void set_parameters(EquationSystems& equation_systems) override;

	// update object due to PBC after check_wall()
	void update_object(std::string stage) override;




private:

	std::string _particle_mesh_type;

	std::vector<std::string> _particle_mesh_file;

	std::vector<Real> _hsize_solid;

	Real _hmins; // surface mesh hmin

	Real _hmaxs;

  	ParticleMesh<3>* _particle_mesh;

  	std::vector<MeshSpringNetwork*> _mesh_spring_network;

  	std::size_t _num_particles;

  	const std::string smesh_file_name = "particle_surface_mesh.e";

  	std::ostringstream smesh_file_name_i;

  	bool _write_particle_mesh;


};


}