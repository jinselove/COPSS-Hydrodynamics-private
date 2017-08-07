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




#include "libmesh/libmesh_logging.h"
#include "libmesh/equation_systems.h"

#include "particle_mesh.h"
#include "point_mesh.h"
#include "pm_periodic_boundary.h"
#include "force_field.h"


using namespace libMesh;


// ======================================================================
ForceField::ForceField(PMLinearImplicitSystem& pm_sys,
                       ElasticitySystem& el_sys)
: _pm_system(&pm_sys), _elastic_system(&el_sys)
{
  // initialize the private memebers
  _dim      = _pm_system->get_mesh().mesh_dimension();
  
  _particle_type = _pm_system->get_equation_systems().parameters.get<std::string>("particle_type");
  
  _particle_mesh = _pm_system->particle_mesh();
  
  _point_mesh = _pm_system->point_mesh();
  
  _num_particles = _particle_mesh -> num_particles();

  _num_points = _point_mesh -> num_particles();

  _wall_type = _pm_system->get_equation_systems().parameters.get<std::string>("wall_type");

  _wall_params = _pm_system->get_equation_systems().parameters.get<std::vector<Real>>(_wall_type);

  _box_min = _point_mesh->pm_periodic_boundary()->box_min();
  
  _box_max = _point_mesh->pm_periodic_boundary()->box_max();

  _box_len =  _point_mesh->pm_periodic_boundary()->box_length();
  
  _periodic = _point_mesh->pm_periodic_boundary()->periodic_direction();
  
  _inlet = _point_mesh->pm_periodic_boundary()->inlet_direction();

  this -> check_force_field();

}



// ======================================================================
ForceField::ForceField(PMLinearImplicitSystem& pm_sys)
: _pm_system(&pm_sys)
{
  // initialize the private memebersi

  _point_mesh = _pm_system->point_mesh();
  
  _dim      = pm_sys.get_mesh().mesh_dimension();
  
  _particle_type = _pm_system->get_equation_systems().parameters.get<std::string>("particle_type");

  _wall_type = _pm_system->get_equation_systems().parameters.get<std::string>("wall_type");

  _wall_params = _pm_system->get_equation_systems().parameters.get<std::vector<Real>>(_wall_type);

  _bead_r   = pm_sys.get_equation_systems().parameters.get<Real>("bead radius");

  _kBT      = pm_sys.get_equation_systems().parameters.get<Real>("kBT");

  if(_particle_type == "rigid_particle"){
    _particle_mesh = _pm_system->particle_mesh();
    _num_particles = _particle_mesh -> num_particles();
  }
  else{
    _point_particle_model = _pm_system->get_equation_systems().parameters.get<std::string>("point_particle_model");
    if(_point_particle_model == "polymer_chain"){
        _Ss2      = pm_sys.get_equation_systems().parameters.get<Real>("Ss2");
        _bk       = pm_sys.get_equation_systems().parameters.get<Real>("bk");
        _Nks      = pm_sys.get_equation_systems().parameters.get<Real>("Nks");
    }  
  }

  _num_points = _point_mesh -> num_particles();
  
  //box info
  _box_min = _point_mesh->pm_periodic_boundary()->box_min();
  
  _box_max = _point_mesh->pm_periodic_boundary()->box_max();

  _box_len =  _point_mesh->pm_periodic_boundary()->box_length();
  
  _periodic = _point_mesh->pm_periodic_boundary()->periodic_direction();
  
  _inlet = _point_mesh->pm_periodic_boundary()->inlet_direction();

  this -> check_force_field();
}



// ======================================================================
ForceField::~ForceField()
{
  // do nothing
}




// ======================================================================
void ForceField::check_force_field()
{
START_LOG("check_force_field()", "ForceField");
_pp_force_types = _pm_system->get_equation_systems().parameters.get<std::vector<std::string>>("pp_force_types");

_num_pp_force_types = _pp_force_types.size();

_pw_force_types = _pm_system->get_equation_systems().parameters.get<std::vector<std::string>>("pw_force_types");

_num_pw_force_types = _pw_force_types.size();

for (unsigned int i = 0; i < _num_pp_force_types + _num_pw_force_types; i++ ){
 (i < _num_pp_force_types) ? (forceTypeMap[_pp_force_types[i]]= -1) : (forceTypeMap[_pw_force_types[i-_num_pp_force_types]] = -1) ;
}
// particle-particle force
forceTypeMap["pp_ev_gaussian"]=pp_ev_gaussian;

forceTypeMap["pp_ev_gaussian_polymerChain"]=pp_ev_gaussian_polymerChain;

forceTypeMap["pp_ev_lj_cut"]=pp_ev_lj_cut;

forceTypeMap["pp_ev_lj_repulsive"]= pp_ev_lj_repulsive;

forceTypeMap["pp_ev_harmonic_repulsive"]=pp_ev_harmonic_repulsive;

forceTypeMap["pp_wormLike_spring"]=pp_wormLike_spring;

forceTypeMap["p_constant"]=p_constant;

forceTypeMap["p_gravity"]=p_gravity;

forceTypeMap["p_surface_constraint"]=p_surface_constraint;

// particle-wall force
forceTypeMap["pw_ev_empirical_polymerChain"]=pw_ev_empirical_polymerChain;

forceTypeMap["pw_ev_lj_cut"]=pw_ev_lj_cut;

forceTypeMap["pw_ev_lj_repulsive"]=pw_ev_lj_repulsive;

forceTypeMap["pw_ev_harmonic_repulsive"]=pw_ev_harmonic_repulsive;  

// Attach particle-particle forces
for(int i = 0; i < _num_pp_force_types; i++){
  std::vector<Real> params = _pm_system->get_equation_systems().parameters.get<std::vector<Real>> (_pp_force_types[i]);
  switch (forceTypeMap[_pp_force_types[i]]){
    case pp_ev_gaussian:
      if(_particle_type !="point_particle" or _point_particle_model != "bead"){
          std::cout << std::endl << "*******************Error message*********************" << std::endl
                  << "The force field 'pw_ev_gaussian' is for 'bead' models, but not for 'polymer_chain' model" <<std::endl
                  << "Try force field 'pw_ev_gaussian_polymerChain' instead!!!!" <<std::endl
                  << "****************************************" << std::endl;    
          libmesh_error();
      }
      if(params.size() != 2){
        std::cout << std::endl << "********************Error message********************" << std::endl
                  << "---------------> The force type 'pp_ev_gaussian' requires 2 parameter (c1,c2)" << std::endl
                  << "****************************************" << std::endl;
        libmesh_error(); 
      }
      break;
    case pp_ev_gaussian_polymerChain:
      if(_particle_type != "point_particle" or _point_particle_model != "polymer_chain"){
          std::cout << std::endl << "*******************Error message*********************" << std::endl
                  << "The force field 'pw_ev_gaussian_polymerChain' is for 'polymer_chain' models, but not for 'bead' models" <<std::endl
                  << "Try force field 'pw_ev_gaussian' instead!!!!" <<std::endl
                  << "****************************************" << std::endl;    
          libmesh_error();
      }  
      Real ev, c1, c2;
      if (params.size()!=1)
      {
        std::cout << std::endl << "********************Error message********************" << std::endl
                  << "---------------> The force type 'pp_ev_gaussian_polymerChain' requires 1 parameter (ev)" << std::endl
                  << "****************************************" << std::endl;
        libmesh_error();        
      }
      break;
    case pp_ev_lj_cut:
      if(params.size()!=3){
        std::cout << std::endl << "********************Error message********************" << std::endl
                  << "---------------> The force type 'pp_ev_lj_cut' requires 3 parameter (epsilon, sigma, rcut) (dimensionless form)" << std::endl
                  << "****************************************" << std::endl;
        libmesh_error();    
      }
      break;
    case pp_ev_lj_repulsive:
      if(params.size()!=2){
        std::cout << std::endl << "********************Error message********************" << std::endl
                  << "---------------> The force type 'pp_ev_lj_cut_repulsive' requires 2 parameter (epsilon, sigma, rcut) (dimensionless form)" << std::endl
                  << "****************************************" << std::endl;
        libmesh_error();    
      }
      break;
    case pp_ev_harmonic_repulsive:
      if(params.size() != 2){
        std::cout << std::endl << "********************Error message********************" << std::endl
                  << "---------------> The force type 'pp_ev_harmonic_repulsive' requires 2 parameter (k, r0) (dimensionless form)" << std::endl
                  << "****************************************" << std::endl;
        libmesh_error();    
      } 
      break;
    case pp_wormLike_spring:
      if(_particle_type != "point_particle" or _point_particle_model != "polymer_chain"){
          std::cout << std::endl << "*******************Error message*********************" << std::endl
          << "The force field 'pp_wormLike_spring' is for 'polymer_chain' models, but not for 'bead' models" <<std::endl
          << "*********************************************" << std::endl;    
          libmesh_error();
      }  
      break;
    case p_constant:
      if(params.size() != 3){
        std::cout << std::endl << "********************Error message********************" << std::endl
          << "---------------> The force type 'p_constant' requires 3 parameter (fx,fy,fz)" << std::endl
          << "****************************************" << std::endl;
        libmesh_error();    
      }

    //  attach_p_constant(_pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("p_constant"));
    //  std::cout << "attached p_constant force" << std::endl; 
      break;
    case p_gravity:
      if(params.size() != 3)
      {
        std::cout << std::endl << "********************Error message********************" << std::endl
                  << "---------------> The force type 'p_gravity' requires 3 parameter (body force density) (fx,fy,fz)" << std::endl
                  << "****************************************" << std::endl;
        libmesh_error();    
      }
      break;
    case p_surface_constraint:
      if(_particle_type != "rigid_particle"){
        std::cout << std::endl << "*******************Error message*********************" << std::endl
            << "p_surface_constraint (surface spring force) only applies to rigid particle" << std::endl
            << "****************************************" << std::endl;    
        libmesh_error();
      }
      if(params.size() != 1){
        std::cout << std::endl << "********************Error message********************" << std::endl
            << "---------------> The force type 'p_surface_constraint' requires 1 parameter (k0)" << std::endl
            << "****************************************" << std::endl;
        libmesh_error();   
      }
      break;
    default:
     std::cout << std::endl << "********************Error message********************" << std::endl
            << "---------------> The force type '"<< _pp_force_types[i] << "' cannot be found!!!" << std::endl
            << "****************************************" << std::endl;
     libmesh_error();
     break;
  }
}
// Attach particle-wall forces
for (int i = 0; i < _num_pw_force_types; i++){
  std::vector<Real> params = _pm_system->get_equation_systems().parameters.get<std::vector<Real>> (_pw_force_types[i]);
  switch(forceTypeMap[_pw_force_types[i]]){
    case pw_ev_empirical_polymerChain:
      if(_particle_type != "point_particle" or _point_particle_model != "polymer_chain"){
        std::cout << std::endl << "*******************Error message*********************" << std::endl
                << "The force field 'pp_wormLike_spring' is for 'polymer_chain' models, but not for 'bead' models" <<std::endl
                << "*********************************************" << std::endl;    
        libmesh_error();
      }    
      break;
    case pw_ev_lj_cut:
      if(params.size() != 3){
        std::cout << std::endl << "********************Error message********************" << std::endl
                  << "---------------> The force type 'pw_ev_lj_cut' requires 3 parameter (epsilon, sigma, rcut) (dimensionless form)" << std::endl
                  << "****************************************" << std::endl;
        libmesh_error();    
      }
      break;
    case pw_ev_lj_repulsive:
      if(params.size() != 2){
        std::cout << std::endl << "********************Error message********************" << std::endl
                  << "---------------> The force type 'pw_ev_lj_repulsive' requires 3 parameter (epsilon, sigma) (dimensionless form)" << std::endl
                  << "****************************************" << std::endl;
        libmesh_error();    
      }        
      break;
    case pw_ev_harmonic_repulsive:
      if (params.size() != 2){
        std::cout << std::endl << "********************Error message********************" << std::endl
                  << "---------------> The force type 'pw_ev_harmonic_repulsive' requires 2 parameter (k, r0) (dimensionless form)" << std::endl
                  << "****************************************" << std::endl;
        libmesh_error();    
      }     
      break;
    default:
      std::cout << std::endl << "********************Error message********************" << std::endl
            << "---------------> The force type '"<< _pw_force_types[i] << "' cannot be found!!!" << std::endl
            << "****************************************" << std::endl;
      libmesh_error();
      break;
  }
}

STOP_LOG("check_force_field()", "ForceField");
}


// ======================================================================
void ForceField::reinit_force_field()
{
  START_LOG ("reinit_force_field()", "ForceField");
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   0. zero the force vectors on all the points.
   This has been done in point_mesh->reinit()
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //this->reinit_point_force_zeros();
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   1. Apply the stiff spring force for nodes on each rigid particle.
   Each node on the surface is connected with its neighboring nodes by
   stiff springs.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    // Attach particle-particle forces
    for(int i = 0; i < _num_pp_force_types; i++){
      switch (forceTypeMap[_pp_force_types[i]]){
        case pp_ev_gaussian:
          attach_pp_ev_gaussian(_pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("pp_ev_gaussian"));
         // std::cout << "attached pp_ev_gaussian force" << std::endl;
          break;
        case pp_ev_gaussian_polymerChain:
          attach_pp_ev_gaussian_polymerChain(_pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("pp_ev_gaussian_polymerChain"));
         // std::cout << "attached pp_ev_gaussian_polymerChain force" << std::endl;
          break;
        case pp_ev_lj_cut:
          attach_pp_ev_lj_cut(_pm_system ->get_equation_systems().parameters.get<std::vector<Real>>("pp_ev_lj_cut"));
         // std::cout << "attached pp_ev_lj_cut force" << std::endl;
          break;
        case pp_ev_lj_repulsive:
          attach_pp_ev_lj_repulsive(_pm_system ->get_equation_systems().parameters.get<std::vector<Real>>("pp_ev_lj_repulsive"));
         // std::cout << "attached pp_ev_lj_repulsive force" << std::endl;
          break;
        case pp_ev_harmonic_repulsive:
          attach_pp_ev_harmonic_repulsive(_pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("pp_ev_harmonic_repulsive"));
        //  std::cout << "attached pp_ev_harmonic_repulsive force" << std::endl;  
          break;
        case pp_wormLike_spring:
          attach_pp_wormLike_spring();
        //  std::cout << "attached pp_wormLike_spring force" << std::endl;
          break;
        case p_constant:
          attach_p_constant(_pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("p_constant"));
        //  std::cout << "attached p_constant force" << std::endl; 
          break;
        case p_gravity:
          attach_p_gravity(_pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("p_gravity"));
          break;
        case p_surface_constraint:
          attach_p_surface_constraint(_pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("p_surface_constraint"));
          break;
        default:
         std::cout << std::endl << "********************Error message********************" << std::endl
                << "---------------> The force type '"<< _pp_force_types[i] << "' cannot be found!!!" << std::endl
                << "****************************************" << std::endl;
         libmesh_error();
         break;
      }
    }
    // Attach particle-wall forces
    for (int i = 0; i < _num_pw_force_types; i++){
      switch(forceTypeMap[_pw_force_types[i]]){
        case pw_ev_empirical_polymerChain:
          attach_pw_ev_empirical_polymerChain();
        //  std::cout << "attached pw_ev_empirical_polymerChain force" << std::endl;         
          break;
        case pw_ev_lj_cut:
	        attach_pw_ev_lj_cut(_pm_system ->get_equation_systems().parameters.get<std::vector<Real>>("pw_ev_lj_cut"));
        //  std::cout << "attached pw_ev_lj_cut force" << std::endl;         
          break;
        case pw_ev_lj_repulsive:
          attach_pw_ev_lj_repulsive(_pm_system ->get_equation_systems().parameters.get<std::vector<Real>>("pw_ev_lj_repulsive"));
        //  std::cout << "attached pw_ev_lj_repulsive force" << std::endl;         
          break;
        case pw_ev_harmonic_repulsive:
          attach_pw_ev_harmonic_repulsive(_pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("pw_ev_harmonic_repulsive"));
        //  std::cout << "attached pw_ev_harmonic_repulsive force" << std::endl;        
          break;
        default:
          std::cout << std::endl << "********************Error message********************" << std::endl
                << "---------------> The force type '"<< _pw_force_types[i] << "' cannot be found!!!" << std::endl
                << "****************************************" << std::endl;
          libmesh_error();
          break;
      }
    }
  STOP_LOG ("reinit_force_field()", "ForceField");
}

// ======================================================================
void ForceField::attach_pp_ev_gaussian(const std::vector<Real>& params)
{
  START_LOG("attach_pp_ev_gaussian(&params)", "ForceField");
  Real c1 = params[0];
  Real c2 = params[1];
  for(std::size_t i=0; i<_num_points; ++i)
  {  
    // apply the excluded volume force to each particle i
    std::vector<Real> evforce(_dim);
    this->compute_pp_ev_gaussian(c1,c2,i,evforce);
    _point_mesh->particles()[i]->add_particle_force(evforce);
  }
  STOP_LOG ("attach_pp_ev_gaussian(&params)", "ForceField");
}

// ======================================================================
void ForceField::attach_pp_ev_gaussian_polymerChain(const std::vector<Real>& params)
{
  START_LOG("attach_pp_ev_gaussian_polymerChain(&params)", "ForceField");
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop each point and apply forces
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
  Real ev = params[0] * _bead_r * _bead_r * _bead_r;
  Real c1   = ev*_Nks*_Nks*std::pow( 3./(4.*PI*_Ss2),1.5 ); //changed the sign of c1
  Real c2   = 3.*_bead_r*_bead_r/(4.*_Ss2);   

  for(std::size_t i=0; i<_num_points; ++i)
  {  
    // apply the excluded volume force to each particle i
    std::vector<Real> pforce(_dim);
    this->compute_pp_ev_gaussian(c1,c2,i,pforce);
    _point_mesh->particles()[i]->add_particle_force(pforce);
  }

  STOP_LOG("attach_pp_ev_gaussian_polymerChain(&params)", "ForceField");
}


// ======================================================================
void ForceField::compute_pp_ev_gaussian(Real& c1,
                                        Real& c2,
                                        const std::size_t&  p_id,
                                        std::vector<Real>& pforce )
{
  START_LOG ("compute_pp_ev_gaussian(c1,c2,p_id, &p_force)", "ForceField");
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Compute excluded volume force: force between particle and its neighboring beads
   Note: the neighbor_list returns p_id and the distance, which has already considered
   the periodic BCs with corrected bead-bead distance, see PointMesh()/ParticleMesh()
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  std::vector<std::pair<std::size_t,Real> > n_list = _point_mesh->particles()[p_id]->neighbor_list();
  const Point pti   = _point_mesh->particles()[p_id]->point();
  
  // Loop over each neigbhor
  for (std::size_t i=0; i<n_list.size(); ++i)
  {
    const std::size_t n_id  = n_list[i].first;
    std::vector<Real> f_ij;
    if(p_id != n_id)  // make sure this bead and the neighboring bead are not the same bead.
    {
      const Point ptj   = _point_mesh->particles()[n_id]->point();
      const Point r_ij  = _point_mesh->pm_periodic_boundary()->point_vector(pti,ptj);
     
      f_ij = this->gaussian_force(r_ij, c1, c2);
      for (std::size_t j=0; j<_dim; ++j) pforce[j] += f_ij[j];
    } // end if
  } // end for i-loop
  
  STOP_LOG ("compute_excluded_volume_force(c1, c2, p_id, &pforce)", "ForceField");
}

// ======================================================================
void ForceField::attach_pp_ev_lj_cut(const std::vector<Real>& params)
{
  START_LOG("attach_pp_ev_lj_cut(&params)", "ForceField");

  Real   epsilon = params[0];
  Real   sigma = params[1];
  Real   rcut = params[2];
  for(std::size_t i=0; i<_num_points; ++i)
  {  
    // apply lj force on particle i
    std::vector<Real> pforce(_dim);
    this->compute_pp_ev_lj_cut(epsilon,sigma,rcut,i,pforce);
    _point_mesh->particles()[i]->add_particle_force(pforce);
  }  

  STOP_LOG("attach_pp_ev_lj_cut(&params)", "ForceField");
}

// ======================================================================
void ForceField::attach_pp_ev_lj_repulsive(const std::vector<Real>& params)
{
  START_LOG("attach_pp_ev_lj_repulsive(&params)", "ForceField");
  Real epsilon = params[0];
  Real sigma = params[1];
  Real rcut = sigma * std::pow(2.,1./6.);
  for(std::size_t i=0; i<_num_points; ++i)
  {  
    // apply lj force on particle i
    std::vector<Real> pforce(_dim);
    // attach_pp_ev_lj_repulsive() function also call pp_ev_lj_cut
    // but r_cut is fixed to be sigma / 2^(1/6)
    this->compute_pp_ev_lj_cut(epsilon,sigma,rcut,i,pforce);
    _point_mesh->particles()[i]->add_particle_force(pforce);
  }  

  STOP_LOG("attach_pp_ev_lj_repulsive(&params)", "ForceField");
}

// ======================================================================
void ForceField::compute_pp_ev_lj_cut(const Real& epsilon,
                                      const Real& sigma,
                                      const Real& rcut,
                                      const std::size_t&  p_id,
                                      std::vector<Real>& pforce) const
{
  START_LOG("compute_pp_ev_lj_cut(&epsilon, &sigma, &rcut, &pforce)","ForceField");

  std::vector<std::pair<std::size_t,Real> > n_list = _point_mesh->particles()[p_id]->neighbor_list();
  const Point pti   = _point_mesh->particles()[p_id]->point();
  
  // Loop over each neigbhor
  for (std::size_t i=0; i<n_list.size(); ++i)
  {
    const std::size_t n_id  = n_list[i].first;
    if(p_id != n_id)  // make sure this bead and the neighboring bead are not the same bead.
    {
      const Point ptj   = _point_mesh->particles()[n_id]->point();
      const Point r_ij  = _point_mesh->pm_periodic_boundary()->point_vector(pti,ptj);
      if(r_ij.norm() <= rcut){
        std::vector<Real> f_ij = this->lj_force(r_ij, epsilon, sigma);
        for (std::size_t j=0; j<_dim; ++j) pforce[j] += f_ij[j];
      }
    } // end if
  } // end for i-loop 

  STOP_LOG("compute_pp_ev_lj_cut(&epsilon, &sigma, &rcut, &pforce)","ForceField");
}

// ======================================================================
void ForceField::attach_pp_ev_harmonic_repulsive(const std::vector<Real>& params){
  START_LOG("attach_pp_ev_harmonic_repulsive(&params)", "ForceField");
  // when r > sigma,e.g. sigma = 2*bead radius, pp_ev_harmonic_repulsive force = 0;
  for(std::size_t i=0; i<_num_points; ++i)
  {  
    std::vector<Real> pforce(_dim);
    this->compute_pp_ev_harmonic_repulsive(params[0],params[1],i,pforce);
    _point_mesh->particles()[i]->add_particle_force(pforce);
  }  

  STOP_LOG("attach_pp_ev_harmonic_repulsive(&params)", "ForceField");
}



void ForceField::compute_pp_ev_harmonic_repulsive(const Real& k,
                                                  const Real& r0,
                                                  const std::size_t& p_id,
                                                  std::vector<Real>& pforce) const
{
  START_LOG("compute_pp_ev_harmonic_repulsive(&epsilon,&sigma,&p_id,&pforce)","ForceField");

  std::vector<std::pair<std::size_t,Real> > n_list = _point_mesh->particles()[p_id]->neighbor_list();
  const Point pti   = _point_mesh->particles()[p_id]->point();

  // Loop over each neigbhor
  for (std::size_t i=0; i<n_list.size(); ++i)
  {
    const std::size_t n_id  = n_list[i].first;
    if(p_id != n_id)  // make sure this bead and the neighboring bead are not the same bead.
    {
      const Point ptj   = _point_mesh->particles()[n_id]->point();
      const Point r_ij  = _point_mesh->pm_periodic_boundary()->point_vector(pti,ptj);
      if(r_ij.norm() <= r0){
        std::vector<Real> f_ij = this->harmonic_force(r_ij, k, r0);
        for (std::size_t j=0; j<_dim; ++j) pforce[j] += f_ij[j];
      }
    } // end if
  } // end for i-loop 
  STOP_LOG("compute_pp_ev_harmonic_repulsive(&epsilon,&sigma,&p_id,&pforce)","ForceField");
}

// ======================================================================
void ForceField::attach_pp_wormLike_spring()
{
  START_LOG ("attach_wormLike_spring_force()", "ForceField");
  // if(_particle_type != "point_particle" or _point_particle_model != "polymer_chain"){
  //   std::cout << std::endl << "*******************Error message*********************" << std::endl
  //           << "The force field 'pp_wormLike_spring' is for 'polymer_chain' models, but not for 'bead' models" <<std::endl
  //           << "*********************************************" << std::endl;    
  //   libmesh_error();
  // }  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get the PointMesh objects
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const std::size_t  n_bonds    = _point_mesh->num_bonds();

  PolymerChain*   polymer_chain = _point_mesh->polymer_chain();
  const std::vector<std::vector<std::size_t> >& bonds = polymer_chain->bonds();

  // Worm-like-spring force parameters
  const Real c1  = _bead_r/(2.0*_bk);
  const Real Ls  = _Nks*_bk/_bead_r;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Loop each bond and apply forces
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for(std::size_t i=0; i<n_bonds; ++i)
  {
    // Connecting bead ids
    std::size_t bead_id_1 = bonds[i][1];  // connect bead 1
    std::size_t bead_id_2 = bonds[i][2];  // connect bead 2

    // Force on bead 1
    const Point pti   = _point_mesh->particles()[bead_id_1]->point();
    const Point ptj   = _point_mesh->particles()[bead_id_2]->point();
    const Point R_ij  = _point_mesh->pm_periodic_boundary()->point_vector(pti,ptj);
    const std::vector<Real> F_ij = this->spring_force_wls(R_ij,c1,Ls);

    // Force on bead 2
    std::vector<Real> F_ji (F_ij);
    for (std::size_t j=0; j<_dim; ++j) F_ji[j] = -F_ij[j];

    // Add forces to beads
    _point_mesh->particles()[bead_id_1]->add_particle_force(F_ij);
    _point_mesh->particles()[bead_id_2]->add_particle_force(F_ji);
  }

  STOP_LOG ("attach_wormLike_spring_force()", "ForceField");
}





// ======================================================================
void ForceField::attach_p_constant(const std::vector<Real>& params)
{
  START_LOG ("attach_p_constant()", "ForceField");
  for(std::size_t i=0; i<_num_points; ++i)
  {
   _point_mesh->particles()[i]->add_particle_force(params);   
  } // end for i-loop
  STOP_LOG ("attach_p_constant()", "ForceField");  
}

// ======================================================================
void ForceField::attach_p_gravity(const std::vector<Real>& params)
{
  START_LOG ("attach_p_gravity()", "ForceField");

  for(std::size_t i=0; i<_num_points; ++i)
  {
//     _point_mesh->particles()[i]->add_particle_force(params * );   
  } // end for i-loop

  STOP_LOG ("attach_p_gravity()", "ForceField");  
}

void ForceField::attach_p_surface_constraint(const std::vector<Real>& params)
{
  START_LOG("attach_surface_constraint()", "ForceField");

  Real k0 = params[0];
  std::size_t  point_start_id     = 0;
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Loop over each particle, and add the body force to each node.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  for(std::size_t i=0; i < _num_particles; ++i)
  {
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Check if this particle is on the periodic boundary. If so, move the nodes
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    const bool on_pb = _particle_mesh->particles()[i]->on_the_periodic_boundary();

    if( on_pb )
    {
      if(_pm_system->comm().rank()==0){
        printf("--->TEST in ForceField::attach_surface_constraint() \n");
        printf("         The particle %lu is on the periodic boundary!\n",i);
      }
      _particle_mesh->particles()[i]->rebuild_periodic_mesh();
    }
    
    
    // compute the rigid constraint force on each node
    std::vector<Point> surface_constraint_force;
    this->compute_surface_constraint_force(i,k0,surface_constraint_force);
  
    /*
     * Loop over each node through node iterator
     */
    const std::size_t n_nodes = _particle_mesh->particles()[i]->num_mesh_nodes();
    MeshBase& p_mesh = _particle_mesh->particles()[i]->mesh();
    MeshBase::node_iterator       nd     = p_mesh.active_nodes_begin();
    const MeshBase::node_iterator end_nd = p_mesh.active_nodes_end();
    for ( ; nd != end_nd; ++nd)
    {
      // Store a pointer to the element we are currently working on.
      Node* node = *nd;
      const dof_id_type node_id = node->id();
      
      // get the dof numbers at this node (only for force vector)
      std::vector<Real> gforce(_dim,0.);
      for(std::size_t k=0; k<_dim; ++k){
        gforce[k] = surface_constraint_force[node_id](k);
      } // end k-loop
      
      // ------------------ TEST: print out the surface constraint force ----------------------
      if(_pm_system->comm().rank()==0){
        printf("--->TEST:attach_surface_constraint() gforce = (%f,%f,%f)\n",gforce[0],gforce[1],gforce[2]);
        printf("         surface_constraint_force = (%f,%f,%f)\n",
               surface_constraint_force[node_id](0),surface_constraint_force[node_id](1),surface_constraint_force[node_id](2));
      }
      // -------------------------------------------------------------------------
      _point_mesh->particles()[point_start_id+node_id]->add_particle_force(gforce);
    } // end for nd-loop
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Move back the mesh nodes.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    if( on_pb ) _particle_mesh->particles()[i]->restore_periodic_mesh();
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     update the point start id
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    point_start_id += n_nodes;
    
  } // end for i-loop
  STOP_LOG("attach_surface_constraint()", "ForceField");  
}

// ======================================================================
void ForceField::compute_surface_constraint_force(const std::size_t& i, // the i-th particle
                                        const Real& k0,
                                        std::vector<Point>& nodal_force)
{
  START_LOG ("compute_surface_constraint_force()", "ForceField");
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   find the partilce center and number of tracking points on the surface.
   We need this center position to construct springs in the radial direction.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  MeshSpringNetwork* mesh_spring  = _particle_mesh->particles()[i]->mesh_spring_network();
  const std::size_t n_nodes = _particle_mesh->particles()[i]->num_mesh_nodes();
  const Point      p_center = _particle_mesh->particles()[i]->compute_centroid();
  // if(_pm_system->comm().rank()==0){
  //   printf("--->TEST in rigid_constraint_force() center of the %lu-th particle is (%f,%f,%f)\n",
  //          i,p_center(0),p_center(1),p_center(2));
  // }
  

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   loop over each tracking point(node), and compute spring forces
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  nodal_force.resize(n_nodes);
  for(std::size_t j=0; j<n_nodes; ++j)
  {
    // 1. Get the coords and neighboring nodes of the j-th node
    std::vector< std::pair<std::size_t, Real> > node_neighbors;
    node_neighbors    = mesh_spring->nodes_neighbors(j);
    const Point& pt0  = _particle_mesh->particles()[i]->mesh_point(j);
    
    // 2. Loop over each neighboring nodes and compute the SPRING forces
    std::vector<Real> spring_f(_dim,0.);
    for(std::size_t k=0; k<node_neighbors.size(); ++k)
    {
      const std::size_t& neigh_id = node_neighbors[k].first;
      const Real&        neigh_l0 = node_neighbors[k].second;
      const Point& ptk = _particle_mesh->particles()[i]->mesh_point(neigh_id);
      const Point R_ij = _particle_mesh->pm_periodic_boundary()->point_vector(pt0,ptk);
      
      std::vector<Real> sf  = this->spring_force_lhs(R_ij,neigh_l0,k0);
      for(std::size_t l=0; l<_dim; ++l) spring_f[l] += sf[l];
    } // end for k-loop
    
    // 3. compute the Spring force between the node the the particle center
    const Real  lc0    = mesh_spring->node_center_equilibrium_dist(j);
    const Point Rc_ij  = _particle_mesh->pm_periodic_boundary()->point_vector(pt0,p_center);
    std::vector<Real> sfc  = this->spring_force_lhs(Rc_ij,lc0,k0);
    for(std::size_t k=0; k<_dim; ++k) spring_f[k] += sfc[k];
    
    // 4. convert the point local id to the global id, and apply forces.
    for(std::size_t k=0; k<_dim; ++k) nodal_force[j](k) = spring_f[k];
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     * test: print out the spring force on each node.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//    if(_pm_system->comm().rank()==0)
//    {
//      printf("--->test in ForceField::rigid_constraint_force(), Particle %lu\n",i);
//      printf("--->Spring force (surface+center=total) on node %lu is fs = (%f,%f,%f)\n",
//             j,spring_f[0],spring_f[1],spring_f[2]);
//      printf("--->Spring force (center) on node %lu is               fs = (%f,%f,%f)\n",
//             j,sfc[0],sfc[1],sfc[2]);
//    }
    
  } // end for j-loop
  
  
  STOP_LOG ("compute_surface_constraint_force()", "ForceField");
}

void ForceField::attach_nodal(const std::vector<Real>& nodal_force_density)
{
  START_LOG("attach_nodal()", "ForceField");
  std::size_t  point_start_id     = 0;
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Loop over each particle, and add the body force to each node.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  for(std::size_t i=0; i < _num_particles; ++i)
  {
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Check if this particle is on the periodic boundary. If so, move the nodes
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    const bool on_pb = _particle_mesh->particles()[i]->on_the_periodic_boundary();

    if( on_pb )
    {
      if(_pm_system->comm().rank()==0){
        printf("--->TEST in ForceField::attach_nodal() \n");
        printf("         The particle %lu is on the periodic boundary!\n",i);
      }
      _particle_mesh->particles()[i]->rebuild_periodic_mesh();
    }

    std::vector<Point> nodal_force;
    _particle_mesh->particles()[i]->build_nodal_force(nodal_force_density,nodal_force);
    
    /*
     * Loop over each node through node iterator
     * Compute the gravitational force vector on each node
     */
    const std::size_t n_nodes = _particle_mesh->particles()[i]->num_mesh_nodes();
    MeshBase& p_mesh = _particle_mesh->particles()[i]->mesh();
    MeshBase::node_iterator       nd     = p_mesh.active_nodes_begin();
    const MeshBase::node_iterator end_nd = p_mesh.active_nodes_end();
    for ( ; nd != end_nd; ++nd)
    {
      // Store a pointer to the element we are currently working on.
      Node* node = *nd;
      const dof_id_type node_id = node->id();
      
      // get the dof numbers at this node (only for force vector)
      std::vector<Real> gforce(_dim,0.);
      for(std::size_t k=0; k<_dim; ++k){
        gforce[k] = nodal_force[node_id](k);
      } // end k-loop
      
      // ------------------ TEST: print out the nodal force ----------------------
      // if(_pm_system->comm().rank()==0){
      //   printf("--->TEST:reinit_force_field() gforce = (%f,%f,%f)\n",gforce[0],gforce[1],gforce[2]);
      //   printf("         nodal force = (%f,%f,%f)\n",
      //          nodal_force[node_id](0),nodal_force[node_id](1),nodal_force[node_id](2));
      //   printf("         rigid_nodal_force = (%f,%f,%f)\n",
      //          rigid_nodal_force[node_id](0),rigid_nodal_force[node_id](1),rigid_nodal_force[node_id](2));
      // }
      // -------------------------------------------------------------------------
      
      _point_mesh->particles()[point_start_id+node_id]->add_particle_force(gforce);
      
    } // end for nd-loop
    
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Move back the mesh nodes.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    if( on_pb ) _particle_mesh->particles()[i]->restore_periodic_mesh();
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     update the point start id
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    point_start_id += n_nodes;
  } // end for i-loop
  STOP_LOG("attach_nodal_force()", "ForceField");
}

// // ======================================================================
// void ForceField::attach_pp_friction(const std::vector<Real>& v_beads)
// {
//   START_LOG ("attach_friction_force(v_beads)", "ForceField");
  
//   const Real dmin = _pm_system->get_equation_systems().parameters.get<Real>("solid mesh size");
//   const Real Hf   = 1.0;  // friction coefficient
  
//   Real c1, c2;
//   if(_particle_type=="point_particle")
//   {
//     if(_point_particle_model == "polymer_chain"){
//  	    c1   = -ev*_Nks*_Nks*std::pow( 3./(4.*PI*_Ss2),1.5 );
//   	  c2   = 3.*_bead_r*_bead_r/(4.*_Ss2);
//     }
//     else if(_point_particle_model == "bead"){
//  	// c1   = -ev*_Nks*_Nks*std::pow( 3./(4.*PI*_Ss2),1.5 );
//   	// c2   = 3.*_bead_r*_bead_r/(4.*_Ss2);

// 	 c1 = -100;
//          c2 = +2;
//     }
//     else{std::cout << "\nWarning: excluded volume force constant is not defined between this type of point particles";}
//   }
//   else if(_particle_type=="rigid_particle")
//   {
//     c1 = -100;
//     c2 = +2;
//   }
//   else
//   {
//     libmesh_assert("*** Invalid particle_type in ForceField::modify_force_field()!");
//     libmesh_error();
//   }
//   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    Loop over each point, and compute the friction force
//    Compute excluded volume force: force between particle and its neighboring beads
//    Compute the friction force, which can be evaluated through the exv force
//    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//   for(std::size_t i=0; i<_num_points; ++i)
//   {
//     // point id, coords, and its neighboring points.
//     const std::size_t p_id      = _point_mesh->particles()[i]->id();
//     const std::size_t parent_id = _point_mesh->particles()[p_id]->parent_id();
//     const Point       pti       = _point_mesh->particles()[p_id]->point();
//     const std::vector<std::pair<std::size_t,Real> > n_list = _point_mesh->particles()[p_id]->neighbor_list();
//     const std::vector<Real> veli = _pm_system->point_velocity(v_beads,p_id);
    
//     // Loop over each neigbhor
//     std::vector<Real> friction(_dim,0.0);
//     for (std::size_t j=0; j<n_list.size(); ++j)
//     {
//       const std::size_t n_id  = n_list[j].first;
//       if(p_id != n_id)  // make sure this bead and the neighboring bead are not the same bead.
//       {
//         // compute the excluded volume force
//         const std::size_t parent_idj = _point_mesh->particles()[n_id]->parent_id();
//         const Point ptj   = _point_mesh->particles()[n_id]->point();
//         const Point r_ij  = _point_mesh->pm_periodic_boundary()->point_vector(pti,ptj);
//         const std::vector<Real> fv_ij = this->excluded_volume_force(r_ij, c1, c2);
//        // if the two points are from different parents, compute the friction force.
//         if(parent_id != parent_idj)
//         {
//           const std::vector<Real> velj = _pm_system->point_velocity(v_beads,n_id);
//           const std::vector<Real> ff_ij = this->friction_force(pti,ptj,veli,velj,fv_ij,Hf,dmin);
//           for (std::size_t j=0; j<_dim; ++j) friction[j] += ff_ij[j];
//         } // end if
//       } // end if
//     } // end for j-loop
//     // add friction force to particle i
//     _point_mesh->particles()[p_id]->add_particle_force(friction);
//   } // end for i-loop
//   STOP_LOG ("attach_friction_force(v_beads)", "ForceField");
// }



// ======================================================================
void ForceField::attach_pw_ev_empirical_polymerChain()
{
  START_LOG("attach_pw_ev_empirical_polymerChain()", "ForceField");

  // this particle-wall force field only works for polymer_chain model
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop each point and apply forces
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for(std::size_t i=0; i<_num_points; ++i)
  {
    std::vector<Real> pforce(_dim);
    this->compute_pw_ev_empirical_polymerChain(i,pforce);
    _point_mesh->particles()[i]->add_particle_force(pforce);    
  } // end for i-loop
  
  STOP_LOG ("attach_pw_ev_empirical_polymerChain()", "ForceField");
}

// ======================================================================
void ForceField::compute_pw_ev_empirical_polymerChain(const std::size_t&  p_id,
                                             std::vector<Real>& pforce ) const
{
  START_LOG ("compute_pw_ev_empirical_polymerChain(p_id, &pforce)", "ForceField");

  // Coefficients.
  const Real c1 = _bead_r/_bk;
  const Real c2  = c1/std::sqrt(_Nks);
  const Real d0 = 0.5/c2;
  const Real c0  = 25.0*c1;
  if(_wall_type == "slit"){  
  //  retrieve bead position and the box boundary
  const Point pti     = _point_mesh->particles()[p_id]->point();
  // compute the particle-wall interaction force
  for (std::size_t j=0; j<_dim; ++j){  // loop over each direction
    // compute the distance to the wall if NO periodic boundary and no inlet 
    if ( _periodic[j]==false and _inlet[j]==false) {
      Point r_i_lo, r_i_hi;
      // distance to the parallel walls
      r_i_lo(j) = _box_min(j) - pti(j) ; 
      r_i_hi(j) = _box_max(j) - pti(j);       
      // bead-wall interaction force
      pforce[j] += this->polymer_wall_empirical_force(r_i_lo, c0, d0)[j];
      pforce[j] += this->polymer_wall_empirical_force(r_i_hi, c0, d0)[j];
    } // end if
  } // end for j-loop
  } // end if (wall_type == "slit")
  else if(_wall_type == "sphere"){
    // Retrieve bead position
    Real cavity_radius = _wall_params[0];
    const Point pti      = _point_mesh->particles()[p_id]->point();
    const Point r_i_sphereWall = pti.unit() * (cavity_radius-pti.norm()); 
    // bead-wall interaction force
    pforce = this->polymer_wall_empirical_force(r_i_sphereWall, c0, d0);
  } // end if (wall_type == "sphere")
  else{
  std::cout << "*** Error: wall_type '"<< _wall_type
                  << "' is not defined (check typo first)!! (in compute particle-wall force)" << std::endl;
  libmesh_error();
  }

  STOP_LOG ("compute_pw_ev_empirical_polymerChain(p_id,&pforce)", "ForceField");
}

// ======================================================================
void ForceField::attach_pw_ev_lj_cut(const std::vector<Real>& params)
{
  START_LOG("attach_pw_ev_lj_cut(&params)", "ForceField");
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop each point and apply forces
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Real  epsilon = params[0];
  Real  sigma = params[1];
  Real  rcut = params[2];
  for(std::size_t i=0; i<_num_points; ++i)
  {
    std::vector<Real> pforce(_dim);
    this->compute_pw_ev_lj_cut(epsilon,sigma,rcut,i,pforce);
    _point_mesh->particles()[i]->add_particle_force(pforce);    
  } // end for i-loop
  
  STOP_LOG ("attach_pw_ev_lj_cut(&params)", "ForceField");
}

// ======================================================================
void ForceField::attach_pw_ev_lj_repulsive(const std::vector<Real>& params)
{
  START_LOG("attach_pw_ev_lj_cut(&params)", "ForceField");
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop each point and apply forces
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Real epsilon = params[0];
  Real sigma = params[1];
  Real rcut = sigma * std::pow(2., 1./6.);
  for(std::size_t i=0; i<_num_points; ++i)
  {
    std::vector<Real> pforce(_dim,0.);
    this->compute_pw_ev_lj_cut(epsilon,sigma,rcut,i,pforce);
    _point_mesh->particles()[i]->add_particle_force(pforce);    
  } // end for i-loop
  
  STOP_LOG ("attach_pw_ev_lj_cut(&params)", "ForceField");
}


// ======================================================================
void ForceField::compute_pw_ev_lj_cut(const Real& epsilon,
                                      const Real& sigma,
                                      const Real& rcut,
                                      const std::size_t& p_id,
                                      std::vector<Real>& pforce) const
{
  START_LOG ("compute_pw_ev_lj_cut(&epsilon,&sigma,&rcut,&p_id,&pforce)", "ForceField");

  const Point pti     = _point_mesh->particles()[p_id]->point();

  if(_wall_type == "slit"){  
    // compute the particle-wall interaction force
    for(std::size_t j = 0; j<_dim;++j){
      Point r_i_lo, r_i_hi;
      if(_periodic[j]==false and _inlet[j]==false){
        r_i_lo(j) = _box_min(j) - pti(j);// opposite sign compared to rij = rj - r_i 
        r_i_hi(j) = _box_max(j)- pti(j);// opposite sign compared to rij = rj - r_i 
        if(r_i_lo.norm() <= rcut) pforce[j] += this->lj_force(r_i_lo, epsilon,sigma)[j];
        if(r_i_hi.norm() <= rcut) pforce[j] += this->lj_force(r_i_hi, epsilon,sigma)[j];
      }// end of
    } // end for (i<_dim)
  }
  else if(_wall_type == "sphere"){
    // Retrieve wall_info
    Real cavity_radius = _wall_params[0];
    const Real  pti_norm = pti.norm();
    const Point pti_unit = pti.unit();
    // vector point from particle i to the sphere wall
    const Point r_i_sphereWall = pti_unit * (cavity_radius-pti_norm); 
    // bead-wall interaction force
    if(r_i_sphereWall.norm() <= rcut) pforce = this->lj_force(r_i_sphereWall, epsilon, sigma);
  } // end (else if)
  else{
    std::cout << "*** Error: wall_type '"<< _wall_type
                    << "' is not defined (check typo first)!! (in compute particle-wall force)" << std::endl;
    libmesh_error();
  } // end else

  STOP_LOG ("compute_pw_ev_lj_cut(&epsilon,&sigma,&rcut,&p_id,&pforce)", "ForceField");
}


// ======================================================================
void ForceField::attach_pw_ev_harmonic_repulsive(const std::vector<Real>& params)
{
  START_LOG ("attach_pw_ev_harmonic_repulsive(&params)", "ForceField");
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop each point and apply forces
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for(std::size_t i=0; i<_num_points; ++i)
  {
    std::vector<Real> pforce(_dim,0.);
    this->compute_pw_ev_harmonic_repulsive(params[0],params[1],i,pforce);
    _point_mesh->particles()[i]->add_particle_force(pforce);    
  } // end for i-loop


  STOP_LOG ("attach_pw_ev_harmonic_repulsive(&params)", "ForceField");
}

// ======================================================================d
void ForceField::compute_pw_ev_harmonic_repulsive(const Real& k,
                                      const Real& r0,
                                      const std::size_t& p_id,
                                      std::vector<Real>& pforce) const
{
  START_LOG ("compute_pw_ev_harmonic_repulsive(&k, &r0, &p_id, &pforce)", "ForceField");

  const Point pti     = _point_mesh->particles()[p_id]->point();

  if(_wall_type == "slit"){  
    // compute the particle-wall interaction force
    for(std::size_t j = 0; j<_dim;++j){
      Point r_i_lo, r_i_hi;
      if(_periodic[j]==false and _inlet[j]==false){
        r_i_lo(j) = _box_min(j) - pti(j);
        r_i_hi(j) = _box_max(j) - pti(j);
        if(r_i_lo.norm() <= r0) pforce[j] += this->harmonic_force(r_i_lo, k, r0)[j];
        if(r_i_hi.norm() <= r0) pforce[j] += this->harmonic_force(r_i_hi, k, r0)[j];
      }// end of
    } // end for (i<_dim)
  }
  else if(_wall_type == "sphere"){
    // Retrieve wall_info
    Real cavity_radius = _wall_params[0];
    const Real  pti_norm = pti.norm();
    const Point pti_unit = pti.unit();
    // vector point from particle i to the sphere wall
    Point r_i_sphereWall = pti_unit * (cavity_radius-pti_norm);
    // bead-wall interaction force
    if(r_i_sphereWall.norm() <= r0) pforce = this->harmonic_force(r_i_sphereWall, k, r0);
  } // end (else if)
  else{
    std::cout << "*** Error: wall_type '"<< _wall_type
                    << "' is not defined (check typo first)!! (in compute particle-wall force)" << std::endl;
    libmesh_error();
  } // end else

  STOP_LOG ("compute_pw_ev_harmonic_repulsive(&k, &r0, &p_id, &pforce)", "ForceField");
}



// // ======================================================================
// void ForceField::modify_force_field(const std::vector<Real>& v_beads)
// {
//   START_LOG ("modify_force_field()", "ForceField");
  
//   const Real ev   = _pm_system->get_equation_systems().parameters.get<Real>("ev");
//   const Real dmin = _pm_system->get_equation_systems().parameters.get<Real>("solid mesh size");
//   const Real Hf   = 1.0;  // friction coefficient
  
//   Real c1, c2;
//   if(_particle_type=="point_particle")
//   {
//     c1   = -ev*_Nks*_Nks*std::pow( 3./(4.*PI*_Ss2),1.5 );
//     c2   = 3.*_bead_r*_bead_r/(4.*_Ss2);
//   }
//   else if(_particle_type=="rigid_particle")
//   {
//     c1 = -100;
//     c2 = +2;
//   }
//   else
//   {
//     libmesh_assert("*** Invalid particle_type in ForceField::modify_force_field()!");
//     libmesh_error();
//   }
// //  printf("--->TEST: ev = %f, _Nks = %f, Ss2 = %f, c1 = %f, c2 = %f\n",
// //         ev, _Nks, _Ss2, c1, c2);
  
  
//   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    Loop over each point, and compute the friction force
//    Compute excluded volume force: force between particle and its neighboring beads
//    Compute the friction force, which can be evaluated through the exv force
//    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//   for(std::size_t i=0; i<_num_points; ++i)
//   {
//     // point id, coords, and its neighboring points.
//     const std::size_t p_id      = _point_mesh->particles()[i]->id();
//     const std::size_t parent_id = _point_mesh->particles()[p_id]->parent_id();
//     const Point       pti       = _point_mesh->particles()[p_id]->point();
//     const std::vector<std::pair<std::size_t,Real> > n_list = _point_mesh->particles()[p_id]->neighbor_list();
//     const std::vector<Real> veli = _pm_system->point_velocity(v_beads,p_id);
    
//     // Loop over each neigbhor
//     std::vector<Real> evforce(_dim,0.0);
//     for (std::size_t j=0; j<n_list.size(); ++j)
//     {
//       const std::size_t n_id  = n_list[j].first;
//       if(p_id != n_id)  // make sure this bead and the neighboring bead are not the same bead.
//       {
//         // compute the excluded volume force
//         const std::size_t parent_idj = _point_mesh->particles()[n_id]->parent_id();
//         const Point ptj   = _point_mesh->particles()[n_id]->point();
//         const Point R_ij  = _point_mesh->pm_periodic_boundary()->point_vector(pti,ptj);
//         const std::vector<Real> Fv_ij = this->excluded_volume_force(R_ij, c1, c2);
//         for (std::size_t j=0; j<_dim; ++j) evforce[j] += Fv_ij[j];
// //        printf("--->TEST: Fv(%lu,%lu) = (%f, %f, %f);\n",
// //               p_id, n_id, Fv_ij[0], Fv_ij[1], Fv_ij[2]);
        
//         // if the two points are from different parents, compute the friction force.
//         if(parent_id != parent_idj)
//         {
//           const std::vector<Real> velj = _pm_system->point_velocity(v_beads,n_id);
//           const std::vector<Real> Ff_ij = this->friction_force(pti,ptj,veli,velj,Fv_ij,Hf,dmin);
//           for (std::size_t j=0; j<_dim; ++j) evforce[j] += Ff_ij[j];
// //          printf("--->TEST: Ff(%lu,%lu) = (%f, %f, %f);\n",
// //                 p_id, n_id, Ff_ij[0], Ff_ij[1],Ff_ij[2]);
//         } // end if
        
//       } // end if
//     } // end for i-loop
    
    
//     // add the resultant force: excluded volume + friction
//     _point_mesh->particles()[p_id]->add_particle_force(evforce);
//   } // end for i-loop
  
  
//   STOP_LOG ("modify_force_field()", "ForceField");
// }





// ======================================================================
void ForceField::reinit_point_force_zeros()
{
  START_LOG ("reinit_point_force_zeros()", "ForceField");
  
  for(std::size_t i=0; i<_num_points; ++i)
  {
    _point_mesh->particles()[i]->zero_particle_force();
  }
  
  STOP_LOG ("reinit_point_force_zeros()", "ForceField");
}




// ======================================================================
void ForceField::check_wall(const std::size_t& p_id)
{
  START_LOG ("check_wall()", "ForceField");
 
  //  retrieve bead position
  Point& pt0     = _point_mesh->particles()[p_id]->point();

  std::string wall_type = _pm_system->get_equation_systems().parameters.get<std::string>("wall_type");

  if(_wall_type == "slit"){
	const PMPeriodicBoundary* pbc = _point_mesh->pm_periodic_boundary();
	std::vector<int> & count = _point_mesh->particles()[p_id]->counter();
	
	// check all three directions of the box
	for(std::size_t i=0; i<_dim; ++i)
	{
	  if(_periodic[i]) // periodic boundary
	  {
	    if( pt0(i) <  _box_min(i) ) {pt0(i) += _box_len(i); count[i] -= 1;}
	    if( pt0(i) >= _box_max(i) ) {pt0(i) -= _box_len(i); count[i] += 1;}
	  }
	  else   // non-periodic boundary (inpenetrable wall)
	  {
	    if( pt0(i) <  _box_min(i) ){
              pt0(i) = 2.*_box_min(i) - pt0(i);  _out_domain_counter += 1;
              std::cout << "*** Warning: " << p_id << "-th bead is out of domain. Out of domain occurs "
                        << _out_domain_counter << " times in total." << std::endl;
            }
	    if( pt0(i) >  _box_max(i) ){
              pt0(i) = 2.*_box_max(i) - pt0(i);  _out_domain_counter += 1;
              std::cout << "*** Warning: " << p_id << "-th bead is out of domain. Out of domain occurs "
                        << _out_domain_counter << " times in total." << std::endl;
            }
	  }
	} // end for i-loop
  }
  else if(_wall_type == "sphere"){
	Real  pt0_norm = pt0.norm();
	Point pt0_unit = pt0.unit();
        Real cavity_radius = _wall_params[0];
        // Move particle back into spherical cavity
        if( pt0_norm > cavity_radius ){ 
          pt0 = pt0_unit * (2.*cavity_radius) - pt0;  _out_domain_counter += 1;
          std::cout << "*** Warning: " << p_id << "-th bead ("<< pt0(0) << "," << pt0(1) << "," <<pt0(2) <<") is out of domain. Out of domain occurs "
                    << _out_domain_counter << " times in total." << std::endl;
        }
  }
  else if(_wall_type == "not_defined"){
	std::cout << "*** Error: wall_type '"<< _wall_type << "' !!! (in check_wall)" << std::endl;
	libmesh_error();
  }
  else{
	std::cout << "*** Error: wall_type '"<< _wall_type
		  << "' is defined, but name is wrong or not supported!!! (in check_wall)" << std::endl;
	libmesh_error();
  }

  STOP_LOG ("check_wall()", "ForceField");
}



// ======================================================================
void ForceField::check_walls()
{
  START_LOG ("check_walls()", "ForceField");
  
  for(std::size_t i=0; i<_num_points; ++i)
  {
    const std::size_t p_id = _point_mesh->particles()[i]->id();
    this->check_wall(p_id);
  }
  
  STOP_LOG ("check_walls()", "ForceField");
}


// // ======================================================================
// void ForceField::compute_pp_wormLike_spring(const std::size_t  p_id,
//   std::vector<Real>& pforce ) const
//   {
//   START_LOG ("compute_wormLike_spring_force(p_id)", "ForceField");

//   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    c1 = a/(2bk), and c2 = q0/a = Nks*bk/a
//    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//   const Real c1  = _bead_r/(2.0*_bk);
//   const Real c2  = _Nks*_bk/_bead_r;


//   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    Wormlike spring model: force between connected beads
//    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//   pforce.resize(_dim,0.0);
//   const Point pti   = _point_mesh->particles()[p_id]->point();
//   if (p_id==0)  // the head bead with only one connection
//   {
//     const Point ptj   = _point_mesh->particles()[p_id+1]->point();
//     const Point R_ij  = _point_mesh->pm_periodic_boundary()->point_vector(pti,ptj);
//     const std::vector<Real> F_ij = this->spring_force_wls(R_ij,c1,c2);
//     for (std::size_t j=0; j<_dim; ++j) pforce[j] = F_ij[j];
//   }
//   else if ( p_id==(_num_points-1) ) // the end bead with one connection
//   {
//     const Point ptj   = _point_mesh->particles()[p_id-1]->point();
//     const Point R_ij  = _point_mesh->pm_periodic_boundary()->point_vector(pti,ptj);
//     const std::vector<Real> F_ij = this->spring_force_wls(R_ij,c1,c2);
//     for (std::size_t j=0; j<_dim; ++j) pforce[j] = F_ij[j];
//   }
//   else          // normal bead with two connected neighboring springs
//   {
//     for (std::size_t i=0; i<2; ++i)
//     {
//       // i=0, j=2*i-1 = -1;   i=1,j=2*i-1 = +1
//       const Point ptj   = _point_mesh->particles()[p_id+2*i-1]->point();
//       const Point R_ij  = _point_mesh->pm_periodic_boundary()->point_vector(pti,ptj);
//       const std::vector<Real> F_ij = this->spring_force_wls(R_ij,c1,c2);
//       for (std::size_t j=0; j<_dim; ++j) pforce[j] += F_ij[j];
//     } // end for i-loop
//   } // end if-else

//   STOP_LOG ("compute_wormLike_spring_force(p_id)", "ForceField");
// }

// // ======================================================================
// void ForceField::compute_pp_fene_spring(const std::size_t  p_id,
//                                            std::vector<Real>& pforce ) const
// {
//   START_LOG ("compute_fene_spring_force()", "ForceField");
  
//   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    In normalize formulation: c1 = 3a/(bk), and c2 = q0/a = Nks*bk/a
//    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//   const Real c1  =  3.0*_bead_r/_bk;
//   const Real c2  = _Nks*_bk/_bead_r;  // Ls: maximum spring length
  
  
//   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    FENE spring model: force between connected beads
//    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//   pforce.resize(_dim,0.0);
//   const Point pti   = _point_mesh->particles()[p_id]->point();
//   if (p_id==0)  // the head bead with only one connection
//   {
//     const Point ptj   = _point_mesh->particles()[p_id+1]->point();
//     const Point R_ij  = _point_mesh->pm_periodic_boundary()->point_vector(pti,ptj);
//     const std::vector<Real> F_ij = this->spring_force_fene(R_ij,c1,c2);
//     for (std::size_t j=0; j<_dim; ++j) pforce[j] = F_ij[j];
//   }
//   else if ( p_id==(_num_points-1) ) // the end bead with one connection
//   {
//     const Point ptj   = _point_mesh->particles()[p_id-1]->point();
//     const Point R_ij  = _point_mesh->pm_periodic_boundary()->point_vector(pti,ptj);
//     const std::vector<Real> F_ij = this->spring_force_fene(R_ij,c1,c2);
//     for (std::size_t j=0; j<_dim; ++j) pforce[j] = F_ij[j];
//   }
//   else          // normal bead with two connected neighboring springs
//   {
//     for (std::size_t i=0; i<2; ++i)
//     {
//       // i=0, j=2*i-1 = -1;   i=1,j=2*i-1 = +1
//       const Point ptj   = _point_mesh->particles()[p_id+2*i-1]->point();
//       const Point R_ij  = _point_mesh->pm_periodic_boundary()->point_vector(pti,ptj);
//       const std::vector<Real> F_ij = this->spring_force_fene(R_ij,c1,c2);
//       for (std::size_t j=0; j<_dim; ++j) pforce[j] = F_ij[j];
//     } // end for i-loop
//   } // end if-else
  
  
//   STOP_LOG ("compute_fene_spring_force()", "ForceField");
// }

// // ======================================================================
// void ForceField::constant_force(std::vector<Real>& pforce ) const
// {
//   START_LOG ("constant_force()", "ForceField");
  
//   // give a constant value along i direction only
//   pforce.resize(_dim,0.0);
//   pforce[0] = 0.0098;
//   pforce[1] = 0.00;
//   pforce[2] = 0.0098;
  
//   STOP_LOG ("constant_force()", "ForceField");
// }

