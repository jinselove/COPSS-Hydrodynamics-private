#include "fix_rigid_analytic_polarization.h"

FixRigidAnalyticPolarization::FixRigidAnalyticPolarization(PMLinearImplicitSystem& pm_sys_)
:FixRigid(pm_sys_)
{
  force_type = "analytic_polarization";
  this -> initParams();
}

void FixRigidAnalyticPolarization::initParams()
{
  force_params = pm_system->get_equation_systems().parameters.get<std::vector<Real>> ("analytic_polarization");
  if(force_params.size()!=2){
      std::cout << std::endl << "********************Error message********************" << std::endl
                << "---------------> The force type 'analytic_polarization' requires 2 parameter (epsilon_out, order),"
                << "epsilon_out: relative dielectric permittivity of medium with respect to vacuum;" 
                <<"order: if order = 0, only compute columbic force, if order = 1, compute both columbic force and polarization using first order image method"
                << std::endl
                << "****************************************" << std::endl;
    libmesh_error();    
  }
  eouter = force_params[0];
  order = force_params[1];
  radius = rigid_particles[0]->radius();
  einner = rigid_particles[0]->epsilon_in();
  _e = (1.0 - einner/eouter) / (1.0 + einner/eouter);
  _g = 1.0 / (1.0 + einner/eouter);
}

void FixRigidAnalyticPolarization::print_fix()
{
  START_LOG("FixRigidAnalyticPolarization::print_fix()", "FixRigidAnalyticPolarization");
  std::cout <<"this is FixRigidAnalyticPolarization" << std::endl;
  STOP_LOG("FixRigidAnalyticPolarization::print_fix()", "FixRigidAnalyticPolarization");
}

void FixRigidAnalyticPolarization::compute()
{
 START_LOG("FixRigidAnalyticPolarization::compute()", "FixRigidAnalyticPolarization");
 std::vector<Point> body_force(num_rigid_particles);
 // compute body columbic force for all particles
 // std::cout<<"compute columbic force" << std::endl;
 this->compute_coulombic_force(body_force);
 // Add first order image method if specified
 if (order==1){
     // std::cout<<"compute polarization force" << std::endl;
     this->compute_first_order_polarization(body_force);
 }
 // std::cout<<"convert unit of electrostatic force" << std::endl;
 for (int i=0; i<body_force.size(); i++) body_force[i] = body_force[i] * (radius * radius * eouter);
 /*
 * notice: the calculated electrostatic force is of unit 
 * e^2/(4*pi*epsilon_0*particle_radius^2), where 
 * e = 1.60 * 10^(-19) C
 * epsilon_0 = 8.854 * 10^(-12) C^2 N^(-1) m^(-2)
 * particle radius is in the unit of bead radius, i.e., for example particle radius = 2
 * which is set in the data file.
 */

 // distribute body force to nodes
 for (int i=0; i < num_rigid_particles; i++){
  // std::cout<<"Electrostatic force on particle "<<i <<" is: "<<body_force[i](0) <<", "<<body_force[i](1) <<", "<<body_force[i](2)<<std::endl;
  rigid_particles[i]->build_nodal_force(body_force[i]);
 }
 STOP_LOG("FixRigidAnalyticPolarization::compute()", "FixRigidAnalyticPolarization");  
}

void FixRigidAnalyticPolarization::compute_coulombic_force(std::vector<Point>& body_force)
{
    START_LOG("FixRigidAnalyticPolarization::compute_coulombic_force", "FixRigidAnalyticPolarization");
    for (int i = 0; i < body_force.size(); i++){
        for (int j = i+1; j < body_force.size(); j++){
            Point rij = rigid_particles[j]->get_centroid() - rigid_particles[i]->get_centroid(); // pointing from particle i to particle j
            Real d2 = rij.norm_sq();
            Point force = rigid_particles[i]->charge() * rigid_particles[j]->charge() / d2 * rij.unit();
            body_force[i] -= force;
            body_force[j] += force;
        }
    }
    STOP_LOG("FixRigidAnalyticPolarization::compute_coulombic_force", "FixRigidAnalyticPolarization");
} 

void FixRigidAnalyticPolarization::compute_first_order_polarization(std::vector<Point>& body_force)
{
    START_LOG("FixRigidAnalyticPolarization::compute_first_order_polarization", "FixRigidAnalyticPolarization");
    for (int j=0; j<body_force.size(); j++){
        for (int i=0; i<body_force.size(); i++) if (i!=j){
            for (int k=0; k<body_force.size(); k++) if (k!=j){
                if (i==k){
                    this->force_pol(body_force, i, j, k);
                }
                else if(i<k){
                    this->force_pol(body_force, i, j, k);
                    this->force_pol(body_force, k, j, i);
                }
            }
        }
    }
    STOP_LOG("FixRigidAnalyticPolarization::compute_first_order_polarization", "FixRigidAnalyticPolarization");
}

void FixRigidAnalyticPolarization::force_pol(std::vector<Point>& body_force, int ith, int jth, int kth)
{
    START_LOG("FixRigidAnalyticPolarization::force_pol", "FixRigidAnalyticPolarization");
    // std::cout<<"debug in force pol:" <<std::endl;
    // for(int i=0; i<rigid_particles.size(); i++) {
    //     Point tmp = rigid_particles[i]->get_centroid();
    //     tmp.print();
    //     std::cout << std::endl;
    // }
    Point Ri = rigid_particles[ith]->get_centroid();
    Point Rj = rigid_particles[jth]->get_centroid();
    Point Rk = rigid_particles[kth]->get_centroid();
    // Rkj: vector points from j-th to k-th
    Point Rkj = Rk - Rj;    
    rkj = Rkj.norm();
    Rxkj = Rkj(0);
    Rykj = Rkj(1);
    Rzkj = Rkj(2);
    ukj = Rxkj / rkj;
    vkj = Rykj / rkj;
    wkj = Rzkj / rkj;
    // Rij: vector points from jth to ith particle
    Point Rij = Ri - Rj;
    rij = Rij.norm();
    Rxij = Rij(0);
    Ryij = Rij(1);
    Rzij = Rij(2);
    // Auxiliary variables
    aux1 = _e * radius / rkj;
    aux2 = radius * radius / rkj;
    // std::cout<<"_e = "<<_e << ", radius = "<<radius << ", rkj="<<rkj<<std::endl;
    // std::cout<<"aux1 = "<<aux1 <<", axu2 = "<<aux2<<std::endl;
    // tmp storage for polarization force
    Point force_pol(0.);
    
    // 1) delta_s, 1 term of equation.29 in method paper
    auxv_x_delta = Rxij - aux2 * ukj;
    auxv_y_delta = Ryij - aux2 * vkj;
    auxv_z_delta = Rzij - aux2 * wkj;
    aux3Sqrt_delta = sqrt(auxv_x_delta * auxv_x_delta + auxv_y_delta * auxv_y_delta + auxv_z_delta * auxv_z_delta);
    aux3_delta = aux3Sqrt_delta * aux3Sqrt_delta * aux3Sqrt_delta;
       
    force_pol(0) = auxv_x_delta / aux3_delta;
    force_pol(1) = auxv_y_delta / aux3_delta;
    force_pol(2) = auxv_z_delta / aux3_delta;
    
    // 2) integral term of equation.29 in method paper 
    std::vector<double> _xg(5);
    for (int ig = 0; ig < ngauss; ++ig)
    {
      _xg[ig] = (xhi - xlo) / 2.0 * _xg0[ig] + (xhi + xlo) / 2.0; 
      auxv_x_integ = Rxij - ( pow(_xg[ig] , 1.0/_g) * aux2 * ukj );
      auxv_y_integ = Ryij - ( pow(_xg[ig] , 1.0/_g) * aux2 * vkj );
      auxv_z_integ = Rzij - ( pow(_xg[ig] , 1.0/_g) * aux2 * wkj );
      aux3Sqrt_integ = sqrt(auxv_x_integ * auxv_x_integ + auxv_y_integ * auxv_y_integ + auxv_z_integ * auxv_z_integ);
      aux3_integ = aux3Sqrt_integ * aux3Sqrt_integ * aux3Sqrt_integ;
      force_pol(0) += (xhi-xlo)/2 * (- auxv_x_integ / aux3_integ * _wg0[ig]); 
      force_pol(1) += (xhi-xlo)/2 * (- auxv_y_integ /aux3_integ * _wg0[ig]);
      force_pol(2) += (xhi-xlo)/2 * (- auxv_z_integ /aux3_integ * _wg0[ig]);
    }
    
    // final polarization force for term {i;j;k}
    Real zi = rigid_particles[ith]->charge();
    Real zk = rigid_particles[kth]->charge();
    force_pol(0) *= 1.0 / 2.0 * aux1 * zi * zk;
    force_pol(1) *= 1.0 / 2.0 * aux1 * zi * zk;
    force_pol(2) *= 1.0 / 2.0 * aux1 * zi * zk;
    // update total force
    // std::cout<<"force pol for (" <<ith <<","<<jth <<","<<kth<<std::endl;
    // force_pol.print();
    // std::cout<<std::endl;
    body_force[ith] += 2 * force_pol;
    body_force[jth] -= 2 * force_pol;
    STOP_LOG("FixRigidAnalyticPolarization::force_pol", "FixRigidAnalyticPolarization");
}
