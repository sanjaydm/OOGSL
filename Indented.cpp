#include "Indented.h"
#include <iomanip>

void Indented :: fdf(){
  _f = 0; // reset energy
  _df.setZero(); //reset residue to zero
  int numEle = _conn.size(); 
  Vector wt = _quad.weights();
  vector<Vector> q = _quad.nodes();
  int numQuad = q[0].size();
  Vector est(2);
  Vector kst(2);
  //Reference surface
  double rho = 1;
  double rho_s = 0; double rho_ss = 0;
  double z_s = 1; double z_ss = 0;
  double fr = 0; double fz = 0;
  double rhofr = rho*fr; double rhofz = rho*fz;
  
  Shape hermite(1, 0);
  int modelOption = 2;
  for (int e=0; e< numEle; e++){
    for (int j=0; j<numQuad; j++) {
      
      Vector ele (2);
      ele(0) = _conn[e](0); ele(1) = _conn[e](1) ;
      double u, v, u_s, v_s, u_ss, v_ss;
      // Quadrature
      double qj = q[0](j);
      double wj = wt(j);
      hermite.compute(qj);

      Vector N = hermite.getN();
      Vector DN = hermite.getDN();
      Vector D2N = hermite.getD2N();

      double slope = (_nodes[ele(1)]-_nodes[ele(0)])/2;
      // Compute u, v, u', v', u'', v'' using shape functions
      // For each node u,v,u' and v' are stored in this order
      double T0 =  _x(2*ele(0));
      double T1 = _x(2*ele(0)+1);
      //double F0 = _x(3*ele(0)+2);
      
      double T2 = _x(2*ele(1));
      double T3 = _x(2*ele(1)+1);
      //double F1 = _x(3*ele(1)+2);

      double S = _nodes[ele(0)]*N(0) + slope*N(1) + _nodes[ele(1)]*N(2) + slope*N(3); 
      double Jac = _nodes[ele(0)]*DN(0) + slope*DN(1) + _nodes[ele(1)]*DN(2) + slope*DN(3); 
      double Jac_xi = _nodes[ele(0)]*D2N(0) + slope*D2N(1) + _nodes[ele(1)]*D2N(2) + slope*D2N(3); // 0
      double Jac2 = (Jac * Jac);
      double Jac1 = Jac;
      double Jac_xiJac2 = (Jac_xi/Jac2);
      
      z_s = 1; //(_nodes[ele(0)]*DN(0) + slope*DN(1) + _nodes[ele(1)]*DN(2) + slope*DN(3))/Jac1;
      z_ss = 0;// (_nodes[ele(0)]*D2N(0) + slope*D2N(1) + _nodes[ele(1)]*D2N(2) + slope*D2N(3))/Jac2;

      //cout << "\033[32m z_ss = " << z_ss << "\033[0m\n";

      
      auto L0 = [](double xi) {return (1-xi)/2.0;};
      auto L1 = [](double xi) {return (1+xi)/2.0;};
      

      
      double DL0 = -0.5;
      double DL1 = 0.5;

      double t = T0*N(0)+T1*N(1) + T2*N(2) + T3*N(3);

      double t_s = (T0*DN(0)+T1*DN(1) + T2*DN(2) + T3*DN(3))/Jac1;
      double t_ss = (T0*D2N(0)+T1*D2N(1) + T2*D2N(2) + T3*D2N(3))/Jac2- t_s*Jac_xiJac2;
      
      //double fn = F0; //F0*L0(qj) + F1*L1(qj);
      
      /*
      fz = fn*(-sin(t) + _mu*cos(t));
      fr = fn*(cos(t) + _mu*sin(t));
      rhofr = rho*fr;
      rhofz = rho*fz;

      u = (_d - _r*cos(t)-rho);
      v = _r*sin(t)-S;

      u_s = _r*sin(t)*t_s;
      v_s = _r*cos(t)*t_s-1;

      u_ss = _r*cos(t)*t_s*t_s + _r*sin(t)*t_ss;
      v_ss = -_r*sin(t)*t_s*t_s + _r*cos(t)*t_ss;
      */


      u = (_d - _r - rho);
      v = t;

      u_s = 0;
      v_s = t_s;

      u_ss = 0;
      v_ss = t_ss;

      fz = -(_C+_D)*v_ss;
      fr = 0; fz/_mu;
      rhofr = rho*fr;
      rhofz = rho*fz;


      double Jacwj = Jac * wj;
      // Strain measures for cylindrical indenter
      double es_lin = u_s*rho_s+z_s*v_s; //u' rho' + z' v' = v'
      double et_lin = u/rho; //u/rho = u
      double phi = -z_s*u_s+rho_s*v_s; //-z' u' + rho' v' = 0
      double phi_s = -z_ss*u_s -z_s*u_ss + rho_ss*v_s + rho_s*v_ss; //=0
      double k_s = phi_s; //phi' = 0
      double k_t = v_s/rho; //v'/rho = v'

      est(0) = es_lin + 0.5*phi*phi; //v'
      est(1) = et_lin; //u
      kst(0) = k_s; //=0
      kst(1) = k_t; //v'

      // Call material law
      NeoHookean(est, kst);

      // If Energy flag is on
      if (_fFlag){
	_f +=  (_lclEnergyDensity*rho - rho*fr*u - rho*fz*v) * Jacwj;
	//_f +=  (- rho*fr*u - rho*fz*v) * Jacwj;
      }
      if (_dfFlag){

	double cos_gamma = rho_s;
	double sin_gamma = z_s;
	double Tr = cos_gamma - phi*sin_gamma;
	double Tz = sin_gamma + phi*cos_gamma;

	// double delta_es = Tr*(0) + Tz*(0); // Tr (du)' + Tz (dv)'
	// double delta_et = 0/rho; //du/rho
	// double delta_ks = -sin_gamma * (0) + cos_gamma*(0); // -sin g*(du)'+cos(g)*(dv)'
	// double delta_kt = 0/rho; //(dv)'/rho

	double Ns = _lclResidue(0); double Nt = _lclResidue(1);
	double Ms = _lclResidue(2); double Mt = _lclResidue(3);
	double rhoNs = rho*Ns;

	// Variations
	double dt0, dt1, dt2, dt3;
	double dfn0, dfn1;

	double dt0_s, dt1_s, dt2_s, dt3_s;
	double dt0_ss, dt1_ss, dt2_ss, dt3_ss;
	
	double du0, du1, du2, du3;
	double dv0, dv1, dv2, dv3;

	double du0_s, du1_s, du2_s, du3_s;
	double dv0_s, dv1_s, dv2_s, dv3_s;

	double du0_ss, du1_ss, du2_ss, du3_ss;
	double dv0_ss, dv1_ss, dv2_ss, dv3_ss;

	double dphi0U_s, dphi1U_s, dphi2U_s, dphi3U_s;
	double dphi0V_s, dphi1V_s, dphi2V_s, dphi3V_s;

	double dfr0, dfr1, dfr2, dfr3;
	double dfz0, dfz1, dfz2, dfz3;

	//if (e!=0 || e!=numEle-1) {
	dt0 = N(0); dt1 = N(1); dt2 = N(2); dt3 = N(3);
	dfn0 = 1;//L0(qj); //dfn1 = L1(qj);

	dt0_s = DN(0)/Jac1; dt1_s = DN(1)/Jac1; dt2_s = DN(2)/Jac1; dt3_s = DN(3)/Jac1;

	dt0_ss = D2N(0)/Jac2 - dt0_s*Jac_xiJac2;
	dt1_ss = D2N(1)/Jac2 - dt1_s*Jac_xiJac2;
	dt2_ss = D2N(2)/Jac2 - dt2_s*Jac_xiJac2;
	dt3_ss = D2N(3)/Jac2 - dt3_s*Jac_xiJac2;


	/*
	// Dependence on du , dv s
	du0 = _r*sin(t)*dt0; du1 = _r*sin(t)*dt1; du2 = _r*sin(t)*dt2; du3 = _r*sin(t)*dt3;
	dv0 = _r*cos(t)*dt0; dv1 = _r*cos(t)*dt1; dv2 = _r*cos(t)*dt2; dv3 = _r*cos(t)*dt3;

	du0_s = _r*cos(t)*t_s*dt0 + _r*sin(t)*dt0_s;
	du1_s = _r*cos(t)*t_s*dt1 + _r*sin(t)*dt1_s;
	du2_s = _r*cos(t)*t_s*dt2 + _r*sin(t)*dt2_s;
	du3_s = _r*cos(t)*t_s*dt3 + _r*sin(t)*dt3_s;

	dv0_s = -_r*sin(t)*t_s*dt0 + _r*cos(t)*dt0_s;
	dv1_s = -_r*sin(t)*t_s*dt1 + _r*cos(t)*dt1_s;
	dv2_s = -_r*sin(t)*t_s*dt2 + _r*cos(t)*dt2_s;
	dv3_s = -_r*sin(t)*t_s*dt3 + _r*cos(t)*dt3_s;

	du0_ss =  -_r*sin(t)*t_s*t_s*dt0 + _r*cos(t)*t_ss*dt0 + 2*_r*cos(t)*t_s*dt0_s + _r*sin(t)*dt0_ss;
	du1_ss =  -_r*sin(t)*t_s*t_s*dt1 + _r*cos(t)*t_ss*dt1 + 2*_r*cos(t)*t_s*dt1_s + _r*sin(t)*dt1_ss;
	du2_ss =  -_r*sin(t)*t_s*t_s*dt2 + _r*cos(t)*t_ss*dt2 + 2*_r*cos(t)*t_s*dt2_s + _r*sin(t)*dt2_ss;
	du3_ss =  -_r*sin(t)*t_s*t_s*dt3 + _r*cos(t)*t_ss*dt3 + 2*_r*cos(t)*t_s*dt3_s + _r*sin(t)*dt3_ss;
	  
	dv0_ss = -_r*cos(t)*t_s*t_s*dt0 -_r*sin(t)*t_ss*dt0 -2*_r*sin(t)*t_s*dt0_s + _r*cos(t)*dt0_ss;
	dv1_ss = -_r*cos(t)*t_s*t_s*dt1 -_r*sin(t)*t_ss*dt1 -2*_r*sin(t)*t_s*dt1_s + _r*cos(t)*dt1_ss;
	dv2_ss = -_r*cos(t)*t_s*t_s*dt2 -_r*sin(t)*t_ss*dt2 -2*_r*sin(t)*t_s*dt2_s + _r*cos(t)*dt2_ss;
	dv3_ss = -_r*cos(t)*t_s*t_s*dt3 -_r*sin(t)*t_ss*dt3 -2*_r*sin(t)*t_s*dt3_s + _r*cos(t)*dt3_ss;

	dfr0 = fn*(-sin(t) + _mu*cos(t))*dt0;
	dfr1 = fn*(-sin(t) + _mu*cos(t))*dt1;
	dfr2 = fn*(-sin(t) + _mu*cos(t))*dt2;
	dfr3 = fn*(-sin(t) + _mu*cos(t))*dt3;

	dfz0 = -fn*(cos(t) + _mu*sin(t))*dt0;
	dfz1 = -fn*(cos(t) + _mu*sin(t))*dt1;
	dfz2 = -fn*(cos(t) + _mu*sin(t))*dt2;
	dfz3 = -fn*(cos(t) + _mu*sin(t))*dt3;

	*/

	
	// Dependence on du , dv s
	du0 = 0;   du1 = 0;   du2 = 0;   du3 = 0;
	dv0 = dt0; dv1 = dt1; dv2 = dt2; dv3 = dt3;

	du0_s = 0;
	du1_s = 0;
	du2_s = 0;
	du3_s = 0;

	dv0_s = dt0_s;
	dv1_s = dt1_s;
	dv2_s = dt2_s;
	dv3_s = dt3_s;

	du0_ss =  0;
	du1_ss =  0;
	du2_ss =  0;
	du3_ss =  0;
	  
	dv0_ss = dt0_ss;
	dv1_ss = dt1_ss;
	dv2_ss = dt2_ss;
	dv3_ss = dt3_ss;


	dfz0 = 0; //-(_C + _D)*dv0_ss;
	dfz1 = 0;//-(_C + _D)*dv1_ss;
	dfz2 = 0;//-(_C + _D)*dv2_ss;
	dfz3 = 0;//-(_C + _D)*dv3_ss;

	dfr0 = 0;//dfz0/_mu;
	dfr1 = 0;//dfz1/_mu;
	dfr2 = 0;//dfz2/_mu;
	dfr3 = 0;//dfz3/_mu;


	// ------------ END -------------------
		
	dphi0U_s = -z_ss* du0_s - z_s* du0_ss; 
	dphi1U_s = -z_ss* du1_s - z_s* du1_ss;
	dphi2U_s = -z_ss* du2_s - z_s* du2_ss;
	dphi3U_s = -z_ss* du3_s - z_s* du3_ss; 

	dphi0V_s = rho_ss* dv0_s + rho_s* dv0_ss; 
	dphi1V_s = rho_ss* dv1_s + rho_s* dv1_ss;
	dphi2V_s = rho_ss* dv2_s + rho_s* dv2_ss;
	dphi3V_s = rho_ss* dv3_s + rho_s* dv3_ss; 

	

	
	_df(2*ele(0)     ) +=  ((rhoNs * Tr)* du0_s  +
				rho*Ms*dphi0U_s + 
				(Nt - rhofr)* du0 - rho*u*dfr0) *Jacwj +
	                        ((rhoNs * Tz)* dv0_s +
				  rho*Ms*dphi0V_s + Mt*dv0_s
				  - rhofz* dv0 - rho*v*dfz0) *Jacwj ;
	
	_df(2*ele(0) + 1 ) +=  ((rhoNs * Tr)* du1_s  +
				rho*Ms*dphi1U_s +
				(Nt - rhofr)* du1 -rho*u*dfr1) *Jacwj + 	
	                        ((rhoNs * Tz)* dv1_s +
				  rho*Ms*dphi1V_s + Mt*dv1_s
				  - rhofz* dv1 - rho*v*dfz1) *Jacwj ;

	/*
	_df(3*ele(0) + 2   ) +=   -((cos(t) + _mu*sin(t))*u + 
				   (-sin(t) + _mu*cos(t))*v )*dfn0*rho*Jacwj ;
	*/
	/*
	_df(3*ele(0) + 2   ) +=   -(u + 
				   _mu*v )*dfn0*rho*Jacwj ;
	*/

	_df(2*ele(1)     ) +=  ((rhoNs * Tr)* du2_s  +
				rho*Ms*dphi2U_s + 
				(Nt - rhofr)* du2 - rho*u*dfr2) *Jacwj +
	                        ((rhoNs * Tz)* dv2_s +
				  rho*Ms*dphi2V_s + Mt*dv2_s
				  - rhofz* dv2 - rho*v*dfz2) *Jacwj ;
	
	_df(2*ele(1) + 1 ) +=  ((rhoNs * Tr)* du3_s  +
				rho*Ms*dphi3U_s +
				(Nt - rhofr)* du3 - rho*u*dfr3) *Jacwj +
	                        ((rhoNs * Tz)* dv3_s +
				  rho*Ms*dphi3V_s + Mt*dv3_s
				  - rhofz* dv3 - rho*v*dfz3) *Jacwj ;

	/*
	_df(3*ele(1) + 2   ) +=  -((cos(t) + _mu*sin(t))*u +
				   (-sin(t) + _mu*cos(t))*v )*dfn1*rho*Jacwj ;

	*/
	/*----
	_df(3*ele(1) + 2   ) +=   -(u + 
				   _mu*v )*dfn1*rho*Jacwj ;
	-----*/

	/*
	_df(3*ele(0)     ) +=  (- rhofr* du0 - rho*u*dfr0
				  - rhofz* dv0 - rho*v*dfz0) *Jacwj ;
	
	_df(3*ele(0) + 1 ) +=  ( - rhofr* du1 -rho*u*dfr1
				  - rhofz* dv1 - rho*v*dfz1) *Jacwj ;

	_df(3*ele(0) + 2   ) +=   -((cos(t) + _mu*sin(t))*u + 
				   (-sin(t) + _mu*cos(t))*v )*dfn0*rho*Jacwj ;


	_df(3*ele(1)     ) +=  (- rhofr* du2 - rho*u*dfr2
				  - rhofz* dv2 - rho*v*dfz2) *Jacwj ;
	
	_df(3*ele(1) + 1 ) +=  (- rhofr* du3 - rho*u*dfr3
				  - rhofz* dv3 - rho*v*dfz3) *Jacwj ;
	
	_df(3*ele(1) + 2   ) +=  -((cos(t) + _mu*sin(t))*u +
				   (-sin(t) + _mu*cos(t))*v )*dfn1*rho*Jacwj ;
	
	*/
      }

    } //end: for quadrature
  } // end: for element

}
void Indented :: NeoHookean (Vector& est, Vector& kst){
  // Material law.
  // Inputs: strain measures
  // Output stored in _lclEnergy, _lclResidue
  double es = est(0); double et = est(1);
  double ks = kst(0); double kt = kst(1);
  if (_fFlag){
    _lclEnergyDensity = 0.5*( _C* (es*es + et*et + 2*_nu*es*et) + 
			      _D* (ks*ks + kt*kt + 2*_nu*ks*kt) );
  }
  if (_dfFlag) {
    
    _lclResidue(0) = _C*(es + _nu*et); //Ns = C (v' + nu u)
    _lclResidue(1) = _C*(et + _nu*es); //Nt = C(u + nu v')
    _lclResidue(2) = _D*(ks + _nu*kt); //Ms = D * nu * v'
    _lclResidue(3) = _D*(kt + _nu*ks); //Mt = D * v'
    
  }
}

void Indented::printU(){
  int NDOF = _x.size();
  cout <<"U=[ ";
  for (int i=0; i < NDOF; i=i+4){
    cout << _x(i) << ", ";
  }
  cout << "];\n";
}

void Indented::writeMesh(string filename){
  ofstream myfile;
  myfile.open (filename.c_str());
  cout << "Writing mesh to a file.\n";
  for (auto it = _nodes.begin(); it != _nodes.end(); ++it) {
    myfile << *it << endl;
}
  myfile.close();
}

// void Indented::writeSolution(string filename){
//   ofstream myfile;
//   myfile.open (filename.c_str());
//   cout << "Writing solution to a file.\n";
//   myfile << setprecision(15);
//   for (auto i = 1; i <= _x.size(); i++) {
//     myfile << _x(i-1);
//     if (i % 2 == 0)
//       myfile << "\n";
//     else
//       myfile << "\t";
// }
//   myfile.close();
// }

void Indented::writeSolution(string filename){
  ofstream myfile;
  _x.print();
  myfile.open (filename.c_str());
  cout << "Writing solution to a file.\n";
  myfile << setprecision(15);
  myfile << "x = [";
  for (auto i = 0; i < _x.size(); i++) {
    if (i % 2 == 0){
      myfile << _nodes[i/2] << ",";
      //myfile <<_d - _r*cos(_x(i)) << ",";
    }
  }
  myfile << "]\n";
  myfile << "y = [";
  for (auto i = 0; i < _x.size(); i++) {
    if (i%2 == 0){
      myfile << _x(i)<< ",";
      //myfile << _r*sin(_x(i)) << ",";
    }
  }
  myfile << "]\n";
  myfile << "from matplotlib import pyplot as plt\n";
  myfile << "from numpy import *\n";
  myfile << "x = matrix(x)\n";
  myfile << "y = matrix(y)\n";
  myfile << "xmin = -matrix(x)\n";
  myfile << "y = matrix(y)\n";

  myfile << "plt.plot(x.T, y.T,'*-')\n";
  //myfile << "plt.plot(xmin.T, y.T,'-')\n";
  myfile << "plt.axis('tight')\n";
  //myfile << "plt.axis('equal')\n";
  myfile << "plt.show()\n";
  myfile.close();
}
