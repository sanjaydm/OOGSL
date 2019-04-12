#include "VirusShell.h"
#include <iomanip>

void VirusShell :: fdf(){
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
  // // //Boundary conditions
  // _x(0) = 0; //u(0) = 0
  // _x(1) = 0; //u'(0) = 0
  // _x(2) = 0; //v(0) = 0
 
  // _x(4*numEle) = 0;   //u(L) = 0;
  // _x(4*numEle+1) = 0; //u'(L) = 0;
  // _x(4*numEle+2) = 0; //v(L) = 0
  int modelOption = 2;
  for (int e=0; e< numEle; e++){
    for (int j=0; j<numQuad; j++) {
      if (e == 3) {
      	fr = 100;
      }
      else
      	fr=0;

      rhofr = rho*fr;
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
      double U0 =  _x(4*ele(0));
      double U1 = _x(4*ele(0)+1);
      double V0 =  _x(4*ele(0)+2);
      double V1 = _x(4*ele(0)+3);

      double U2 =  _x(4*ele(1));
      double U3 = _x(4*ele(1)+1);
      double V2 =  _x(4*ele(1)+2);
      double V3 = _x(4*ele(1)+3);

      // cout << "\033[31m " << U0 << " " << U1 << " " << U2 << " " << U3 << "\033[0m\n";
      // cout << "\033[32m " << V0 << " " << V1 << " " << V2 << " " << V3 << "\033[0m\n";
      
      double Jac = _nodes[ele(0)]*DN(0) + slope*DN(1) + _nodes[ele(1)]*DN(2) + slope*DN(3); 
      double Jac_xi = _nodes[ele(0)]*D2N(0) + slope*D2N(1) + _nodes[ele(1)]*D2N(2) + slope*D2N(3); // 0
      double Jac2 = (Jac * Jac);
      double Jac1 = Jac;
      double Jac_xiJac2 = (Jac_xi/Jac2);
      // U1 = U1/slope;
      // U3 = U3/slope;
      // V1 = V1/slope;
      // V3 = V3/slope;
      
      z_s = 1; //(_nodes[ele(0)]*DN(0) + slope*DN(1) + _nodes[ele(1)]*DN(2) + slope*DN(3))/Jac1;
      z_ss = 0;// (_nodes[ele(0)]*D2N(0) + slope*D2N(1) + _nodes[ele(1)]*D2N(2) + slope*D2N(3))/Jac2;

      //cout << "\033[32m z_ss = " << z_ss << "\033[0m\n";
      /*
      N(1) *= slope;
      N(3) *= slope;
      DN(1) *= slope;
      DN(3) *= slope;
      D2N(1) *= slope;
      D2N(3) *= slope;
      */

      u = U0*N(0)+U1*N(1) + U2*N(2) + U3*N(3);
      v = V0*N(0)+V1*N(1) + V2*N(2) + V3*N(3);

      u_s = (U0*DN(0)+U1*DN(1) + U2*DN(2) + U3*DN(3))/Jac1;
      v_s = (V0*DN(0)+V1*DN(1) + V2*DN(2) + V3*DN(3))/Jac1;

      u_ss = (U0*D2N(0)+U1*D2N(1) + U2*D2N(2) + U3*D2N(3))/Jac2- u_s*Jac_xiJac2;
      v_ss = (V0*D2N(0)+V1*D2N(1) + V2*D2N(2) + V3*D2N(3))/Jac2 - v_s*Jac_xiJac2;
      
      //cout << "U0 = " << U0 << ", U1 = " << U1 << ", U2 = " << U2 << ", U3 = " << U3 << endl;
      //cout << "u = " << u << ", u_s = " << u_s << ", u_ss = " << u_ss << endl;
      //cout << Jac << endl;
      //cout << 0.05*D2N(0)+.01*D2N(1) + 0.07*D2N(2) + .01*D2N(3) << endl;
      //cout << D2N(0) << " " << D2N(1) << " " << D2N(2) << " " << D2N(3) << endl;

      double Jacwj = Jac * wj;
      // Strain measures
      double es_lin = u_s*rho_s+z_s*v_s; //u' rho' + z' v'
      double et_lin = u/rho; //u/rho
      double phi = -z_s*u_s+rho_s*v_s; //-z' u' + rho' v'
      double phi_s = -z_ss*u_s -z_s*u_ss + rho_ss*v_s + rho_s*v_ss;
      double k_s = phi_s; //phi'
      double k_t = v_s/rho; //v'/rho

      est(0) = es_lin + 0.5*phi*phi;
      est(1) = et_lin;
      kst(0) = k_s;
      kst(1) = k_t;

      // Call material law
      NeoHookean(est, kst);
      // If Energy flag is on
      if (_fFlag){
	switch(modelOption){
	case 0:
	  _f += (v_s*v_s -fr*v)* Jacwj;
	  break;
	case 1:
	  _f += (u_ss*u_ss -fr*u)* Jacwj;
	  break;
	case 2:
	  _f +=  (_lclEnergyDensity*rho - rhofr*u - rhofz*v) * Jacwj;
	  break;
	}
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
	double du0, du1, du2, du3;
	double dv0, dv1, dv2, dv3;

	double du0_s, du1_s, du2_s, du3_s;
	double dv0_s, dv1_s, dv2_s, dv3_s;

	double du0_ss, du1_ss, du2_ss, du3_ss;
	double dv0_ss, dv1_ss, dv2_ss, dv3_ss;

	double dphi0U_s, dphi1U_s, dphi2U_s, dphi3U_s;
	double dphi0V_s, dphi1V_s, dphi2V_s, dphi3V_s;

	//if (e!=0 || e!=numEle-1) {
	  du0 = N(0); du1 = N(1); du2 = N(2); du3 = N(3);
	  dv0 = N(0); dv1 = N(1); dv2 = N(2); dv3 = N(3);
	  
	  du0_s = DN(0)/Jac1; du1_s = DN(1)/Jac1; du2_s = DN(2)/Jac1; du3_s = DN(3)/Jac1;
	  dv0_s = DN(0)/Jac1; dv1_s = DN(1)/Jac1; dv2_s = DN(2)/Jac1; dv3_s = DN(3)/Jac1;

	  du0_ss = D2N(0)/Jac2 - du0_s*Jac_xiJac2;
	  du1_ss = D2N(1)/Jac2 - du1_s*Jac_xiJac2;
	  du2_ss = D2N(2)/Jac2 - du2_s*Jac_xiJac2;
	  du3_ss = D2N(3)/Jac2 - du3_s*Jac_xiJac2;
	  
	  dv0_ss = D2N(0)/Jac2 - dv0_s*Jac_xiJac2;
	  dv1_ss = D2N(1)/Jac2 - dv1_s*Jac_xiJac2;
	  dv2_ss = D2N(2)/Jac2 - dv2_s*Jac_xiJac2;
	  dv3_ss = D2N(3)/Jac2 - dv3_s*Jac_xiJac2;

	  //}
	// else if (e==0){
	//   du0 = 0; du1 = 0; du2 = N(2); du3 = N(3);
	//   dv0 = 0; dv1 = N(1); dv2 = N(2); dv3 = N(3);

	//   du0_s = 0; du1_s = 0; du2_s = DN(2); du3_s = DN(3);
	//   dv0_s = 0; dv1_s = DN(1); dv2_s = DN(2); dv3_s = DN(3);

	//   du0_ss = 0; du1_ss = 0; du2_ss = D2N(2); du3_ss = D2N(3);
	//   dv0_ss = 0; dv1_ss = D2N(1); dv2_ss = D2N(2); dv3_ss = D2N(3);
	// }
	// else if (e==numEle-1){
	//   du0 = N(0); du1 = N(1); du2 = 0; du3 = 0;
	//   dv0 = N(0); dv1 = N(1); dv2 = 0; dv3 = N(3);
	//   du0_s = DN(0); du1_s = DN(1); du2_s = 0; du3_s = 0;
	//   dv0_s = DN(0); dv1_s = DN(1); dv2_s = 0; dv3_s = DN(3);
	//   du0_ss = D2N(0); du1_ss = D2N(1); du2_ss = 0; du3_ss = 0;
	//   dv0_ss = D2N(0); dv1_ss = D2N(1); dv2_ss = 0; dv3_ss = D2N(3);
	// }
		
	dphi0U_s = -z_ss* du0_s - z_s* du0_ss; 
	dphi1U_s = -z_ss* du1_s - z_s* du1_ss;
	dphi2U_s = -z_ss* du2_s - z_s* du2_ss;
	dphi3U_s = -z_ss* du3_s - z_s* du3_ss; 

	dphi0V_s = rho_ss* dv0_s + rho_s* dv0_ss; 
	dphi1V_s = rho_ss* dv1_s + rho_s* dv1_ss;
	dphi2V_s = rho_ss* dv2_s + rho_s* dv2_ss;
	dphi3V_s = rho_ss* dv3_s + rho_s* dv3_ss; 

	switch (modelOption) {
	case 0: 
	  _df(4*ele(0)     ) +=  0;
	  _df(4*ele(0) + 1 ) +=  0;
	
	  _df(4*ele(0) + 2   ) += (  2*v_s*dv0_s - fr * dv0 ) *Jacwj ;	  
	  _df(4*ele(0) + 3   ) += (  2*v_s*dv1_s- fr * dv1 ) *Jacwj ;	
	
	  _df(4*ele(1)     ) +=  0;
	  _df(4*ele(1) + 1 ) +=  0;
	
	  _df(4*ele(1) + 2   ) += (  2*v_s*dv2_s - fr*dv2  ) *Jacwj ;		
	  _df(4*ele(1) + 3   ) += (  2*v_s*dv3_s- fr*dv3 ) *Jacwj ;
	  break;

	case 1:
	    _df(4*ele(0)     ) +=  (  2*u_ss*du0_ss - fr * du0 ) *Jacwj;
	    _df(4*ele(0) + 1 ) +=  (  2*u_ss*du1_ss- fr * du1 ) *Jacwj ;
	
	    _df(4*ele(0) + 2   ) += 0;	  
	    _df(4*ele(0) + 3   ) += 0 ;	
	
	    _df(4*ele(1)     ) +=  (  2*u_ss*du2_ss - fr*du2  ) *Jacwj ;
	    _df(4*ele(1) + 1 ) +=  (  2*u_ss*du3_ss- fr*du3 ) *Jacwj ;
	
	    _df(4*ele(1) + 2   ) += 0;		
	    _df(4*ele(1) + 3   ) += 0;		
	    break;

	case 2:
	    _df(4*ele(0)     ) +=  ((rhoNs * Tr)* du0_s  +
	    rho*Ms*dphi0U_s + 
	    (Nt - rhofr)* du0) *Jacwj ;	
	    _df(4*ele(0) + 1 ) +=  ((rhoNs * Tr)* du1_s  +
	    rho*Ms*dphi1U_s +
	    (Nt - rhofr)* du1) *Jacwj ;	
	
	    _df(4*ele(0) + 2   ) +=  ((rhoNs * Tz)* dv0_s +
	    rho*Ms*dphi0V_s + Mt*dv0_s
	    - rhofz* dv0) *Jacwj ;
	    _df(4*ele(0) + 3   ) +=  ((rhoNs * Tz)* dv1_s +
	    rho*Ms*dphi1V_s + Mt*dv1_s
	    - rhofz* dv1) *Jacwj ;

	    _df(4*ele(1)     ) +=  ((rhoNs * Tr)* du2_s  +
	    rho*Ms*dphi2U_s + 
	    (Nt - rhofr)* du2) *Jacwj ;	
	    _df(4*ele(1) + 1 ) +=  ((rhoNs * Tr)* du3_s  +
	    rho*Ms*dphi3U_s +
	    (Nt - rhofr)* du3) *Jacwj ;	
	
	    _df(4*ele(1) + 2   ) +=  ((rhoNs * Tz)* dv2_s +
	    rho*Ms*dphi2V_s + Mt*dv2_s
	    - rhofz* dv2) *Jacwj ;
	    _df(4*ele(1) + 3   ) +=  ((rhoNs * Tz)* dv3_s +
	    rho*Ms*dphi3V_s + Mt*dv3_s
	    - rhofz* dv3) *Jacwj ;
	    break;

	}
	
      }

    } //end: for quadrature
  } // end: for element

}
void VirusShell :: NeoHookean (Vector& est, Vector& kst){
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
    
    _lclResidue(0) = _C*(es + _nu*et); //Ns
    _lclResidue(1) = _C*(et + _nu*es); //Nt
    _lclResidue(2) = _D*(ks + _nu*kt); //Ms
    _lclResidue(3) = _D*(kt + _nu*ks); //Mt
    
  }
}

void VirusShell::printU(){
  int NDOF = _x.size();
  cout <<"U=[ ";
  for (int i=0; i < NDOF; i=i+4){
    cout << _x(i) << ", ";
  }
  cout << "];\n";
}

void VirusShell::writeMesh(string filename){
  ofstream myfile;
  myfile.open (filename.c_str());
  cout << "Writing mesh to a file.\n";
  for (auto it = _nodes.begin(); it != _nodes.end(); ++it) {
    myfile << *it << endl;
}
  myfile.close();
}

// void VirusShell::writeSolution(string filename){
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
void VirusShell::writeSolution(string filename){
  ofstream myfile;
  myfile.open (filename.c_str());
  cout << "Writing solution to a file.\n";
  myfile << setprecision(15);
  myfile << "x = [";
  for (auto i = 1; i <= _x.size(); i++) {
    if (i % 2 == 1)
      myfile << "\n[" << _x(i-1) << ",\t";
    else
      myfile << _x(i-1) << "],";
}
  myfile << "]\n";
  myfile << "from matplotlib import pyplot as plt\n";
  myfile << "from numpy import *\n";
  myfile << "x = matrix(x)\n";
  myfile << "plt.plot(x[0:len(x):2],'*-')\n";
  myfile << "plt.show()\n";
  myfile.close();
}
