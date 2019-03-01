#include "VirusShell.h"

void VirusShell :: fdf(){
  _f = 0; // reset energy
  _df.setZero(); //reset residue to zero
  int numEle = _conn.size(); 
  Vector wt = _quad.weights();
  vector<Vector> q = _quad.nodes();
  int numQuad = q.size();
  Vector est(2);
  Vector kst(2);
  //Reference surface
  double rho = 1;
  double rho_s = 0; double rho_ss = 0;
  double z_s = 1; double z_ss = 0;
  double fr = 1; double fz = 0;
  for (int e=0; e< numEle; e++){
    for (int j=0; j<numQuad; j++) {
      
      Vector ele (2);
      ele << _conn[e](0); ele << _conn[e](1) ;
      cout << "\033[31m" << numEle << " " << numQuad << " here \033[0m\n";
      double u, v, u_s, v_s, u_ss, v_ss;
      double Jac;
      // Quadrature
      Vector qj = q[j];
      double wj = wt(j);
      
      // Compute u, v, u', v', u'', v'' using shape functions
      // For each node u,v,u' and v' are stored in this order

      // Strain measures
      double es_lin = 0; //u' rho' + z' v'
      double et_lin = u/rho; //u/rho
      double phi = 0; //-z' u' + rho' v'

      double k_s = 0; //phi'
      double k_t = v_s/rho; //v'/rho
      
      est(0) = es_lin + 0.5*phi*phi;
      est(1) = et_lin;
      kst(0) = k_s;
      kst(1) = k_t;

      // Call material law
      NeoHookean(est, kst);
      
      // If Energy flag is on
      if (_fFlag){
	_f += _lclEnergyDensity * Jac * wj;
      }
      // If residue flag is on
      if (_dfFlag){

	double cos_gamma = rho_s;
	double sin_gamma = z_s;
	double Tr = cos_gamma - phi*sin_gamma;
	double Tz = sin_gamma + phi*cos_gamma;

	double delta_es = Tr*(0) + Tz*(0); // Tr (du)' + Tz (dv)'
	double delta_et = 0/rho; //du/rho
	double delta_ks = -sin_gamma * (0) + cos_gamma*(0); // -sin g*(du)'+cos(g)*(dv)'
	double delta_kt = 0/rho; //(dv)'/rho

	double Ns = _lclResidue(0); double Nt = _lclResidue(1);
	double Ms = _lclResidue(2); double Mt = _lclResidue(2);
	double Jacwj = Jac * wj;
	double rhoNs = rho*Ns;

	// Variations
	double du0, du1, du2, du3;
	double dv0, dv1, dv2, dv3;

	double du0_s, du1_s, du2_s, du3_s;
	double dv0_s, dv1_s, dv2_s, dv3_s;

	double du0_ss, du1_ss, du2_ss, du3_ss;
	double dv0_ss, dv1_ss, dv2_ss, dv3_ss;

	double dphi0_s = -z_s* du0_s - z_ss* du0_ss + rho_s* dv0_s + rho_ss* dv0_ss; 
	double dphi1_s = -z_s* du1_s - z_ss* du1_ss + rho_s* dv1_s + rho_ss* dv1_ss;
	double dphi2_s = -z_s* du2_s - z_ss* du2_ss + rho_s* dv2_s + rho_ss* dv2_ss;
	double dphi3_s = -z_s* du3_s - z_ss* du3_ss + rho_s* dv3_s + rho_ss* dv3_ss; 
	
	_df(2*ele(0)     ) +=  ((rhoNs * Tr)* du0_s + (rhoNs * Tz)* dv0_s +
			       rho*Ms*dphi0_s + Mt*dv0_s + 
			       (Nt - fr)* du0 - fz* dv0) *Jacwj ;
	_df(2*ele(0) + 1 ) +=  ((rhoNs * Tr)* du1_s + (rhoNs * Tz)* dv1_s +
			       rho*Ms*dphi1_s + Mt*dv1_s + 
			       (Nt - fr)* du1 - fz* dv1) *Jacwj ;
	_df(2*ele(1)     ) += ((rhoNs * Tr)* du2_s + (rhoNs * Tz)* dv2_s +
			       rho*Ms*dphi2_s + Mt*dv2_s + 
			       (Nt - fr)* du2 - fz* dv2) *Jacwj ;
	_df(2*ele(1) + 1 ) += ((rhoNs * Tr)* du3_s + (rhoNs * Tz)* dv3_s +
			       rho*Ms*dphi3_s + Mt*dv3_s + 
			       (Nt - fr)* du3 - fz* dv3) *Jacwj ;
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
			      _C* (ks*ks + kt*kt + 2*_nu*ks*kt) );
  }
  if (_dfFlag) {
    
    _lclResidue(0) = _C*(es + _nu*et); //Ns
    _lclResidue(1) = _C*(et + _nu*es); //Nt
    _lclResidue(2) = _D*(ks + _nu*kt); //Ms
    _lclResidue(3) = _D*(kt + _nu*ks); //Mt
    
  }
}

