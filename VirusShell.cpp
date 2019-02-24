#include "VirusShell.h"

void VirusShell :: fdf(){
  _f = 0; // reset energy
  _df.setZero(); //reset residue to zero
  int numEle = 0; // SET NUM ELEMENTS
  Vector wt = _quad.weights();
  vector<Vector> q = _quad.nodes();
  int numQuad = q.size();
  Vector est(2);
  Vector kst(2);

  //Reference surface
  double rho = 1;
  double rho_s = 0;
  double z_s = 1;
  
  for (int e=0; e< numEle; e++){
    for (int j=0; j<numQuad; j++) {
      double u, v, u_s, v_s, u_ss, v_ss;
      double Jac;

      // Quadrature
      Vector qj = q[j];
      double wj = wt(j);
      
      // Compute u, v, u', v', u'', v'' using shape functions


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

	_df(0) = 0;
	_df(1) = 0;

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

