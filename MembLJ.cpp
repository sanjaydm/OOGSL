#include "MembLJ.h"
#include<gsl/gsl_sf.h>

using namespace std;

// Fixed particle coordinates
#define tN 1.0172220997017924 //PI/2
#define pN 0.00

void MembLJ::fdf(){

  _f = 0; // reset energy
  _df.setZero(); //reset residue to zero

  //Initialize variables
  Vector uP(_NP);      //  Particle positions
  Vector uP_t(_NP);    //  Particle gradients
  Vector uP_p(_NP);    //      - do - 
  Vector tP(_NP);      //  Particle 'theta' coords
  Vector pP(_NP);      //  Particle 'phi' coords

  Vector ai(_Ntot);  //  Area penalty residue contribs 
  Vector xi(_Ntot); //  X-cm penalty residue contribs 
  Vector yi(_Ntot);  //  Y-cm penalty residue contribs 
  Vector zi(_Ntot);  //  Z-cm penalty residue contribs 
  
  // Size of the system
  int Total_len = _Ntot;       // Global variable, number of modes
  double eps = 1e4;
  // START: Lookup table Sph Harm at Particles
  Matrix YlmP(_NP, _Ntot);

  for (int pi=0; pi < _NP; pi++){
    double ctP, pPi;
    if (pi == (_NP-1)) {
      ctP = cos(tN);
      pPi = pN;
    }
    else{
      ctP = cos(_x(_Ntot+2*pi));
      pPi = _x(_Ntot+2*pi + 1);
    }
    double* plmP;
    int arr_size = gsl_sf_legendre_array_n(_Lmax);
    plmP = new double[arr_size];
    gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM,_Lmax, ctP,
				       CSPHASE,  plmP);
    int absm, idx;
    int i=0;
    //------------Legendre Poly----------------------
    for(int l=0; l<=_Lmax; l++){
      for (int m=-l; m<=l; m++){
	absm = abs(m);
	idx = gsl_sf_legendre_array_index(l,absm);
	if (m<0){
	  YlmP(pi, i) = S2*plmP[idx]*sin(absm*pPi);
	}
	if (m==0){
	  YlmP(pi, i) = plmP[idx];
	}
	if (m>0){
	  YlmP(pi, i) = S2*plmP[idx]*cos(m*pPi);
	}
	i ++;
      }
    }
    delete[] plmP;
  }

  // END: Lookup table Sph Harm at Particles
  // Particle coordinates
  int arr_size = gsl_sf_legendre_array_n(_Lmax);
  double* plmP = new double[arr_size]; 
  double* plmP_1 = new double[arr_size];
  for (int pid=0; pid < _NP; pid++) {
    // The particle position follows modal coeffs
    if (pid == _NP-1 ){
      tP(pid) = tN;
      pP(pid) = pN;
    }
    else{
      tP(pid) = _x(_Ntot+2*pid);
      pP(pid) = _x(_Ntot+2*pid+1);
    }
    // Cos and sin values at particle pos.
    double ctP = cos(tP(pid));
    double stP = sin(tP(pid));
    double cpP = cos(pP(pid));
    double spP = sin(pP(pid));

    gsl_sf_legendre_deriv_alt_array_e(GSL_SF_LEGENDRE_SPHARM, _Lmax, ctP,
                                      CSPHASE,  plmP, plmP_1);
    // Compute particle positions
    uP(pid) = 1; 
    uP_t(pid) = 0;
    uP_p(pid) = 0;
    int idx;
    int i = 0;
    for(int l=0; l<= _Lmax; l++){
      for (int m=-l; m<=l; m++){
        int absm = abs(m);
        idx = gsl_sf_legendre_array_index(l,absm);
        if (m<0){
          uP  (pid) += S2*plmP[idx]*sin(absm*pP(pid))*_x(i);
          uP_t(pid) += S2*(plmP_1[idx])*sin(absm*pP(pid))*_x(i);
          uP_p(pid) += absm*S2*(plmP[idx])*cos(absm*pP(pid))*_x(i);
        }
        if (m==0){
          uP  (pid) +=   plmP[idx]*_x(i);
          uP_t(pid) += (plmP_1[idx])*_x(i);
        }
        if (m>0){
          uP  (pid) +=   S2*plmP[idx]*cos(m*pP(pid))*_x(i);
          uP_t(pid) += S2*(plmP_1[idx])*cos(m*pP(pid))*_x(i);
          uP_p(pid) += -m*S2*(plmP[idx])*sin(m*pP(pid))*_x(i);
        }
        i ++;
      }// end loop m 
    } // end loop l
  } // end loop pid

  // Computing particle energy
  for (int pidB=0; pidB < _NP; pidB++) {
    for (int pidA=pidB+1; pidA < _NP; pidA++) {
      if (pidA != pidB) {
        double enBA   = 0;
        Vector fBA(3);
        Vector PBA(3);

        double ctA = cos(tP(pidA));
        double stA = sin(tP(pidA));
        double cpA = cos(pP(pidA));
        double spA = sin(pP(pidA));

        double ctB = cos(tP(pidB));
        double stB = sin(tP(pidB));
        double cpB = cos(pP(pidB));
        double spB = sin(pP(pidB));

        double uA = uP(pidA); double uB = uP(pidB);
        fBA(0) = uB*stB*cpB - uA*stA*cpA;
        fBA(1) = uB*stB*spB - uA*stA*spA;
        fBA(2) = uB*ctB     - uA*ctA;
        LJ(fBA, _eps);
        _f += _LJEnergy;
      }
    }
  }
  // // Bypass and return particle energy (stored in energyPtr) if option=PARTICLE
  // if (option == PARTICLE) {return;}

  //Initialize penatly related quantities
  double area = 0;
  double volume = 0;
  double X0 = 0;
  double Y0 = 0;
  double Z0 = 0 ;
  double refArea = 4*PI;
  double refVolume = 4*PI/3.0;
  
  // START: Compute surface pos, derivatives at all quad points as a vector
  Vector clm = _x.view (0, _Ntot);
  //Vector u_vec(_qThetas.size());
  // Ylm must be defined outside
  Vector u_vec = _Ylm*clm;
  
  //Vector ut_vec(qThetas.size());
  Vector ut_vec = _Ylm_t*clm;

  //Vector utt_vec(qThetas.size());
  Vector utt_vec = _Ylm_tt*clm;

  //Vector up_vec(qThetas.size());
  Vector up_vec = _Ylm_p*clm;

  //Vector upp_vec(qThetas.size());
  Vector upp_vec = _Ylm_pp*clm;

  //Vector utp_vec(qThetas.size());
  Vector utp_vec = _Ylm_tp*clm;

  // END: Compute surface pos, derivatives at all quad points as a vector

  double t, p, wt;
  bool quadFlag = true; // Particle interaction is outside the integral (bypass flag)
  for(int qi = 0; qi < _qThetas.size(); qi++){ // quadrature loop
    t = _qThetas[qi];
    p = _qPhis[qi];
    wt = _qWts[qi];
    
    double ct=cos(t); double st=sin(t);
    double cp=cos(p); double sp=sin(p);
    // gsl_sf_legendre_deriv2_alt_array_e(GSL_SF_LEGENDRE_SPHARM, _Lmax, ct,
    //     			       CSPHASE,  plm, plm_1, plm_2);
    Vector& plm = _plms[qi];
    Vector& plm_1 = _plms_1[qi];
    Vector& plm_2 = _plms_2[qi];
    
    double u=1; 
    double u_t=0;
    double u_p=0.;
    double u_tt=0;
    double u_tp=0.;
    double u_pp=0.;

    int i = 0;
    int absm, idx;

    //START: Using vectorized version
    u = u_vec(qi)+1;
    u_t = ut_vec(qi);
    u_tt = utt_vec(qi);
    u_p = up_vec(qi);
    u_pp = upp_vec(qi);
    u_tp = utp_vec(qi);
    //END: Using vectorized version
    
    // Metric Tensor
    double At = (u_t); double Bt = (u); double Ct = 0;
    double At_T = (u_tt); double Bt_T = (u_t); double Ct_T = 0;
    double At_P = (u_tp); double Bt_P = (u_p); double Ct_P = 0;
  
    double Ap = (u_p); double Bp = (0); double Cp = (u*st);
    double Ap_T = (u_tp); double Bp_T = (0); double Cp_T = (u_t*st+u*ct);
    double Ap_P = (u_pp); double Bp_P = (0); double Cp_P = (u_p*st);
    
    double g_tt = At*At+Bt*Bt+Ct*Ct;
    double g_tp = At*Ap+Bt*Bp+Ct*Cp;
    double g_pp = Ap*Ap+Bp*Bp+Cp*Cp;
    double g = g_tt*g_pp-g_tp*g_tp;
    double Sg = pow(g,.5); 
    double gtt = g_pp/g; double gtp = -g_tp/g; double gpp = g_tt/g;
	
    double g_tt_T = 2*(At*At_T+Bt*Bt_T+Ct*Ct_T);
    double g_tp_T = At_T*Ap+At*Ap_T + Bt_T*Bp+Bt*Bp_T + Ct_T*Cp+Ct*Cp_T;
    double g_pp_T = 2*(Ap*Ap_T+Bp*Bp_T+Cp*Cp_T);
		
    double g_tt_P = 2*(At*At_P+Bt*Bt_P+Ct*Ct_P);
    double g_tp_P = At_P*Ap+At*Ap_P + Bt_P*Bp+Bt*Bp_P + Ct_P*Cp+Ct*Cp_P;
    double g_pp_P = 2*(Ap*Ap_P+Bp*Bp_P+Cp*Cp_P);
	
    double g_T = g_tt_T*g_pp+g_tt*g_pp_T -2*g_tp*g_tp_T;
    double g_P = g_tt_P*g_pp+g_tt*g_pp_P -2*g_tp*g_tp_P;
    double Sg_T = (g_T)/(2*Sg);
	
    double gtt_T = g_pp_T/g - g_pp*g_T/(g*g);
    double gtt_P = g_pp_P/g - g_pp*g_P/(g*g);
    double gtp_T = -(g_tp_T/g - g_tp*g_T/(g*g));
    double gtp_P = -(g_tp_P/g - g_tp*g_P/(g*g));
    double gpp_T = g_tt_T/g - g_tt*g_T/(g*g);
    double gpp_P = g_tt_P/g - g_tt*g_P/(g*g);
	
    // Second Derivative of Position Vector
    double Rr_tt = At_T-Bt; 
    double Tr_tt = Bt_T+At; 
    double Pr_tt = Ct_T;
	
	
    double Rr_tp = At_P-st*Ct; 
    double Tr_tp = Bt_P-ct*Ct; 
    double Pr_tp = At*st+Bt*ct+Ct_P;
	
    double Rr_pp = Ap_P-st*Cp;
    double Tr_pp = Bp_P-ct*Cp;
    double Pr_pp = Ap*st+Bp*ct+Cp_P;
  
    // Normal Vector 
    double Rn = (Bt*Cp-Bp*Ct)/Sg;
    double Tn = (Ct*Ap-Cp*At)/Sg;
    double Pn = (At*Bp-Ap*Bt)/Sg;
	
  
    double Rn_T = (Bt_T*Cp+Bt*Cp_T-Bp_T*Ct-Bp*Ct_T)/Sg-0.5*Rn*g_T/g;
    double Tn_T = (Ct_T*Ap+Ct*Ap_T-Cp_T*At-Cp*At_T)/Sg-0.5*Tn*g_T/g;
    double Pn_T = (At_T*Bp+At*Bp_T-Ap_T*Bt-Ap*Bt_T)/Sg-0.5*Pn*g_T/g;
	
    double Rn_P = (Bt_P*Cp+Bt*Cp_P-Bp_P*Ct-Bp*Ct_P)/Sg-0.5*Rn*g_P/g;
    double Tn_P = (Ct_P*Ap+Ct*Ap_P-Cp_P*At-Cp*At_P)/Sg-0.5*Tn*g_P/g;
    double Pn_P = (At_P*Bp+At*Bp_P-Ap_P*Bt-Ap*Bt_P)/Sg-0.5*Pn*g_P/g;
  
    // Second Fundamental Form
    double b_tt = Rn*Rr_tt+Tn*Tr_tt+Pn*Pr_tt;
    double b_tp = Rn*Rr_tp+Tn*Tr_tp+Pn*Pr_tp;
    double b_pp = Rn*Rr_pp+Tn*Tr_pp+Pn*Pr_pp;
	
    //Curvature Tensor
    double bt_t = gtt*b_tt+gtp*b_tp;
    double bt_p = gtt*b_tp+gtp*b_pp;
    double bp_t = gtp*b_tt+gpp*b_tp;
    double bp_p = gtp*b_tp+gpp*b_pp;
	
    double btt = gtt*bt_t+gtp*bt_p; 
    double btp = gtt*bp_t+gtp*bp_p; 
    double bpp = gtp*bp_t+gpp*bp_p;

    //Other quantities
    double dTR = Rn_T-Tn; double dPR = Rn_P - st*Pn;
    double dTRT = Tn_T + Rn; double dPRT = Tn_P - ct* Pn;
    double dTRP = st*Pn_T;   double dPRP = st*( ct*Tn + st*Rn + Pn_P);
    
    //----------------------------------------------------------------
    //Curvatures
    double H = 0.5*(bt_t+bp_p); double K = (bt_t*bp_p - bt_p*bp_t);
    
    double ui= 0;//r_i(i,t,p);
    double upi= 0; //rp_i(i,t,p);
    double uti= 0;
    double uppi= 0; //rpp_i(i,t,p);
    double utpi= 0; //rzp_i(i,t,p);
    double utti= 0;
    double uiA = 0;
    double uiB = 0;
    i = 0;

    area += Sg*wt;
    //volume += u/3*(Rn)*Sg*wt;
    X0 += u*st*cp*Sg*wt;
    Y0 += u*st*sp*Sg*wt;
    Z0 += u*ct*Sg*wt;

    // Energy
    // Volume is computed as a surface integral 1/3. <f,n>
    if (_fFlag) {
     _f += ( _k*H*H)*Sg*wt; 
    //*energyPtr += ( _k*H*H + gamma  +  pressure/3.*u*Rn )*Sg/st*wt-gamma*wt; 
    }
    
    if (_dfFlag) {
      
      // START: Vectorization of residue
      // Variation of Sg
      Vector Ylmqi(_Ylm.row(qi));
      Vector Ylm_tqi( _Ylm_t.row(qi));
      Vector Ylm_pqi( _Ylm_p.row(qi));
      Vector Ylm_ttqi( _Ylm_tt.row(qi));
      Vector Ylm_tpqi( _Ylm_tp.row(qi));
      Vector Ylm_ppqi( _Ylm_pp.row(qi));
      
      //Vector dSg(_Ntot);
      gsl_vector* dSg = gsl_vector_calloc(_Ntot);
      //dSg += (gtt*u_t + gtp*u_p)*Ylm_tqi;
      gsl_blas_daxpy (gtt*u_t + gtp*u_p, Ylm_tqi._gsl_vec, dSg); //gtt*u_t*uti
      //dSg += (gtp*u_t + gpp*u_p)*Ylm_pqi;
      gsl_blas_daxpy (gtp*u_t + gpp*u_p, Ylm_pqi._gsl_vec, dSg); //+gtp*u_t*upi
      //dSg += (u*(gtt+st*st*gpp))*Ylmqi;
      gsl_blas_daxpy (u*(gtt+st*st*gpp), Ylmqi._gsl_vec, dSg); //+u*(gtt+st*st*gpp)*ui

      // Variation of H
      //Vector dH(_Ntot);
      gsl_vector* dH = gsl_vector_calloc(_Ntot);
      
      double temp_a = 0.5*((gtt*g_T+gtp*g_P)/(2*g) + (gtt_T+gtp_P));
      double temp_b = 0.5*((gtp*g_T+gpp*g_P)/(2*g) + (gtp_T+gpp_P));

      //dH += (gtt*dTR+gtp*dPR + temp_a*Rn + gtt*Tn + gtp*st*Pn)*Ylm_tqi;
      gsl_blas_daxpy( gtt*dTR+gtp*dPR + temp_a*Rn + gtt*Tn + gtp*st*Pn, Ylm_tqi._gsl_vec, dH); 
      //dH += (gtp*dTR+gpp*dPR + temp_b*Rn + gtp*Tn + gpp*st*Pn)*Ylm_pqi;
      gsl_blas_daxpy( gtp*dTR+gpp*dPR + temp_b*Rn + gtp*Tn + gpp*st*Pn, Ylm_pqi._gsl_vec, dH);
      //dH += (gtt*dTRT + gtp*(dTRP+dPRT) + gpp*dPRP + temp_a*Tn + temp_b*st*Pn +
      //0.5*( -gtt*Rn + 2*gtp*ct*Pn  - gpp*st*(st*Rn + ct*Tn)))*Ylmqi;
      gsl_blas_daxpy( gtt*dTRT + gtp*(dTRP+dPRT) + gpp*dPRP + temp_a*Tn + temp_b*st*Pn +
		      0.5*( -gtt*Rn + 2*gtp*ct*Pn  - gpp*st*(st*Rn + ct*Tn)), Ylmqi._gsl_vec, dH);
      //dH += ( 0.5*gtt*Rn)* Ylm_ttqi;
      gsl_blas_daxpy( 0.5*gtt*Rn, Ylm_ttqi._gsl_vec, dH);
      //dH += (gtp*Rn)* Ylm_tpqi;
      gsl_blas_daxpy( gtp*Rn , Ylm_tpqi._gsl_vec, dH);
      //dH += (0.5*gpp*Rn)* Ylm_ppqi;
      gsl_blas_daxpy( 0.5*gpp*Rn , Ylm_ppqi._gsl_vec, dH); 

      Vector outview = _df.view(0, _Ntot);
      gsl_blas_daxpy(2*wt*_k*H*Sg, dH, outview._gsl_vec); gsl_blas_daxpy(wt*_k*H*H*Sg, dSg, outview._gsl_vec); //out->data(i) += wt*( 2*_k*H*dH + _k*H*H * dSg )*Sg;     
      gsl_blas_daxpy(wt*Sg, dSg, ai._gsl_vec); //ai->data(i) += wt*dSg*Sg;
      gsl_blas_daxpy(wt*u*st*cp*Sg, dSg, xi._gsl_vec); gsl_blas_daxpy(wt*st*cp*Sg, Ylmqi._gsl_vec, xi._gsl_vec); //xi->data(i) += wt*(u*st*cp*dSg + st*cp*ui)*Sg;
      gsl_blas_daxpy(wt*Sg*u*st*sp, dSg, yi._gsl_vec); gsl_blas_daxpy(wt*Sg*st*sp, Ylmqi._gsl_vec, yi._gsl_vec);  //yi->data(i) += wt*(u*st*sp*dSg + st*sp*ui)*Sg;
      gsl_blas_daxpy(wt*Sg*u*ct, dSg, zi._gsl_vec); gsl_blas_daxpy(wt*ct*Sg, Ylmqi._gsl_vec, zi._gsl_vec); //zi->data(i) += wt*(u*ct*dSg    + ct*ui)*Sg;

      // Particle interactions
      // Efficetively "Outside" quadrature loop, as it should be
      if (quadFlag) {
        for (int l=0;l<= _Lmax;l++){
          for (int m=-l; m<=l; m++){
            int absm = abs(m);
            idx = gsl_sf_legendre_array_index(l,absm);

            for (int pidB=0; pidB < _NP; pidB++) {
              for (int pidA=pidB+1; pidA < _NP; pidA++) {
                if (pidA != pidB) {
                  double enBA   = 0;
                  Vector fBA(3);
                 

                  double ctA = cos(tP(pidA));
                  double stA = sin(tP(pidA));
                  double cpA = cos(pP(pidA));
                  double spA = sin(pP(pidA));

                  double ctB = cos(tP(pidB));
                  double stB = sin(tP(pidB));
                  double cpB = cos(pP(pidB));
                  double spB = sin(pP(pidB));

                  double uA = uP(pidA); double uB = uP(pidB);
                  fBA(0) = uB*stB*cpB - uA*stA*cpA;
                  fBA(1) = uB*stB*spB - uA*stA*spA;
                  fBA(2) = uB*ctB     - uA*ctA;
                  LJ(fBA,  _eps);

		  Vector PBA = _LJForce;
		  
                  uiA = YlmP(pidA, i);
                  uiB = YlmP(pidB, i);
                  _df(i) += (PBA(0)*(stB*cpB*uiB-stA*cpA*uiA) +
                                     PBA(1)*(stB*spB*uiB-stA*spA*uiA) +
                                     PBA(2)*(ctB*uiB-ctA*uiA));
                  if (i==0) {
                    // Particle-Particle Residue
                    // Only once per residue computation
                    double uA_t = uP_t(pidA); double uB_t = uP_t(pidB);
                    double uA_p = uP_p(pidA); double uB_p = uP_p(pidB);
                    double At_A = (uA_t); double Bt_A = (uA); double Ct_A = 0;
                    double Ap_A = (uA_p); double Bp_A = (0); double Cp_A = (uA*stA);
                  
                    double At_B = (uB_t); double Bt_B = (uB); double Ct_B = 0;
                    double Ap_B = (uB_p); double Bp_B = (0); double Cp_B = (uB*stB);

                    double fBt_x = (uB_t*stB + uB*ctB)*cpB;
                    double fBt_y = (uB_t*stB + uB*ctB)*spB;
                    double fBt_z = (uB_t*ctB - uB*stB);

                    double fAt_x = (uA_t*stA + uA*ctA)*cpA;
                    double fAt_y = (uA_t*stA + uA*ctA)*spA;
                    double fAt_z = (uA_t*ctA - uA*stA);

                    double fBp_x = (uB_p*stB*cpB - uB*stB*spB);
                    double fBp_y = (uB_p*stB*spB + uB*stB*cpB);
                    double fBp_z =  uB_p*ctB;

                    double fAp_x = (uA_p*stA*cpA - uA*stA*spA);
                    double fAp_y = (uA_p*stA*spA + uA*stA*cpA);
                    double fAp_z =  uA_p*ctA;

                    if (pidB != (_NP-1)){
                      _df(_Ntot+pidB*2)   +=   (PBA(0)*fBt_x + PBA(1)*fBt_y + PBA(2)*fBt_z);
                      _df(_Ntot+pidB*2+1) += (PBA(0)*fBp_x + PBA(1)*fBp_y + PBA(2)*fBp_z);
                    }
                    if (pidA!=(_NP-1)){
		      _df(_Ntot+pidA*2)   +=   -(PBA(0)*fAt_x + PBA(1)*fAt_y + PBA(2)*fAt_z);
		      _df(_Ntot+pidA*2+1) += -(PBA(0)*fAp_x + PBA(1)*fAp_y + PBA(2)*fAp_z);
                      }
                
                  } // end if i==0
                }
              } // end loop pidA
            } // end loop pidB
            i++;
          } // end loop m
        } // end loop l
      } // end if quadFlag
      gsl_vector_free(dH);
      gsl_vector_free(dSg);
      // END: Vectorization of residue
    } // end if RESIDUE
    quadFlag = false;
  } // end quad loop
  
  // Energy contribution from penalty
  _f += pow(area - refArea, 2)*0.5*eps;
  //_f += 0*pow(volume - refVolume,2)*0.5*eps;
  _f += (pow(X0,2)+pow(Y0,2) + pow(Z0,2) )*0.5*eps;
  //_f += (pow(tP[_NP-1]-PI/2,2) + pow(pP[_NP-1],2) )*0.5*eps;
  //_f += (pow(_x[1],2) + pow(_x[2],2) + pow(_x[3],2) )*0.5*eps; // setting l=1 modes to zero
  // Rotation constraint
  
  //Compute constraint forces
  double gamma= (area-refArea)*eps;
  double lambX = (X0)*eps;
  double lambY = (Y0)*eps;
  double lambZ = (Z0)*eps;

  //double lambTheta = (tP[_NP-1]-PI/2)*eps;
  //double lambPhi = (pP[_NP-1])*eps;
  if (_dfFlag) {
    // START: Vectorization
    //gsl_vector_view outview = gsl_vector_subvector(out, 0, _Ntot);
    Vector outview = _df.view(0, _Ntot);
    gsl_blas_daxpy(gamma, ai._gsl_vec, outview._gsl_vec);
    gsl_blas_daxpy(lambX, xi._gsl_vec, outview._gsl_vec);
    gsl_blas_daxpy(lambY, yi._gsl_vec, outview._gsl_vec);
    gsl_blas_daxpy(lambZ, zi._gsl_vec, outview._gsl_vec);
    // END: Vectorizatoin      
  }

}


void MembLJ::LJ(Vector f12, double e){
  // //Harmonic 
  // double r = sqrt( f12(0)*f12(0) + f12(1)*f12(1) + f12(2)*f12(2) );
  // if(option==ENERGY)
  //   *energy = 0.5*e*pow(r-rm, 2);
  // else if (option==RESIDUE) {
  //   force(0) =  e*(r-rm)*f12(0)/r;
  //   force(1) =  e*(r-rm)*f12(1)/r;
  //   force(2) =  e*(r-rm)*f12(2)/r;
  // }

    // --------------- Morse -----------------------
  double r = sqrt( f12(0)*f12(0) + f12(1)*f12(1) + f12(2)*f12(2) );
  //De = e
  double a = 6/_re; //At this value r = 2re and exponent is e-1
  double fact = (1-exp(-a*(r-_re)));
  _LJEnergy = e* pow(fact,2);
  //  cout << "\033[32m " << r << "\033[0m\n";
  if (_dfFlag) {
    double fact2 = 2*e*fact*exp(-a*(r-_re))*a/r;
    _LJForce(0) =  fact2 * f12(0);
    _LJForce(1) =  fact2 * f12(1);
    _LJForce(2) =  fact2 * f12(2);
  }


  // // --------------- LJ -----------------------
  // double sigma = _re*pow(2,-1./6.);
  // double r = sqrt( f12(0)*f12(0) + f12(1)*f12(1) + f12(2)*f12(2) );
  // double sigma_r6 = pow(sigma/r,6);
  // *energy = 4*e*( pow(sigma_r6,2) - sigma_r6 );
  // //  cout << "\033[32m " << r << "\033[0m\n";
  // if (option==RESIDUE) {
  //   double fact = 24 *e/r * sigma_r6 * (1-2*sigma_r6)/r ;
  //   force(0) =  fact * f12(0);
  //   force(1) =  fact * f12(1);
  //   force(2) =  fact * f12(2);
  // }

  // //---------------- Truncated and Shifted LJ ----------------
  // double sigma = _re*pow(2,-1./6.);
  // double rc = 2.5*sigma;
  // double r = sqrt( f12(0)*f12(0) + f12(1)*f12(1) + f12(2)*f12(2) );
  // double sigma_r6 = pow(sigma/r,6);
  // double sigma_rc6 = pow(sigma/rc,6);
  // if (r<=rc) {
  //   *energy = 4*e*(( pow(sigma_r6,2) - sigma_r6 ) -
  //                  ( pow(sigma_rc6,2) - sigma_rc6 ));
  // }
  // else {
  //   *energy = 0;
  // }
  // if (option==RESIDUE) {
  //   double fact = 24 *e/r * sigma_r6 * (1-2*sigma_r6)/r ;
  //   if (r<=rc) {
  //     force(0) =  fact * f12(0);
  //     force(1) =  fact * f12(1);
  //     force(2) =  fact * f12(2);
  //   }
  //   else {
  //     force(0) =  0;
  //     force(1) =  0;
  //     force(2) =  0;
  //   }
    
  // }

} 

// --------- Helper function: PRINT DATA TO VTK FILE ---------------
void MembLJ::printToVTK(string filename){
  /*
  ifstream inp("sphere_2562_nodes.dat");
  ifstream inpConn("sphere_2562_conn.dat");
  string filenameExt = filename + ".vtk";
  string filenamePartExt = filename + "_Ps.vtk";
  ofstream out(filenameExt.c_str());
  ofstream outAB(filenamePartExt.c_str());
  double* arg = uvec->data;
  
  int NNodes =0;
  int dim = 0;
  inp >> NNodes >> dim;
  
  // Header
  char spchar = '#';
  out << spchar << " vtk DataFile Version 3.1" << endl;
  out << "vtk output" << endl;
  out << "ASCII" << endl;
  out << "DATASET POLYDATA" << endl;
  out << "POINTS " << NNodes << " FLOAT" << endl;

  double x, y, z;
  double t, p;
  vector <double> uP(_NP);
  vector <double> tP(_NP);
  vector <double> pP(_NP);

  // Read Quadrature Rule
  //string lebfile("quad.dat");

  // Particle coordinates
  for (int pid=0; pid < _NP; pid++) {
    if (pid == (_NP-1) ){
      tP(pid) = tN;
      pP(pid) = pN;
    }
    else{
      tP(pid) = arg[_Ntot+2*pid];
      pP(pid) = arg[_Ntot+2*pid+1];
    }
    double ctP = cos(tP(pid));
    double stP = sin(tP(pid));
    double cpP = cos(pP(pid));
    double spP = sin(pP(pid));
    gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM, _Lmax, ctP,
				       CSPHASE,  plmP);

    // gsl_sf_legendre_deriv_alt_array_e(GSL_SF_LEGENDRE_SPHARM,__Lmax, ctP,
    //                                   CSPHASE,  plmP, plmP_1);
    uP(pid) = 1; 
    int idx;
    int i = 0;
    //------------Legendre Poly----------------------
    for(int l=0; l<= _Lmax; l++){
      for (int m=-l; m<=l; m++){
        int absm = abs(m);
        idx = gsl_sf_legendre_array_index(l,absm);
        if (m<0){

          uP(pid) +=   S2*plmP[idx]*sin(absm*pP(pid))*arg(i);
        }
        if (m==0){
          uP(pid) +=   plmP[idx]*arg(i);
        }
        if (m>0){
          uP(pid) +=   S2*plmP[idx]*cos(m*pP(pid))*arg(i);
        }
        i ++;
      }// end loop m 
    } // end loop l
  } // end loop pid

  int three, v0, v1, v2;
  for (int i=0; i < NNodes; i++){

    inp >> x >> y >> z;
    double nrm = sqrt(x*x+y*y+z*z);
    x = x/nrm; y = y/nrm; z= z/nrm;
    t = acos(z/nrm);
    p = atan2(y,x);
    //int Total_len = _Ntot;
    int idx;
    int md_i = 0;
    double ct=cos(t); double st=sin(t);
    double cp=cos(p); double sp=sin(p);

    gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM, _Lmax, ct,
				       CSPHASE,  plm);

    double u  = 1.0;
    for(int l=0; l<= _Lmax; l++){
      for (int m=-l; m<=l; m++){
	int absm = abs(m);
	idx = gsl_sf_legendre_array_index(l,absm);
	if (m<0){
	  u +=   S2*plm[idx]*sin(absm*p)*arg[md_i];
	  
	}
	if (m==0){
	  u +=   plm[idx]*arg[md_i];
	}
	if (m>0){

	  u +=   S2*plm[idx]*cos(m*p)*arg[md_i];
	}
	md_i++;
      }
      
    }
    //u=u+1; // In the following (u,v,w) represents current configuration
    out << u*st*cp << " " << u*st*sp << " " << u*ct<< endl;
        
  }

  int NFaces =0;
  int temp=0;
  inpConn >> NFaces >> temp;

  out << "POLYGONS " << NFaces << " " 
      << NFaces * 4 << endl;

  for(int i=0; i<NFaces; i++){
    inpConn >> three >> v0 >> v1 >> v2;
    out << three << " " << v0 << " " << v1 << " "<< v2 << endl;
  }

  // Write Particle Post
  outAB << spchar << " vtk DataFile Version 3.1" << endl;
  outAB << "vtk output" << endl;
  outAB << "ASCII" << endl;
  outAB << "DATASET POLYDATA" << endl;
  outAB << "POINTS " << _NP << " FLOAT" << endl;
  for (int pid=0; pid < _NP; pid++) {
    outAB << uP(pid)*sin(tP(pid))*cos(pP(pid)) << " " << uP(pid)*sin(tP(pid))*sin(pP(pid)) << " " << uP(pid)*cos(tP(pid))<< endl;
  }
  
  
  // Close file
  out.close();
  outAB.close();
  inp.close();
  inpConn.close();
  */
}

void MembLJ::printModes(string filename){
  /*
  string filenameExt = filename + "_modes.txt";
  ofstream out(filenameExt.c_str());
  double* para = (double*) paraVoid;
  out << setprecision(15);
  out << para[0] << endl; //kappa
  out << para[1] << endl; //eps
  out << para[2] << endl; //re
  
  for (int i=0; i<_Ntot; i++){
    out << uvec->data[i] << "\n";
  }
  */
}
void MembLJ::printParticleCoords(string filename){
  /*
  string filenameExt = filename + "_Particle_theta_phi.txt";
  ofstream out(filenameExt.c_str());
  double* para = (double*) paraVoid;
  out << setprecision(15);
  out << para[0] << endl; //kappa
  out << para[1] << endl; //eps
  out << para[2] << endl; //re

  for (int i=0; i<_NP-1; i++){
    out << uvec->data[_Ntot+2*i] <<"\n";
    out << uvec->data[_Ntot+2*i+1] <<"\n";

  }
  */
}
// void generateGausLegendreQuad(const char* filename) {
//   ofstream fid(filename);
//   // Use only even order 
//   int orderX = 10;
//   int orderY = 10;
//   double* X = x10;
//   double* Y = x10;
//   double* WX = w10;
//   double* WY = w10;
//   for (int ix=0; ix<orderX/2; ix++){
//     for (int iy=0; iy<orderY/2; iy++){
//       double tP = PI/2*(X[ix] + 1);
//       double pP = PI*(Y[iy] + 1);
//       double tM = PI/2*(-X[ix] + 1);
//       double pM = PI*(-Y[iy] + 1);
//       fid << setprecision(15) ;
//       fid << pP << " " << tP << " " << WX[ix]*WY[iy]*PI*PI/2 << endl;
//       fid << pP << " " << tM << " " << WX[ix]*WY[iy]*PI*PI/2 << endl;
//       fid << pM << " " << tP << " " << WX[ix]*WY[iy]*PI*PI/2 << endl;
//       fid << pM << " " << tM << " " << WX[ix]*WY[iy]*PI*PI/2 << endl;
//     }
//   }
//   fid.close();
// }

void MembLJ::generateCeliaReinaQuad(int orderX, int orderY) {
  // Celia Reina's Quad
  //ofstream fid(filename);
  double* X ;
  double* WX;
  switch (orderX) {
  case 16:
    X = x16;
    WX = w16;
    break;
  case 20:
    X = x20;
    WX = w20;
    break;
  case 32:
    X = x32;
    WX = w32;
    break;
  case 64:
    X = x64;
    WX = w64;
    break;
  }
  cout << "Creating quadrature nodes and weights ...";
  for (int ix=0; ix<orderX/2; ix++){
    for (int iy=0; iy<orderY; iy++){
      double tP = acos(X[ix]);
      double tM = acos(-X[ix]);
      double p = 2*PI/(orderY)*iy;

      _qThetas.push_back(tP);
      _qPhis.push_back(p);
      _qWts.push_back(WX[ix]/sin(tP)*(2*PI)/orderY);

      _qThetas.push_back(tM);
      _qPhis.push_back(p);
      _qWts.push_back(WX[ix]/sin(tM)*(2*PI)/orderY);
    }
  }
  cout << " done" << endl;
}

void MembLJ::computeYlmLookUp(){
  // Generate Spherical Harmonic Lookup
  cout << "Generating spherical harmonic look up table ... ";
  _Ylm.setDim(_qThetas.size(), _Ntot);
  _Ylm_t.setDim(_qThetas.size(), _Ntot);
  _Ylm_tt.setDim(_qThetas.size(), _Ntot);
  
  _Ylm_p.setDim(_qThetas.size(), _Ntot);
  _Ylm_pp.setDim(_qThetas.size(), _Ntot);
  _Ylm_tp.setDim(_qThetas.size(), _Ntot);

  for (int qi=0; qi < _qThetas.size(); qi++){
    double ct = cos(_qThetas[qi]);
    double p = _qPhis[qi];
    double* plm;
    double* plm_1;
    double* plm_2;
    int arr_size = gsl_sf_legendre_array_n(_Lmax);
    plm = new double[arr_size];
    plm_1 = new double[arr_size];
    plm_2 = new double[arr_size];
    gsl_sf_legendre_deriv2_alt_array_e(GSL_SF_LEGENDRE_SPHARM,_Lmax, ct,
				       CSPHASE,  plm, plm_1, plm_2);
    int absm, idx;
    int i=0;
    //------------Legendre Poly----------------------
    for(int l=0; l<=_Lmax; l++){
      for (int m=-l; m<=l; m++){
	absm = abs(m);
	idx = gsl_sf_legendre_array_index(l,absm);
	if (m<0){
	  _Ylm   (qi, i) = S2*plm[idx]*sin(absm*p);
	  _Ylm_t (qi, i) = S2*(plm_1[idx])*sin(absm*p);
	  _Ylm_tt(qi, i) = S2*(plm_2[idx])*sin(absm*p);

	  _Ylm_p (qi, i) = absm*S2*(plm[idx])*cos(absm*p);
	  _Ylm_pp(qi, i) = -m*m*S2*(plm[idx])*sin(absm*p);

	  _Ylm_tp(qi, i) = absm*S2*(plm_1[idx])*cos(absm*p);

	}
	if (m==0){
	  _Ylm    (qi, i) = plm[idx];
	  _Ylm_t  (qi, i) = plm_1[idx];
	  _Ylm_tt (qi, i) = plm_2[idx];	  
	}
	if (m>0){

	  _Ylm    (qi, i) =  S2*plm[idx]*cos(m*p)     ;
	  _Ylm_t  (qi, i) = S2*(plm_1[idx])*cos(m*p) ;
	  _Ylm_tt (qi, i) = S2*(plm_2[idx])*cos(m*p) ;

	  _Ylm_p  (qi, i) = -m*S2*(plm[idx])*sin(m*p);
	  _Ylm_pp (qi, i) = -m*m*S2*(plm[idx])*cos(m*p);

          _Ylm_tp (qi, i) = -m*S2*(plm_1[idx])*sin(m*p);
	}
	i ++;
      }

    }
    delete[] plm;
    delete[] plm_1;
    delete[] plm_2;
    
  }
  cout << "done\n";
  // END: Generating spherical harmonic lookup table

}

void MembLJ::generateAssociatedLegendreLookUp(){
  
  //START: Generating associated legendre lookup
  cout << "Generating associated legendre look up table ... ";
  int arr_size = gsl_sf_legendre_array_n(_Lmax);
  for (int i=0; i < _qThetas.size(); i++){
    double ct = cos(_qThetas[i]);
    Vector temp(arr_size);
    Vector temp_1(arr_size);
    Vector temp_2(arr_size);
    gsl_sf_legendre_deriv2_alt_array_e(GSL_SF_LEGENDRE_SPHARM, _Lmax, ct,
				       CSPHASE,  (temp._gsl_vec)->data, (temp_1._gsl_vec)->data, (temp_2._gsl_vec)->data);
    _plms.push_back(temp);
    _plms_1.push_back(temp_1);
    _plms_2.push_back(temp_2);
  }
  cout << "done\n";
  // END: Lengendre lookup
  
}
