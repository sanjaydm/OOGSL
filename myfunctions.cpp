#include"myfunctions.h"

int Plm_all(int l,int m,double x, double* plm){
	// To save time in computing the legendre polynomials (which takes up to bulk of computational
	//time) this function returns an array containing Pm_0 .... Pm_(l-1) evaluated at an x.
	//PLM_COUNT++; 
	
	double fact,oldfact,pll,pmm,pmmp1,omx2;
	int i,ll;
	
	if(m<0||m>l||fabs(x)>1.0){
		//cout << "Bad arguments in routine plgndr";
		return 0.0;
	}
	pmm=1.0;
	if(m>0){
		omx2=((1.0-x)*(1.0+x));
		fact=1.0;
		for(i=1;i<=m;i++){
			pmm = pmm*omx2*fact/(fact+1); 
			fact +=2.0;
		}
	}
	pmm=sqrt((2.0*m+1.0)*pmm/(4*PI));
	if(m % 2==1) pmm = -pmm;
	oldfact=sqrt(2.0*m+3.0);
	for(ll=0;ll<=l;ll++){
		if (ll<m) plm[ll]=0.0;
		else if(ll==m) plm[ll]=pmm;
		else if(ll==m+1) plm[ll]=x*sqrt(2.0*m+3.0)*pmm;
		else{
			fact = sqrt((4.0*ll*ll-1.0)/(ll*ll-m*m));
			plm[ll]=(x*plm[ll-1]-plm[ll-2]/oldfact)*fact;
			oldfact=fact;
		}
	}
	return 0.0;		
		
}

double W(double Phi,double beta,double m1,double m2){
	return beta*pow((Phi-m1)*(Phi-m2),2);
}
double Wprime(double Phi,double beta,double m1,double m2){
	return 4*beta*(Phi-m1)*(Phi-m2)*(Phi-(m1+m2)/2.);
}
double W2prime(double Phi,double beta,double m1,double m2){
	return 4*beta*((Phi-m2)*(Phi-(m1+m2)/2.)+(Phi-m1)*(Phi-(m1+m2)/2.)+(Phi-m1)*(Phi-m2));
}
double stiffness(double Phi,double m1,double m2,double k1,double k2){
	return (( k2 - k1 )/( exp( (2* Phi-(m1+m2))/(2*0.1)) + 1 )+k1)*0.5;
}
double stiffness_prime(double Phi,double m1,double m2,double k1,double k2){
	return (-(k2-k1)/pow(exp(1/2*(2*Phi-m1-m2)/.1)+1,2)/.1*exp(1/2*(2*Phi-m1-m2)/.1))*0.5;
}
double sqfac(double l,int m){
	double ret=l+m;
	for(int i=m-1;i>-m;i--)
		ret=ret*(l+i);
	return sqrt(ret);
		
}
