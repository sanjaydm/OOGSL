#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "Problem.h"
#include "Group.h"
#include "MultiMin.h"
#include "MultiRoot.h"
#include "wigner.cpp"
#include "compute_D_matrix.cpp"
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include "SymmReduced.h"
#include "MembLJ.h"
#include <iomanip>
#include<string.h>
#include <sstream>
#include <fstream>
#include "kbhit.h"
#include"continuation.h"

using namespace std;
int sz, sz_para;
SymmReduced symmP;

int f_wrapper(double* indvar, double* invar, gsl_vector* out){
  /* Wrapper for the fdf function so that it can be used with continuation code */
  Vector in(invar, sz);
  Vector par(sz_para);
  for (int i=0; i < sz_para; i++)
    par(i) = invar[sz+i];

  symmP._para = par;
  symmP._x = in;
  symmP.fdf();

  // TODO: This for loop must be replaced by an efficient alternative
  for (int i=0; i<sz; i++)
    out->data[i] = symmP._df(i);
  
  return 0;
};

int main(int argc, char** argv){
  // Number of modes and particles  
  int N = 12;
  int NP = 12;
  int Ntot= (N+1)*(N+1);
  Vector in(Ntot+3*NP); 
  Vector para(3);
  Vector discPara(2); //Discretization parameters
  discPara(0) = N; //Lmax
  discPara(1) = NP; //Number of particles
  

  // Parameters of the model
  double re = 2*sqrt(4.0/NP)-0.1; 
  double k = 5;
  double eps = 1;
  para(0) = k;   //bending stiffness kappa
  para(1) = eps;  //eps
  para(2) = re; //rm = equilibirum distance

  
  // Initialize the input vector
  srand(time(NULL));
  for (int i=0; i<Ntot; i++){
    in(i) = 0.0;//*rand()/RAND_MAX;
  }
  for (int j=0; j< NP; j++) {
    in(Ntot + 3*j) =  double(rand())/RAND_MAX;
    in(Ntot + 3*j+1) = double(rand())/RAND_MAX;
    in(Ntot + 3*j+2) = double(rand())/RAND_MAX;
  }

  // Construct the model
  MembLJ prob(in, para, discPara);
  prob.checkConsistency();

  //return 0;

  // Generators of Icosahedral Group
    Matrix g2(3,3);
   g2 <<  -0.809016994374947 <<  0.500000000000000 << -0.309016994374947
   << 0.500000000000000 <<  0.309016994374947 << -0.809016994374947
   << -0.309016994374948 << -0.809016994374947 << -0.500000000000000;

   Matrix g3(3,3);
   g3 <<   -0.0000  <<  -0.0000  <<  -1.0000
      <<    1.0000  <<   0.0000  <<  -0.0000
      <<      0     <<  -1.0000  <<   0.0000;

   // Construct Icosahedral group-rep on spherical harmonics
   Matrix r2 = compute_R(0,N,g2);
   Matrix r3 = compute_R(0,N,g3);

   // Construct permutation representation
   vector<int> g3p = {8,5,6,9,7,12,2,10,11,1,4,3};
   vector<int> g2p = {7,5,8,12,2,10,1,3,11,6,9,4};

   Matrix g2_Mat = constructMat(g2p,g2);
   Matrix g3_Mat = constructMat(g3p,g3);


   // Icosahedral group-rep on spherical harmonics+particles
   int GN = r2.size()[0] +g2_Mat.size()[0];
   Matrix G2(GN, GN);
   Matrix G3(GN, GN);
   for (int i=0; i <GN; i++){
     for (int j=0; j<GN; j++){
       if(i<r2.size()[0] && j < r2.size()[0]){
	 G2(i,j) = r2(i,j);
	 G3(i,j) = r3(i,j);
       }
       if(i>=r2.size()[0] && j >= r2.size()[0]){
	 G2(i,j) = g2_Mat(i-r2.size()[0],j-r2.size()[0]);
	 G3(i,j) = g3_Mat(i-r2.size()[0],j-r2.size()[0]);
       }
     }
   }

   Icosahedral im(G2,G3);

   // Construct projection operators
   Projection pm(im);
   Matrix rg = pm._P.range();
   //rg.print();
   
   int nn = 1;
   
   // Construct symmetry reduced problem
   Vector x0 = rg.T()*in;
   //x0.print();
   
   symmP = SymmReduced(x0,para,rg,&prob);
   symmP.checkConsistency();
   sz = x0.size();
   sz_para = para.size();
   

   MultiMin M("lbfgs", &symmP);
   M._LBFGSB_Initialize();
  
   //Print Energy to output
   ofstream outFileAp("Energy.txt", ofstream::app);
   outFileAp << "kappa\t rm\t Total Energy\n";
   outFileAp.close();

   // One step LBFGS minimization
   for (int i=0; i<1; i++) {
     // Update parameter re
     symmP._para(2) = symmP._para(2) + 0.05*(i/30.0);
     cout << "re = " << prob._re << endl;

     // Solve
     M.LBFGSB_Solve();
     //symmP._x.print();

     //Print Energy to output
     ofstream outFileAp("Energy.txt", ofstream::app);

     outFileAp.setf(ios_base::scientific);
     double projg = M._dsave[12];
     double enTot = symmP._f;
     outFileAp << symmP._para(0) << " \t " << symmP._para(2)<< "\t" << enTot << "\t" <<  projg << endl;
     outFileAp.close();

     // Print to VTK file
     ostringstream toStringNP, toStringCntr; 
     toStringNP << NP;
     toStringCntr << i;
     string temp("Icosa");
     temp = temp + "_N_" + toStringNP.str() +  "_" + toStringCntr.str();
     prob.printToVTK(temp);

   }

   // Numerical Continuation Starts
   double* trivial = new double[sz+sz_para];
   for (int i=0; i< sz; i++)
     trivial[i] = symmP._x(i);
   for (int i=0; i< sz_para; i++)
     trivial[i+sz] = para(i);

   Continuer c;
   c.ptr_F = f_wrapper;

   c.setNoOfUnknownVar(sz);
   c.setNoOfInArgs(sz+sz_para);
   c.noOfIndVar = 0;
   c.noOfPara=sz_para;
   c.performQuad = false;
   c.posOfPara=sz+2;  //Corresponds to re
   c.correctorType="quasi-newton";
	
		
   c.setInitialSolution(trivial);
   int key;
   int i=1;
   do{
     if(kbhit())
       {
	 key = fgetc(stdin);
	 if (key == 'i'){
	   c.stepSize += .01;
	   cout<< "\n \033[1;35m Step Size increased to "<< c.stepSize<<"\033[m\n\n";
	 }
	 else if(key=='d'){
	   c.stepSize -= .01;
	   cout<< "\n \033[1;36m Step Size decreased to "<< c.stepSize<<"\033[m\n\n";
	 }
	 else if(key=='r'){
	   c.stepSize=-c.stepSize;
	   cout<< "\n \033[1;36m Step Size reversed \033[m\n\n";
	 }
	 else if(key =='n'){
	   c.stepSize = 2*c.stepSize;
	   cout<< "\n \033[1;36m Step Size doubled to"<<c.stepSize<<" \033[m\n\n";
	 }
	 else if(key =='h'){
	   c.stepSize = c.stepSize/2.0;
	   cout<< "\n \033[1;33m Step Size halved to"<<c.stepSize<<" \033[m\n\n";
	 }
	 fflush(stdin);
       } 
     else{
       c.Continue();
       for (int j=0; j<sz; j++){
	 symmP._x(j) = c.solution[j];
       }
       for (int j=0; j<sz_para; j++)
	 symmP._para(j) = c.solution[sz+j];

       //Print Energy to output
       ofstream outFileAp("Energy.txt", ofstream::app);

       outFileAp.setf(ios_base::scientific);
       double enTot = symmP.f();
       outFileAp << symmP._para(0) << " \t " << symmP._para(2)<< "\t" << enTot << "\t" << endl;
       outFileAp.close();

       // Print to VTK file
       ostringstream toStringNP, toStringCntr; 
       toStringNP << NP;
       toStringCntr << i;
       string temp("Icosa");
       temp = temp + "_N_" + toStringNP.str() +  "_" + toStringCntr.str();
     
       prob.printToVTK(temp);
       i++;
     }

   }while(key!='x');
   
   //MultiMin M("lbfgs", &prob);
   //M._LBFGSB_Initialize();
  // // // Set bounds
  //   for (int j=0; j< NP-1; j++) {
  //     M._l(Ntot + 2*j) = 0.00001;
  //     M._u(Ntot + 2*j) = PI-0.00001;
  //     M._l(Ntot + 2*j + 1) = -PI; //0
  //     M._u(Ntot + 2*j + 1) = PI;//2*PI;
  //     M._nbd[Ntot + 2*j] = 2;
  //     M._nbd[Ntot + 2*j+1] = 2;        
  //   }
   //M.LBFGSB_Solve();
   
  return 0;
}
