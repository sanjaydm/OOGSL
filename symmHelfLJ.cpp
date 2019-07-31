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

int main(int argc, char** argv){
  // Parameters of the model
  double re =   2.0; //2*sqrt(4.0/NP)*(1-0.5); //1.303648e-01; // Start r_e //pow(2.0, 1./6.)*0.08757413;//
  double k = 1.0;
  double eps = 1;

  int N = 6;
  int NP = 2;
  int Ntot= (N+1)*(N+1);
  
  Vector in(Ntot+2*NP-2);     // -2 because one particle is fixed
  Vector para(3);
  para(0) = k;   //bending stiffness kappa
  para(1) = eps;  //eps
  para(2) = re; //rm = equilibirum distance

  Vector discPara(2); //Discretization parameters
  discPara(0) = N; //Lmax
  discPara(1) = NP; //Number of particles
  
  
  srand(time(NULL));
  for (int i=0; i<Ntot; i++){
    in(i) = 0.00*rand()/RAND_MAX;
  }
  for (int j=0; j< NP-1; j++) {
    in(Ntot + 2*j) = double(rand())/RAND_MAX*PI;
    in(Ntot + 2*j+1) = (0.5-double(rand())/RAND_MAX)*2*PI;
    //cout << in->data[Ntot + 2*j]*180/PI << " " << in->data[Ntot + 2*j + 1]*180/PI << endl;
  }

  MembLJ prob(in, para, discPara);
  prob.checkConsistency();
  //prob.fdf();
  //prob._df.print();

   // Matrix R(6,6);
   //  R <<
   //  0     <<    0  << cos(2*M_PI/3)  << -sin(2*M_PI/3)    <<    0    <<     0 <<
   //  0    <<     0  <<  sin(2*M_PI/3)  << cos(2*M_PI/3)     <<    0    <<     0 <<
   //  0    <<     0  <<       0   <<      0  << cos(2*M_PI/3)  << -sin(2*M_PI/3)  <<
   //  0    <<     0   <<      0     <<    0  <<  sin(2*M_PI/3)  << cos(2*M_PI/3) <<
   // cos(2*M_PI/3)  << -sin(2*M_PI/3)   <<      0    <<     0    <<     0   <<      0 <<
   //  sin(2*M_PI/3)  << cos(2*M_PI/3)  <<      0    <<     0     <<    0   <<      0;
    
   //  Cyclic g(3,R);
   //  Projection pj(g);
   //  cout << pj._P.rank() << endl;;
    
      
  //icosahedral

   int l = 6;

   Matrix g2(3,3);
   g2 <<  -0.809016994374947 <<  0.500000000000000 << -0.309016994374947
   << 0.500000000000000 <<  0.309016994374947 << -0.809016994374947
   << -0.309016994374948 << -0.809016994374947 << -0.500000000000000;

   Matrix g3(3,3);
   g3 <<   -0.0000  <<  -0.0000  <<  -1.0000
      <<    1.0000  <<   0.0000  <<  -0.0000
      <<      0     <<  -1.0000  <<   0.0000;

   Matrix r2 = compute_R(l,l,g2);
   Matrix r3 = compute_R(l,l,g3);

   Icosahedral ico(r2,r3);
   Projection p(ico);

   cout << p._P.rank() << endl;
   //p._P.range().print();
   // cout << ";" << endl;

   // vector<int> g3p = {8,5,6,9,7,12,2,10,11,1,4,3};
   // vector<int> g2p = {7,5,8,12,2,10,1,3,11,6,9,4};

   // Matrix g2_Mat = constructMat(g2p,g2);
   // Matrix g3_Mat = constructMat(g3p,g3);

   // Icosahedral im(g2_Mat,g3_Mat);
   // Projection pm(im);
   // cout << pm._P.rank() << endl;
   //Matrix bas = pj._P.range();
   //SymmReduced symmP(bas.T()*x,para,bas,p);


  // MultiMin M("lbfgs", &prob);
  // M._LBFGSB_Initialize();
  // M.LBFGSB_Solve();
  // return 0;
}
