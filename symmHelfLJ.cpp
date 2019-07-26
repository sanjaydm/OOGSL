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
  int NP = 1;
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
    in(Ntot + 2*j) = /*PI/(NP-1)*j+0.3;*/double(rand())/RAND_MAX*PI;
    in(Ntot + 2*j+1) = /*PI/(NP-1)*j;*/(0.5-double(rand())/RAND_MAX)*2*PI;
    //cout << in->data[Ntot + 2*j]*180/PI << " " << in->data[Ntot + 2*j + 1]*180/PI << endl;
  }

  MembLJ prob(in, para, discPara);
  prob.fdf();
  prob._df.print();
  // Matrix mat(2,3);
  // mat << 1 << 2 << 3 <<4 << 5 << 6;
  // Vector v = mat.row(1);
  // v.print();
  // cout << v(1) << endl;
  return 0;
}
