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


using namespace std;
int main(int argc, char** argv){
  // Number of modes and particles  
  int N = 3;
  int NP = 10;
  int Ntot= (N+1)*(N+1);
  Vector in(Ntot+2*NP); 
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
    in(Ntot + 2*j) =  double(rand())/RAND_MAX*PI;
    in(Ntot + 2*j+1) = (0.5-double(rand())/RAND_MAX)*2*PI;
  }

  // Construct the model
  MembLJ prob(in, para, discPara);
  //prob.checkConsistency();
  
//
//  // Generators of Icosahedral Group
//    Matrix g2(3,3);
//   g2 <<  -0.809016994374947 <<  0.500000000000000 << -0.309016994374947
//   << 0.500000000000000 <<  0.309016994374947 << -0.809016994374947
//   << -0.309016994374948 << -0.809016994374947 << -0.500000000000000;
//
//   Matrix g3(3,3);
//   g3 <<   -0.0000  <<  -0.0000  <<  -1.0000
//      <<    1.0000  <<   0.0000  <<  -0.0000
//      <<      0     <<  -1.0000  <<   0.0000;
//
//   // Construct Icosahedral group-rep on spherical harmonics
//   Matrix r2 = compute_R(0,N,g2);
//   Matrix r3 = compute_R(0,N,g3);
//
//   // Construct permutation representation
//   vector<int> g3p = {8,5,6,9,7,12,2,10,11,1,4,3};
//   vector<int> g2p = {7,5,8,12,2,10,1,3,11,6,9,4};
//
//   Matrix g2_Mat = constructMat(g2p,g2);
//   Matrix g3_Mat = constructMat(g3p,g3);

    
//    // Generators of C3
//    Matrix r(3,3);
//    r << cos(2*PI/3) << -sin(2*PI/3) << 0
//      << sin(2*PI/3) << cos(2*PI/3) << 0
//      << 0 << 0 << 1;
//
//    // Construct cyclic group-rep on spherical harmonics
//    Matrix rr = compute_R(0,N,r);
//
//
//    // Construct permutation representation
//    vector<int> rp = {2,3,1,5,6,4};
//
//    Matrix rp_Mat = constructMat(rp,r);
//
//
//    Cyclic C_sh(3,rr);
//    Cyclic C_p(3,rp_Mat);
//
//    Projection p_sp(C_sh);
//    Projection p_p(C_p);
//
//    Matrix rg_sh = p_sp._P.range(); // bases for spherical coordinates
//    Matrix rg_p = p_p._P.range(); // bases for particles in cartesian
//
//    int n_sh = rg_sh.rank(); // number of bases for spherical coord.
//    int rg_prk = rg_p.rank(); // number of bases for particles in cartesian
//
//
//    cout << rg_sh.rank() << endl;
//    rg_sh.print();
//
//
//    cout << rg_p.rank() << endl;
//    rg_p.print();
//
//    cout << " " << endl;
//
//
    
    // Generators of D5
    Matrix r(3,3);
    r << cos(2*PI/5) << -sin(2*PI/5) << 0
      << sin(2*PI/5) << cos(2*PI/5) << 0
      << 0 << 0 << 1;
    
    Matrix s(3,3);
    s << 1 << 0 << 0
      << 0 <<-1 << 0
      << 0 << 0 << -1;
    
    // Construct D5 group-rep on spherical harmonics
    Matrix rr = compute_R(0,N,r);
    Matrix ss = compute_R(0,N,s);
    
    vector<int> rp = {2,3,4,5,1,7,8,9,10,6};
    vector<int> sp = {6,10,9,8,7,1,5,4,3,2};
    
    // Construct D5 group-rep on particles
    Matrix rp_Mat = constructMat(rp,r);
    Matrix sp_Mat = constructMat(sp,s);
    
    // group construction
    D5 D5_p(sp_Mat,rp_Mat);
    D5 D5_sh(ss,rr);
    
    // Projection operator
    Projection p_sh(D5_sh);
    Projection p_p(D5_p);
    
    // bases in cartesian
    Matrix rg_sh = p_sh._P.range()
    Matrix rg_p = p_p._P.range()
    
    
    
    
    // bruteforce a bases matrix for particles
    Matrix rgp_bf(NP*3,rg_prk);
    for (int j=0; j<rg_prk/2; j++){
        for (int i=0; i<NP*3; i++){
            rgp_bf(i,2*j) = rg_p(i,2*j) + rg_p(i,2*j+1);
            rgp_bf(i,2*j+1) = rg_p(i,2*j) - rg_p(i,2*j+1);
        }
    }
    
    return 0;
    
    
    
    
    // Convert rg_p to spherical coordinates
    
    Matrix rg_sc(NP*2,rg_prk);
    for (int j=0; j<rg_prk; j++){
        for (int ip = 0; ip < NP; ip++){
            double z = rgp_bf(ip*3+2,j);
            double x = rgp_bf(ip*3,j);
            double y = rgp_bf(ip*3+1,j);
            double r = sqrt(x*x+y*y+z*z);
            cout << ip << " " << r << endl;
            rg_sc(ip*2,j) = acos(z/r);
            rg_sc(ip*2+1,j) = atan2(y, x);
        }
    }
    
    rg_sc.print();

    int n_p = rg_sc.rank();
    cout << n_p << endl;
    cout << rg_sc.size()[0] << " " << rg_sc.size()[1]<< endl;
    
    
    Matrix bs(rg_sh.size()[0]+rg_sc.size()[0],rg_sh.size()[1]+rg_sc.size()[1]);
    for (int i=0; i<rg_sh.size()[0]+rg_sc.size()[0]; i++){
        for (int j=0; j<n_sh+n_p; j++){
            if (i<rg_sh.size()[0] && j<rg_sh.size()[1]){
                bs(i,j) = rg_sh(i,j);
            }
            if (i>=rg_sh.size()[0] && j>= rg_sh.size()[1]){
                bs(i,j) = rg_sc(i-rg_sh.size()[0],j-rg_sh.size()[1]);
            }
        }
    }
    
    bs.print();
    cout << bs.size()[0] << " " << bs.size()[1] << endl;
 
    

//
//   // Icosahedral group-rep on spherical harmonics+particles
//   int GN = r2.size()[0] +g2_Mat.size()[0];
//   Matrix G2(GN, GN);
//   Matrix G3(GN, GN);
//   for (int i=0; i <GN; i++){
//     for (int j=0; j<GN; j++){
//       if(i<r2.size()[0] && j < r2.size()[0]){
//	 G2(i,j) = r2(i,j);
//	 G3(i,j) = r3(i,j);
//       }
//       if(i>=r2.size()[0] && j >= r2.size()[0]){
//	 G2(i,j) = g2_Mat(i-r2.size()[0],j-r2.size()[0]);
//	 G3(i,j) = g3_Mat(i-r2.size()[0],j-r2.size()[0]);
//       }
//     }
//   }
//
//   Icosahedral im(G2,G3);
    

//    // Cyclic group-rep on spherical harmonics+particles
//    int GN = rr.size()[0] +rp_Mat.size()[0];
//    Matrix RR(GN, GN);
//    for (int i=0; i <GN; i++){
//      for (int j=0; j<GN; j++){
//        if(i<rr.size()[0] && j < rr.size()[0]){
//      RR(i,j) = rr(i,j);
//        }
//        if(i>=rr.size()[0] && j >= rr.size()[0]){
//      RR(i,j) = rp_Mat(i-rr.size()[0],j-rr.size()[0]);
//        }
//      }
//    }
//
//    Cyclic C3(3,RR);


//   // Construct projection operators
//   Projection pm(C3);
//   Matrix rg = pm._P.range();

    
   // cout << "projection size" << rg.size()[0]<< "  " << rg.size()[1] << endl;
//   int nn = 7;


    
    
//   // Basis vectors adapted for spherical coordinates
//   Matrix bas((N+1)*(N+1)+NP*2, nn);
//   for (int j=0; j<=nn-1; j++){
//     for (int i=0; i < (N+1)*(N+1)+NP*2; i ++) {
//       if (i < (N+1)*(N+1) && j<nn-6){
//	 bas(i,j)=rg(i,j);
//       }
//       if (j>=nn-6 && i >= (N+1)*(N+1)) {
//	 double z = rg(3*((i-(N+1)*(N+1))/2) + 2 + (N+1)*(N+1),j);
//	 double x = rg(3*((i-(N+1)*(N+1))/2) + (N+1)*(N+1),j);
//	 double y = rg(3*((i-(N+1)*(N+1))/2) + 1 + (N+1)*(N+1),j);
//	 double r = sqrt(x*x+y*y+z*z);
//	 //cout << 3*((i-(N+1)*(N+1))/2) <<" :" << x << ";" << y << ";" << z <<endl;
//
//	 if ((i-(N+1)*(N+1))%2==0)
//	   bas(i,j) = acos(z/r);
//
//	 else
//	   bas(i,j) = atan2(y, x);
//       }
//     }
//   }

    
   // cout << bas.size()[0] << bas.size()[1]<< endl;

   // Construct symmetry reduced problem
   Vector x0 = bs.T()*in;
    bs.print();



   // cout << "initial vector size" << x0.size() << endl;
    

//   x0(nn-1) = 1;
//   x0(0) = 1;
    x0(0) = 0.0;
    x0(1) = 0.0;
    x0(2) = 0.0;
    x0(3) = 0.0;
    x0(4) = 0.0;
    x0(5) = 0.0;

    for (int i = 0; i<rg_sc.size()[1]; i++){
        x0(i+rg_sh.size()[1]) = 0.1;
    }
    
    x0.print();

    
   SymmReduced symmP(x0,para,bs,&prob);
   // symmP.checkConsistency();
    symmP.fdf();
    symmP._df.print();
    


   MultiMin M("lbfgs", &symmP);
   M._LBFGSB_Initialize();
   // Optionally constrain the particle coordinates to lie within spherical coordinate domain
   // M._l(2) = 1;
   // M._u(2) = 1;
   // M._l(2) = 1;
   // M._u(2) = 1;
   // M._nbd[2] = 2;
   // M._nbd[2] = 2;  

   //Print Energy to output
   ofstream outFileAp("Energy.txt", ofstream::app);
   outFileAp << "kappa\t rm\t Total Energy\t Projg\n";
   outFileAp.close();
   for (int i=0; i<30; i++) {
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
     string temp("C3");
     temp = temp + "_N_" + toStringNP.str() +  "_" + toStringCntr.str();
     prob.printToVTK(temp);

   }
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
