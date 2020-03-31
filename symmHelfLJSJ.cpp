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
#include <string>

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
  int N = 5;
  int NP = 8;
  int Ntot= (N+1)*(N+1);
  int x = 4;
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

//   // Construct Icosahedral group-rep on spherical harmonics
//   Matrix r2 = compute_R(0,N,g2);
//   Matrix r3 = compute_R(0,N,g3);

//   // Construct permutation representation
//   vector<int> g3p = {8,5,6,9,7,12,2,10,11,1,4,3};
//   vector<int> g2p = {7,5,8,12,2,10,1,3,11,6,9,4};
//
//   Matrix g2_Mat = constructMat(g2p,g2);
//   Matrix g3_Mat = constructMat(g3p,g3);

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

//   Icosahedral im(G2,G3);
//
//   // Construct projection operators
//   Projection pm(im);
//   Matrix rg = pm._P.range();
   //rg.print();
   
    
    
    
//    // Generators of D9
//    Matrix r(3,3);
//    r << cos(2*PI/9) << -sin(2*PI/9) << 0
//      << sin(2*PI/9) << cos(2*PI/9) << 0
//      << 0 << 0 << 1;

    
//    // Generators of D5
//    Matrix r(3,3);
//      r << cos(2*PI/5) << -sin(2*PI/5) << 0
//        << sin(2*PI/5) << cos(2*PI/5) << 0
//        << 0 << 0 << 1;
    
    // Generators of Dx
    Matrix r(3,3);
      r << cos(2*PI/x) << -sin(2*PI/x) << 0
        << sin(2*PI/x) << cos(2*PI/x) << 0
        << 0 << 0 << 1;

    
    Matrix s(3,3);
      s << 1 << 0 << 0
        << 0 <<-1 << 0
        << 0 << 0 << -1;
    
    
    // Construct group-rep on spherical harmonics
    Matrix rr = compute_R(0,N,r);
    Matrix ss = compute_R(0,N,s);
    
    
    //candidate for D5
    vector<int> rp;
    vector<int> sp;
    
    
    
    ifstream myfiler;
    if (NP == 6){
        
        if (x == 3){
            myfiler.open("/Users/clarec/Documents/MATLAB/clare-4/D3_NP6_r_perm.txt");//, ios_base::in
            int y=0;
            for (int i=0; i< NP; i++){
                myfiler >> y;
                cout << y << endl;
                rp.push_back(y);
            }
            cout << "---------------"<< endl;
            
            ifstream myfiles;
            myfiles.open("/Users/clarec/Documents/MATLAB/clare-4/D3_NP6_s_perm.txt");//, ios_base::in
            for (int i=0; i< NP; i++){
                myfiles >> y;
                cout << y << endl;
                sp.push_back(y);
            }
        }
        if (x == 4){
            myfiler.open("/Users/clarec/Documents/MATLAB/clare-4/D4_NP6_r_perm.txt");//, ios_base::in
            int y=0;
            for (int i=0; i< NP; i++){
                myfiler >> y;
                cout << y << endl;
                rp.push_back(y);
            }
            cout << "---------------"<< endl;
            
            ifstream myfiles;
            myfiles.open("/Users/clarec/Documents/MATLAB/clare-4/D4_NP6_s_perm.txt");//, ios_base::in
            for (int i=0; i< NP; i++){
                myfiles >> y;
                cout << y << endl;
                sp.push_back(y);
            }
        }
    
            if (x == 6){
                myfiler.open("/Users/clarec/Documents/MATLAB/clare-4/D6_NP6_r_perm.txt");//, ios_base::in
                int y=0;
                for (int i=0; i< NP; i++){
                    myfiler >> y;
                    cout << y << endl;
                    rp.push_back(y);
                }
                cout << "---------------"<< endl;
                
                ifstream myfiles;
                myfiles.open("/Users/clarec/Documents/MATLAB/clare-4/D3_NP6_s_perm.txt");//, ios_base::in
                for (int i=0; i< NP; i++){
                    myfiles >> y;
                    cout << y << endl;
                    sp.push_back(y);
                }
            }
    }
    if (NP == 8){
        
        if (x == 3){
            myfiler.open("/Users/clarec/Documents/MATLAB/clare-4/D3_NP8_r_perm.txt");//, ios_base::in
            int y=0;
            for (int i=0; i< NP; i++){
                myfiler >> y;
                cout << y << endl;
                rp.push_back(y);
            }
            cout << "---------------"<< endl;
            
            ifstream myfiles;
            myfiles.open("/Users/clarec/Documents/MATLAB/clare-4/D3_NP8_s_perm.txt");//, ios_base::in
            for (int i=0; i< NP; i++){
                myfiles >> y;
                cout << y << endl;
                sp.push_back(y);
            }
        }
        if (x == 4){
            myfiler.open("/Users/clarec/Documents/MATLAB/clare-4/D4_NP8_r_perm.txt");//, ios_base::in
            int y=0;
            for (int i=0; i< NP; i++){
                myfiler >> y;
                cout << y << endl;
                rp.push_back(y);
            }
            cout << "---------------"<< endl;
            
            ifstream myfiles;
            myfiles.open("/Users/clarec/Documents/MATLAB/clare-4/D4_NP8_s_perm.txt");//, ios_base::in
            for (int i=0; i< NP; i++){
                myfiles >> y;
                cout << y << endl;
                sp.push_back(y);
            }
        }
        if (x == 6){
            myfiler.open("/Users/clarec/Documents/MATLAB/clare-4/D6_NP8_r_perm.txt");//, ios_base::in
            int y=0;
            for (int i=0; i< NP; i++){
                myfiler >> y;
                cout << y << endl;
                rp.push_back(y);
            }
            cout << "---------------"<< endl;
            
            ifstream myfiles;
            myfiles.open("/Users/clarec/Documents/MATLAB/clare-4/D6_NP8_s_perm.txt");//, ios_base::in
            for (int i=0; i< NP; i++){
                myfiles >> y;
                cout << y << endl;
                sp.push_back(y);
            }
        }
        if (x == 8){
            myfiler.open("/Users/clarec/Documents/MATLAB/clare-4/D8_NP8_r_perm.txt");//, ios_base::in
            int y=0;
            for (int i=0; i< NP; i++){
                myfiler >> y;
                cout << y << endl;
                rp.push_back(y);
            }
            cout << "---------------"<< endl;
            
            ifstream myfiles;
            myfiles.open("/Users/clarec/Documents/MATLAB/clare-4/D8_NP8_s_perm.txt");//, ios_base::in
            for (int i=0; i< NP; i++){
                myfiles >> y;
                cout << y << endl;
                sp.push_back(y);
            }
        }
        
    }
        
        



//    if (x == 2) {
//        rp = {2,1,4,3,5}; //
//        sp = {3,4,1,2,5};
//    }
//    if (x == 3){
//        rp = {2,3,1,4,5};
//        sp = {3,2,1,5,4};
//    }
//    if (x == 4){
//        rp = {2,3,4,1,5};
//        sp = {2,1,4,3,5};
//    }
//    if (x == 5){
//        rp = {2,3,4,5,1};
//        sp = {4,3,2,1,5};
//    }
//    if (x == 6){
//        fstream myfiler("~/Documents/MATLAB/clare-4/r_perm.txt", ios_base::in);
//        int y;
//        for (int i=0; i< NP; i++){
//            myfiler >> y;
//            rp.push_back(y);
//        }
//
//        fstream myfiles("~/Documents/MATLAB/clare-4/s_perm.txt", ios_base::in);
//        for (int i=0; i< NP; i++){
//            myfiles >> y;
//            sp.push_back(y);
//        }
//       }
    
    

    
//    fstream myfiler("~/Documents/MATLAB/clare-4/r_perm.txt", ios_base::in);
//    int x;
//    for (int i=0; i< NP; i++){
//        myfiler >> x;
//        rp.push_back(x);
//    }
//
//
//    fstream myfiles("~/Documents/MATLAB/clare-4/s_perm.txt", ios_base::in);
//    for (int i=0; i< NP; i++){
//        myfiles >> x;
//        sp.push_back(x);
//    }

    Projection p_sh;
    Projection p_p;
    
    
    // Construct group-rep on particles
    Matrix rp_Mat = constructMat(rp,r);
    Matrix sp_Mat = constructMat(sp,s);
    
    // group construction
    if (x == 2) {
        D2 D2_p(sp_Mat,rp_Mat);
        D2 D2_sh(ss,rr);
        
        p_sh = Projection(D2_sh);
        p_p = Projection(D2_p);
    }
    if (x == 3){
        D3 D3_p(sp_Mat,rp_Mat);
        D3 D3_sh(ss,rr);
         
        p_sh = Projection(D3_sh);
        p_p = Projection(D3_p);
    }
    if (x == 4){
        D4 D4_p(sp_Mat,rp_Mat);
        D4 D4_sh(ss,rr);
        
        p_sh = Projection(D4_sh);
        p_p = Projection(D4_p);
    }
    if (x == 5){
        D5 D5_p(sp_Mat,rp_Mat);
        D5 D5_sh(ss,rr);
        
        p_sh = Projection(D5_sh);
        p_p = Projection(D5_p);
    }
    if (x == 6){
        D6 D6_p(sp_Mat,rp_Mat);
        D6 D6_sh(ss,rr);
        
        p_sh = Projection(D6_sh);
        p_p = Projection(D6_p);
    }

    
    // bases in cartesian
    Matrix rg_sh = p_sh._P.range();
    Matrix rg_p = p_p._P.range();
    
    int n_sh = rg_sh.rank(); // number of bases for spherical coord.
    int rg_prk = rg_p.rank(); // number of bases for particles in cartesian
    int n_p = rg_p.rank();

    // Combine the bases for particles and spherical harmonics into one
    Matrix bs(rg_sh.size()[0]+rg_p.size()[0],rg_sh.size()[1]+rg_p.size()[1]);
    for (int i=0; i<rg_sh.size()[0]+rg_p.size()[0]; i++){
        for (int j=0; j<n_sh+n_p; j++){
            if (i<rg_sh.size()[0] && j<rg_sh.size()[1]){
                bs(i,j) = rg_sh(i,j);
            }
            if (i>=rg_sh.size()[0] && j>= rg_sh.size()[1]){
                bs(i,j) = rg_p(i-rg_sh.size()[0],j-rg_sh.size()[1]);
            }
        }
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
//   int nn = 1;
   
   // Construct symmetry reduced problem
//   Vector x0 = rg.T()*in;
   //x0.print();
    bs.print();
    
   Vector x0 = bs.T()*in;
    int nn = x0.size();
    
    srand(time(NULL));
    for (int j=0; j< nn; j++) {
        x0(j) =  double(rand())/RAND_MAX;
    }
   
    
   symmP = SymmReduced(x0,para,bs,&prob);
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



     // Print to VTK file
     ostringstream toStringNP, toStringN, toStringCntr;
     toStringNP << NP;
     toStringN << N;
     toStringCntr << i;
     string temp(to_string(x));
     temp = "/Users/clarec/Desktop/Energy_Data/D" + temp + "_N_" + toStringNP.str() + "_L_"+ toStringN.str() ;
     string tempp = temp + "_" + toStringCntr.str();
     prob.printToVTK(tempp);
       
       
       //Print Energy to output
       ofstream outFileAp(temp + "_Energy.txt", ofstream::out);
       outFileAp.setf(ios_base::scientific);
       double projg = M._dsave[12];
       double enTot = symmP._f;
       outFileAp << symmP._para(0) << " \t " << symmP._para(2)<< "\t" << enTot << "\t" <<  projg << endl;
       outFileAp.close();


   }
    

    
    
    
    
    
    // Print to Data.txt
    ofstream outFileApp("Data.txt", ofstream::out);
    outFileApp.setf(ios_base::scientific);
    for (int i=0; i< symmP._fullX.size(); i++){
        outFileApp << symmP._fullX(i) << " \t ";
    }
    for (int i=0; i< symmP._para.size(); i++){
        outFileApp << symmP._para(i) << " \t ";
    }
    outFileApp << "\n";
    // outFileApp << symmP._f << endl;
    outFileApp.close();

    
    
    
    
    
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


         

       // Print to VTK file
       ostringstream toStringNP, toStringN, toStringCntr;
       toStringNP << NP;
       toStringN << N;
       toStringCntr << i;
       string temp(to_string(x));
       temp = "/Users/clarec/Desktop/Energy_Data/D" + temp + "_N_" + toStringNP.str() + "_L_"+ toStringN.str() ;
       string tempp = temp + "_" + toStringCntr.str();
         
        prob.printToVTK(tempp);
        i++;
         
         
         
         //Print Energy to output
         ofstream outFileAp(temp + "_Energy.txt", ofstream::app);

         outFileAp.setf(ios_base::scientific);
         double enTot = symmP.f();
         outFileAp << symmP._para(0) << " \t " << symmP._para(2)<< "\t" << enTot << "\t" << endl;
         outFileAp.close();
         
   
         
         
         
         ofstream outFileApp("Data.txt", ofstream::app);
          outFileApp.setf(ios_base::scientific);
         for (int i=0; i< symmP._fullX.size(); i++){
             outFileApp << symmP._fullX(i) << " \t ";
         }
         for (int i=0; i< symmP._para.size(); i++){
             outFileApp << symmP._para(i) << " \t ";
         }
         // outFileApp << symmP._f << endl;
         outFileApp << endl;
         outFileApp.close();
         
         
         
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
