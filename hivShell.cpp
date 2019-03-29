#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "VirusShellBC.h"
#include "VirusShell.h"
#include "MultiMin.h"


int main(int argc, char** argv){

  int N = 100; //Num of elements
  double a = 0; double b = 2.0; // end points of domain
  vector<double> nodes;
  vector<Vector> conn;

  // Create mesh
  cout << "\033[32m" << "Creating mesh and connectivity table ...";
  for (int i=0; i <= N; i++) {
    nodes.push_back((b-a)/N * i);
    if (i != N) {
      Vector ele(2);
      ele(0) = i; ele(1) = i+1;
      conn.push_back(ele);
      //cout << "\n" << ele(0) << " " << ele(1);
    }
  }
  cout << "done. \033[0m\n";

  
  // Create initial guess
  Vector x0(4*N+4); //u,u',v,v'
  //Vector x0(4*N+4-6); //u,u',v,v'
  Vector para(1);
  para(0) = 1;
  cout << "\033[32m" << "Initializing x0 (guess) ...";
  // for (int i=0; i< x0.size(); i++){
  //   x0(i) = 0.0;
  // }
  for (int i=0; i<=N; i++){
    x0(4*i) = 0.1;
    x0(4*i+1) = 0; //u'
    x0(4*i+2) = 0.1; //v
    x0(4*i+3) = 0.; //v'
    
  }
  cout << "done. \033[0m\n";

  VirusShell* p = new  VirusShell(x0, para, nodes, conn);
  //VirusShellBC* p = new  VirusShellBC(x0, para, nodes, conn);
  //  p->fdf();
  // p->_df.print();

  // cout << "energy = " << p->_f << endl;
  p->checkConsistency();

  // MultiMin M("lbfgs", p);
  // M._LBFGSB_Initialize();


  // // Set boundary conditions
  // M._l(0) = 0; //u_0  = 0
  // M._l(1) = 0; //u'_0 = 0
  // M._l(2) = 0; //v_0  = 0
  // //M._l(3) = 0; //v'_0 = 0
  // M._u(0) = 0; //u_0  = 0
  // M._u(1) = 0; //u'_0 = 0
  // M._u(2) = 0; //v_0  = 0
  // //M._u(3) = 0; //v'_0 = 0
  // M._nbd[0] = 2;
  // M._nbd[1] = 2; //both lower and upper bounds
  // M._nbd[2] = 2;
  // // M._nbd[3] = 2; //both lower and upper bounds

  // //M._l(4*N+3) = 0; //v'_0  = 0
  // M._l(4*N+2) = 0; //v_0 = 0
  // M._l(4*N+1) = 0; //u'_0  = 0
  // M._l(4*N  ) = 0; //u_0 = 0
  // //M._u(4*N+3) = 0; 
  // M._u(4*N+2) = 0; 
  // M._u(4*N+1) = 0; 
  // M._u(4*N  ) = 0; 
  // //M._nbd[4*N+3] = 2;
  // M._nbd[4*N+2] = 2; //both lower and upper bounds
  // M._nbd[4*N+1] = 2;
  // M._nbd[4*N  ] = 2; //both lower and upper bounds

  
  // M.LBFGSB_Solve();
  
  // // p->_df.print();
  // cout << "x = ";
  // p->_x.print();
  // p->printU();
  return 0;
}
