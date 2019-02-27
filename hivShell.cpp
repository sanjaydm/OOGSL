#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "VirusShell.h"
#include "MultiMin.h"


int main(int argc, char** argv){

  int N = 3; //Num of dof
  double a = 0; double b = 5; // end points of domain
  vector<double> nodes;
  vector<Vector> conn;

  // Create mesh
  double ba = (b-a);
  cout << "\033[32m" << "Creating mesh and connectivity table ...";
  for (int i=0; i <= N; i++) {
    nodes.push_back((b-a)/N * i);
    if (i != N) {
      Vector ele(2);
      ele(0) = i; ele(1) = i+1;
    }
  }
  cout << "done. \033[0m\n";

  
  // Create initial guess
  Vector x0(N);
  Vector para(1);
  para(0) = 1;
  cout << "\033[32m" << "Initializing x0 (guess) ...";
  for (int i=0; i<=N; i++){
    x0(i) = nodes[i];
  }
  cout << "done. \033[0m\n";
  
  VirusShell* p = new  VirusShell(x0, para, nodes, conn);

  MultiMin M("lbfgs", p);
  M._LBFGSB_Initialize();
  M.LBFGSB_Solve();
  // // M._GSLMin_Initialize();
  // // M.GSLMin_Solve();
  // p->_x.print();

  return 0;
}
