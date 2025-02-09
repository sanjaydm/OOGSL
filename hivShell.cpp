#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "VirusShellBC.h"
#include "VirusShell.h"
#include "MultiMin.h"
#include <cstdlib>
#include <ctime>

int main(int argc, char** argv){

  int N = 30; //Num of elements
  double a = 5; double b = 15.0; // end points of domain
  vector<double> nodes;
  vector<Vector> conn;

  // Create mesh
  cout << "\033[32m" << "Creating mesh and connectivity table ...";
  for (int i=0; i <= N; i++) {
    nodes.push_back((b-a)/N * i + a );
    if (i != N) {
      Vector ele(2);
      ele(0) = i; ele(1) = i+1;
      conn.push_back(ele);
      //cout << "\n" << ele(0) << " " << ele(1);
    }
  }
  cout << "done. \033[0m\n";

  
  // Create initial guess
  //Vector x0(4*N+4); //u,u',v,v'
  //Vector x0(4*N+4-6); //u,u',v,v' for cantilever bc

  Vector x0(4*N-4); //u,u',v,v' for indent-free edge
  Vector para(8);
  double C = 1; double D = 1; double nu = 0.3;
  double R = 1;
  double rho = 1; double d = 1.9; double alpha = 0.52;
  double u0 = d - rho*cos(alpha)-R;
  double v0 = rho*sin(alpha)-a;
  para(0) = C; 
  para(1) = D;
  para(2) = nu;
  para(3) = rho;
  para(4) = d;
  para(5) = alpha;
  para(6) = u0;
  para(7) = v0;
  cout << "\033[32m" << "Initializing x0 (guess) ...";

  /* initialize random seed: */
  srand (time(NULL));
  for (int i=0; i< x0.size(); i++){
    double rnd = double(rand())/RAND_MAX; 
    x0(i) = 0.0+1*rnd;
  }
  // double rnd;
  // for (int i=0; i<=N; i++){
  //   rnd = double(rand())/RAND_MAX; 
  //   x0(4*i) = 0 + 0.1*rnd ;
  //   rnd = double(rand())/RAND_MAX; 
  //   x0(4*i + 1) = 0 + 0.1*rnd ; 
  //   rnd = double(rand())/RAND_MAX;
  //   x0(4*i + 2) = 0 + 0.1*rnd ; 
  //   rnd = double(rand())/RAND_MAX;
  //   x0(4*i + 3) = 0 + 0.1*rnd ; 
    
  // }
  cout << "done. \033[0m\n";

  //x0.print();
  //VirusShell* p = new  VirusShell(x0, para, nodes, conn);
  VirusShellBC* p = new  VirusShellBC(x0, para, nodes, conn);
  //  p->fdf();
  //p->_df.print();

  //cout << "energy = " << p->_f << endl;
  p->checkConsistency();

  MultiMin M("lbfgs", p);
  M._LBFGSB_Initialize();


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

  
  M.LBFGSB_Solve();
  
  // // p->_df.print();
  // cout << "x = ";
  p->_x.print();
  p->writeMesh("mesh_run.txt");
  p->writeSolution("solution_run.py");
  // p->printU();
  return 0;
}
