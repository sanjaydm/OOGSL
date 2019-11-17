#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "Indented.h"
#include "IndentedBC.h"
#include "MultiMin.h"
#include "MultiRoot.h"
#include <cstdlib>
#include <ctime>

int main(int argc, char** argv){

  int N =30; //Num of elements
  double a = 0; double b =1; // end points of domain
  vector<double> nodes;
  vector<Vector> conn;

  // Create mesh
  cout << "\033[32m" << "Creating mesh and connectivity table ...";
  for (int i=0; i <= N; i++) {
    nodes.push_back((b-a)/N * i + a );
    //cout << (b-a)/N * i + a << endl;
    if (i != N) {
      Vector ele(2);
      ele(0) = i; ele(1) = i+1;
      conn.push_back(ele);
      //cout << "\n" << ele(0) << " " << ele(1);
    }
  }
  cout << "done. \033[0m\n";

  
  // Create initial guess
  //Vector x0(5*N+5); //u,u',v,v'
  Vector x0(5*N-2); //u,u',v,v'

  Vector para(8);
  double C = 1.0; double D = 1.0; double nu = 0.3;
  double R = 1.0;
  double rho = 1; double d = 0.7+rho; double alpha = 0.52;
  para(0) = C; 
  para(1) = D;
  para(2) = nu;
  para(3) = rho;
  para(4) = d;
  para(5) = alpha;
  cout << "\033[32m" << "Initializing x0 (guess) ...";

  /* initialize random seed: */
  srand (time(NULL));
  for (int i=0; i< x0.size(); i++){
    double rnd = double(rand())/RAND_MAX; 
    if (i%5==0)
      x0(i) = .1*rnd;;
    if(i%5 == 1)
      x0(i) = .1*rnd;
    else
      x0(i) = .1*rnd;
  }
  cout << "done. \033[0m\n";

  //x0.print();
  //Indented* p = new  Indented(x0, para, nodes, conn);
  IndentedBC* p = new  IndentedBC(x0, para, nodes, conn);

  /*
  p->f();
  cout << "energy = " << p->_f << endl;
  (p)->df();
  (p)->_df.print();

  (p->vs)->df();
  (p->vs)->_df.print();
  return 0;
  */
  cout << "energy = " << p->_f << endl;
  p->checkConsistency();

  //return 0;

  
  MultiRoot rt("generic", p);
  rt._GSLRoot_Initialize();
  rt.GSLRoot_Solve();
  (p->vs)->_x.print();
  //cout << (p->vs)->_x.size() << endl;
  //p->writeMesh("mesh_run.txt");
  (p->vs)->writeSolution("solution_run.py");
  cout << "Norm = " << (p->_df).norm()<< endl;
  
  
  /*
  MultiMin M("lbfgs", p);
  M._tol = 1e-1;
  M._LBFGSB_Initialize();
  M.LBFGSB_Solve();
  
  (p->vs)->_x.print();
  p->writeMesh("mesh_run.txt");
  (p->vs)->writeSolution("solution_run.py");
  cout << "size of vs = " << (p->vs)->_x.size() << endl;
  cout << "size of s = " << (p->_x.size()) << endl;
  */
  return 0;
  
}
