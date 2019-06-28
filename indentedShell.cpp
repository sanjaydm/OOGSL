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

  int N =20; //Num of elements
  double a = 0; double b =1; // end points of domain
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
  //Vector x0(3*N+3); //u,u',v,v'
   Vector x0(3*N); 

  Vector para(8);
  double C = 1; double D = 1; double nu = 0.3;
  double R = 0.5;
  double rho = 1; double d = 1.9; double alpha = 0.52;
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
    if (i%3==2)
      x0(i) = rnd;
    else
      x0(i) = rnd;
  }
  cout << "done. \033[0m\n";

  //x0.print();
  //Indented* p = new  Indented(x0, para, nodes, conn);
  IndentedBC* p = new  IndentedBC(x0, para, nodes, conn);

  p->fdf();
  (p->vs)->_df.print();

  
  cout << "energy = " << p->_f << endl;
  p->checkConsistency();

  
  MultiRoot rt("generic", p);
  rt._GSLRoot_Initialize();
  rt.GSLRoot_Solve();
  (p->vs)->_x.print();
  p->writeMesh("mesh_run.txt");
  (p->vs)->writeSolution("solution_run.py");
  
    
  /*
  MultiMin M("lbfgs", p);
  M._tol = 1e-4;
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
