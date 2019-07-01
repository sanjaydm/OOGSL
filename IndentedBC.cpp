#include "IndentedBC.h"
#include "Indented.h"

void IndentedBC :: fdf(){
  // Prescribed BC: Cantilever ends
  int N = _x.size();
  _f = 0; // reset energy
  _df.setZero(); //reset residue to zero
  int numEle = vs->_conn.size();
  double slope = (vs->_nodes[1]-vs->_nodes[0])/2;
  /// x = [x0'f0 x1 x1' f1 x2 x2' f2 .... xn-1 xn-1' fn-1 xn   fn]
  //       0  1  2   3  4  5   6  7 ....  3n-4 3n-3  3n-2 3n-1 3n

  //vsx = [x0 x0' f0 x1 x1' f1 .... xn   xn'   fn]
  //        0  1   2  3  4   5  .... 3n 3n+1   3n+2
  // //Boundary conditions applied to "Indented dofs"
  vs->_x(0) = 0; //v(0) = 0
  vs->_x(1) = 0.75*_x(1)-0.5*_x(2); // 
  vs->_x(4) = _x(2); // v3 = v3      v''(0)=0
  vs->_x(2) = _x(1); //f0
  // ------------------
  vs->_x(3*numEle) = _x(3*numEle-1);
  double u = (vs->_d - vs->_r - 1);
  double T = 1;
  vs->_x(3*numEle+1) = slope*(T - vs->_nu*u);  
  vs->_x(3*numEle+2) = _x(3*numEle);
  for (int i=2; i< N-1; i++){
    // pass dofs to Indented's
    vs->_x(i+1) = _x(i);
  }
  
  vs->fdf();
  if (_fFlag) {
    // retrieve energy from virusShell's
    _f = vs->_f;
  }
  if (_dfFlag){
    // retrieve residue from virusShell's
    for (int i=2; i< N-1; i++){
      _df(i) = vs->_df(i+1);
    }
    _df(0) = vs->_df(2);
    _df(1) = vs->_df(3)+0.75*vs->_df(1);
    _df(2) = -0.5*vs->_df(1)+vs->_df(4);
    _df(3*numEle-1) = vs->_df(3*numEle) ;
    _df(3*numEle) = vs->_df(3*numEle+2) ;
  }

}

/*
void IndentedBC :: fdf(){
  // Prescribed BC: Cantilever ends
  int N = _x.size();
  _f = 0; // reset energy
  _df.setZero(); //reset residue to zero
  int numEle = vs->_conn.size(); 

  // //Boundary conditions applied to "Indented dofs"
  vs->_x(0) = 0; //u(0) = 0
  vs->_x(1) = 0; //u'(0) = 0
  vs->_x(2) = 0; //v(0) = 0
 
  vs->_x(4*numEle) = 0;   //u(L) = 0;
  vs->_x(4*numEle+1) = 0; //u'(L) = 0;
  vs->_x(4*numEle+2) = 0; //v(L) = 0
  
  for (int i=0; i<N-1; i++){
    // pass dofs to Indented's
    vs->_x(i+3) = _x(i);
  }
  vs->_x(N+5) = _x(N-1);
  
  vs->fdf();
  if (_fFlag) {
    // retrieve energy from virusShell's
    _f = vs->_f;
  }
  if (_dfFlag){
    // retrieve residue from virusShell's
    for (int i=0; i<N-1; i++){
      _df(i) = vs->_df(i+3);
    }
    _df(N-1) = vs->_df(N+5);
  }
}
*/

// void Indented::writeSolution(string filename){
//   ofstream myfile;
//   myfile.open (filename.c_str());
//   cout << "Writing solution to a file.\n";
//   //myfile << setprecision(15);
//   for (auto i = 1; i <= _x.size(); i++) {
//     myfile << _x(i-1);
//     if (i % 2 == 0)
//       myfile << "\n";
//     else
//       myfile << "\t";
// }
//   myfile.close();
// }
// void IndentedBC::writeSolution(string filename){
//   ofstream myfile;
//   myfile.open (filename.c_str());
//   cout << "Writing solution to a file.\n";
//   myfile << setprecision(15);
//   myfile << "x = [";
//   for (auto i = 1; i <= _x.size(); i++) {
//     if (i % 2 == 1)
//       myfile << "\n[" << _x(i-1) << ",\t";
//     else
//       myfile << _x(i-1) << "],";
// }
//   myfile << "]\n";
//   myfile << "from matplotlib import pyplot as plt\n";
//   myfile << "from numpy import *\n";
//   myfile << "x = matrix(x)\n";
//   myfile << "plt.plot(x[0:len(x):2],'*-')\n";
//   myfile << "plt.show()\n";
//   myfile.close();
// }

