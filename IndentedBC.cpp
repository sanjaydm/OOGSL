#include "IndentedBC.h"
#include "Indented.h"

void IndentedBC :: fdf(){
  // Prescribed BC: Cantilever ends
  int N = _x.size();
  _f = 0; // reset energy
  _df.setZero(); //reset residue to zero
  int numEle = vs->_conn.size();
  double slope = (vs->_nodes[1]-vs->_nodes[0])/2;
  /// x = [X X   X   X  f0 u1 u1' v1  v1' f1 .... fn-1  X   un'   X  X  fn]
  //                    0  1  2   3   4   5 ....  5n-5     5n-4        5n-3

  //vsx = [u0 u0' v0 v0' f0 u1 u1' v1 v1' f1 ... un   un'   vn   vn'   fn]
  //       0  1   2  3   4  5  6   7  8   9      5n   5n+1  5n+2  5n+3 5n+4
  //      
  // //Boundary conditions applied to "Indented dofs"
  double delta = (vs->_d -1);
  double theta = 1;
  double u = (delta - vs->_r);
  double T = 0; // Cannot be nonzero, there must be a balancing force on the other end

  vs->_x(0) = u;
  vs->_x(1) = 0;
  vs->_x(2) = 0; 
  vs->_x(3) = (0.75*_x(3)-0.5*_x(4)); //
  // ------------------
  for (int i=4; i<= 5*numEle-1; i++){ 
    // pass dofs to Indented's
    vs->_x(i) = _x(i-4);
  }
  double u1 = delta + vs->_r*cos(theta); //Set u1
  double v1 = -1 + vs->_r*sin(theta); //Set v1
  double fac = 1/(vs->_nu*slope);
  
  vs->_x(5*numEle) = u1;
  vs->_x(5*numEle+1) = _x(5*numEle-4);

  vs->_x(5*numEle+2) = v1;
  vs->_x(5*numEle+3) = fac*(1.5*_x(5*numEle-9) + _x(5*numEle-8) -1.5*u1 + 2*_x(5*numEle-4)); //Comes from u'' = nu/R v';
  vs->_x(5*numEle+4) = _x(5*numEle-3);

  // Compute energy
  vs->fdf();
  if (_fFlag) {
    // retrieve energy from virusShell's
    _f = vs->_f;
  }
  if (_dfFlag){
    _df(0) = vs->_df(4);
    _df(1) = vs->_df(5);
    _df(2) = vs->_df(6);
    _df(3) = vs->_df(7)+0.75*vs->_df(3);
    _df(4) = -0.5*vs->_df(3)+vs->_df(8);
    
    for (int i=5; i<=5*numEle-5; i++){ //N-3 for linear
      _df(i) = vs->_df(i+4);
    }

    _df(5*numEle-9) += fac*1.5*vs->_df(5*numEle+3);
    _df(5*numEle-8) += fac*vs->_df(5*numEle+3);
    _df(5*numEle-4) = 2*fac*vs->_df(5*numEle+3) + vs->_df(5*numEle+1); 
    _df(5*numEle-3) = vs->_df(5*numEle+4) ;
  }

}
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

