#include "VirusShellBC.h"
#include "VirusShell.h"

void VirusShellBC :: fdf(){
  // Prescribed BC: Cantilever ends
  int N = _x.size();
  _f = 0; // reset energy
  _df.setZero(); //reset residue to zero
  int numEle = vs->_conn.size();
  double slope = (vs->_nodes[1]-vs->_nodes[0])/2;
  // //Boundary conditions applied to "VirusShell dofs"
  vs->_x(0) = _u0; //u(0) = 0
  vs->_x(1) = cos(_alpha)*slope; //u'(0) = 0
  vs->_x(2) = _v0; //v(0) = 0
  vs->_x(3) = (1-sin(_alpha))*slope; //v(0) = 0

  vs->_x(4*numEle-3) = -1.5*_x(N-4);
  vs->_x(4*numEle) = 0;   //u(L) = 0;
  vs->_x(4*numEle+1) = 0; //u'(L) = 0;
  vs->_x(4*numEle+2) = 0; //v(L) = 0
  
  for (int i=0; i<= N-4; i++){
    // pass dofs to VirusShell's
    vs->_x(i+4) = _x(i);
  }
  vs->_x(4*numEle-2) = _x(N-3); //v
  vs->_x(4*numEle-1) = _x(N-2); //v'
  vs->_x(4*numEle+3) = _x(N-1); //v'
  
  vs->fdf();
  if (_fFlag) {
    // retrieve energy from virusShell's
    _f = vs->_f;
  }
  if (_dfFlag){
    // retrieve residue from virusShell's
    for (int i=0; i<= N-4; i++){
      _df(i) = vs->_df(i+4);
    }
    _df(N-4) = _df(N-4) - 1.5*vs->_df(4*numEle-3);
    _df(N-3) = vs->_df(4*numEle-2);
    _df(N-2) = vs->_df(4*numEle-1);
    _df(N-1) = vs->_df(4*numEle+3);
  }
}

/*
void VirusShellBC :: fdf(){
  // Prescribed BC: Cantilever ends
  int N = _x.size();
  _f = 0; // reset energy
  _df.setZero(); //reset residue to zero
  int numEle = vs->_conn.size(); 

  // //Boundary conditions applied to "VirusShell dofs"
  vs->_x(0) = 0; //u(0) = 0
  vs->_x(1) = 0; //u'(0) = 0
  vs->_x(2) = 0; //v(0) = 0
 
  vs->_x(4*numEle) = 0;   //u(L) = 0;
  vs->_x(4*numEle+1) = 0; //u'(L) = 0;
  vs->_x(4*numEle+2) = 0; //v(L) = 0
  
  for (int i=0; i<N-1; i++){
    // pass dofs to VirusShell's
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

// void VirusShell::writeSolution(string filename){
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
// void VirusShellBC::writeSolution(string filename){
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

