#include "VirusShellBC.h"

#include "VirusShell.h"
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
