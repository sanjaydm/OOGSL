#include "VirusShellBC.h"

#include "VirusShell.h"

void VirusShellBC :: fdf(){
  int N = _x.size();
  for (int i=0; i<N-1; i++){
    vs->_x(i+3) = _x(i);
  }
  vs->_x(N+5) = _x(N-1);
  
  vs->fdf();
  if (_fFlag) {
    _f = vs->_f;
  }
  if (_dfFlag){
    for (int i=0; i<N-1; i++){
      _df(i) = vs->_df(i+3);
    }
    _df(N-1) = vs->_df(N+5);
  }
}
