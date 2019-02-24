#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "VirusShell.h"
#include "MultiMin.h"


int main(int argc, char** argv){

  Vector x0(2);
  x0(0) = .3223; x0(1) = .1;
  Vector para(1);
  VirusShell* p = new  VirusShell(x0, para);

  // MultiMin M("lbfgs", p);
  // M._LBFGSB_Initialize();
  // M.LBFGSB_Solve();
  // // M._GSLMin_Initialize();
  // // M.GSLMin_Solve();
  // p->_x.print();

  return 0;
}
