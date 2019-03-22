#ifndef VIRUSSHELLBC_H
#define VIRUSSHELLBC_H
#include <iostream>
#include "VirusShell.h"


class VirusShellBC: public Problem {
public:
 VirusShellBC(Vector x, Vector para, vector<double> nodes, vector<Vector> conn): Problem(x,para) {
    Vector xvs(x.size()+6);
    vs = new VirusShell(xvs, para, nodes, conn);
    
  }

  ~VirusShellBC(){
    delete vs;
  }
  void fdf();
  VirusShell* vs;
};
#endif
