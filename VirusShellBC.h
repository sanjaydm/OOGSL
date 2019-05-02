#ifndef VIRUSSHELLBC_H
#define VIRUSSHELLBC_H
#include <iostream>
#include "VirusShell.h"


class VirusShellBC: public Problem {
public:
  double _alpha; //angle 
  double _u0; // rad displacement at s=0
  double _v0; // vertical displacement at s=0
  double _d; // offset distance: center of indentor and center of cylinder
  double _rho; // radius of the indenter

  
 VirusShellBC(Vector x, Vector para, vector<double> nodes, vector<Vector> conn): Problem(x,para) {
    Vector xvs(x.size()+8);
    //Vector xvs(x.size()+6); For cantilever
    _rho = 1;
    _d = 2-0.1;
    _alpha = 0.52;

    _u0 = _d - _rho*cos(_alpha)-1; //1, here represents the radius of the cylinder
    _v0 = _rho*sin(_alpha);
    
    vs = new VirusShell(xvs, para, nodes, conn);
    
  }

  ~VirusShellBC(){
    delete vs;
  }
  void fdf();
  void writeMesh(string filename) {
    vs->writeMesh(filename);
  }
  void writeSolution(string filename) {
    // int N = _x.size();
    // int numEle = vs->_conn.size();
    // double slope = (vs->_nodes[1]-vs->_nodes[0])/2;
    // // //Boundary conditions applied to "VirusShell dofs"
    // vs->_x(0) = _u0; //u(0) = 0
    // vs->_x(1) = cos(_alpha)*slope; //u'(0) = 0
    // vs->_x(2) = _v0; //v(0) = 0
    // vs->_x(3) = (1-sin(_alpha))*slope; //v(0) = 0

    // vs->_x(4*numEle-3) = -1.5*_x(N-4);
    // vs->_x(4*numEle) = 0;   //u(L) = 0;
    // vs->_x(4*numEle+1) = 0; //u'(L) = 0;
    // vs->_x(4*numEle+2) = 0; //v(L) = 0
    vs->writeSolution(filename);
  }
  VirusShell* vs;
};
#endif
