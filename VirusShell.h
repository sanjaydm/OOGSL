#ifndef VIRUSSHELL_H
#define VIRUSSHELL_H
#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "Problem.h"
#include "MultiMin.h"
#include "Quadrature.h"
#include "Shape.h"
#include <fstream>

class VirusShell: public Problem {
public:
  VirusShell(Vector x, Vector para, vector<double> nodes, vector<Vector> conn) : _nodes(nodes), _conn(conn), Problem(x,para) {
    _quad.setOrder(4);
    _lclResidue.setDim(4); 
    _D = 1.0; //reference from para;
    _C = 1.0; //reference from para;
    _nu = 0.3;
    
  }

  ~VirusShell(){
    
  }
  void NeoHookean(Vector& est, Vector& kst);
  void printU();
  void writeMesh(string filename);
  void writeSolution(string filename);
  
  void fdf();

  // Data specific to this problem
  Quadrature _quad;
  double _lclEnergyDensity;
  Vector _lclResidue;
  double _D;
  double _C;
  double _nu;
  vector<double> _nodes;
  vector<Vector> _conn;
};
#endif
