#ifndef Indented_H
#define Indented_H
#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "Problem.h"
#include "MultiMin.h"
#include "Quadrature.h"
#include "Shape.h"
#include <fstream>

class Indented: public Problem {
public:
  Indented(Vector x, Vector para, vector<double> nodes, vector<Vector> conn) : _nodes(nodes), _conn(conn), Problem(x,para) {
    _quad.setOrder(4);
    _lclResidue.setDim(4); 
    _D = para(0); //reference from para;
    _C = para(1); //reference from para;
    _nu = para(2);
    _mu = 0.0;
    _r = 1;
    _d = para(4);

    cout << "\033[32m " 
	 << "D = " << _D << "\n"
      	 << "C = " << _C << "\n" 
      	 << "nu = " << _nu << "\n"
      	 << "mu = " << _mu << "\n"
      	 << "r = " << _r << "\n"
      	 << "d = " << _d << "\n" << "\033[0m";
    
  }

  ~Indented(){
    
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
  // Indenter parameters
  double _d;
  double _r;
  double _mu;
  
  vector<double> _nodes;
  vector<Vector> _conn;
};
#endif
