#ifndef SymmReduced_H
#define SymmReduced_H
#include <iostream>
#include <string.h>
#include "Vector.h"
#include "Matrix.h"
#include "Problem.h"
#include "MultiMin.h"
#define PI (3.141592653589793)

class SymmReduced: public Problem {
public:
  SymmReduced(){;}
  SymmReduced(Vector x0, Vector para, Matrix& bas, Problem* prob) {
    _bas = bas;
    _basT = _bas.T();
    array<int,2> rowCol = _bas.size();
    _fullX.setDim(rowCol[0]); //Initialize _fullX to the original dimension of the problem
    _x.setDim(rowCol[1]);
    _x = x0;
    _para = para;
    _prob = prob;
    _fFlag = true;
    _dfFlag = true;
    _f = 0.0;
    _Ndof = _x.size();
    _df.setDim(_Ndof);
    _d2f.setDim(_Ndof, _Ndof);
  }
  //~SymmReduced();
  void convertBasis2Full();
  void fdf();
  
  // Data 
  Vector _fullX; //"Full" dofs of the originial (non symmetry reduced problem)
  Matrix _bas; //Symmetry adapted basis
  //Matrix _P; // Projection matrix
  Matrix _basT;
  Problem* _prob;

};
#endif
