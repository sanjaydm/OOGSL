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
 SymmReduced(Vector x0, Vector para, Matrix& bas, Matrix& Proj, Problem* prob) : Problem(x0,para) {
    _bas = bas;
    _P = Proj;
    array<int,2> rowCol = _bas.size();
    Vector _fullX(rowCol[0]); //Initialize _fullX to the original dimension of the problem
  }
  ~SymmReduced();
  void convertBasis2Full();
  void fdf();
  
  // Data 
  Vector _fullX; //"Full" dofs of the originial (non symmetry reduced problem)
  Matrix _bas; //Symmetry adapted basis
  Matrix _P; // Projection matrix
  Problem* _prob;

};
#endif
