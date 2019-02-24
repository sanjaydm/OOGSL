#ifndef QUADRATURE_H
#define QUADRATURE_H
#include <vector>
#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include "gausslegendre.h"

class Quadrature {
public:
  Quadrature();
  Quadrature(int order, int dim);   
  ~Quadrature();  

  void createNodesWeights(); 
  
  vector<Vector>& nodes(); 
  Vector & weights(); 

  void setOrder(int order);
  
private:

  vector<Vector> _nodes; 
  Vector       _weights; 
  int          _order;
  int          _dim;

};

#endif
