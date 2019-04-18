#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <array>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "Vector.h"
using namespace std;
class Matrix {
private:
  int _dim1, _dim2;
  unsigned int _idx1, _idx2;
public:
  gsl_matrix* _gsl_mat;
  Vector _vec;
  gsl_vector_view _vview;

  Matrix() ;
  Matrix(int n);
  Matrix(int n1, int n2);
  Matrix(const Matrix& v); //Copy constructor
  ~Matrix();
  array<int,2> size();
  void setDim(int n1, int n2);
  void setDim(int n);
  double norm();
  // ------------------------------------------------------------
  double& operator()(int i, int j);
  Vector& operator()(int i);
  Matrix& operator+= (Matrix& v2);
  Matrix& operator-= (Matrix& v2);
  Matrix& elementMult (Matrix& v2);
  Matrix& elementDiv (Matrix& v2);
  Matrix& operator*= (double alpha);
  Matrix& operator^= (int n);
  Matrix operator^ (int n);
  Matrix& operator*= (Matrix& m);
  Matrix& operator+= (double alpha);
  Matrix operator+ (Matrix& v2);
  Matrix operator* (Matrix& v2);
  Matrix inv();
  Matrix operator- (Matrix& v2);    
  Matrix& operator= (Matrix& v2);
  Matrix& operator<< (double val);
  Vector row(int i);
  Vector col(int i);
  Vector operator* (Vector& v);
  void print();
};
#endif
