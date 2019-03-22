#ifndef VECTOR_H
#define VECTOR_H
#include <iostream>
#include <array>
#include <cassert>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

using namespace std;
#define RED (stuff) cout << "\033[36m " << stuff << cout << "\033[0m";
class Vector {
public:
  int _dim;
  int _idx;

public:
  gsl_vector* _gsl_vec;
  gsl_vector_view _view;

  Vector();
  Vector(int n);
  Vector(const Vector& v); //Copy constructor
  Vector(gsl_vector* v); 
  Vector(double* v, int size);
  Vector(void* v, int size);
  ~Vector();
  int size();
  void setDim(int n);
  double norm();
  double& operator()(int i);
  double operator|(Vector& v2);
  Vector view(int i, int j);
  Vector& operator+= (Vector& v2);
  Vector& operator-= (Vector& v2);
  Vector& operator*= (Vector& v2);
  Vector& operator/= (Vector& v2);
  Vector& operator*= (double alpha);
  Vector& operator/= (double alpha);  
  Vector& operator+= (double alpha);
  Vector operator+ (Vector& v2);
  Vector operator- (Vector& v2);
  Vector& operator= (Vector& v2);
  Vector& operator= (gsl_vector* v2);
  Vector& operator<< (double val);
  void setZero();
  double* data();
  void print();
  void write();
};
#endif
