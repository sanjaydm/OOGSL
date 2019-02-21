#ifndef PROBLEM_H
#define PROBLEM_H
#include "Vector.h"
#include "Matrix.h"
#include <gsl/gsl_vector.h>
class Problem;
struct CCallbackHolder{
  Problem* cls;
  void* params;
};

class Problem {
public:
  int _Ndof;
  Vector _x;
  Vector _para;
  double _f;
  Vector _df;
  Matrix _d2f;
  bool _fFlag;
  bool _dfFlag;

public:  
  // double _gsl_f(const gsl_vector * x, void * params);
  // void _gsl_df(const gsl_vector * x, void * params, gsl_vector * g);
  // void _gsl_fdf(const gsl_vector * x, void * params, double * f, gsl_vector * g);

  Problem (){;}
  Problem (Vector x, Vector para);
  ~Problem(){
  };
  double f();
  void df();
  virtual void fdf() {;};
  void d2f();
  static double _C_gsl_f(const gsl_vector* x,void* classNParams){
    CCallbackHolder* temp = static_cast<CCallbackHolder*>(classNParams);
    return temp->cls->_gsl_f(x, temp->params);
  }
  static void _C_gsl_df(const gsl_vector* x,void* classNParams, gsl_vector* g){
    CCallbackHolder* temp = static_cast<CCallbackHolder*>(classNParams);
    temp->cls->_gsl_df(x, temp->params, g);
  }
  static void _C_gsl_fdf(const gsl_vector* x,void* classNParams, double* f, gsl_vector* g){
    CCallbackHolder* temp = static_cast<CCallbackHolder*>(classNParams);
    temp->cls->_gsl_fdf(x, temp->params, f, g);
  }

  double _gsl_f(const gsl_vector* x, void* params){
    _x._gsl_vec = const_cast<gsl_vector*>(x);
    (_para._gsl_vec)->data = (double*) params; 
    _f = f();
    return _f;  
  }

  void _gsl_df(const gsl_vector* x, void* params, gsl_vector* g){
    _x._gsl_vec = const_cast<gsl_vector*>(x);
    (_para._gsl_vec)->data = (double*) params;
    _df._gsl_vec = g;
    df();
    //gsl_vector_memcpy(g, _df._gsl_vec);
  }

  void _gsl_fdf(const gsl_vector * x, void * params, double * f, gsl_vector * g){

    _x._gsl_vec = const_cast<gsl_vector*>(x);
    (_para._gsl_vec)->data = (double*) params;
    _df._gsl_vec = g;
    fdf();
    //gsl_vector_memcpy(g, _df._gsl_vec);
    *f = _f;
  }
};

#endif
