#ifndef MultiRoot_H
#define MultiRoot_H
#include <gsl/gsl_multiroots.h>
#include <string>
#include "Problem.h"
#include <vector>
#include "Vector.h"
#include <cstdlib>
#include <string.h>

class MultiRoot{
public:
  gsl_multiroot_fsolver* _gsl_fsolver;
  gsl_multiroot_function _f;
  const gsl_multiroot_fsolver_type* _T;

  int _dim; //dimension of the system
  string _type;
  double _tol = 1e-6;
  void _setMinimizerType(string type);
  Problem* _problem;

public:
  MultiRoot();
  MultiRoot(string type, Problem* prob);

  void _GSLRoot_Initialize();
  int GSLRoot_Iterate();
  void GSLRoot_Solve();
};
#endif
