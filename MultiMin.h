#ifndef MULTIMIN_H
#define MULTIMIN_H
#include <gsl/gsl_multimin.h>
#include <string>
#include "Problem.h"
#include <vector>
#include "Vector.h"
#include <cstdlib>
#include <string.h>

extern "C" void setulb_(int * n, int *m, 
                        double * x, double * l, double * u, 
                        int * nbd, 
                        double * f, double * g, 
                        double * factr, double * pgtol,
                        double * wa, 
                        int * iwa,
                        char * task, 
                        int * iprint,  
                        char * csave, 
                        int * lsave, 
                        int * isave, 
                        double * dsave );
class MultiMin{
public:
  gsl_multimin_fdfminimizer* _gsl_fdfminimizer;
  gsl_multimin_function_fdf _fdf;
  gsl_multimin_function* _f;
  const gsl_multimin_fdfminimizer_type* _T;

  int _dim; //dimension of the system
  string _type;
  double _stepSize = 0.1;
  double _tol = 1e-8;
  double _lineSearchTol = 0.1;
  void _setMinimizerType(string type);
  Problem* _problem;

  //LBFGS Related variables
  Vector _l;
  Vector _u;
  vector<int> _nbd;
  vector<double> _wa;
  vector<int> _iwa;
  int _m = 5;
  double _factr = 0.1;
  char _csave[60];
  char _task[60];
  int _lsave[4];
  double _dsave[29];
  int _isave[44];
  int _iprint=1;
    
public:
  MultiMin();
  MultiMin(string type, Problem* prob);
  void _LBFGSB_Initialize();
  void LBFGSB_Iterate();
  void LBFGSB_Solve();

  void _GSLMin_Initialize();
  int GSLMin_Iterate();
  void GSLMin_Solve();
};
#endif
