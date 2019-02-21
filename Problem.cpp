#include "Problem.h"
#include "gsl/gsl_vector.h"
#include "Vector.h"
Problem::Problem(Vector x, Vector para) : _x(x), _para(para) {
  _f = 0.0;
  _fFlag = true;
  _dfFlag = true;
  _Ndof = _x.size();
  _df.setDim(_Ndof);
  _d2f.setDim(_Ndof,_Ndof);
};

double Problem::f(){
  _fFlag = true;
  _dfFlag = false;
  fdf();
  _dfFlag = true;
  return _f;
}
void Problem::df(){
  _fFlag = false;
  _dfFlag = true;
  fdf();
  _fFlag = true;
}
// virtual void Problem::fdf(){
//   if (_fFlag)
//     _f = (_x(0)-3.0)*(_x(0)-3.0)+ (_x(1)-1.0)*(_x(1)-1.0);
//   if (_dfFlag){
//     _df(0) = 2*(_x(0)-3.0);
//     _df(1) = 2*(_x(1)-1.0);
//   }
// //cout << "in fdf" << endl;
// //cout << _df(0) <<" " << _df(1) << endl;
// }
void Problem::d2f(){
  Vector df0 = _df;
  double eps = 1e-4;
  for (int i=0; i<_Ndof; i++){
    _x(i) += eps;
    fdf();
    Vector deltaf = (_df - df0);
    deltaf *= 1/eps;
    _d2f(i) = deltaf;
    _x(i) -= eps;
  }
  _df = df0;
}
