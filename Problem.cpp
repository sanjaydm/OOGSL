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

int Problem::checkConsistency(double eps, double tol){

  Vector dfappx(_Ndof);
  dfappx.setZero();
  double f0 = f();
  //cout << "f0 = " << f0 << endl;
  for (int i=0; i<_Ndof; i++){
    _x(i) += eps;
    dfappx(i) = (f()-f0)/eps;
    //cout << "\033[32m" << f() << " " << f0 << " " << (f()-f0) << "\033[0m\n";
    _x(i) -= eps;
  }
  df();
  Vector absError = (dfappx- _df);
  double absErrorNrm = absError.norm();
  dfappx.print();
  _df.print();
  absError.print();
  cout << "energy = " << f0 << endl;
  //_df.print();
  //absError.print();
  cout << "\033[7;36m\n";
  cout << "absolute error = " << absErrorNrm << endl;
  cout << "relative error = " << absErrorNrm/_df.norm() << endl;
  cout << "\033[0m\n";
  if (absErrorNrm < tol ){
    cout << "\033[92m Consistency check passed \033[0m\n";
    return 1;
  }
  else {
    cout << "\033[31m Consistency check failed \033[0m\n";
    return -1;
  }
    
}
