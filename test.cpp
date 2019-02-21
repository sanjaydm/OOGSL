#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "Problem.h"
#include "Group.h"
#include "MultiMin.h"

class P2: public Problem {
public:
  P2(Vector x, Vector para) : Problem(x,para) {
  }

  void fdf(){
  double a = 1.12312; double b = 12.123132;
  if (_fFlag)
    _f = (_x(0)-a)*(_x(0)-a)+ (_x(1)-b)*(_x(1)-b);
  if (_dfFlag){
    _df(0) = 2*(_x(0)-a);
    _df(1) = 2*(_x(1)-b);
  }
//cout << "in fdf" << endl;
//cout << _df(0) <<" " << _df(1) << endl;
  }
};
int main(int argc, char** argv){

  // v(0) = 1.1231;
  // v(1) = 2;
  // Vector v2(3);
  // v *= 2.5;
  // v = v+v2;
  // v.print();
  // v << 1 << 2 << 3;
  // v.print();
  // cout << "-------\n";
  // cout << v.norm() << endl;
  // cout << "-------\n";
  // double s = (v|v);
  // cout << s << endl;

  // cout << "----- Matrix -----\n";
  // Matrix m(3);
  // m << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9 << 10 << 11 <<12 << 13;
  // cout << "-------" << endl;
  // m.print();
  // cout << "-------\n";
  // m(0).print();
  // m(1).print();
  // cout << "-------\n";
  // m(0)(0) = 5.5;
  // m.print();

  // gsl_vector_view v = gsl_matrix_column(m._gsl_mat,1);
  // cout << v(0) << endl;
  // cout << "owner of m(0) " << m(0)._gsl_vec->owner << endl;
  // Vector v(3);
  // cout << "owner of v " << v._gsl_vec->owner << endl;
  // gsl_vector* v = gsl_vector_calloc(5);
  // Vector a(v);
  // a.print();

  // gsl_matrix* m = gsl_matrix_alloc(3,3);
  // Matrix M(m);
  // M.print();

  Vector x0(2);
  x0(0) = .3223; x0(1) = .1;
  Vector para(1);
  P2* p = new  P2(x0, para);
  // x0.print();
  // (p._x).print();
  // p.f();
  // p->d2f();
  // p->_d2f.print();

  MultiMin M("lbfgs", p);
  M._LBFGSB_Initialize();
  M.LBFGSB_Solve();
  // M._GSLMin_Initialize();
  // M.GSLMin_Solve();
  p->_x.print();

  //p.d2f();
  //p._d2f.print();
  // Group I;
  // cout << "Size of I = " << I._g.size() << endl;

  // Constructing projection operator

  return 0;
}
