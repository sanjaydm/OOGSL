#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "Problem.h"
#include "Group.h"
#include "MultiMin.h"
#include "MultiRoot.h"
#include "wigner.cpp"
#include "compute_D_matrix.cpp"
#include <cmath>

class P2: public Problem {
public:
  P2(Vector x, Vector para) : Problem(x,para) {
  }

  void fdf(){
  double a = 11.2312; double b = 12.123132;
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

  Matrix M(3,3);
  double alpha = 0.5;
  M << 1 << 0 << 0 << 
  0 << cos(alpha) << sin(alpha) <<
  0 << -sin(alpha) << cos(alpha);


  Vector abc = computeEulerAngles(M);
  abc.print();
<<<<<<< HEAD

  double beta = M_PI/6;
  Wigner_d d(1, beta);
  //cout << d.get_d(1, 1, 0) << " " << -sin(beta)/sqrt(2) << endl;
  //cout << d.get_d(2,-1,0) << endl;

  //Matrix D = compute_D (1,M_PI/2,0.1,-M_PI/2);
  // D.print();
=======
    double alphaa = abc(0);
    double betaa  = abc(1);
    double gammaa = abc(2);
  
//  Wigner_d d(1, M_PI/4);
//  d._d[1].print();
   // cout << d.get_d(2,-1,0) << endl;
//
//    Matrix D = compute_D (1,M_PI/2,0.1,-M_PI/2);
//    D.print();
//
//
//    cout << -sin(0.1)/sqrt(2) << endl;

>>>>>>> 1a050af6297e23bf0eeb18015191ba6a0ad34103
    
    int ub = 2;
    Wigner_d d(ub,betaa);
    for (int i = 0; i <= ub; i++){
        Matrix D = compute_D (i,alphaa,betaa,gammaa,d);
        D.print();
    }
    
<<<<<<< HEAD
  //cout << -sin(0.1)/sqrt(2) << endl;

  Matrix N(3);
  N << 1 << 1 << 1
  << 2 << 2 << 2
  << 4 << 5 << 6;
  Matrix rg = N.range();
  rg.print();
    
=======
    Matrix R = compute_R(0,2,M);
    R.print();
>>>>>>> 1a050af6297e23bf0eeb18015191ba6a0ad34103


    
    
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

  // Vector x0(2);
  // x0(0) = .3223; x0(1) = .1;
  // Vector para(1);
  // P2* p = new  P2(x0, para);
  
  // MultiMin M("lbfgs", p);
  // M._LBFGSB_Initialize();
  // M.LBFGSB_Solve();

  // MultiRoot rt("generic", p);
  // rt._GSLRoot_Initialize();
  // rt.GSLRoot_Solve();

  // p->_x.print();

  //p.d2f();
  //p._d2f.print();
  // Group I;
  // cout << "Size of I = " << I._g.size() << endl;

  // Constructing projection operator
  
  // Vector v(5);
  // v << 1 << 2 << 3 << 4 << 5;
  // v.print();
  // Vector w = v.view(0,3);
  // cout <<"asfasdf\n";
  // w.print();

  // Matrix m(5,5);
  // m << 1 << 1 << 1 << 0 << 0
  //   << 0 << 1 << 1 << 0 << 0
  //   << 0 << 0 << 1 << 1 << 0
  //   << 0 << 0 << 0 << 0 << 1
  //   << 1 << 1 << 1 << 1 <<1;
  // m.print();
  // Vector x(5);
  
  // (m*v).print();

  // m.col(0).print();
  return 0;
}
