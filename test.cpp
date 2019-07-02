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

//  Matrix M(3,3);
//  double alpha = 0.5;
//  M << 1 << 0 << 0 <<
//  0 << cos(alpha) << sin(alpha) <<
//  0 << -sin(alpha) << cos(alpha);


//  Vector abc = computeEulerAngles(M);
//  abc.print();
//    double alphaa = abc(0);
//    double betaa  = abc(1);
//    double gammaa = abc(2);
//
//  Wigner_d d(1, M_PI/4);
//  d._d[1].print();
   // cout << d.get_d(2,-1,0) << endl;
//
//    Matrix D = compute_D (1,M_PI/2,0.1,-M_PI/2,d);
//    D.print();
//
//
//    cout << -sin(0.1)/sqrt(2) << endl;

  
  // M1 << 
  //   1 << 0 << 0 <<
  //   0 << cos(alpha) << sin(alpha) << 
  //   0 << -sin(alpha)<< cos(alpha) ;
   

//    int ub = 2;
//    Wigner_d d(ub,betaa);
//    for (int i = 0; i <= ub; i++){
//        Matrix D = compute_D (i,alphaa,betaa,gammaa,d);
//        D.print();
//    }
//
//    Matrix R = compute_R(0,2,M);
//    R.print();

    
    

    
    
    // cyclic
    
//    Matrix Q(3,3);
//    double alpha = 2*M_PI/3;
//    Q << cos(alpha) << sin(alpha) << 0 <<
//        -sin(alpha) << cos(alpha) << 0 <<
//        0 << 0 << 1;
//    Q.print();
//    cout << " " << endl;
//
//    Matrix R = compute_R(3, 3, Q.T());
//
//
//
//    Cyclic c(3,R);
//    // c.listElements();
//    Projection p(c);
//    p._P.print();
//

    
    //icosahedral

//    Matrix I(3,3);
//    I << 1 << 0 << 0 <<
//         0 << 1 << 0 <<
//         0 << 0 << 1;
//
//    Matrix A(3,3);
//    A << 0 << (1+sqrt(5))/(2*sqrt(2)*sqrt(3+sqrt(5))) << sqrt(3+sqrt(5))/(2*sqrt(2)) <<
//
//        -(1+sqrt(5))/(2*sqrt(2)*sqrt(3+sqrt(5))) << 0 << -1/(sqrt(2)*sqrt(3+sqrt(5))) <<
//        -sqrt(3+sqrt(5))/(2*sqrt(2)) << 1/(sqrt(2)*sqrt(3+sqrt(5))) << 0;
//
//
//// R = I + A*sin(M_PI) + (A^2)*(1-cos(M_PI));
//    A.print();
//    cout << " " << endl;
//
//    Matrix temp1 = A;
//    temp1 *= sin(M_PI);
//
//    Matrix temp2 = A^2;
//    temp2 *= (1-cos(M_PI));
//    Matrix g2 = temp1 + temp2 + I ;
//
//
//    Matrix C(3,3);
//    C << 0 << -(1+sqrt(5))/(sqrt(2)*sqrt(5+sqrt(5))) << sqrt((sqrt(5)-1)/(2*sqrt(5))) <<
//        (1+sqrt(5))/(sqrt(2)*sqrt(5+sqrt(5))) << 0 << 0 <<
//        -sqrt((sqrt(5)-1)/(2*sqrt(5))) << 0 << 0;
//
//
//// S = I + sin(2*pi/3)*C + (1-cos(2*pi/3))*C^2;
//    temp1 = C;
//    Matrix temp3 = C^2;
//    temp1 *= sin(2*M_PI/3);
//    temp3 *= (1-cos(2*M_PI/3));
//    Matrix g5 = I + temp1 + temp3;
//


    
// T = R * inv(S) * R * S * R * inv(S);
//    Matrix temp4 = S.inv();
//    Matrix T = g2 * temp4 * g2 * g5 * g2 * temp4;
//
//    Matrix S2 = (g5^2);
//    Matrix g3 = (g5^2) * g2 * (g5^3) * T; //g3
//    Matrix g2 = R;
    
    
    Matrix g2(3,3);
    g2 <<  -0.809016994374947 <<  0.500000000000000 << -0.309016994374947
    << 0.500000000000000  <<  0.309016994374947 << -0.809016994374947
    << -0.309016994374948 <<  -0.809016994374947 << -0.500000000000000;
    
    Matrix g3(3,3);
    g3 <<   -0.499999999999997 <<  0.736685209782634 << -0.455296498655014
   << -0.736685209782634 << -0.085410196624967  << 0.670820393249936
    << 0.455296498655014 <<  0.670820393249936  << 0.585410196624968;
    
    

    
    
    Matrix r2 = compute_R(3,3,g2);
    Matrix r3 = compute_R(3,3,g3);
    
    Matrix temp = (r2^2);
    temp.print();
    
    cout << fmod(6*2*M_PI+M_PI/4,2*M_PI) << endl;
//    Icosahedral ico(r2,r3);
//    Projection p(ico);
//    p._P.print();
    

// T = g2d
// S = g5
// R = g2
    
    
    
    
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
