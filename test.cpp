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
#include <stdlib.h>
#include <time.h>
#include "SymmReduced.h"

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


class P3: public Problem {
public:
    P3(Vector x, Vector para) : Problem(x,para) {
    }
    
    void fdf(){

        if (_fFlag)
            _f = (pow((pow(_x(0),2)+pow(_x(1),2)-1),2)+pow((pow(_x(2),2)+pow(_x(3),2)-1),2)+pow((pow(_x(4),2)+pow(_x(5),2)-1),2))/2;
        if (_dfFlag){
            _df(0) = 2*_x(0)*(pow(_x(0),2)+pow(_x(1),2)-1);
            _df(1) = 2*_x(1)*(pow(_x(0),2)+pow(_x(1),2)-1);
            _df(2) = 2*_x(2)*(pow(_x(2),2)+pow(_x(3),2)-1);
            _df(3) = 2*_x(3)*(pow(_x(2),2)+pow(_x(3),2)-1);
            _df(4) = 2*_x(4)*(pow(_x(4),2)+pow(_x(5),2)-1);
            _df(5) = 2*_x(5)*(pow(_x(4),2)+pow(_x(5),2)-1);
        }
        //cout << "in fdf" << endl;
        //cout << _df(0) <<" " << _df(1) << endl;
    }
};


int main(int argc, char** argv){

    
    //    Icosahedral im(g2_Mat,g3_Mat);
    //    Projection pm(im);
    //    cout << pm._P.rank() << endl;
    //    pm._P.range().print();
    
    
    Matrix R(6,6);
    R <<
    0     <<    0  << cos(2*M_PI/3)  << -sin(2*M_PI/3)    <<    0    <<     0 <<
    0    <<     0  <<  sin(2*M_PI/3)  << cos(2*M_PI/3)     <<    0    <<     0 <<
    0    <<     0  <<       0   <<      0  << cos(2*M_PI/3)  << -sin(2*M_PI/3)  <<
    0    <<     0   <<      0     <<    0  <<  sin(2*M_PI/3)  << cos(2*M_PI/3) <<
   cos(2*M_PI/3)  << -sin(2*M_PI/3)   <<      0    <<     0    <<     0   <<      0 <<
    sin(2*M_PI/3)  << cos(2*M_PI/3)  <<      0    <<     0     <<    0   <<      0;
    
    R.print();
    Cyclic g(3,R);
    g.listElements();
    Projection pj(g);
    pj._P.range().print();
    
    
    
    
    srand (time(NULL));


    Vector x(6);
    for (int i = 0; i < 6; i++){
        x(i) =  0* rand() % 100 ;
    }

//    x(0) = 2.3; x(1) = 5;
//    x(2) = 2.4; x(3) = 4.1;
//    x(4) = 10; x(5) = 8;
    Vector para(1);
    P3* p = new  P3(x, para);
    //p -> df();
    
    
    Matrix bas = pj._P.range();
    SymmReduced symmP(bas.T()*x,para,bas,p);
   //symmP.df();
    MultiMin M("lbfgs", &symmP);
    M._LBFGSB_Initialize();
    M.LBFGSB_Solve();
    symmP._x.print();
    
    Vector check = bas*symmP._x;
    check.print();
    cout << pow(check(0),2)+ pow(check(1),2) << endl;
    cout << pow(check(2),2)+ pow(check(3),2) << endl;
    cout << pow(check(4),2)+ pow(check(5),2) << endl;
    
    
    
//
//    MultiMin M("lbfgs", p);
//    M._LBFGSB_Initialize();
//    M.LBFGSB_Solve();
//    cout << "Starting point: " << endl;
//    x.print();
//    cout << "Solution: " << endl;
//    p -> _x.print();
//    cout << "Distance: " << endl;
//    cout << pow(p ->_x(0),2) + pow(p ->_x(1),2) << endl;
//    cout << pow(p ->_x(2),2) + pow(p ->_x(3),2) << endl;
//    cout << pow(p ->_x(4),2) + pow(p ->_x(5),2) << endl;
//
    
    
    

//
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
//
//    double ord = 5;
//    int l = 6;
//
//    double alpha = 2*M_PI/ord;
//    Q << cos(alpha) << sin(alpha) << 0 <<
//        -sin(alpha) << cos(alpha) << 0 <<
//        0 << 0 << 1;
//    Q.print();
//    cout << " " << endl;
//
//    Matrix R = compute_R(l, l, Q.T());
//
//
//
//    Cyclic c(ord,R);
//    // c.listElements();
//    Projection p(c);
//    p._P.print();


    
    //icosahedral
//
//    int l = 6;
//
//    cout << "inc = 0.01;" << endl;
//    cout << "theta1 = [0:inc/2:pi];" << endl;
//    cout << "phi1 = [0:inc:2*pi];" << endl;
//    cout << "[theta,phi]=meshgrid(theta1,phi1);" << endl;
//    cout << "l = " << l << ";" << endl;
//    cout << "leg=legendre(l,cos(theta));" << endl;
//    cout << "[whatever,row,col]=size(leg);" << endl;
//    cout << "sph_har=reshape(0*leg(1,:,:),row,col);" << endl;
//    cout << "Y = zeros(2*l+1,row,col);" << endl;
//    cout << "s = size(leg);" << endl;
//
//
//
//    Matrix g2(3,3);
//    g2 <<  -0.809016994374947 <<  0.500000000000000 << -0.309016994374947
//    << 0.500000000000000 <<  0.309016994374947 << -0.809016994374947
//    << -0.309016994374948 << -0.809016994374947 << -0.500000000000000;
//
//
//
//
//
//    Matrix g3(3,3);
//    g3 <<   -0.0000  <<  -0.0000  <<  -1.0000
//       <<    1.0000  <<   0.0000  <<  -0.0000
//       <<      0     <<  -1.0000  <<   0.0000;



//
//
//
//    Matrix r2 = compute_R(l,l,g2);
//    Matrix r3 = compute_R(l,l,g3);
//
////    g3.print();
////    Matrix temp = g3.inv();
////    temp.print();
//
////    Matrix temp = (r2*r3)^5;
////    temp.print();
//
////    cout << fmod(6*2*M_PI+M_PI/4,2*M_PI) << endl;
//
//
//    Icosahedral ico(r2,r3);
////    ico._g[0].print();
//    Projection p(ico);
////    if (p._P.rank() == 1){
//    cout << "sp_coeff =" ;
//
//    p._P.range().print();
//    cout << ";" << endl;
//
//    cout << "for m = 0: (s(1)-1)" << endl;
//    cout << "if m==0" << endl;
//    cout << "leg0=reshape(leg(1,:,:),row,col);" << endl;
//    cout << "sph_har=sph_har+sp_coeff(m+1+l)*sqrt((2*l+1)/(4*pi))*leg0;" << endl;
//    cout << "else" << endl;
//    cout << "legm=reshape(leg(m+1,:,:),row,col);" << endl;
//    cout << "sph_har=sph_har+sp_coeff(m+1+l)*(-1)^(m)*sqrt(2)*sqrt((2*l+1)/(4*pi))*sqrt(factorial(l-m)/factorial(l+m))*legm.*cos(m*phi);" << endl;
//    cout << "sph_har=sph_har+sp_coeff(-m+1+l)*(-1)^(m)*sqrt(2)*sqrt((2*l+1)/(4*pi))*sqrt(factorial(l-m)/factorial(l+m))*legm.*sin(m*phi);" << endl;
//    cout << "end" << endl;
//    cout << "end" << endl;
//    cout << "x=cos(phi).*sin(theta);\ny=sin(phi).*sin(theta);\nz=cos(theta); \nsurf(x,y,z,(sph_har)); \nshading interp \naxis equal" << endl;
//
//
    
    
    
//
//    vector<int> g3p = {8,5,6,9,7,12,2,10,11,1,4,3};
//    vector<int> g2p = {7,5,8,12,2,10,1,3,11,6,9,4};
////    cout << inv_perm(g2p,3) << " should be 8" << endl;
////    cout << inv_rep(g3p)[0] << " should be 10"<< endl;
//
//    Matrix g2_Mat = constructMat(g2p,g2);
//    Matrix g3_Mat = constructMat(g3p,g3);
//
//    Icosahedral im(g2_Mat,g3_Mat);
//    Projection pm(im);
//    cout << pm._P.rank() << endl;
//    pm._P.range().print();
//
//
//
//
//
    
    
    
    
    
    
   // cout << p._P.rank() << endl;
   // p._P.print();


    
    
    
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
    // p->_df();
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
