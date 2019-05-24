#include "Group.h"
#include <cmath>
#include <gsl/gsl_math.h>
Group::Group(){
  constructGroup();
}
void Group::listElements(){
  for (int i=0; i < _g.size(); i++) {
    _g[i].print();
    cout << "------\n";
  }
}
/* ---------------- Cyclic -----------*/
Cyclic::Cyclic(int n): _n(n){
  _gen.push_back(Matrix(1,'i'));
  constructGroup();
}
Cyclic::Cyclic(int n, Matrix gen) : _n(n){
  _gen.push_back(gen);
  constructGroup();
}
void Cyclic::constructGroup(){
  _g.push_back(_gen[0]);
  for (int m=1; m<_n; m++){
    _g.push_back(_gen[0]*_g[m-1]);
  }
}



/*
void Group::constructGroup(){
  // double tau = sqrt(5);
  // Matrix g2(3);
  // g2 << -0.75-1/(4*tau) << -sqrt(2/(5+tau))/(tau-1) << 0.1*(-5+tau)
  //    << (5*sqrt(10-2*tau)-2*sqrt(10*(5+tau)) )/(10*(-1+tau)) << 0.5*sqrt(0.5*(3-tau)) << sqrt (0.1*(5+tau))
  //    << .1*(-5+tau) << 2*sqrt(2/(5+tau))/(-1+tau) << -1/tau;
  // //_g.push_back(g2);

  // Matrix g3(3);
  // g3 << -0.5+1/tau-0.5*sqrt(0.5*(3-tau))
  //    << (sqrt(10)-sqrt(5*(3-tau)))/(2*sqrt(5-tau)-2*sqrt(5+tau))
  //    << 0.1*(5+tau)
  //    << -1/sqrt(2*(5+tau))
  //    <<0.25*(2+sqrt(6-2*tau))
  //    <<sqrt(2/(5+tau))
  //    <<-2/tau << 0 << - 1/tau;


  // // g5 and g2d
  // Matrix g5 = (g3*g2*g3);
  // Matrix g5inv = g5.inv();
  // Matrix g2_g5inv = g2*g5inv;
  // Matrix g2d = g2_g5inv * g2*g5* g2_g5inv;

  // // Create Generator
  // _gen.push_back(g2);
  // _gen.push_back(g3);
  // _gen.push_back(g5);

  // cout << "\033[34m Generating group elements ....";
  // for (int mu=0; mu<=4; mu++ ){
  //   for (int sigma=0; sigma<=1; sigma++){
  //     Matrix g2d_sigma = (g2d^sigma);
  //     _g.push_back( (g5^mu)* g2d_sigma);
  //     for (int nu=0; nu<=4; nu++){
  //       Matrix g5_nu = g5^nu;
  //       _g.push_back(g2* g5_nu * g2d_sigma);
	
  //     }
  //   }
  // }
  // cout << "done\033[0m\n";  
    
  }*/
