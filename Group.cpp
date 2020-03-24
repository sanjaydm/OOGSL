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


/* ----------- Icosahedral ------------- */

Icosahedral::Icosahedral(Matrix g2, Matrix g3){
    _g2 = g2;
    _g3 = g3;
    _gen.push_back(_g2);
    _gen.push_back(_g3);
    constructGroup();
}

void Icosahedral::constructGroup(){
    Matrix g5 = (_g3*_g2*_g3);
    Matrix g5inv = g5.inv();
    Matrix g2d = _g2 * g5inv * _g2 * g5 * _g2 * g5inv;
    for (int mu = 0; mu <= 4; mu++){
        for (int sigma = 0; sigma <= 1; sigma ++){
            Matrix g5_mu = g5^mu;
            Matrix g2d_sigma =g2d^sigma;
            _g.push_back( g5_mu * g2d_sigma );
            for (int nu = 0; nu <= 4; nu++){
                Matrix g5_nu = g5^nu;
                Matrix g2d_sigma = g2d^sigma;
                _g.push_back( g5_mu * _g2 * g5_nu * g2d_sigma);
            }
        }
    }
}
/* --------------- D2 -------------- */

D2::D2(Matrix g2, Matrix g22){
    _g2 = g2;
    _g22 = g22;
    _gen.push_back(_g2);
    _gen.push_back(_g22);
    constructGroup();
}

void D2::constructGroup(){
    for (int nu = 0; nu <= 1; nu ++){
        for (int sigma = 0; sigma <= 1; sigma ++){
            Matrix g22_nu = _g22^nu;
            Matrix g2_sigma = _g2^sigma;
            _g.push_back( g22_nu * g2_sigma );
        }
    }
}

/* --------------- D3 -------------- */

D3::D3(Matrix g2, Matrix g3){
    _g2 = g2;
    _g3 = g3;
    _gen.push_back(_g2);
    _gen.push_back(_g3);
    constructGroup();
}

void D3::constructGroup(){
    for (int nu = 0; nu <= 2; nu ++){
        for (int sigma = 0; sigma <= 1; sigma ++){
            Matrix g3_nu = _g3^nu;
            Matrix g2_sigma = _g2^sigma;
            _g.push_back( g3_nu * g2_sigma );
        }
    }
}
/* --------------- D4 -------------- */

D4::D4(Matrix g2, Matrix g4){
    _g2 = g2;
    _g4 = g4;
    _gen.push_back(_g2);
    _gen.push_back(_g4);
    constructGroup();
}

void D4::constructGroup(){
    for (int nu = 0; nu <= 3; nu ++){
        for (int sigma = 0; sigma <= 1; sigma ++){
            Matrix g4_nu = _g4^nu;
            Matrix g2_sigma = _g2^sigma;
            _g.push_back( g4_nu * g2_sigma );
        }
    }
}


/* --------------- D5 -------------- */

D5::D5(Matrix g2d, Matrix g5){
    _g2d = g2d;
    _g5 = g5;
    _gen.push_back(_g2d);
    _gen.push_back(_g5);
    constructGroup();
}

void D5::constructGroup(){
    for (int nu = 0; nu <= 4; nu ++){
        for (int sigma = 0; sigma <= 1; sigma ++){
            Matrix g5_nu = _g5^nu;
            Matrix g2d_sigma = _g2d^sigma;
            _g.push_back( g5_nu * g2d_sigma );
        }
    }
}
/* --------------- D6 -------------- */

D6::D6(Matrix g2, Matrix g6){
    _g2 = g2;
    _g6 = g6;
    _gen.push_back(_g2);
    _gen.push_back(_g6);
    constructGroup();
}

void D6::constructGroup(){
    for (int nu = 0; nu <= 5; nu ++){
        for (int sigma = 0; sigma <= 1; sigma ++){
            Matrix g6_nu = _g6^nu;
            Matrix g2_sigma = _g2^sigma;
            _g.push_back( g6_nu * g2_sigma );
        }
    }
}
/* --------------- D9 -------------- */

D9::D9(Matrix g2, Matrix g9){
    _g2 = g2;
    _g9 = g9;
    _gen.push_back(_g2);
    _gen.push_back(_g9);
    constructGroup();
}

void D9::constructGroup(){
    for (int nu = 0; nu <= 8; nu ++){
        for (int sigma = 0; sigma <= 1; sigma ++){
            Matrix g9_nu = _g9^nu;
            Matrix g2_sigma = _g2^sigma;
            _g.push_back( g9_nu * g2_sigma );
        }
    }
}

/* --------------- Tetrahedral ------------- */

Tetrahedral:: Tetrahedral(Matrix g2, Matrix g3d){
    _g2 = g2;
    _g3d = g3d;
    _gen.push_back(_g2);
    _gen.push_back(_g3d);
    constructGroup();
}

void Tetrahedral::constructGroup(){
    for (int tau = 0; tau <= 2; tau ++){
        Matrix g3d_tau = _g3d^tau;
        _g.push_back(g3d_tau);
        for (int mu = 0; mu <= 2; mu ++){
            Matrix g3d_mu = _g3d^mu;
            _g.push_back( g3d_tau * _g2 * g3d_mu);
        }
        
    }
}




Projection:: Projection(Group group):_group(group){
    computeProjection();
}

void Projection::computeProjection(){
    float n = _group._g.size();
    array<int,2> rc = _group._g[0].size();
    _P.setDim(rc[0],rc[1]);
    _P = _group._g[0];
    for (int i = 1; i < n; i++){
        _P += (_group._g[i]);
    }
    _P *= (1/n);
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
