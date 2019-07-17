//
//  compute_D_matrix.cpp
//  
//
//  Created by Clare Cheng on 6/26/19.
//


#include <iostream>
#include <cmath>
#include "Matrix.h"
#include "Vector.h"
#define i(m, L) m+L

using namespace std;



double compute(int l, int m, int mp, double alpha, double beta, double gamma, Wigner_d d){
    double temp=0;
    if (m < 0){
        if (mp < 0){
            double d1 = d.get_d(l,m,-mp);
            double d2 = d.get_d(l,m,mp);
            temp = -pow(-1,mp)*d1*cos(-mp*gamma+m*alpha) + d2*cos(mp*gamma+m*alpha);
        }
        else if (mp == 0){
            double d1 = d.get_d(l,m,mp);
            temp = -sqrt(2)*d1*sin(m*alpha);
        }
        else if (mp > 0){
            double d1 = d.get_d(l,m,mp);
            double d2 = d.get_d(l,m,-mp);
            temp = -pow(-1,mp)*d1*sin(mp*gamma+m*alpha) - d2*sin(-mp*gamma+m*alpha);
        }
    }
    else if (m == 0){
        if (mp < 0){
            double d1 = d.get_d(l,m,mp);
            temp = sqrt(2)*d1*sin(mp*gamma);
        }
        else if (mp == 0){
            double d1 = d.get_d(l,m,mp);
            temp = d1;
        }
        else if (mp > 0){
            double d1 = d.get_d(l,m,mp);
            temp = pow(-1,mp)*sqrt(2)*d1*cos(mp*gamma);
        }
    }
    else if (m > 0){
        if (mp < 0){
            double d1 = d.get_d(l,m,mp);
            double d2 = d.get_d(l,m,-mp);
            temp = pow(-1,m)*(d1*sin(mp*gamma+m*alpha) - pow(-1,mp)*d2*sin(-mp*gamma+m*alpha));
        }
        else if (mp == 0){
            double d1 = d.get_d(l,m,mp);
            temp = pow(-1,m)*sqrt(2)*d1*cos(m*alpha);
        }
        else if (mp > 0){
            double d1 = d.get_d(l,m,mp);
            double d2 = d.get_d(l,m,-mp);
            temp = pow(-1,m)*(pow(-1,mp)*d1*cos(mp*gamma+m*alpha) + d2*cos(-mp*gamma+m*alpha));
        }
    }
    return temp;
}





Matrix compute_D (int l, double alpha, double beta, double gamma, Wigner_d& d){
    Matrix D(2*l+1);
    for (int m = -l; m <= l; m++ ){
        for (int mp = -l; mp <= l; mp++){
            D(i(m,l),i(mp,l)) = compute(l,m,mp,alpha,beta,gamma,d);
        }
    }
    
    return D;
}

Matrix compute_R(int lmin, int lmax, Matrix Q){
    Vector v = computeEulerAngles(Q);
    double alpha = -v(0);
    double beta = fmod(v(1), 2*M_PI);
    double gamma = -v(2);
//    cout << "alpha " << alpha << endl;
//    cout << "beta " << beta << endl;
//    cout << "gamma " << gamma << endl;
    Wigner_d d(lmax,beta);
    int n = (lmin+lmax+1)*(lmax-lmin+1);
    Matrix R(n);
    int begin = 0;
    for (int l = lmin; l <= lmax; l++){
        Matrix D = compute_D(l,alpha,beta,gamma,d);
        for (int i = begin; i <= begin+2*l; i++){
            for (int j = begin; j <= begin+2*l; j++){
                R(i,j) = D(i-begin,j-begin);
            }
        }
        begin = begin + 2*l + 1;
    }
    
    return R;
}


int inv_perm(vector<int> v, int pos){
    int len = v.size();
    for (int i = 0; i < len; i++){
        if (v[i] == pos){
            return i+1;
        }
    }
}

vector<int> inv_rep(vector<int> v){
    int len = v.size();
    vector<int> ret_v;
    for (int i = 1; i <= len; i++){
        ret_v.push_back(inv_perm(v,i));
    }
    return ret_v;
}


Matrix constructMat(vector<int> v, Matrix Q){
    int size = 3*v.size();
    Matrix R(size);
    vector<int> pos = inv_rep(v);
    int begin = 0;
    for (int i = 0; i < v.size(); i++){
        begin = (pos[i] - 1) *3;
        int j = i*3;
        R(j,begin) = Q(0,0);
        R(j,begin+1) = Q(0,1);
        R(j,begin+2) = Q(0,2);
        R(j+1,begin) = Q(1,0);
        R(j+1,begin+1) = Q(1,1);
        R(j+1,begin+2) = Q(1,2);
        R(j+2,begin) = Q(2,0);
        R(j+2,begin+1) = Q(2,1);
        R(j+2,begin+2) = Q(2,2);
    }
    return R;
}

