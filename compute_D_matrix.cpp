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



double compute(int l, int m, int mp, double alpha, double beta, double gamma){
    Wigner_d d(l,beta);
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
            cout << d1 << endl;
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





Matrix compute_D (int l, double alpha, double beta, double gamma){
    Matrix D(2*l+1);
    for (int m = -l; m <= l; m++ ){
        for (int mp = -l; mp <= l; mp++){
            D(i(m,l),i(mp,l)) = compute(l,m,mp,alpha,beta,gamma);
        }
    }
    
    return D;
}
