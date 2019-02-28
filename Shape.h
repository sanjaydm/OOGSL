#ifndef SHAPE_H
#define SHAPE_H

#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <iostream>
#include <cstdint>
#include <cstring>

class Shape {
    public:
        Shape(int dim, double point);   // This
        ~Shape();  //Decontructor

        void compute_shape_functions();
        
        Vector & get_shape_functions_N();
        Vector & get_shape_functions_DN();
        Vector & get_shape_functions_D2N();


    private:
        int  _dim;
        //Vector _point; //point to find functions on. Is a Vector because this can be 1,2,3dim
        double _point;
        Vector _N; //First set of shape functions (Vector of px1) or (Number of shape functions x 1)
        Vector _DN; //Second eval (number of shape fncts x d)
        Vector _D2N;//vector<Matrix> _D2N; // each matrix is (dxd) and there are # of functions of those.

};

#endif
