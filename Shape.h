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
        Shape(int dim);   // This
        ~Shape();  //Decontructor

        void input_data(double * x_values, double *f_values, double *df_values, int num_points);
        void generate_table();
        //compute_shape(P);

    private:

        int _num_points; 

        double *  _x_values;
        double *  _f_values;
        double *  _df_values;

        int  _dim;

        //Vector _P; //known points to interpolate on 
        //vector _N; //First eval (Vector of px1) or (Number of shape functions x 1)
        //vector<Vector> _DN; //Second eval (number of shape fncts x d)
        //vector<Matrix> _D2N; // each matrix is (dxd) and there are # of functions of those.
};

#endif
