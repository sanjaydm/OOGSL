//Rafael Orozco  
//interpolate a function   
//Uses hermitian matrix

#include "Shape.h"

Shape::Shape(int dim, string method_select) {
    int _dim = dim;
}

Shape::~Shape(void) {
}
    
void Shape::input_data(double * x_values, double *f_values, double *df_values, int num_points) {
    _num_points = num_points;

    _x_values = new double[num_points];
    _f_values = new double[num_points];
    _df_values = new double[num_points];

    memcpy(_x_values, x_values, num_points * sizeof(double));
    memcpy(_f_values, f_values, num_points * sizeof(double));
    memcpy(_df_values, df_values, num_points *sizeof(double));
}


void Shape::generate_table() {
    double main_table[_num_points * 2][_num_points * 2];
    double z_values[_num_points*2];

    int i;
    for (i=0; i<_num_points; i++) {
         z_values[2*i] = _x_values[i];
         z_values [2*i+1] = _x_values[i];
         main_table[2*i+1][0] = main_table[2*i][0];

         main_table[2*i][0]   = _f_values[i];//main_table[2*I][0];
         main_table[2*i+1][0] = _f_values[i];// main_table[2*I][0];
         main_table[2*i+1][1] = _df_values[i];//main_table[2*I][0];

         if (i != 0)
            main_table[2*i][1] = (main_table[2*i][0] - main_table[2*i-1][0]) / (z_values[2*i] - z_values[2*i-1]);
    }

    int k = 2 * _num_points + 1;
    int j;
    for (i=2; i<k; i++) {
        for (j=2; j<i; j++) {
            main_table[i][j] = ( main_table[i][j - 1] - main_table[i - 1][j - 1] ) / ( z_values[i] - z_values[i - j] );
        }
    }


    //Evaluate
    double interp_value;
    double x_to_interp = 0.25;

    interp_value = main_table[k][k] * (x_to_interp - z_values[k-1]);

    for (i=2; i<=k; i++) {
        j = k - i + 1;
        interp_value = (interp_value + main_table[j][j]) * (x_to_interp - z_values[j-1]);
    }
    interp_value = interp_value + main_table[0][0];
    printf("Value: %f\n", interp_value);

}
