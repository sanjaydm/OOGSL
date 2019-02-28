//Rafael Orozco  
//interpolate a function   
//Uses hermitian matrix

#include "Shape.h"

Shape::Shape(int dim, double point) {
    int _dim = dim;
    double _point = point;
}

Shape::~Shape(void) {
}

void Shape :: compute_shape_functions() {
    printf("made it here\n");
    fflush(stdout);
    _N.setDim(4);
    _DN.setDim(4);
    _D2N.setDim(4);

    //First
    _N(0) = (1.0/4.0)*(2.0 + _point)*(1.0 - _point)*(1.0 - _point);
    _N(1) = (1.0/4.0)*(2.0 - _point)*(1.0 + _point)*(1.0 + _point);
    _N(2) = (1.0/4.0)*(1.0 + _point)*(1.0 - _point)*(1.0 -_point);
    _N(3) = (1.0/4.0)*(_point - 1.0)*(1.0 + _point)*(1.0 + _point);

    //Derivative
    _DN(0) = (3.0/4.0)*(_point*_point - 1.0);
    _DN(1) = -(3.0/4.0)*(_point*_point - 1.0);
    _DN(2) = (1.0/4.0)*(3.0*_point*_point - 2.0*_point - 1.0);
    _DN(3) = (1.0/4.0)*(3.0*_point*_point + 2.0*_point - 1.0);

    //Derivative 2
    _D2N(0) = (3.0/2.0)*_point; 
    _D2N(1) = (3.0/2.0)*_point; 
    _D2N(2) = (1.0/2.0)*(3.0*_point - 1.0);
    _D2N(3) = (1.0/2.0)*(3.0*_point + 1.0);
}

Vector& Shape :: get_shape_functions_N() {
    
    static Vector functions_N(4);
    functions_N = _N;
    return functions_N;
}

Vector& Shape :: get_shape_functions_DN () {
    static Vector functions_DN(4);
    functions_DN = _DN;
    return functions_DN;
}

Vector& Shape :: get_shape_functions_D2N () {
    static Vector functions_D2N(4);
    functions_D2N = _D2N;
    return functions_D2N;
}



