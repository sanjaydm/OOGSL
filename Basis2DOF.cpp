#include "Matrix.h"
#include <cassert>
using namespace std;

Vector Basis2DOF(Matrix& bas, Vector& coeffs){
  /* Input: =bas= with columns as basis vectors
     Input: =coeffs= contain the coefficients, the dof, for the problem
     Output: Vector coefficients that will mutliply /all/ spherical harmonics 

     If =n= is the number of basis vectors and the size of coeffs must also be =n=. 
     If =N$ is the length of basis vector, then the returned vector must also be of size =N=. 
     This number corresponds to the total number of spherical harmonics used.
  */
  Vector ret = bas*coeffs;
}
