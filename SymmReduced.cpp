#include "SymmReduced.h"

void SymmReduced::convertBasis2Full(){
  /* Input: =bas= with columns as basis vectors
     Input: =coeffs= contain the coefficients, the dof, for the problem
     Output: Vector coefficients that will mutliply /all/ spherical harmonics 

     If =n= is the number of basis vectors and the size of coeffs must also be =n=. 
     If =N$ is the length of basis vector, then the returned vector must also be of size =N=. 
     This number corresponds to the total number of spherical harmonics used.
  */
  _fullX = _bas*_x;
}
void SymmReduced::fdf(){
  //Convert symmetry reduced dofs to full
  convertBasis2Full();
  _prob->_x = _fullX;

  // Pass flags to full problem
  _prob->_fFlag = _fFlag;
  _prob->_dfFlag = _dfFlag;

  // compute f, df, d2f for full problem
  _prob->fdf();

  // Store values of returned from full into reduced
  _f = _prob->_f;
  _df = _basT*_prob->_df;
  //_d2f = _prob->_d2f; Projection operators P^T d2f P, here?
  
}
