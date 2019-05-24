#include "Matrix.h"
#include <assert.h>

Matrix::Matrix() {
  _dim1 = 0;
  _dim2 = 0;
  _idx1 = 0;
  _idx2 = 0;
}
Matrix::Matrix(int n) : _dim1(n), _dim2(n){
  _gsl_mat = gsl_matrix_calloc(n, n);
  _idx1 = 0; _idx2 = 0;
  _vec.setDim(_dim2);
}
Matrix::Matrix(int n1, int n2) : _dim1(n1), _dim2(n2){
  _gsl_mat = gsl_matrix_calloc(n1, n2);
  _idx1 = 0; _idx2=0;
  _vec.setDim(_dim2);
}
Matrix::Matrix(int n, char c) : _dim1(n), _dim2(n){
  _gsl_mat = gsl_matrix_calloc(n,n);
  _vec.setDim(_dim2);
  _idx1 = n-1; _idx2 = n-1;
  gsl_matrix_set_identity(_gsl_mat);
}
Matrix::Matrix(const Matrix& v) {
  array<int,2> ret = const_cast<Matrix&>(v).size();
  _dim1 = ret[0]; _dim2 = ret[1]; 
  _gsl_mat = gsl_matrix_alloc(_dim1,_dim2);
  gsl_matrix_memcpy(_gsl_mat, v._gsl_mat);
  _idx1 = 0; _idx2 = 0;
  
}
Matrix::~Matrix(){
  if (_dim1 !=0 && _dim2 !=0 && _gsl_mat->owner == 1)
    gsl_matrix_free(_gsl_mat);
}
array<int,2> Matrix::size(){
  int dims[] = {_dim1, _dim2};
  array<int,2> ret;
  ret[0] = _dim1; ret[1] = _dim2;
  return ret;
}
void Matrix::setDim(int n1, int n2){
  array<int,2> num = size();
  if (num[0] ==0 && num[1] ==0 ){
    _gsl_mat = gsl_matrix_calloc(n1, n2);
    _dim1 = n1; _dim2 = n2;
    _vec.setDim(_dim2);
  }
  else {
    cout << "The matrix already created" << endl;
  }
}
void Matrix::setDim(int n){
  setDim(n,n);
}
double Matrix::norm(){
  //return gsl_blas_dnrm2(_gsl_vec);
}
// ------------------------------------------------------------
double& Matrix::operator()(int i, int j){
  return _gsl_mat->data[i*_gsl_mat->tda + j ];
      
}
Vector& Matrix::operator()(int i){
  _vview = gsl_matrix_row(_gsl_mat, i);
  _vec._gsl_vec = &_vview.vector;
  return _vec;
}
// double operator|(Matrix& v2){
//   double ret;
//   gsl_blas_ddot(_gsl_vec, v2._gsl_vec, &ret);
//   return ret;
// }
Matrix& Matrix::operator+= (Matrix& v2){
  // Note that v2 has to be passed by reference. But why?
  gsl_matrix_add(_gsl_mat, v2._gsl_mat);
  return *this;
}
Matrix& Matrix::operator-= (Matrix& v2){
  // Note that v2 has to be passed by reference. But why?
  gsl_matrix_sub(_gsl_mat, v2._gsl_mat);
  return *this;
}
Matrix& Matrix::elementMult (Matrix& v2){
  // Note that v2 has to be passed by reference. But why?
  gsl_matrix_mul_elements(_gsl_mat, v2._gsl_mat);
  return *this;
}
Matrix& Matrix::elementDiv (Matrix& v2){
  // Note that v2 has to be passed by reference. But why?
  gsl_matrix_div_elements(_gsl_mat, v2._gsl_mat);
  return *this;
}
Matrix& Matrix::operator*= (double alpha){
  gsl_matrix_scale(_gsl_mat, alpha);
  return *this;
}
Matrix& Matrix::operator*= (Matrix& m2){
  Matrix temp1 = *this;
  if (m2._gsl_mat == _gsl_mat){
    Matrix temp2 = m2;
     gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp1._gsl_mat, temp2._gsl_mat, 0.0, _gsl_mat); 
     }
  else{
     gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp1._gsl_mat, m2._gsl_mat, 0.0, _gsl_mat); 
  }
  return *this;
}
Matrix& Matrix::operator^= (int n){
  if (n==0){
    gsl_matrix_set_identity (this->_gsl_mat);
  }
  else{
    Matrix tempOriginal = *this;
    for (int i=1; i<n; i++){
      Matrix tempAccumulator = *this;
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tempAccumulator._gsl_mat, tempOriginal._gsl_mat, 0.0, _gsl_mat); 
    }
  }
  return *this;
}
Matrix Matrix::operator^ (int n){
  Matrix ret(_dim1,_dim2);
  if (n==0){
    gsl_matrix_set_identity (this->_gsl_mat);
  }
  else{
    Matrix tempOriginal = *this;
    for (int i=1; i<n; i++){
      Matrix tempAccumulator = *this;
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tempAccumulator._gsl_mat, tempOriginal._gsl_mat, 0.0, ret._gsl_mat); 
    }
  }
  return ret;
}
Matrix& Matrix::operator+= (double alpha){
  gsl_matrix_add_constant(_gsl_mat, alpha);
  return *this;
}
Matrix Matrix::operator+ (Matrix& v2){
  array <int,2> rowCol = v2.size();
  Matrix temp(rowCol[0], rowCol[1]);
  temp = *this;
  temp += v2;
  return temp;
}
Matrix Matrix::operator* (Matrix& v2){
  array <int,2> rowCol = v2.size();
  assert(_dim2 == rowCol[0]);
  Matrix temp(_dim1, rowCol[1]);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, _gsl_mat, v2._gsl_mat, 0.0, temp._gsl_mat);
  return temp;
}
Matrix Matrix::operator -(Matrix& v2){
  array <int,2> rowCol = v2.size();
  Matrix temp(rowCol[0], rowCol[1]);
  temp = *this;
  temp -= v2;
  return temp;

}
Matrix& Matrix::operator = (Matrix& v2){
  array<int,2> dims = v2.size();
  if( size()[0]==0 && size()[1] ==0){
    setDim(dims[0],dims[1]);
  }
  gsl_matrix_memcpy(_gsl_mat, v2._gsl_mat );
  _idx1 = 0; _idx2 = 0;
  return *this;
}
Matrix& Matrix::operator<< (double val){
  (*this)(_idx1,_idx2) = val;
  if (_idx2 < _dim2-1){
    _idx2++;
  }
  else {
    _idx2 = 0;
    if (_idx1 < _dim1-1){
      _idx1++;
    }
    else{
      _idx1 = 0;
    }
  }
  return *this;
}
void Matrix::print(){
  // Print to screen
  for (int i=0; i < _dim1; i++) {
    for (int j=0; j<_dim2; j++){
      cout << (*this)(i,j) << " ";
    }
    cout << endl;
  }

}
Matrix Matrix::inv(){
  assert (_dim1 == _dim2);
  Matrix copy = *this;
  Matrix inv(_dim1,_dim2);
  gsl_permutation* p = gsl_permutation_alloc (_dim1);
  int signum;
  gsl_linalg_LU_decomp (copy._gsl_mat, p, &signum);
  gsl_linalg_LU_invert (copy._gsl_mat, p, inv._gsl_mat);
  gsl_permutation_free(p);
  return inv;
}

Vector Matrix::operator*(Vector & v){
  assert(_dim2 == v._dim);
  Vector ret(_dim1);
  gsl_blas_dgemv(CblasNoTrans, 1.0, _gsl_mat, v._gsl_vec, 0.0, ret._gsl_vec);
  return ret;
}

Vector Matrix::row(int i){
  gsl_vector_view v = gsl_matrix_row(_gsl_mat, i);
  Vector ret(&v.vector);
  return ret;
}

Vector Matrix::col(int i){
  gsl_vector_view v = gsl_matrix_column(_gsl_mat, i);
  Vector ret(&v.vector);
  return ret;
}

int Matrix::rank(double tol){
  // Computes the rank of a matrix using SVD
  // The number of nonzero singular values gives the rank.
  // tol controls what is considered as zero.
  Matrix copy = *this;
  Matrix V(_dim2, _dim2);
  Vector S(_dim2);
  Vector work(_dim2);
  gsl_linalg_SV_decomp(copy._gsl_mat, V._gsl_mat, S._gsl_vec, work._gsl_vec);
  int counter = 0;
  for (int i=0; i<_dim2; i++) {
    if (fabs(S(i)) > tol)
      counter++;
  }
  return counter;
}
