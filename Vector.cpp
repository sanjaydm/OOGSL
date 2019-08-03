#include "Vector.h"
#include <iomanip>

using namespace std;

Vector::Vector() {
  //cout << "here\n";
  _dim = 0;
  _idx = 0;
}
Vector::Vector(int n) : _dim(n){
  _gsl_vec = gsl_vector_calloc(n);
  _idx = 0;
}
Vector::Vector(const Vector& v) {
  //cout << "in copy constructor" << endl;
  _dim = const_cast<Vector&>(v).size();
  _idx = _dim;
  _gsl_vec = gsl_vector_calloc(_dim);
  gsl_vector_memcpy(_gsl_vec, v._gsl_vec);
}
Vector::Vector(gsl_vector* v) {
  _dim = v->size;
  _idx = _dim;
  _gsl_vec = gsl_vector_calloc(_dim);
  _gsl_vec->data = v->data;
  //gsl_vector_memcpy(_gsl_vec, v);
}
Vector::Vector(double* v, int size) : _dim(size) {
  _gsl_vec = gsl_vector_calloc(_dim);
  _gsl_vec->data = v;
  _idx = _dim;
}
Vector::Vector(void* v, int size){
  double* temp = (double*) v;
  Vector(temp, size);
}
Vector:: ~Vector(){
  if (_dim != 0 && _gsl_vec->owner == 1){
    gsl_vector_free(_gsl_vec);
    _dim = 0;
    _idx = 0;
  }
}
int Vector::size(){
  return _dim;
}
void Vector::setDim(int n){
  if (size()== 0 ){
    _gsl_vec = gsl_vector_calloc(n);
    _dim = n;
    //gsl_vector_set_zero(_gsl_vec);
  }
  else {
    cout << "The vector is already created" << endl;
  }
}
double Vector::norm(){
  return gsl_blas_dnrm2(_gsl_vec);
}
double& Vector::operator()(int i){
  // Note the return is a reference. This lets us make assignment
  return _gsl_vec->data[i];
}
double Vector::operator|(Vector& v2){
  double ret;
  gsl_blas_ddot(_gsl_vec, v2._gsl_vec, &ret);
  return ret;
}
Vector& Vector::operator+= (Vector& v2){
  gsl_vector_add(_gsl_vec, v2._gsl_vec);
  return *this;
}
Vector& Vector::operator-= (Vector& v2){
  gsl_vector_sub(_gsl_vec, v2._gsl_vec);
  return *this;
}
Vector& Vector::operator*= (Vector& v2){
  gsl_vector_mul(_gsl_vec, v2._gsl_vec);
  return *this;
}
Vector& Vector::operator/= (Vector& v2){
  gsl_vector_div(_gsl_vec, v2._gsl_vec);
  return *this;
}
Vector& Vector::operator*= (double alpha){
  gsl_vector_scale(_gsl_vec, alpha);
  return *this;
}
Vector& Vector::operator/= (double alpha){
  *this *= (1/alpha);
  return *this;
}
Vector& Vector::operator+= (double alpha){
  gsl_vector_add_constant(_gsl_vec, alpha);
  return *this;
}
Vector Vector::operator+ (Vector& v2){
  int m = min(v2.size(), this->size());
  if (m  == 0 && this->size()==0) this->setDim(v2.size());
  else if (m == 0 && v2.size()== 0) v2.setDim(this->size());
        
  Vector temp(v2.size());
  temp = *this;
  temp += v2;
  return temp;
}
Vector Vector::operator- (Vector& v2){
  int m = min(v2.size(), this->size());
  if (m  == 0 && this->size()==0) this->setDim(v2.size());
  else if (m == 0 && v2.size()== 0) v2.setDim(this->size());
  
  Vector temp(v2.size());
  temp = *this;
  temp -= v2;
  return temp;
}
Vector& Vector::operator= (const Vector& v2){
  //cout << "assignment constructor\n" ;
  if(size()==0){
    setDim(const_cast<Vector&>(v2).size());
  }
  gsl_vector_memcpy(_gsl_vec, v2._gsl_vec );
  _idx = _dim-1;
  return *this;
}
Vector& Vector::operator= (gsl_vector* v2){
  if(size()==0){
    setDim(v2->size);
  }
  gsl_vector_memcpy(_gsl_vec, v2 );
  return *this;
}
Vector& Vector::operator<< (double val){
  (*this)(_idx) = val;
  if (_idx == _dim-1) _idx=0;
  _idx++;
  return *this;
}
double* Vector::data(){
  return _gsl_vec->data;
}
void Vector::print(){
  // Print to screen
  cout << scientific << setprecision(10);
  cout << "[\n";
  for (int i=0; i<_dim; i++){
    cout << (*this)(i) << ", " << endl;
  }
  cout << "]\n";
}
void Vector::write(){
  // Write to file
}

void Vector::setZero(){
  gsl_vector_set_zero(_gsl_vec);
}

Vector Vector::view(int i, int j){
  /* Subvector from i to j*/
  assert (j>i);
  _view = gsl_vector_subvector(_gsl_vec, i, j-i);
  Vector ret(&_view.vector);
  return ret;
}


