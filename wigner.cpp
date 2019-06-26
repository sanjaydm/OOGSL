// Generated on Sun Mar 29 2015 08:35:13 for TomoMiner C++ core by   doxygen 1.8.4
// Courtesy: http://web.cmb.usc.edu/people/alber/Software/tomominer/docs/cpp/wigner_8cpp_source.html
// Adapted to work with our Matrix and Vector classes

#include <vector>
#include <cassert>
#include "Matrix.h"
#include "Vector.h"
using namespace std;
//#define parity(X) ((X % 2) == 0 ? 1 : -1)
#define i(m, L) m+L

inline
int parity(int x)
{
  if(x % 2 == 0) return 1;
  return -1;
}

Vector computeEulerAngles(Matrix& Q){
  // Computed wrt ZYZ convention
  Vector abc(3); //alpha, beta, gamma
  abc(0) = atan2(Q(1,2),Q(0,2));
  abc(2) = atan2(Q(2,1),-Q(2,0));
  double s2 = sqrt(Q(2,0)*Q(2,0)+Q(2,1)*Q(2,1));
  abc(1) = atan2(s2,Q(2,2));
  return abc;
}
class Wigner_d {
public:

  Wigner_d(int Lmax, double theta): _Lmax(Lmax), _theta(theta) {
    compute();
  }
  double get_d(int L, int m, int k) {return _d[L]( i(m,L), i(k,L));}
  void compute() {
    // Based on the work of J. Chem Phys. 124, 144115 (2006)
    int L = _Lmax;
    if (_theta==0) {
      for(int l = 0; l <= L; l++)
	_d.push_back(Matrix(2*l+1,'i'));
    }
  
    else if (_theta > 0 && _theta <= M_PI/2) {
      _d = wigner_d0Pi2(_theta, L);
    }
  
    else if (_theta > M_PI/2 && _theta < M_PI){
      
      for(int l = 0; l <= L; l++) {
	_d.push_back(Matrix(2*l+1,2*l+1));
      }
      vector <Matrix> Dret = wigner_d0Pi2(M_PI - _theta, L);
      for (int l=0; l<= L; l ++) {
	for (int m= -l; m <= l; m++) {
	  for (int k=-l; k<= l; k++) {
	    _d[l]( i(m,l), i(k,l) ) = parity(l+k)* Dret[l]( i(-m,l), i(k,l) );
	  }
	}
      }
    }
  
    else if (_theta == M_PI){
      for(int l = 0; l <= L; l++)
	_d.push_back(Matrix(2*l+1,'i'));

      for (int l=0; l<= L; l ++) {
	for (int m= -l; m <= l; m++) {
	  for (int k=-l; k<= l; k++) {
	    _d[l]( i(m,l), i(k,l) ) = parity(l+k)*(-m == k);
	  }
	}
      }
    }
  
    else if (_theta > M_PI && _theta <= 3*M_PI/2) {
      
      for(int l = 0; l <= L; l++) {
	_d.push_back(Matrix(2*l+1,2*l+1));
      }
      vector <Matrix> Dret = wigner_d0Pi2(_theta - M_PI, L);
      // .....
      for (int l=0; l<= L; l ++) {
	for (int m= -l; m <= l; m++) {
	  for (int k=-l; k<= l; k++) {
	     _d[l]( i(m,l), i(k,l) ) = parity(l+k)* Dret[l]( i(m,l), i(-k,l) );
	  }
	}
      }
    }
  
    else if (_theta >= 3*M_PI/2 && _theta <= 2*M_PI){
      for(int l = 0; l <= L; l++) {
	_d.push_back(Matrix(2*l+1,2*l+1));
      }
      vector <Matrix> Dret = wigner_d0Pi2(2*M_PI - _theta, L);
      // .....
      for (int l=0; l<= L; l ++) {
	for (int m= -l; m <= l; m++) {
	  for (int k=-l; k<= l; k++) {
	    _d[l]( i(m,l), i(k,l) ) = parity(m+k)* Dret[l]( i(m,l), i(k,l) );
	  }
	}
      }
    
    }

  }
  
  vector<Matrix> wigner_d0Pi2(double theta, int L) {

    // The Wigner D-matrix of order l, is a matrix with entries from -l:l in both
    // dimensions.  The code work for 0 < theta <= pi/2
    int l, k, m;
  
    // According to paper, recursion formulas are only valid for 0 < theta <= pi/2.0. 
    // Alternatives are given to adjust if we are out of range (Eqn #30).

    assert( 0 < theta && theta <= M_PI/2.0 );

    vector<Matrix> D;
    // generate matrices to be filled in.
    for(l = 0; l <= L; l++)
      D.push_back(Matrix(2*l+1,2*l+1));
 
    // precompute sin/cos.
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
 
    Matrix g(L+1, L+1);
    g(0,0) = 1;
   
    // g recursion Eqn #29
    for(l = 1; l <= L; l++)
      {
	g(l,0) = sqrt( ((double)(2*l-1))/(2*l) ) * g(l-1,0);
     
	for(m = 1; m <= l; m++)
	  {
	    g(l,m) = sqrt( ((double)(l-m+1))/(l+m) )*g(l,m-1);
	  }
      }
   
    // precompute sin/cos powers.
    Vector c1(L+1);
    Vector c2(L+1);
    c1(0) = 1;
    c2(0) = 1;
    for(int i=1; i<=L;i++)
      {
	c1(i) = c1(i-1) * (1.0 + cos_theta);
	c2(i) = c2(i-1) * sin_theta;
      }
 
 
    // Fill in k=l. Eqn #28
    for(l=0; l <= L; l++)
      for(m=0; m <= l; m++)
	D[l](m+l,l+l) = parity(m+l) * g(l,m)* c1(m) * c2(l-m);
   
    // Precompute sin(theta) / (1+cos(theta)).
    double c3 = sin_theta / (1.0 + cos_theta);
 
    // Precompute 1.0/sqrt( (l*(l+1) - k*(k-1)) ).
    Matrix C4(L+1, 2*L+1);
 
    for(l = 0; l <= L; l++)
      for(k = l; k > -l; k--)
	C4(l,k+L) = 1.0/sqrt( l*(l+1) - k*(k-1) );
 
    // Fill in bottom row: Eqn #26
    for(l=0; l<= L; l++)
      for(k = l; k > -l; k--)
	D[l](l+l,k-1+l) = (l+k) * C4(l,k+L) * c3 * D[l](l+l,k+l);
 
    // Fill in from bottom up. Eqn #25
    for(l = 0; l <= L; l++)
      for(m = l-1; m >= 0; m--)
	for(k = l; k > -l; k--)
	  D[l](m+l,k-1+l) = sqrt(l*(l+1)-m*(m+1)) * C4(l,k+L) * D[l](m+1+l,k+l) + (m+k) * C4(l,k+L) * c3 * D[l](m+l,k+l);
 
    // fill in negative m. Eqn #27
    for(l = 0; l <= L; l++)
      for(m = -l; m < 0; m++)
	for(k = -l; k <= l; k++)
	  D[l](m+l,k+l) = parity(m+k) * D[l](l-m,l-k);
     
    // Adjust the sign pattern.
    for (l=0; l<=L; l++){
      array<int,2> rc = D[l].size();
      for (int r=0; r< rc[0]; r++){
	for (int c=0; c<rc[1]; c++){
	  if ( (r+c)%2 != 0)
	    D[l](r,c) *= -D[l](r,c);
	}
      }
    
    }
    return D;
  }

  
  int _Lmax;
  vector<Matrix> _d;
  double _theta;
  
}; 
