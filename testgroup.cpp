#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "Problem.h"
#include "Group.h"
#include "MultiMin.h"

int main(int argc, char** argv){


  Matrix I(3,'i');
  Cyclic Z3(3,I);
  Z3._gen[0].print();
  cout << "order = " << Z3.order() << endl;
  Z3.listElements();
  // Matrix P = I._g[0] + I._g[1] + I._g[2];
  // P.print();
  // cout << "rank = " << P.rank() << endl;;
  // Constructing projection operator
  
  // Vector v(5);
  // v << 1 << 2 << 3 << 4 << 5;
  // v.print();
  // Vector w = v.view(0,3);
  // cout <<"asfasdf\n";
  // w.print();

  // Matrix m(5,5);
  // m << 1 << 1 << 1 << 0 << 0
  //   << 0 << 1 << 1 << 0 << 0
  //   << 0 << 0 << 1 << 1 << 0
  //   << 0 << 0 << 0 << 0 << 1
  //   << 1 << 1 << 1 << 1 <<1;
  // m.print();
  // Vector x(5);
  
  // (m*v).print();

  // m.col(0).print();
  return 0;
}
