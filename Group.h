#ifndef GROUP_H
#define GROUP_H
#include "Vector.h"
#include <vector>
#include "Matrix.h"
class Group{
public:
  
  vector<Matrix> _g;
  vector<Matrix> _gen;
  Group();
  void constructGroup(){;}
  int order() {return _g.size();}
  int numGenerators() {return _gen.size();}
  void listElements();
};
class Cyclic: public Group{
public:
  int _n;
  Cyclic(int n);
  Cyclic(int n, Matrix gen);
  void constructGroup();
};
#endif
