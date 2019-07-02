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

class Icosahedral: public Group{
public:
    Matrix _g2;
    Matrix _g3;
    Icosahedral(Matrix g2, Matrix g3);
    void constructGroup();
};

class D5: public Group{
public:
    Matrix _g2d;
    Matrix _g5;
    D5(Matrix g2d, Matrix g5);
    void constructGroup();
};

class Tetrahedral: public Group{
public:
    Matrix _g2;
    Matrix _g3d;
    Tetrahedral(Matrix g2, Matrix g3d);
    void constructGroup();
};

class Projection{
public:
    Group _group;
    Projection(Group group);
    Matrix computeProjection();
};





#endif
