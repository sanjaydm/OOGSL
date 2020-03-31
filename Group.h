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

class D2: public Group{
public:
    Matrix _g2;
    Matrix _g22;
    D2(Matrix g2, Matrix g22);
    void constructGroup();
};

class D3: public Group{
public:
    Matrix _g2;
    Matrix _g3;
    D3(Matrix g2, Matrix g3);
    void constructGroup();
};

class D4: public Group{
public:
    Matrix _g2;
    Matrix _g4;
    D4(Matrix g2, Matrix g4);
    void constructGroup();
};

class D5: public Group{
public:
    Matrix _g2d;
    Matrix _g5;
    D5(Matrix g2d, Matrix g5);
    void constructGroup();
};

class D6: public Group{
public:
    Matrix _g2;
    Matrix _g6;
    D6(Matrix g2, Matrix g6);
    void constructGroup();
};

class D8: public Group{
public:
    Matrix _g2;
    Matrix _g8;
    D6(Matrix g2, Matrix g8);
    void constructGroup();
};

class D9: public Group{
public:
    Matrix _g2;
    Matrix _g9;
    D9(Matrix g2, Matrix g9);
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
    Matrix _P;
    Projection() {;};
    Projection(Group group);
    void computeProjection();
};





#endif
