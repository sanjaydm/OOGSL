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
  void constructGroup();
};
#endif
