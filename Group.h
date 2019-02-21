#ifndef GROUP_H
#define GROUP_H
#include "Vector.h"
#include <vector>
#include "Matrix.h"
class Group{
public:
  vector<Matrix> _g;

  Group();
  void constructGroup();
};
#endif
