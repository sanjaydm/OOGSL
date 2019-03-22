#ifndef MembLJ_H
#define MembLJ_H
#include <iostream>
#include <string.h>
#include "Vector.h"
#include "Matrix.h"
#include "Problem.h"
#include "MultiMin.h"
#include "Quadrature.h"
#include "Shape.h"
#include "gausslegendre.h"

#define PI (3.141592653589793)
#define CSPHASE 1
#define S2 1.4142135623730951
#define ENERGY 0
#define RESIDUE 1
#define PARTICLE 2

class MembLJ: public Problem {
public:
 MembLJ(Vector x0, Vector para) : Problem(x0,para) {
    _quad.setOrder(4);
    _lclResidue.setDim(4); 
    _eps1 = 1;
    _eps2 = 1;
    _re1 = 1;
    _re2 = 1.4;
    
  }

  ~MembLJ(){
    
  }
  void LJ(double r);
  
  void fdf();

  void printToVTK(string filename);
  void printModes(string filename);
  void printParticleCoords(string filename);
  void printModes(string filename);

  // Data specific to this problem
  Quadrature _quad;
  double _lclEnergyDensity;
  Vector _lclResidue;
  int _NP; // number of particles
  int _Lmax; // max modes
  int _Ntot;
  
  // Lennard Jones Parameters
  double _re1;  double _eps1;
  double _re2;  double _eps2;

  // Spherical harmonic lookup table
  vector<Vector> _plms;
  vector<Vector> _plms_1;
  vector<Vector> _plms_2;

  Vector _plm;
  Vector _plm_1;
  Vector _plm_2;

  Vector _plmP;
  Vector _plmP_1;
  Vector _plmPB;


  Matrix _Ylm;
  Matrix _Ylm_t;
  Matrix _Ylm_tt;

  Matrix _Ylm_p;
  Matrix _Ylm_pp;
  Matrix _Ylm_tp;

  vector<double> qThetas;
  vector<double> qPhis;
  vector<double> qWts;
};
#endif
