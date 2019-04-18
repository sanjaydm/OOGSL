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
    _LJForce.setDim(3); 
    _eps1 = 1;
    _eps2 = 1;
    _re1 = 1;
    _re2 = 1.4;
    _Lmax = 16;
    _NP = 12;
    _k1 = 1;

    _Ntot = (_Lmax + 1 )*(_Lmax + 1);
    generateCeliaReinaQuad(20, 2*_Lmax + 2);
    computeYlmLookUp();
  }

  ~MembLJ(){
    
  }
  void LJ(Vector f12, double e);
  
  void fdf();

  void printToVTK(string filename);
  void printModes(string filename);
  void printParticleCoords(string filename);
  //void generateGausLegendreQuad(const char* filename);
  void generateAssociatedLegendreLookUp();
  void generateCeliaReinaQuad(int orderX, int orderY);
  void computeYlmLookUp();

  // Data specific to this problem
  int _NP; // number of particles
  int _Lmax; // max modes
  int _Ntot;
  
  // Lennard Jones Parameters
  double _re1;  double _eps1;
  double _re2;  double _eps2;
  double _k1;

  // Spherical harmonic lookup table
  vector<Vector> _plms;
  vector<Vector> _plms_1;
  vector<Vector> _plms_2;

  Vector _plm;
  Vector _plm_1;
  Vector _plm_2;

  Vector _plmP;
  Vector _plmP_1;


  Matrix _Ylm;
  Matrix _Ylm_t;
  Matrix _Ylm_tt;

  Matrix _Ylm_p;
  Matrix _Ylm_pp;
  Matrix _Ylm_tp;

  double _LJEnergy;
  Vector _LJForce;

  vector<double> _qThetas;
  vector<double> _qPhis;
  vector<double> _qWts;
};
#endif
