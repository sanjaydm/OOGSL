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
  MembLJ(Vector x0, Vector para, Vector discPara) : Problem(x0,para) {
    _LJForce.setDim(3); 
    _k = para(0);
    _eps = para(1);
    _re = para(2);
    
    _Lmax = discPara(0);
    _NP = discPara(1);

    _Ntot = (_Lmax + 1 )*(_Lmax + 1);
    generateCeliaReinaQuad(64, 3*_Lmax + 1);
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
  double _re;  double _eps;
  double _k;

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
