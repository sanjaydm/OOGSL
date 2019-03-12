#include "MultiMin.h"

MultiMin::MultiMin(){;}

MultiMin::MultiMin(string type, Problem* prob){
  _dim = prob->_Ndof;
  _type = type;
  _problem = prob;

}
void MultiMin::_LBFGSB_Initialize(){
  int n = _dim;
  _l.setDim(n);
  _u.setDim(n);
  _nbd.assign(n,0);
  _wa.assign(2*_m*n+5*n+11*_m*_m+8*_m,0.0);
  _iwa.assign(3*n, 0);

  for(int i=0; i<60; i++) {_task[i] = _csave[i] = '\0';}
  for(int i=0; i<4; i++) _lsave[i]=0;
  for(int i=0; i<29; i++) _dsave[i]=0.0;
  for(int i=0; i<44; i++) _isave[i]=0;

  sprintf(_task,"START");

}

void MultiMin::LBFGSB_Iterate(){
  _problem->fdf();
  setulb_(&_dim, &_m, (_problem->_x).data(), _l.data(), _u.data(), _nbd.data(), 
          &(_problem->_f), (_problem->_df).data(), &_factr, &_tol, _wa.data(),_iwa.data(), &(_task[0]), &_iprint, &(_csave[0]),
              &(_lsave[0]),&(_isave[0]),&(_dsave[0]));
  }

void MultiMin::LBFGSB_Solve(){
  while (true){
    LBFGSB_Iterate();
    if (strncmp(_task, "FG", 2) == 0) continue;
    else if (strncmp(_task, "CONV", 4) == 0) break;
    else if (strncmp(_task, "ABNORM", 6) == 0) break;
  }
}

void MultiMin::_GSLMin_Initialize(){
  _T = gsl_multimin_fdfminimizer_vector_bfgs;
  _fdf.n = _dim;
  _fdf.f = Problem::_C_gsl_f;
  _fdf.df = Problem::_C_gsl_df;
  _fdf.fdf = Problem::_C_gsl_fdf;
  CCallbackHolder* classNParams = new CCallbackHolder;
  classNParams->params = (void*) (_problem->_para).data();
  classNParams->cls = _problem;
  _fdf.params = (void*) classNParams;
  _T = gsl_multimin_fdfminimizer_vector_bfgs;
  //_T = gsl_multimin_fdfminimizer_steepest_descent;
  //_T = gsl_multimin_fdfminimizer_conjugate_fr;
  _gsl_fdfminimizer = gsl_multimin_fdfminimizer_alloc(_T, _fdf.n);
  gsl_multimin_fdfminimizer_set(_gsl_fdfminimizer,
                                &_fdf, (_problem->_x)._gsl_vec, _stepSize, _tol); 
(_problem->_x)._gsl_vec = gsl_multimin_fdfminimizer_gradient(_gsl_fdfminimizer);
  
}

int MultiMin:: GSLMin_Iterate(){
  // gsl_multimin_fdfminimizer_set(_gsl_fdfminimizer,
  //                               &_fdf, (_problem->_x)._gsl_vec, _stepSize, _tol);
  int status = gsl_multimin_fdfminimizer_iterate(_gsl_fdfminimizer);
    //cout << "\033[36m" << gsl_multimin_fdfminimizer_gradient(_gsl_fdfminimizer)->data[0] << " " << "\033[0m" << endl;

  //gsl_multimin_fdfminimizer_restart(_gsl_fdfminimizer);
  //(_problem->_x) = _gsl_fdfminimizer->x; //casting gsl_vector to Vector
  return status;
}

void MultiMin:: GSLMin_Solve(){
  int status = 1;
  do{
    status = GSLMin_Iterate();
    
    if (status == GSL_ENOPROG) cout << status << " Eno prog\n";
    status = gsl_multimin_test_gradient (_gsl_fdfminimizer->gradient, _tol);
    // printf ("%5d %.5f %.5f %10.5f\n", GSL_SUCCESS,
    //         gsl_vector_get (_gsl_fdfminimizer->x, 0),
    //           gsl_vector_get (_gsl_fdfminimizer->x, 1),
    //           _gsl_fdfminimizer->f);
    if (status == GSL_SUCCESS)
        printf ("Minimum found.\n");
    
  }while (status == GSL_CONTINUE);
}
