#include "MultiRoot.h"

MultiRoot::MultiRoot(){;}

MultiRoot::MultiRoot(string type, Problem* prob){
  _dim = prob->_Ndof;
  _type = type;
  _problem = prob;

}
void MultiRoot::_GSLRoot_Initialize(){
  _T = gsl_multiroot_fsolver_hybrids;
  _f = {}
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

int MultiRoot:: GSLRoot_Iterate(){
  // gsl_multimin_fdfminimizer_set(_gsl_fdfminimizer,
  //                               &_fdf, (_problem->_x)._gsl_vec, _stepSize, _tol);
  int status = gsl_multimin_fdfminimizer_iterate(_gsl_fdfminimizer);
    //cout << "\033[36m" << gsl_multimin_fdfminimizer_gradient(_gsl_fdfminimizer)->data[0] << " " << "\033[0m" << endl;

  //gsl_multimin_fdfminimizer_restart(_gsl_fdfminimizer);
  //(_problem->_x) = _gsl_fdfminimizer->x; //casting gsl_vector to Vector
  return status;
}

void MultiRoot:: GSLRoot_Solve(){
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
