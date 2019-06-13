#include "MultiRoot.h"

MultiRoot::MultiRoot(){;}

MultiRoot::MultiRoot(string type, Problem* prob){
  _dim = prob->_Ndof;
  _type = type;
  _problem = prob;

}
void MultiRoot::_GSLRoot_Initialize(){
  _T = gsl_multiroot_fsolver_hybrids;
  //_f = {};
  _f.n = _dim;
  _f.f = Problem::_C_gsl_df_solver;
  //_f.df = Problem::_C_gsl_df;
  //_f.f = Problem::_C_gsl_fdf;
  CCallbackHolder* classNParams = new CCallbackHolder;
  classNParams->params = (void*) (_problem->_para).data();
  classNParams->cls = _problem;
  _f.params = (void*) classNParams;


  _gsl_fsolver = gsl_multiroot_fsolver_alloc(_T, _f.n);
  gsl_multiroot_fsolver_set(_gsl_fsolver,
                                &_f, (_problem->_x)._gsl_vec); 
  //(_problem->_x)._gsl_vec = gsl_multimin_fdfminimizer_gradient(_gsl_fdfminimizer);
  
}

int MultiRoot:: GSLRoot_Iterate(){
  // gsl_multimin_fdfminimizer_set(_gsl_fdfminimizer,
  //                               &_fdf, (_problem->_x)._gsl_vec, _stepSize, _tol);
  int status = gsl_multiroot_fsolver_iterate(_gsl_fsolver);
    //cout << "\033[36m" << gsl_multimin_fdfminimizer_gradient(_gsl_fdfminimizer)->data[0] << " " << "\033[0m" << endl;

  //gsl_multimin_fdfminimizer_restart(_gsl_fdfminimizer);
  return status;
}

void MultiRoot:: GSLRoot_Solve(){
  int status = 1;
  do{
    status = GSLRoot_Iterate();
    
    if (status == GSL_ENOPROG) cout << status << " Eno prog\n";
    status = gsl_multiroot_test_residual (_gsl_fsolver->f, _tol);
    cout << (_problem->_df.norm()) << endl;
    if (_problem->_df.norm() < _tol) return;
    // printf ("%5d %.5f %.5f %10.5f\n", GSL_SUCCESS,
    //         gsl_vector_get (_gsl_fdfminimizer->x, 0),
    //           gsl_vector_get (_gsl_fdfminimizer->x, 1),
    //           _gsl_fdfminimizer->f);
    if (status == GSL_SUCCESS)
        printf ("Solution found.\n");
    
  }while (status == GSL_CONTINUE);
}
