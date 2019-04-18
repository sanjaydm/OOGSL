#include<iostream>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_cblas.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_sf.h>
#include <unistd.h>
#include <cstdlib>
#include<string.h>
#include <sstream>
#include "myfdf.h"
#include <fstream>
#include <iomanip>
#include <gsl/gsl_multimin.h>
#include <vector>
#include <gsl/gsl_multiroots.h>
#define ENERGY 0
#define RESIDUE 1
#define PARTICLE 2

extern "C" void setulb_(int * n, int *m, 
			double * x, double * l, double * u, 
			int * nbd, 
			double * f, double * g, 
			double * factr, double * pgtol,
			double * wa, 
			int * iwa,
			char * task, 
			int * iprint,  
			char * csave, 
			int * lsave, 
			int * isave, 
			double * dsave );
using namespace std;
#define PI (3.141592653589793)
#define CSPHASE 1
#define S2 1.4142135623730951
int N = 20;  // No of modes for r
int NP= 254;//12 + 10*(49-1);//12 + 10*(49 - 1);
int arr_size = gsl_sf_legendre_array_n(N);
int Ntot = (N+1)*(N+1);

vector <double*> plms;
vector <double*> plms_1;
vector <double*> plms_2;
  
double* plm = new double[arr_size]; 
double* plm_1 = new double[arr_size];
double* plm_2 =new double[arr_size];
double* plmP = new double[arr_size]; 
double* plmP_1 = new double[arr_size];
double* plmPB = new double[arr_size]; 

vector<double> qThetas;
vector<double> qPhis;
vector<double> qWts;

gsl_matrix* Ylm;
gsl_matrix* Ylm_t;
gsl_matrix* Ylm_tt;

gsl_matrix* Ylm_p;
gsl_matrix* Ylm_pp;
gsl_matrix* Ylm_tp;

int main(int argc,char** argv){
  // Parameters of the model
  double ri =   2*sqrt(4.0/NP)*(1-0.5); //1.303648e-01; // Start r_e //pow(2.0, 1./6.)*0.08757413;//
  double rmax = 2*sqrt(4.0/NP)*(1+0.05);                  // End r_e
  double drm = 2*sqrt(4.0/NP)*0.01;                  // Increment
  double ki = 5.0;
  double kmin = 1.0;
  double kmax = 5.0;
  double dk = -5;

  cout << "number of particles = " << NP << endl;
  // Generate Quadrature
  int ordX = 32; 
  int ordY = 2*N+2;

  generateCeliaReinaQuad("quad.dat", ordX, ordY);
  //generateGausLegendreQuad("quad.dat");
  string lebfile("quad.dat");
  FILE* pfile;
  pfile=fopen(lebfile.c_str(),"r");
  double t, p, wt;
  cout << "Reading quadrature file ...";
  while(fscanf(pfile,"%lf %lf %lf",&p,&t,&wt)!=EOF){
    qThetas.push_back(t);
    qPhis.push_back(p);
    qWts.push_back(wt);
  }
  fclose(pfile);
  cout << " done" << endl;


  //START: Generating associated legendre lookup
  cout << "Generating associated legendre look up table ... ";
  for (int i=0; i < qThetas.size(); i++){
    double ct = cos(qThetas[i]);
    double* temp;
    double* temp_1;
    double* temp_2;
    temp = new double[arr_size];
    temp_1 = new double[arr_size];
    temp_2 = new double[arr_size];
    gsl_sf_legendre_deriv2_alt_array_e(GSL_SF_LEGENDRE_SPHARM, N, ct,
				       CSPHASE,  temp, temp_1, temp_2);
    plms.push_back(temp);
    plms_1.push_back(temp_1);
    plms_2.push_back(temp_2);
  }
  cout << "done\n";
  // END: Lengendre lookup
  
  int Tot_len=Ntot; //Total size of the discretization 
  gsl_vector* in=gsl_vector_calloc(Ntot+2*NP-2);     // -2 because one particle is fixed
  gsl_vector* out=gsl_vector_calloc(Ntot+2*NP-2);
  gsl_vector* out_FD=gsl_vector_calloc(Ntot+2*NP-2);
  string filename = "Particle";
  gsl_vector_set_zero(in);
  double para[3];
  para[0] = ki;   //bending stiffness kappa
  para[1] = 200.; //penatly
  para[2] = ri; //rm = equilibirum distance
  srand(time(NULL));
  for (int i=0; i<Tot_len; i++){
    in->data[i] = 0.00*rand()/RAND_MAX;
  }
  for (int j=0; j< NP-1; j++) {
    in->data[Tot_len + 2*j] = /*PI/(NP-1)*j+0.3;*/double(rand())/RAND_MAX*PI;
    in->data[Tot_len + 2*j+1 ] = /*PI/(NP-1)*j;*/(0.5-double(rand())/RAND_MAX)*2*PI;
    //cout << in->data[Tot_len + 2*j]*180/PI << " " << in->data[Tot_len + 2*j + 1]*180/PI << endl;
  }

  // Restart feature
  // Read modes information
  // string modeFile("Particle_N_502_42_modes.txt");
  // FILE* modeFile_ptr;
  // modeFile_ptr=fopen(modeFile.c_str(),"r");
  // double clm;
  // cout << "Reading mode file ...";
  // int modeCtr = 0;
  // fscanf(modeFile_ptr,"%lf\n",&para[0]);
  // fscanf(modeFile_ptr,"%lf\n",&para[1]);
  // fscanf(modeFile_ptr,"%lf\n",&para[2]);
  // while(fscanf(modeFile_ptr,"%lf\n",&clm)!=EOF){
  //   in->data[modeCtr] = clm;
  //   modeCtr ++;
  // }
  // fclose(modeFile_ptr);
  // cout << " done" << endl;

  //Read particle positions
  string particleFile("ThirdCovered.dat");
  FILE* particleFile_ptr;
  particleFile_ptr=fopen(particleFile.c_str(),"r");
  cout << "Reading particle positions ...\n";
  int pj = 0;
  double thetaPhi;
  fscanf(particleFile_ptr,"%lf\n",&para[0]);
  fscanf(particleFile_ptr,"%lf\n",&para[1]);
  fscanf(particleFile_ptr,"%lf\n",&para[2]);

  while(fscanf(particleFile_ptr,"%lf\n",&thetaPhi )!=EOF){
    if (pj% 2 ==1) thetaPhi -= PI;
    in->data[Tot_len+pj] = thetaPhi;// + thetaPhi * 0.1*(rand())/RAND_MAX;;
    //if (fabs(thetaPhi)<=1e-6 && (pj%2)==0) thetaPhi = 1e-2;
    //if (fabs(thetaPhi-PI)<=1e-6 && (pj%2)==0) thetaPhi = PI-1e-2;
    cout << thetaPhi*180/PI << " ";
    if (pj%2 == 1 ) cout << "\n";
    pj++;
  }
  fclose(particleFile_ptr);
  cout << " done" << endl;

  cout << "\033[32m" << para[0] << "," << para[1] <<", " << para[2] <<"\033[0m\n";
  // //Check Performance 
  // time_t start, end;
  // time(&start);    
  // for (int i=0; i<400; i++){
  //   double energy = 0;
  //   //energy = myf(in, para);
  //   myfdf(in, para, &energy, out );
  // }
  // time (&end);
  // cout << endl << "\033[31m solved in " << difftime(end,start) << " s\033[0m" << endl;

  
  // //Check consistency
  // time_t start, end;
  // time(&start);    
  // double energy = 0;
  // myfdf(in, para, &energy, out );
  // //gsl_vector_fprintf(stdout, out, "%le");
  // cout << "energy = "<< energy << endl;
  // cout << "-----------\n" ;
  // cout << "Residue norm = " << gsl_blas_dnrm2(out) << endl;
  
  // double eps = 1e-8;
  // //double en0 = myf(in, para);
  // for (int i=0; i< Tot_len+2*NP-2; i++){
  //   in->data[i] += eps;
  //   double enp = myf(in, para);
  //   in->data[i] -= 2*eps;
  //   double enm = myf(in, para);
  //   out_FD->data[i] = (enp-enm)/(2*eps);
  //   //cout << "in = " << in->data[i] << ", en=" << en << endl;
  //   in->data[i] += eps;
  //   //cout << i << " of " << Tot_len+2*NP << endl;
  // }
  // gsl_vector_sub(out_FD, out);
  // gsl_vector_fprintf(stdout, out_FD, "%le");
  // cout << scientific << setprecision(6);
  // //cout << "|F_FD - F_ac| = " <<  gsl_blas_dnrm2(out_FD) << endl;
  // cout << "|F_FD - F_ac|/|F_ac| = " <<  gsl_blas_dnrm2(out_FD)/gsl_blas_dnrm2(out) << endl;
  // time (&end);
  // cout << endl << "\033[31m solved in " << difftime(end,start) << " s\033[0m" << endl;
  // return 0;

  // LBFGS Minmization
  ofstream outFile("Energy.txt", ofstream::out);
  outFile << "k" <<"\t" <<"r" << "\t"<<"Par-energy" <<"\t" <<"Tot-energy"<<"\t" << "Projg \n";
  outFile.close();
  int cntr = 0;

  // cout << "\033[34m " << myf(in, para) << "\033[0m" << endl ;

  // return 0;

 loop:
  while ( para[2] >= ri && para[2] <= rmax) {
    double f;
    gsl_vector* g=gsl_vector_calloc(Ntot+2*NP-2);
    vector<double > l;
    vector<double > u;
    vector<int > nbd;
    vector<double > wa;
    vector<int > iwa;
    char task[60];
    int n, m;
    int maxIterations=-1;
    double factr, pgtol, projg;

    projg = 0.0;
    pgtol = 1e-5;
    m = 30;
    factr = 1;
  
    int iprint = 50;

    int iterNo;

    n = Ntot+2*NP-2;
    f = 0.0;
    //x.assign(n, 0.0);		
    //g.assign(n, 0.0);
    l.assign(n, 0.0);
    u.assign(n, 0.0);
    nbd.assign(n, 0);
    wa.assign(2*m*n+5*n+11*m*m+8*m, 0.0); // for LBFGS 3.0
    //wa.assign(2*m*n+4*n+12*m*m+12*m, 0.0); // for LBFGS 2.0
    iwa.assign(3*n, 0);
    
    char csave[60];
    for(int i=0; i<60; i++) task[i] = csave[i] = '\0';
    
    int lsave[4];
    for(int i=0; i<4; i++) lsave[i]=0;
    
    double dsave[29];
    for(int i=0; i<29; i++) dsave[i]=0.0;
    
    int isave[44];
    for(int i=0; i<44; i++) isave[i]=0;

    // Set bounds
    for (int j=0; j< NP-1; j++) {
      l[Tot_len + 2*j] = 0.00001;
      u[Tot_len + 2*j] = PI-0.00001;
      l[Tot_len + 2*j + 1] = -PI; //0
      u[Tot_len + 2*j + 1] = PI;//2*PI;
      nbd[Tot_len + 2*j] = 2;
      nbd[Tot_len + 2*j+1] = 2;        
    }
    // We start the iteration by initializing task.
    sprintf(task,"START");

    //this->setFieldAndCompute();
    compute(in, para, &f, g, RESIDUE);

    // ------- the beginning of the loop ----------
    while(true) {

      // This is the call to the L-BFGS-B code.
      setulb_(&n, &m, in->data, l.data(), u.data(), nbd.data(), 
              &f, g->data, &factr, &pgtol, wa.data(),iwa.data(), 
              &(task[0]), &iprint, &(csave[0]),
              &(lsave[0]),&(isave[0]),&(dsave[0])); 

      if(strncmp(task,"FG",2) == 0) {
        // The minimization routine has returned to request the
        // function f and gradient g values at the current x.

        //this->setFieldAndCompute();
        compute(in, para, &f, g, RESIDUE);
        // Go back to the minimization routine.
        continue;
      }
      else if( strncmp(task,"NEW_X",5) == 0 ) {
        // stop if maximum number of iterations has been reached
        if ( maxIterations > 0 && isave[29] > maxIterations ) {
          break;
        }
        continue;
      }
      else if( strncmp(task,"CONV",4) == 0 ) {
        //this->setFieldAndCompute();
        compute(in, para, &f, g, RESIDUE);
        cout << task << endl;
        break;
      }
      else if( strncmp(task,"ABNORM",6) == 0 ) {
        //this->setFieldAndCompute();
        compute(in, para, &f, g, RESIDUE);
        cout << "Abnormal line search : " << task << endl;
        for (int i=0; i<Tot_len; i++){
          cout << in->data[i] << endl;
        }
        for (int j=0; j< NP-1; j++) {
          //cout << in->data[Tot_len + 2*j]*180/PI << " " << in->data[Tot_len + 2*j + 1]*180/PI << endl;
        }
        break;
      }

      // If task is neither FG nor NEWX we terminate execution.
      else {
        cout << task << endl;
        break;
      }

      // ---------- the end of the loop -------------
    }

    // // --------- Newton Ralphs minimization :START
    // cout << "\033[1;31mSolving the root\033[0m" << endl;
    // gsl_multimin_function_fdf myfunc;
    // const gsl_multiroot_fsolver_type* Troot;
    // gsl_multiroot_fsolver *sroot;
          

    // myfunc.n = Tot_len+2*NP;
    // myfunc.f = &myf;
    // myfunc.df = &mydf;
    // myfunc.fdf = &myfdf;
    // myfunc.params = (void*) para;
    // //struct rparams p = para;
    // gsl_multiroot_function ff = {&F_root, myfunc.n, para };
  
    // int iter_root = 0;
    // int status_root;
    // Troot = gsl_multiroot_fsolver_hybrids;
    // sroot = gsl_multiroot_fsolver_alloc(Troot, myfunc.n);
    // gsl_multiroot_fsolver_set(sroot, &ff, in);
    // //print_state(iter_root, sroot);
    // do{
    //   iter_root++;
    //   status_root = gsl_multiroot_fsolver_iterate(sroot);
    //   status_root = gsl_multiroot_test_residual(sroot->f, 1e-4);
    //   cout << "\033[1;32m norm = "<< gsl_blas_dnrm2(sroot->f) <<"\033[0m"<< endl;
    //   //if(status_root)
    //   // break;
    // }while(status_root==GSL_CONTINUE && iter_root < 100);
    // cout << "\033[1;34m"<<gsl_strerror(status_root) <<" !!\033[0m"<< endl;
    // // END: NR Minimization
    
    cout.unsetf(ios_base::scientific);


    // ---------- the end of solve() -------------
    projg = dsave[12];
    iterNo = isave[29];

    //this->setFieldAndCompute();
    cntr ++; 
    compute(in, para, &f, g, RESIDUE);
    char* NPstr;
    ostringstream toStringNP, toStringCntr; //Hoffman does not have c++11
    toStringNP << NP;
    toStringCntr << cntr;
    string temp = filename + "_N_" + toStringNP.str() +  "_" + toStringCntr.str();
    printToVTK(in, temp);
    printModes(in, para, temp);
    printParticleCoords(in, para, temp);
    double enParticles =  myfParticle(in, para);
    double enTot = myf(in, para);
    gsl_vector* vSphere = gsl_vector_calloc(Ntot+2*NP);
    // Particle on rigid sphere
    for (int i=1; i<Tot_len; i++)
      vSphere->data[i] = 0;
    for (int j=0; j< NP; j++) {
      vSphere->data[Tot_len + 2*j] = in->data[Tot_len + 2*j]; 
      vSphere->data[Tot_len + 2*j+1] = in->data[Tot_len + 2*j+1 ];
    }
    vSphere->data[0] = (-1+0.9510565163 *para[2])*sqrt(4*PI);
    double enParticlesRigid =  myfParticle(vSphere, para);
    cout << "\033[31m Particle Energy = "<< enParticles << "\033[0m\n";
    cout << "\033[32m Total Energy = "<< enTot  << "\033[0m\n";
    cout << "\033[32m Rigid Energy = "<< enParticlesRigid  << "\033[0m\n";
    ofstream outFileAp("Energy.txt", ofstream::app);
    outFileAp.setf(ios_base::scientific);
    outFileAp << para[0] << " \t " << para[2]<< "\t" << enParticles << " \t " << enTot << " \t " << projg << endl;
    outFileAp.close();
    para[2] += drm;
    //para[0] += dk;
    cout << "\033[35m rm = " << para[2] << "\033[0m" << endl;
    cout << "\033[35m k = " << para[0] << "\033[0m" << endl;
  }
  // if (dk < 0){
  //   dk = -dk;
  //   para[0] += dk;
  //   goto loop;
  // }
  // //Start minimization (------- GSL ----------)
  // const gsl_multimin_fdfminimizer_type *T;
  // gsl_multimin_fdfminimizer *s;
  // gsl_multimin_function_fdf myfunc;
  // const gsl_multiroot_fsolver_type* Troot;
  // gsl_multiroot_fsolver *sroot;
          

  // myfunc.n = Tot_len+2*NP;
  // myfunc.f = &myf;
  // myfunc.df = &mydf;
  // myfunc.fdf = &myfdf;
  // myfunc.params = (void*) para;

  // T = gsl_multimin_fdfminimizer_vector_bfgs;
  // //T = gsl_multimin_fdfminimizer_steepest_descent;
  // //T = gsl_multimin_fdfminimizer_conjugate_fr;
  // s = gsl_multimin_fdfminimizer_alloc(T, myfunc.n);

  // gsl_multimin_fdfminimizer_set(s, &myfunc, in, 1, .1);
  // int itMax =5;
  // cout << setprecision(15);
  // for (int j=0; j<1500; j++) {
  //   para[1] = 1.0;
  //   for (int i=0; i< itMax; i++){
  //     size_t iter = 0;
  //     int status;
  //     double oldnorm = 1;
  //     double norm;
  //     do{
  //       iter++;
  //       status = gsl_multimin_fdfminimizer_iterate(s);
  //       norm = gsl_blas_dnrm2(s->gradient);
  //       if (iter%2 == 0)
  //         cout << "iter=" << iter << ", Gradient norm = " << norm << ", energy = " << (s->f) << endl;
  //       if (fabs(oldnorm-norm)<1e-8){
  //         if (status == GSL_ENOPROG){
  //           cout << "\033[1;35m ENOPRG \033[0m" << endl;
  //         }
  //         cout << "\033[1;31mSaturation\033[0m\n";
  //         //printToVTK(s->x, "nonconv");
  //         //exit(0);
  //         cout << "\033[1;31mSolving the root\033[0m" << endl;
  //         //struct rparams p = para;
  //         gsl_multiroot_function f = {&F_root, myfunc.n, para };
  
  //         int iter_root = 0;
  //         int status_root;
  //         Troot = gsl_multiroot_fsolver_hybrids;
  //         sroot = gsl_multiroot_fsolver_alloc(Troot, myfunc.n);
  //         gsl_multiroot_fsolver_set(sroot, &f, s->x);
  //         //print_state(iter_root, sroot);
  //         do{
  //           iter_root++;
  //           status_root = gsl_multiroot_fsolver_iterate(sroot);
  //           status_root = gsl_multiroot_test_residual(sroot->f, 1e-4);
  //           cout << "\033[1;32m norm = "<< gsl_blas_dnrm2(sroot->f) <<"\033[0m"<< endl;
  //           //if(status_root)
  //           // break;
  //         }while(status_root==GSL_CONTINUE && iter_root < 100);
  //         cout << "\033[1;34m"<<gsl_strerror(status_root) <<" !!\033[0m"<< endl;
  //         s->x = sroot->x;
  //         continue;
  //       }
  //       status = gsl_multimin_test_gradient(s->gradient, 1e-4);
  //       if (status == GSL_SUCCESS){
  //         norm = gsl_blas_dnrm2(s->gradient);
  //         cout << "Minimimum found\n";
  //         cout << "Gradient norm = " << norm << ", energy = " << s->f << endl;
  //       }
  //       oldnorm = norm;
  //     }while(status == GSL_CONTINUE && iter < 1000);
  //     para[1] = para[1]/5.0;
  //     cout << "eps  set to " << para[1] << endl;
  //     myfunc.params = (void*) para;
  //     gsl_multimin_fdfminimizer_set(s, &myfunc, s->x, .01, .1);
  
  //   }
  //   string temp = filename + "_" + to_string(j);
  //   printToVTK(s->x, temp);
  //   para[2] = para[2] + 0.1;
  //   cout << "\033[1;31mrm  set to \033[0m" << para[2] << endl;
  //   myfunc.params = (void*) para;
  //   gsl_multimin_fdfminimizer_set(s, &myfunc, s->x, 1, .1);
  
  // }

}

