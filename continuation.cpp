#include<iostream>
#include<gsl/gsl_linalg.h>
#include"continuation.h"
#include<gsl/gsl_cblas.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf.h>
#include<gsl/gsl_eigen.h>
#include<string.h>
#include<sstream>
#include"kbhit.h"
#include"gausslegendre.h"

using namespace std;
Quadr::Quadr(){
	order=23;
}
void Quadr::lebedevQuad(int (*fun)(double*,double*,gsl_vector*), double* args,gsl_vector* out){
	string lebfile("Lebedev/OurLebedev");
	stringstream ss;
	ss<<order;
	lebfile=lebfile+ss.str()+".dat";
	
	FILE* pfile;
	double t,p,weight;
	double* tp=new double[2];
	gsl_vector* temp=gsl_vector_calloc(out->size);
	gsl_vector_set_zero(out);
	pfile=fopen(lebfile.c_str(),"r");
	//in OurLebedev files the first element is phi .in. [-180,180]
	// second element is theta .in. [0,180] and the third element is the weight
	while(fscanf(pfile,"%lf %lf %lf",&p,&t,&weight)!=EOF){
		if(t==0.0) t=1e-4;
		if(t==180.0) t=180.0-1e-4;
		tp[0]=t;
		tp[1]=p;
		fun(tp,args,temp);
		//out<--out+weight*temp
		gsl_blas_daxpy(weight,temp,out);
	}
	gsl_blas_dscal(4*M_PI,out);
	gsl_vector_free(temp);
	delete[] tp;
	fclose(pfile);
};
void Quadr::gaussLegendreQuad(int (*fun)(double*, double*, gsl_vector*),double* args, gsl_vector* out){
	gsl_vector* temp = gsl_vector_calloc(out->size); //temp =[0,0,...,0]_size
	gsl_vector_set_zero(out);
	
	int n = 1024/2; //THIS NEEDS TO BE MODIFIED
	double tp[2]={0.,0.};
	for (int i=0; i< n; i++){
		tp[0] = (1+ x1024[i])*180.0/2; // in degrees
		fun(tp,args,temp);
		//out <-- out+weight*temp
		gsl_blas_daxpy(w1024[i]*M_PI/2,temp,out);
		
		//Now for PI-theta. The weights from PI-theta and theta are the same
		tp[0]=(1-x1024[i])*180.0/2;
		fun(tp,args,temp);
		//out <-- out+weight*temp
		gsl_blas_daxpy(w1024[i]*M_PI/2,temp,out);
		
	}
	/*
	int n = 128/2; //THIS NEEDS TO BE MODIFIED
	double tp[2]={0.,0.};
	for (int i=0; i< n; i++){
		tp[0] = (1+ x128[i])*180.0/2; // in degrees
		fun(tp,args,temp);
		//out <-- out+weight*temp
		gsl_blas_daxpy(w128[i]*M_PI/2,temp,out);
		
		//Now for PI-theta. The weights from PI-theta and theta are the same
		tp[0]=(1-x128[i])*180.0/2;
		fun(tp,args,temp);
		//out <-- out+weight*temp
		gsl_blas_daxpy(w128[i]*M_PI/2,temp,out);
		
	}
	*/
	gsl_blas_dscal(2*M_PI,out); //this step accounts for integration in phi. Very simple for axisymmetric version
	gsl_vector_free(temp);
};
void Quadr::adaptiveQuad(int (*fun)(double*, double*, gsl_vector*),double* args, double a, double b, double tau, gsl_vector* out){
	gsl_vector* temp = gsl_vector_calloc(out->size); //temp =[0,0,...,0]_size
	gsl_vector_set_zero(out);
	// out1 - vector to store lower order quad
	// out - vector stores higher order quad
	gsl_vector* out1 = gsl_vector_calloc(out->size); //temp =[0,0,...,0]_size
	// Lower order quad 
	int n = 12/2; 
	double tp[2]={0.,0.};
	gsl_vector_set_zero(temp);
	for (int i=0; i< n; i++){
		tp[0] = ((a+b)/2.0 + (b-a)/2.0*x12[i]); // in degrees
		fun(tp,args,temp);
		//out <-- out+weight*temp
		gsl_blas_daxpy(w12[i]*(b-a)*M_PI/(2*180.0),temp,out1);
		
		//Now for PI-theta. The weights from PI-theta and theta are the same
		tp[0]=((a+b)/2-x12[i]*(b-a)/2);
		fun(tp,args,temp);
		//out <-- out+weight*temp
		gsl_blas_daxpy(w12[i]*(b-a)*M_PI/(2*180),temp,out1);
		
	}
	gsl_blas_dscal(2*M_PI,out1); //this step accounts for integration in phi. Very simple for axisymmetric version
	// Higher order quad 
	n = 16/2; 
	gsl_vector_set_zero(temp);
	for (int i=0; i< n; i++){
		tp[0] = ((a+b)/2+ x16[i]*(b-a)/2); // in degrees
		fun(tp,args,temp);
		//out <-- out+weight*temp
		gsl_blas_daxpy(w16[i]*(b-a)*M_PI/(2*180),temp,out);
		
		//Now for PI-theta. The weights from PI-theta and theta are the same
		tp[0]=((a+b)/2-x16[i]*(b-a)/2);
		fun(tp,args,temp);
		//out <-- out+weight*temp
		gsl_blas_daxpy(w16[i]*(b-a)*M_PI/(2*180),temp,out);
		
	}
	gsl_blas_dscal(2*M_PI,out); //this step accounts for integration in phi. Very simple for axisymmetric version
	gsl_vector_memcpy(temp,out1);
	gsl_blas_daxpy(-1,out,temp);
	if(gsl_blas_dnrm2(temp) > tau && (b-a)> 1e-3 ){
		gsl_vector_free(temp);
		adaptiveQuad(fun, args, a, (a+b)/2, tau/2, out1);
		adaptiveQuad(fun, args, (a+b)/2, b, tau/2, out);
		//out<--out1+out
		gsl_blas_daxpy(1.0,out1,out);
	}
	gsl_vector_free(out1);
	return;

};
//* SOLVER*/
double LinearSolver::luSolve(gsl_matrix* A, gsl_vector* b,gsl_vector* x){
	
	//I don't want to change the matrix A when I leave this fun.
	gsl_matrix* temp=gsl_matrix_alloc(A->size1,A->size2);
	gsl_matrix_memcpy(temp,A);
		
	gsl_permutation* p=gsl_permutation_alloc(b->size);
	int s;
	gsl_linalg_LU_decomp(temp,p,&s); //After this line A is an "LU-ed" matrix
	gsl_linalg_LU_solve(temp,p,b,x);
	
	/* ****************SV SOLVE*****************************/
	/*
	gsl_matrix* V=gsl_matrix_alloc(A->size2,A->size2);
	gsl_vector* work=gsl_vector_alloc(A->size2);
	gsl_vector* S=gsl_vector_alloc(A->size2);
	gsl_linalg_SV_decomp(A,V,S,work);
	gsl_linalg_SV_solve(A,V,S,b,x);
	*/
	//double determinant=gsl_linalg_LU_det(A,s);
	//gsl_matrix_memcpy(A,temp);
	gsl_matrix_free(temp);
	return 0.;
}
double LinearSolver::linSolve(gsl_matrix* A,gsl_vector* b,gsl_vector* x){
	if(strcmp(solver,"LU")==0) return luSolve(A,b,x);
	else {
		cout<<"Unkown Solver: "<<solver<<endl;
		return NULL;
	}
}
NonlinearSolver::NonlinearSolver(){
	predictorType="generic";
	correctorType="quasi-newton";
	solver="LU";
	maxIterations=50;
	noOfIterations=0;
	tolerance=1e-7;
	eps=1e-5;
	indVar=NULL;
	
}
int NonlinearSolver::setNoOfUnknownVar(int n){
	noOfUnknownVar=n;
	jacobian=gsl_matrix_calloc(n,n);
	gsl_matrix_set_identity(jacobian);
	tangent=gsl_vector_calloc(n);
	return 0;
}
int NonlinearSolver::setNoOfInArgs(int n){
	noOfInArgs=n;
	solution=new double[n];
	return 0;
}
NonlinearSolver::~NonlinearSolver(){
	//delete[] solution;
	//gsl_vector_free(tangent);
	//gsl_matrix_free(jacobian);
	}
/*
void NonlinearSolver::setNoOfIndVar(int n){
	noOfIndVar=2;
	indVar[0]=0.;
	indVar[1]=0.;
}
*/
int NonlinearSolver::F(double* args, gsl_vector* out){
	//cout<<"In Nonlinear::F "<<endl;
	//copying args into a new variable in
	double* in=new double[noOfInArgs];
	memcpy(in,args,sizeof(double)*noOfInArgs);
	ptr_F(indVar,in,out);
	delete[] in;
	return 0;
}
int NonlinearSolver::Jac(double* args, gsl_matrix* out){
	//copying args into a new variable in
	//cout<<"In Nonlinear:: Jac "<<endl;
	double* in=new double[noOfInArgs];
	memcpy(in,args,sizeof(double)*noOfInArgs);
	
	int row=(int)out->size1;
	int col=(int)out->size2;
	gsl_vector* temp_v=gsl_vector_calloc(row);
	gsl_vector* f0=gsl_vector_calloc(row);
	double temp_x;
	//problem->myFun(problem->indVar,args,f0);
	//cout<<"In Nonlinear::Jac\n";
	F(args,f0);
	for(int i=0;i<col;i++){
		temp_x=in[i];
		in[i]=in[i]+eps;
		F(in,temp_v);
		//temp_v<---temp_v-f0
		gsl_blas_daxpy(-1,f0,temp_v);
		//temp_v<--(temp_v)/eps;
		gsl_blas_dscal(1/eps,temp_v);
		//set the ith column as (Fpert-F0)/eps
		gsl_matrix_set_col(out,i,temp_v);
		//reset the perturbation
		in[i]=temp_x;
		
	}
	gsl_vector_free(temp_v);
	gsl_vector_free(f0);
	//gsl_matrix_memcpy(out,out_temp);
	//gsl_matrix_free(out_temp);
	delete[] in;
	return 0;
}
void NonlinearSolver::newtonMethod(double* args){
	//cout<<"Entering Newton Method\n";
	double* x_np1=new double[noOfInArgs];
	gsl_vector* v_np1=gsl_vector_calloc(tangent->size);
	gsl_vector* F_np1=gsl_vector_calloc(tangent->size);
	gsl_vector* Dv=gsl_vector_calloc(tangent->size);
	
	//gsl_matrix* jac=gsl_matrix_calloc(jacobian->size1,jacobian->size2);
	memcpy(x_np1,args,sizeof(double)*noOfInArgs);
	array2vector(args,v_np1);
	F(args,F_np1);
	while(gsl_blas_dnrm2(F_np1)>tolerance){
		Jac(x_np1,jacobian);
		linSolve(jacobian,F_np1,Dv);
		//v_np1=v_np1-Dv
		//the -ve sign is because "Dv"= - inv(jac)*F_np1
		gsl_blas_daxpy(-1.0,Dv,v_np1);
		vector2array(v_np1,x_np1);
		F(x_np1,F_np1);
		cout<<"Iterating Norm= "<<gsl_blas_dnrm2(F_np1)<<endl;
	}
	memcpy(solution,x_np1,sizeof(double)*noOfInArgs);
	
	gsl_vector_free(v_np1);
	gsl_vector_free(F_np1);
	gsl_vector_free(Dv);
	delete[] x_np1;
	//gsl_matrix_free(jac);
}

/*
 * Implementing the quasi newton scheme using the SR1 method. 
 * See wikipedia entry on quasi newton scheme for the method.
 * Actually, see the wiki entry on Broyden. There seems to be 
 * a boo-boo in Quasi-Newton entry of Broyden. Kellers book
 * on Numerical Continution is also a good reference.
*/

void NonlinearSolver::quasiNewtonMethod(double* Args){
	//First, I don't want to change Args:
	double* args=new double[noOfInArgs];
	memcpy(args,Args,sizeof(double)*noOfInArgs);
	noOfIterations=0;
	gsl_vector* x_k=gsl_vector_alloc(jacobian->size1);
	gsl_vector* F_k=gsl_vector_alloc(jacobian->size1);
	gsl_vector* y_k=gsl_vector_alloc(jacobian->size1);
	gsl_vector* Dx_k=gsl_vector_alloc(jacobian->size1);
	gsl_matrix* H_k=gsl_matrix_alloc(jacobian->size1,jacobian->size2);
	gsl_matrix* H_temp=gsl_matrix_alloc(jacobian->size1,jacobian->size2);
	gsl_permutation* p=gsl_permutation_alloc(jacobian->size1);
	int s;
	
	double denominator;
	double numerator;
	gsl_vector* temp_x=gsl_vector_alloc(jacobian->size1);
	gsl_vector* temp_F=gsl_vector_alloc(jacobian->size1);

	Jac(args,H_temp); 
	gsl_matrix_memcpy(jacobian,H_temp);
	//gsl_matrix_fprintf(stdout, H_temp, "%f");
	//H_k from now is actually the inverse of the jacobian
	printf("---here---2\n");
	gsl_linalg_LU_decomp(H_temp,p,&s);
	gsl_linalg_LU_invert(H_temp,p,H_k);
	printf("---here---2^\n");
	//gsl_matrix_set_identity(H_k);
	
	array2vector(args,x_k); F(args,F_k);
	//gsl_vector_fprintf(stdout,F_k,"%f");cout<<endl;
	while(gsl_blas_dnrm2(F_k)>tolerance){
		//Dx_k=-alpha_k*H_k*F_k, here alpha_k=1
		gsl_blas_dgemv(CblasNoTrans,-1.0,H_k,F_k,0.,Dx_k);
		//store x_k in temp_k
		gsl_vector_memcpy(temp_x,x_k);
		//x_k<--x_k+Dx_k
		gsl_vector_add(x_k,Dx_k); //x_k+1
		/***  y_k=F_k+1-F_k ***/
		//temp_k<---F_k
		gsl_vector_memcpy(temp_F,F_k);
		//Compute F_k+1
		vector2array(x_k,args); F(args,F_k); //F_k+1, args="x_k+1"
		//F_k<--F_k-F_k+1
		gsl_blas_daxpy(-1.0,F_k,temp_F);
		//F_k<-- -F_k
		gsl_blas_dscal(-1.0,temp_F);
		//y_k<-- F_k
		gsl_vector_memcpy(y_k,temp_F);
		//Updating H i.e, computing H_k+1
		//temp_x<--H_k*y_k 
		gsl_blas_dgemv(CblasNoTrans,1.0,H_k,y_k,0.,temp_x);
		//denominator<--Dx_k^T * H_k *y_k
		gsl_blas_ddot(Dx_k,temp_x,&denominator);
		//First saving Dx_k in temp_x
		gsl_vector_memcpy(temp_x,Dx_k);
		//Dx_k<--Dx_k-H_k*y_k
		gsl_blas_dgemv(CblasNoTrans,-1.0,H_k,y_k,1.0,Dx_k);
		//numerator<--(Dx_k-H_k*y_k)*Dx_k^T
		gsl_blas_ddot(Dx_k,temp_x,&numerator);
		//H_k<---(1+numerator/denominator)*H_k
		gsl_matrix_scale(H_k,1+numerator/denominator);
		cout<<"Iterating Norm= "<<gsl_blas_dnrm2(F_k)<<endl;
	}
	memcpy(solution,args,sizeof(double)*noOfInArgs);
	gsl_vector_free(x_k);
	gsl_vector_free(F_k);
	gsl_vector_free(y_k);
	gsl_vector_free(Dx_k);
	gsl_vector_free(temp_x);
	gsl_vector_free(temp_F);
	gsl_matrix_free(H_k);
	gsl_matrix_free(H_temp);
	gsl_permutation_free(p);
	delete[] args;
}
void NonlinearSolver::correct(double* args){
	if(strcmp(correctorType,"quasi-newton")==0)
		quasiNewtonMethod(args);
	else if(strcmp(correctorType,"newton")==0)
		newtonMethod(args);
	else
		cout<<"Unknown Corrector!"<<endl;
}
void NonlinearSolver::predict(double* args){
	if(strcmp(predictorType,"generic")==0){
		
		//input: args, outputs: jacobian
		Jac(args, jacobian); 
		//gsl_matrix_fprintf(stdout, jacobian, "%f");
		//jacobian->size1=no. of rows
		gsl_vector* vecPredict=gsl_vector_calloc(jacobian->size1);
		
		//vecPredict=[0,0,0,...,1]
		
		gsl_vector_set(vecPredict,(int)jacobian->size1-1, 1.0);
		memcpy(solution,args,sizeof(double)*noOfInArgs);
		//gsl_matrix_fprintf(stdout, jacobian, "%le");
		printf("---here---\n");
		linSolve(jacobian,vecPredict,tangent);
		printf("---here---^\n");
		/* If the tangent points in the opposite direction of increasign arc
		length, then reverse the direction
		*/
		/*
		if(gsl_vector_get(tangent,tangent->size-1)<0)
			gsl_blas_dscal(-1,tangent);
		*/
		//Normalize the tangent. Since the pseudo-arc equation in F
		//is tangent_(n-1).tangent_n = 1, tangent_n will not
		// have a norm =1.  
		gsl_blas_dscal(1/gsl_blas_dnrm2(tangent),tangent);
		gsl_vector_free(vecPredict);
	}
	else{
		cout<<"Predictor Type : "<<predictorType<<" not available"<<endl;
	}
}
void NonlinearSolver::array2vector(double* args,gsl_vector* out){
	//cout<<"In Nonlinear::array2vector\n";
	if(out->size!=noOfUnknownVar) 
		cout<<"No. Of variables not equal to size of vector!"<<endl;
	for(int i=0;i<noOfUnknownVar;i++)
		gsl_vector_set(out,i,args[i]);
	//gsl_vector_set(out,noOfUnknownVar,args[posOfPara]);
}
void NonlinearSolver::vector2array(gsl_vector* vecIn, double* out){
	if(vecIn->size!=noOfUnknownVar)
		cout<<"Size of vector not equal to no. Of variables!"<<endl;
	for(int i=0;i<noOfUnknownVar;i++)
		out[i]=gsl_vector_get(vecIn,i);
	//out[posOfPara]=gsl_vector_get(vecIn,noOfUnknownVar);
}


/* CONTINUER*/
Continuer::Continuer(){
	performQuad=true;
	stepSize=0.01;
	minStepSize=1e-4;
	maxStepSize=10.;
	saveBPData=true;
	BPFileName="BP_Filename";
	DataFileName="Data_Filename.txt";
	continuerType="pseudo-arc";
	LP=0;
	BP=0;
	
}
void Continuer::array2vector(double* args,gsl_vector* out){
	//cout<<"In Continuer::array2vector\n";
	if(out->size!=noOfUnknownVar+1) 
		cout<<"No. Of variables not equal to size of vector!"<<endl;
	for(int i=0;i<noOfUnknownVar;i++)
		gsl_vector_set(out,i,args[i]);
	gsl_vector_set(out,noOfUnknownVar,args[posOfPara]);
}
void Continuer::vector2array(gsl_vector* vecIn, double* out){
	if(vecIn->size!=noOfUnknownVar+1)
		cout<<"Size of vector not equal to no. Of variables!"<<endl;
	for(int i=0;i<noOfUnknownVar;i++)
		out[i]=gsl_vector_get(vecIn,i);
	out[posOfPara]=gsl_vector_get(vecIn,noOfUnknownVar);
}
int Continuer::F(double* args,gsl_vector* out){
		//cout<< "In Continuer::F\n";
		//if(continuerType=="pseudo-arc" && performQuad==true){
		gsl_vector* v1=gsl_vector_calloc(noOfUnknownVar+1);
		gsl_vector* v2=gsl_vector_calloc(noOfUnknownVar+1);
		array2vector(args,v1);
		array2vector(solution,v2);
		//v1<--v1-v2
		gsl_blas_daxpy(-1.,v2,v1);
		//Pseudo-arc equation
		double v1DotTgt;
		gsl_blas_ddot(v1,tangent,&v1DotTgt);
		double N=v1DotTgt-stepSize;
		//temporarily set last equation(value) to zero so that
		//while doin the quadrature in the next step we are not
		//accumulating anything
		out->data[noOfUnknownVar]=0.;
		// NOW I NEED TO PERFROM QUADRATURE
		if(performQuad==true){
			//Quad.lebedevQuad(ptr_F,args,out);
			Quad.gaussLegendreQuad(ptr_F,args,out);
			//Quad.adaptiveQuad(ptr_F, args, 0.0, 180.0, 1e-6, out);

		}
		else{
			ptr_F(NULL,args,out);
			//gsl_vector_fprintf(stdout, out, "%lf");
			
		}
		//set the last equation(value) to be the pseudo-arc equation
		out->data[noOfUnknownVar]=N;
		gsl_vector_free(v1);
		gsl_vector_free(v2);
	
}

int Continuer::Jac(double* args,gsl_matrix* out){

	//cout<<"In Continuer::Jac\n"<<noOfUnknownVar<<endl;
	int row=(int)out->size1;
	int col=(int)out->size2;
	gsl_vector* temp_v=gsl_vector_calloc(row);
	gsl_vector* f0=gsl_vector_calloc(row);
	double temp_x;
	F(args,f0);
	//cout << "here2\n";

	//Part of the Jacobian by varying unknowns
	for(int i=0;i<noOfUnknownVar;i++){
		temp_x=args[i];
		args[i]=args[i]+eps;
		F(args,temp_v);
		//temp_v<---temp_v-f0
		gsl_blas_daxpy(-1,f0,temp_v);
		//temp_v<--(temp_v)/eps;
		gsl_blas_dscal(1/eps,temp_v);
		//set the ith column as (Fpert-F0)/eps
		gsl_matrix_set_col(out,i,temp_v);
		//reset the perturbation
		args[i]=temp_x;
	}
	//Part of the jacobian by varying parameter
	temp_x=args[posOfPara];
	args[posOfPara]=args[posOfPara]+eps;
	F(args,temp_v);
	//temp_v<---temp_v-f0
	gsl_blas_daxpy(-1,f0,temp_v);
	//temp_v<--(temp_v)/eps;
	gsl_blas_dscal(1/eps,temp_v);
	//set the next, i.e, "noOfUnknownVar"th position
	gsl_matrix_set_col(out,noOfUnknownVar,temp_v);
	//reset the perturbation
	args[posOfPara]=temp_x;
	
	//The last row of th Jacobian contains the tangent vector
	gsl_matrix_set_row(out,out->size1-1,tangent);
	gsl_vector_free(temp_v);
	gsl_vector_free(f0);
	return 0;
	
	//return NonlinearSolver::Jac(args,out);
}
int Continuer::setNoOfUnknownVar(int n){
	//+1 to account for the continuation parameter
	double ret=NonlinearSolver::setNoOfUnknownVar(n+1);
	//resets noOfUnknown variables fron n+1 to n
	noOfUnknownVar=n;
	return ret;
}
int Continuer::setNoOfInArgs(int n){
	return NonlinearSolver::setNoOfInArgs(n);
}
char Continuer::checkSingularPt(double* solution_toSave,gsl_vector* tgt_BP){
	gsl_matrix* jac_b=gsl_matrix_calloc(jacobian->size1,jacobian->size2);
	gsl_matrix_memcpy(jac_b,jacobian);
	gsl_matrix_view temp=gsl_matrix_submatrix(jac_b,0,0,jac_b->size1-1,jac_b->size2-1);
	gsl_matrix* jac_l=gsl_matrix_calloc(jac_b->size1-1,jac_b->size2-1);
	
	gsl_matrix_memcpy(jac_l,&temp.matrix);
	
	gsl_permutation* p_b=gsl_permutation_alloc(jac_b->size1);
	int s;
	gsl_linalg_LU_decomp(jac_b,p_b,&s);
	double BP_next=gsl_linalg_LU_det(jac_b,s);
	
	gsl_permutation* p_l=gsl_permutation_alloc(jac_l->size1);
	gsl_linalg_LU_decomp(jac_l,p_l,&s);
	double LP_next=gsl_linalg_LU_det(jac_l,s);
	printf("\033[0;32m %le %le\n \033[m",LP_next,BP_next);
	if(LP==0 && BP==0){
		//Which is the case for the very first step
		LP=LP_next;
		BP=BP_next;
		
		gsl_matrix_free(jac_l);
		gsl_matrix_free(jac_b);
		gsl_permutation_free(p_b);
		gsl_permutation_free(p_l);
		
		return 'f';
	}
	else if(LP*LP_next<0 && BP*BP_next<0){
		//Which is the case when we move past a BP
		LP=LP_next;
		BP=BP_next;
		
		//Matrix and Vector definitions needed for SVD
		gsl_matrix* V=gsl_matrix_alloc(jacobian->size2,jacobian->size2);
		gsl_vector* S=gsl_vector_alloc(jacobian->size2);
		gsl_vector* work=gsl_vector_alloc(jacobian->size2);
		gsl_matrix_memcpy(jac_b,jacobian);
		//jac_b=USV^t; U is stored in jac_b
		gsl_linalg_SV_decomp(jac_b,V,S,work);
		/*Places the last row of the matrix V^T in tgt_BP (gives 
		the bifurcation direction*/
		
		gsl_matrix_transpose(V);
		gsl_matrix_get_row(tgt_BP,V,V->size1-1);
		
		//Now, we find the solution
		gsl_vector* sol_temp=gsl_vector_alloc(tgt_BP->size);
		array2vector(solution,sol_temp);
		//sol<--sol+stepSize*tgt_BP
		//gsl_blas_daxpy(stepSize,tgt_BP,sol_temp);
		vector2array(sol_temp,solution_toSave);
		
			//save_BP(BPFileName,solution_toSave,tgt_BP,noOfInArgs);
			
		gsl_vector_free(sol_temp);
		gsl_vector_free(work);
		gsl_vector_free(S);
		gsl_matrix_free(V);
		gsl_matrix_free(jac_l);
		gsl_matrix_free(jac_b);
		gsl_permutation_free(p_b);
		gsl_permutation_free(p_l);
		
		return 'B';
	}
	else if(LP*LP_next<0 && BP*BP_next>0){
		//Which happens at a LP
		LP=LP_next;
		BP=BP_next;
		gsl_matrix_free(jac_l);
		gsl_matrix_free(jac_b);
		gsl_permutation_free(p_b);
		gsl_permutation_free(p_l);
		return 'L';
	}
	else{
		LP=LP_next;
		BP=BP_next;
		gsl_matrix_free(jac_b);
		gsl_permutation_free(p_b);
		gsl_permutation_free(p_l);
		return '-';
	}

}
void Continuer::setInitialSolution(double* args){
	//arraycpy(args,solution,noOfInArgs);
	
	memcpy(solution,args,sizeof(double)*noOfInArgs);
	//gsl_vector_set_zero(tangent);
	Jac(args,jacobian);
	/*The jacobian will, at this stage, always have 0's in the last row because
	a trivial Augmented Pseudo arc equation. We will remove the last row to 
	find the tangential direction to move
	*/
	//gsl_matrix_transpose(jacobian);
	
	gsl_matrix_view jac_clip_view=gsl_matrix_submatrix(jacobian,0,0,jacobian->size1,jacobian->size2);
	gsl_matrix* jac_clip=gsl_matrix_alloc(jacobian->size1,jacobian->size2);
	gsl_matrix_memcpy(jac_clip,&jac_clip_view.matrix);
	/*Now, performing an SVD to find the tangential direction*/
	/*------------ Matrices, vectors needed for svd------------*/
	
	gsl_matrix* U=gsl_matrix_alloc(jac_clip->size2,jac_clip->size2);
	gsl_vector* S=gsl_vector_alloc(jac_clip->size2);
	gsl_vector* work=gsl_vector_alloc(jac_clip->size2);
	gsl_linalg_SV_decomp(jac_clip,U,S,work);
	//gsl_matrix_transpose(&jac_clip.matrix);
	gsl_matrix_get_row(tangent,jac_clip,jac_clip->size2-1);
	gsl_matrix_free(U);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_matrix_free(jac_clip);
	//gsl_matrix_transpose(jacobian);
	
}
int Continuer::Continue(){
	gsl_vector* tgt_BP=gsl_vector_calloc(tangent->size);
	double* solution_toSave=new double[noOfInArgs];
	memcpy(solution_toSave,solution,sizeof(double)*noOfInArgs);
	int BPnum=0;
	int key;
	do{
		if(kbhit())
		{
			key = fgetc(stdin);
			if (key == 'i'){
				stepSize += .01;
				cout<< "\n \033[1;35m Step Size increased to "<< stepSize<<"\033[m\n\n";
				}
			else if(key=='d'){
			 	stepSize -= .01;
				cout<< "\n \033[1;36m Step Size decreased to "<< stepSize<<"\033[m\n\n";
				}
			else if(key=='r'){
				stepSize=-stepSize;
				cout<< "\n \033[1;36m Step Size reversed \033[m\n\n";
			}
			else if(key =='n'){
				stepSize = 2*stepSize;
				cout<< "\n \033[1;36m Step Size doubled to"<<stepSize<<" \033[m\n\n";
			}
			else if(key =='h'){
				stepSize = stepSize/2.0;
				cout<< "\n \033[1;33m Step Size halved to"<<stepSize<<" \033[m\n\n";
			}
			fflush(stdin);
		} 
		else{
			predict(solution);
			correct(solution);
			gsl_vector * eval = gsl_vector_alloc(jacobian->size1);
			gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (jacobian->size1);
			gsl_matrix * evec = gsl_matrix_alloc(jacobian->size1, jacobian->size2);
			gsl_matrix* Jac =gsl_matrix_alloc(jacobian->size1, jacobian->size2);
			
			gsl_matrix_memcpy(Jac, jacobian);
			gsl_eigen_symmv(Jac, eval, evec, w);
			gsl_eigen_gensymmv_sort (eval, evec,GSL_EIGEN_SORT_ABS_DESC);
			//printf(" -------------------------Lowest Eigen Value = %lf \n", eval->data[eval->size-1]);
			//printf(" -------------------------Highest Eigen Value = %lf \n", eval->data[0]);
			printf("Condition Number: %le\n\n", fabs(eval->data[0]/eval->data[eval->size-1]));
			//gsl_vector_fprintf(stdout, eval, "%f");
			gsl_eigen_symmv_free(w);
			gsl_vector_free(eval);
			gsl_matrix_free(evec);
			gsl_matrix_free(Jac);
			//------------------------------------------------------------
			char sing_pt=checkSingularPt(solution_toSave,tgt_BP);
			cout<<solution[posOfPara]<<" "<<sing_pt<<endl;
			
			if(sing_pt=='L'){ 
				cout<<"\033[1;34m Limit Pt\033[m Detected"<<endl;
			}
			else if(sing_pt=='B'){
				BPnum=BPnum+1;
				cout<<"\033[1;33m Branch Pt\033[m Detected"<<endl;
				save_BP(BPFileName,solution_toSave,tgt_BP,noOfInArgs,BPnum);
			}
		}
		
		
		savearray(DataFileName,solution,noOfInArgs);
	}while(key!='x');
};
int Continuer::load_BP(const char* filename,int BPnum){
	string fname(filename);
	stringstream out;
	out<<BPnum;
	fname=fname+"_"+out.str()+".txt";
	string fnameT(filename);
	fnameT=fnameT+"_"+out.str()+"_Tgt"+".txt";
	loadarray(fname.c_str(),solution,noOfInArgs);
	
	cout<<endl;
	FILE* pfile;
	pfile=fopen(fnameT.c_str(),"r");
	gsl_vector_fscanf(pfile,tangent);
	
	fclose(pfile);
	gsl_vector* sol_temp=gsl_vector_alloc(tangent->size);
	//solution=solution+stepSize*tangent
	array2vector(solution,sol_temp);
	gsl_blas_daxpy(stepSize,tangent,sol_temp);
	vector2array(sol_temp,solution);
	gsl_vector_free(sol_temp);
	return 0;
}
Continuer::~Continuer(){
}
int save_BP(const char* filename,double* array,gsl_vector* tgt,int dim,int BPnum){
	string fname(filename);
	stringstream out;
	out<<BPnum;
	fname=fname+"_"+out.str()+".txt";
	savearray(fname.c_str(),array,dim);
	FILE* pfile;

	string fnameT(filename);
	fnameT=fnameT+"_"+out.str()+"_Tgt"+".txt";
	pfile=fopen(fnameT.c_str(),"a");
	gsl_vector_fprintf(pfile,tgt,"%le");
	fclose(pfile);
}
int savearray(const char* filename,double* array,int dim){
	FILE* pfile;
	pfile=fopen(filename,"a");
	for(int i=0;i<dim;i++)
		fprintf(pfile,"%1.20f ", array[i]);
	fprintf(pfile,"\n");
	fclose(pfile);
	return 0;
}
int loadarray(const char* filename,double* out,int dim){
	//double* array=new double[dim];
	FILE* pfile;
	pfile=fopen(filename,"r");
	if (pfile!=NULL){
		for(int i=0;i<dim;i++){
			double temp;
			fscanf(pfile,"%le ", &temp);
			out[i]=temp;
		}
	fclose(pfile);
	//memcpy(out,array,dim);
	//printarray(out,dim,"%f|~~");
	
	}
	else{
		cout<<"File: "<<filename<<" does not exist"<<endl;
	}
	//delete[] array;
	return 0;
}
void pprintMat(gsl_matrix* M){
	cout<<endl;
	double val;
	for(int i=0;i<M->size1;i++){
		for(int j=0;j<M->size2;j++){
			val=gsl_matrix_get(M,i,j);
			if(fabs(val)<1e-10){ 
				printf("%1.5f ",val);
			}
			else{
			  //printf("\033[1;33m%1.2f\033[m ",val);
			  printf("%1.5f ",val);
			}
		}
		cout<<endl;
	}
}
void printarray(double* A,int dim,char* format){
	for(int i=0;i<dim;i++)
		printf(format,A[i]);
	cout<<endl;
}

