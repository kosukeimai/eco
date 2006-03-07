/******************************************************************
  This file is a part of eco: R Package for Estimating Fitting
  Bayesian Models of Ecological Inference for 2X2 tables
  by Ying Lu and Kosuke Imai
  Copyright: GPL version 2 or later.
*******************************************************************/

#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <Rinterface.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Memory.h>
#include <R_ext/Random.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "sample.h"
#include "bayes.h"
#include "macros.h"
#include "fintegrate.h"


void readData(Param* params, int n_dim, double* pdX, double* sur_W, double* x1_W1, double* x0_W2,
                int n_samp, int s_samp, int x1_samp, int x0_samp);
void ecoSEM(double* optTheta, double* pdTheta, Param* params, double Rmat_old[7][7], double Rmat[7][7]);
void ecoEStep(Param* params, double* suff);
void ecoMStep(double* Suff, double* pdTheta, Param* params);
void ecoMStepNCAR(double* Suff, double* pdTheta, Param* params);
void initNCAR(Param* params,double mu0, double mu1);
void setHistory(double* t_pdTheta, double loglik, int iter,int len,double history_full[][8]);
int closeEnough(double* pdTheta, double* pdTheta_old, int len, double maxerr);
int semDoneCheck(setParam* setP);
void gridEStep(Param* params, int n_samp, int s_samp, int x1_samp, int x0_samp, double* suff, int verbose, double minW1, double maxW1);
void transformTheta(double* pdTheta, double* t_pdTheta, int len);
void untransformTheta(double* t_pdTheta,double* pdTheta, int len);

void cEMeco(
	    /*data input */
	    double *pdX,         /* data (X, Y) */
	    double *pdTheta_in,  /* Theta^ t
				    mu1, mu2, var1, var2, rho */
	    int *pin_samp,       /* sample size */

	    /* loop vairables */
	    int *iteration_max,          /* number of maximum interations */
	    double *convergence,          /* abs value limit before stopping */

	    /*incorporating survey data */
	    int *survey,         /*1 if survey data available(W_1, W_2)
				   0 not*/
	    int *sur_samp,       /*sample size of survey data*/
	    double *sur_W,       /*set of known W_1, W_2 */

	    /*incorporating homeogenous areas */
	    int *x1,       /* 1 if X=1 type areas available W_1 known,
			      W_2 unknown */
	    int *sampx1,   /* number X=1 type areas */
	    double *x1_W1, /* values of W_1 for X1 type areas */

	    int *x0,       /* 1 if X=0 type areas available W_2 known,
			      W_1 unknown */
	    int *sampx0,   /* number X=0 type areas */
	    double *x0_W2, /* values of W_2 for X0 type areas */

	    /* bounds of W1 */
	    double *minW1, double *maxW1,

	    /* options */
	    int *flag,    /*0th (rightmost) bit: 1 = NCAR, 0=normal; 1st bit: 1 = fixed rho, 0 = not fixed rho*/
	    int *verbosiosity,    /*How much to print out, 0=silent, 1=cycle, 2=data*/
      int *calcLoglik,    /*How much to print out, 0=silent, 1=cycle, 2=data*/
	    double *optTheta,  /*optimal theta obtained from previous EM result; if set, then we're doing SEM*/

	    /* storage */
	    double *pdTheta,  /*EM result for Theta^(t+1) */
	    double *Suff,      /*out put suffucient statistics (E(W_1i|Y_i),
				E(E_1i*W_1i|Y_i..) when  conveges */
      double *inSample, /* In Sample info */
      double *DMmatrix,  /* DM matrix for SEM*/
      int *itersUsed, /* number of iterations used */
      double *history /* history of param (transformed) as well as logliklihood*/
	    ){

  int n_samp  = *pin_samp;    /* sample size */
  int s_samp  = *survey ? *sur_samp : 0;     /* sample size of survey data */
  int x1_samp = *x1 ? *sampx1 : 0;       /* sample size for X=1 */
  int x0_samp = *x0 ? *sampx0 : 0;       /* sample size for X=0 */
  int t_samp=n_samp+s_samp+x1_samp+x0_samp;  /* total sample size*/
  int n_dim=2;        /* dimensions */

  setParam setP;
  //set options
  setP.ncar=bit(*flag,0);
  setP.fixedRho=bit(*flag,1);
  setP.sem=bit(*flag,2) & (optTheta[2]>0);
  if (setP.verbose>=1) Rprintf("OPTIONS (flag: %d)   Ncar: %s; Fixed Rho: %s; SEM: %s\n",*flag,setP.ncar==1 ? "Yes" : "No",
   setP.fixedRho==1 ? "Yes" : "No",setP.sem==1 ? "Second run" : (bit(*flag,2)==1 ? "First run" : "No"));
  setP.verbose=*verbosiosity;
  setP.calcLoglik=*calcLoglik;
  setP.convergence=*convergence;
  setP.t_samp=t_samp; setP.n_samp=n_samp; setP.s_samp=s_samp; setP.x1_samp=x1_samp; setP.x0_samp=x0_samp;
  int param_len=(setP.ncar ? 7 : 5) - setP.fixedRho;

  /* model parameters */
  //double **Sigma=doubleMatrix(n_dim,n_dim);/* inverse covariance matrix*/
  //double **InvSigma=doubleMatrix(n_dim,n_dim);/* inverse covariance matrix*/

  double *pdTheta_old=doubleArray(param_len);
  double *t_pdTheta=doubleArray(param_len); //transformed theta
  double *t_pdTheta_old=doubleArray(param_len);
  double Rmat_old[7][7];
  double Rmat[7][7];
  double history_full[*iteration_max+1][8];

  /* misc variables */
  int i, j,main_loop, start;   /* used for various loops */

  /* get random seed */
  GetRNGstate();

  //assign param
  Param* params=(Param*) R_alloc(t_samp,sizeof(Param));

  for(i=0;i<t_samp;i++) params[i].setP=&setP;
  readData(params, n_dim, pdX, sur_W, x1_W1, x0_W2, n_samp, s_samp, x1_samp, x0_samp);



/***Begin main loop ***/
main_loop=1;start=1;
while (main_loop<=*iteration_max && (start==1 ||
        (setP.sem==0 && !closeEnough(t_pdTheta,t_pdTheta_old,param_len,*convergence)) ||
        (setP.sem==1 && !semDoneCheck((setParam*)&setP)))) {
//while (main_loop<=*iteration_max && (start==1 || !closeEnough(transformTheta(pdTheta),transformTheta(pdTheta_old),param_len,*convergence))) {

  setP.iter=main_loop;
  if (start) {
    for(i=0;i<param_len;i++) pdTheta[i]=pdTheta_in[i];
    transformTheta(pdTheta,t_pdTheta,param_len);
    setHistory(t_pdTheta,0,0,param_len,history_full);
    for(i=0;i<t_samp;i++){
      params[i].caseP.mu[0] = pdTheta[0];
      params[i].caseP.mu[1] = pdTheta[1];
    }
    setP.Sigma[0][0] = pdTheta[2];
    setP.Sigma[1][1] = pdTheta[3];
    setP.Sigma[0][1] = pdTheta[4]*sqrt(pdTheta[2]*pdTheta[3]);
    setP.Sigma[1][0] = setP.Sigma[0][1];
    dinv2D((double*)&setP.Sigma[0][0], 2, (double*)&setP.InvSigma[0][0], "Start of main loop");
    //for SEM
    for(i=0;i<param_len;i++) setP.semDone[i]=0;
    if(setP.fixedRho) setP.semDone[4]=1; //no need to worry about last row
    start=0;
  }


  if (setP.verbose>=1) {
    Rprintf("cycle %d/%d: %5g %5g %5g %5g rho: %5g",main_loop,*iteration_max,pdTheta[0],pdTheta[1],pdTheta[2],pdTheta[3],pdTheta[4]);
    if (setP.calcLoglik==1 && main_loop>2)
      Rprintf(" Prev LL: %5g",Suff[5]);
    Rprintf("\n");
  }
  //keep the old theta around for comaprison
  for(i=0;i<param_len;i++) pdTheta_old[i]=pdTheta[i];
  transformTheta(pdTheta_old,t_pdTheta_old,param_len);



  ecoEStep(params, Suff);
  if (!setP.ncar)
    ecoMStep(Suff,pdTheta,params);
  else
    ecoMStepNCAR(Suff,pdTheta,params);
  transformTheta(pdTheta,t_pdTheta,param_len);
  //char ch;
  //scanf(" %c", &ch );

  //if we're in the second run through of SEM
  if (setP.sem==1) {
    ecoSEM(optTheta, pdTheta, params, Rmat_old, Rmat);
  }
  else {
    setHistory(t_pdTheta,(main_loop<=1) ? 0 : Suff[5],main_loop,param_len,history_full);
//    Rprintf("hist %d/%d: %5g %5g %5g %5g rho: %5g ll: %5g",main_loop,*iteration_max,history_full[main_loop][0],history_full[main_loop][1],history_full[main_loop][2],history_full[main_loop][3],history_full[main_loop][4],history_full[main_loop][5]);
  }


  if (setP.verbose>=2) {
    Rprintf("theta and suff\n");
    Rprintf("%10g%10g%10g%10g%10g (%10g)\n",pdTheta[0],pdTheta[1],pdTheta[2],pdTheta[3],pdTheta[4],pdTheta[4]*sqrt(pdTheta[2]*pdTheta[3]));
    Rprintf("%10g%10g%10g%10g%10g\n",Suff[0],Suff[1],Suff[2],Suff[3],Suff[4]);
  }
  main_loop++;
  R_FlushConsole();
  R_CheckUserInterrupt();
}

/***End main loop ***/

if (setP.verbose>=1) printf("End loop PDT:%5g %5g %5g %5g %5g\n",pdTheta[0],pdTheta[1],pdTheta[2],pdTheta[3],pdTheta[4]);

//finish up: record results and loglik
Param* param;
Suff[5]=0.0;
for(i=0;i<t_samp;i++) {
  if(i<n_samp) {
   param=&(params[i]);
   for(j=0;j<2;j++)
     inSample[i*2+j]=param->caseP.W[j];
    setBounds(param);
    setNormConst(param);
    Suff[5]+=getLogLikelihood(param);
  }
}
//set the DM matrix (only matters for SEM)
if (setP.sem==1) {
  int DMlen=param_len-setP.fixedRho;
  for(i=0;i<DMlen;i++)
    for(j=0;j<DMlen;j++)
      DMmatrix[i*param_len+j]=Rmat[i][j];
}

*itersUsed=main_loop;
for(i=0;i<(*itersUsed);i++) {
  for(j=0;j<(param_len+1);j++)
    history[i*6+j]=history_full[i][j];
}

/* write out the random seed */
PutRNGstate();

/* Freeing the memory */
Free(pdTheta_old);
//FreeMatrix(Rmat_old,5);
//FreeMatrix(Rmat,5);
}

/*
 * input: optTheta,pdTheta,params,Rmat
 * output: param_lenxparam_len matrices Rmat and Rmat_old
 * optTheta is optimal theta
 * pdTheta is current theta
 * Rmat_old contains the input Rmat
 */
 void ecoSEM(double* optTheta, double* pdTheta, Param* params, double Rmat_old[7][7], double Rmat[7][7]) {
    //assume we have optTheta, ie \hat{phi}
    //pdTheta is phi^{t+1}
    int i,j,verbose,len;
    double SuffSem[5]; //sufficient stats
    double phiTI[7]; //phi^t_i
    double phiTp1I[7]; //phi^{t+1}_i
    double t_optTheta[7]; //transformed optimal
    double t_phiTI[7]; //transformed phi^t_i
    double t_phiTp1I[7]; //transformed phi^{t+1}_i
    Param* params_sem=(Param*) Calloc(params->setP->t_samp,Param);
    setParam setP_sem=*(params[0].setP);
    verbose=setP_sem.verbose;
    len=setP_sem.param_len-setP_sem.fixedRho;
    //first, save old Rmat
    for(i=0;i<len;i++)
      for(j=0;j<len;j++)
        Rmat_old[i][j]=Rmat[i][j];

    for(i=0;i<len;i++) {
      if (!setP_sem.semDone[i]) { //we're not done with this row
        //step 1: set phi^t_i
        if (verbose>=2) Rprintf("Theta: ");
        for(j=0;j<len;j++) {
          if (i==j)
            phiTI[j]=pdTheta[j]; //current value
          else phiTI[j]=optTheta[j]; //optimal value
          if (verbose>=2) Rprintf(" %5g ", phiTI[j]);
        }
        if (setP_sem.fixedRho) {
          phiTI[len-1]=pdTheta[len-1];
          phiTp1I[len-1]=pdTheta[len-1];
          if (verbose>=2) Rprintf(" %5g ", phiTI[4]);
        }
        if (verbose>=2) Rprintf("\n");


        //step 2: run an E-step and an M-step with phi^t_i
        //initialize params
        for(j=0;j<setP_sem.t_samp;j++) {
          params_sem[j].setP=&setP_sem;
          params_sem[j].caseP=params[j].caseP;
          params_sem[j].caseP.mu[0] = phiTI[0];
          params_sem[j].caseP.mu[1] = phiTI[1];
        }
        setP_sem.Sigma[0][0] = phiTI[2];
        setP_sem.Sigma[1][1] = phiTI[3];
        setP_sem.Sigma[0][1] = phiTI[4]*sqrt(phiTI[2]*phiTI[3]);
        setP_sem.Sigma[1][0] = setP_sem.Sigma[0][1];
        if (setP_sem.ncar) {
          setP_sem.Sigma3[0][0] = phiTI[2];
          setP_sem.Sigma3[1][1] = phiTI[3];
          setP_sem.Sigma3[0][1] = phiTI[4];
          setP_sem.Sigma3[0][2] = phiTI[5];
          setP_sem.Sigma3[1][2] = phiTI[6];
          setP_sem.Sigma3[1][0] = setP_sem.Sigma3[0][1];
          setP_sem.Sigma3[2][0] = setP_sem.Sigma3[0][2];
          setP_sem.Sigma3[2][1] = setP_sem.Sigma3[1][2];
          setP_sem.Sigma3[2][2] = setP_sem.Sigma3[2][2]; //constant
          initNCAR(params,phiTI[0],phiTI[1]);
        }

        if (verbose>=2) {
          Rprintf("Sigma: %5g %5g %5g %5g\n",setP_sem.Sigma[0][0],setP_sem.Sigma[0][1],setP_sem.Sigma[1][0],setP_sem.Sigma[1][1]);
        }
        dinv2D((double*)(&(setP_sem.Sigma[0][0])), 2, (double*)(&(setP_sem.InvSigma[0][0])), "ecoSem");

        ecoEStep(params_sem, SuffSem);
        if (!params[0].setP->ncar)
          ecoMStep(SuffSem,phiTp1I,params_sem);
        else
          ecoMStepNCAR(SuffSem,phiTp1I,params_sem);

        //step 3: create new R matrix row
        transformTheta(phiTp1I,t_phiTp1I,setP_sem.param_len);
        transformTheta(optTheta,t_optTheta,setP_sem.param_len);
        transformTheta(phiTI,t_phiTI,setP_sem.param_len);
        for(j = 0; j<len; j++)
          Rmat[i][j]=(t_phiTp1I[j]-t_optTheta[j])/(t_phiTI[i]-t_optTheta[i]);

        //step 4: check for difference
        params[0].setP->semDone[i]=closeEnough((double*)Rmat[i],(double*)Rmat_old[i],len,sqrt(params[0].setP->convergence));

      }
      else { //keep row the same
        for(j = 0; j<len; j++)
          Rmat[i][j]=Rmat_old[i][j];
      }
    }
    if(verbose>=1) {
      for(i=0;i<len;i++) {
        Rprintf("\nR Matrix row %d (%s): ", (i+1), (params[0].setP->semDone[i]) ? "Done" : "Not done");
        for(j=0;j<len;j++) {
          Rprintf(" %6.4g ",Rmat[i][j]);
        }
      }
      Rprintf("\n");
    }
    Free(params_sem);
 }


/*
 * The E-step for parametric ecological inference
 * Takes in a Param array of length n_samp + t_samp + x0_samp + x1_samp
 * Suff should be an array of length 5
 * On exit: suff holds the sufficient statistics as follows
 * suff[0]=E[W1*]
 * suff[1]=E[W2*]
 * suff[2]=E[W1*^2]
 * suff[3]=E[W1*W2*]
 * suff[4]=E[W2*^2]
 * suff[5]=log liklihood
 */

void ecoEStep(Param* params, double* suff) {

int t_samp,n_samp,s_samp,x1_samp,x0_samp,i,j,temp0,temp1, verbose;
double loglik,testdens;
Param* param; setParam* setP; caseParam* caseP;
setP=params[0].setP;
verbose=setP->verbose;

t_samp=setP->t_samp;
n_samp=setP->n_samp;
x1_samp=setP->x1_samp;
x0_samp=setP->x0_samp;
s_samp=setP->s_samp;

  double **Wstar=doubleMatrix(t_samp,5);     /* pseudo data(transformed)*/
loglik=0;
if (verbose>=2) Rprintf("E-step start\n");
  for (i = 0; i<n_samp; i++) {
    param = &(params[i]);
    caseP=&(param->caseP);
    if (caseP->Y>=.99 || caseP->Y<=.01) { //if Y is near the edge, then W1 and W2 are very constrained
      Wstar[i][0]=logit(caseP->Y,"Y maxmin W1");
      Wstar[i][1]=logit(caseP->Y,"Y maxmin W2");
      Wstar[i][2]=Wstar[i][0]*Wstar[i][0];
      Wstar[i][3]=Wstar[i][0]*Wstar[i][1];
      Wstar[i][4]=Wstar[i][1]*Wstar[i][1];
      caseP->Wstar[0]=Wstar[i][0];
      caseP->Wstar[1]=Wstar[i][1];
      caseP->W[0]=caseP->Y;
      caseP->W[1]=caseP->Y;
    }
    else {
      setBounds(param); //I think you only have to do this once...check later
//if (verbose>=3) Rprintf("past bounds\n");
      setNormConst(param);
      for (j=0;j<5;j++) {
        caseP->suff=j;
        Wstar[i][j]=paramIntegration(&SuffExp,param);
        if (j<2)
          caseP->Wstar[j]=Wstar[i][j];
      }
      caseP->suff=5;
      caseP->W[0]=paramIntegration(&SuffExp,param);;
      caseP->suff=6;
      caseP->W[1]=paramIntegration(&SuffExp,param);;
      caseP->suff=-1;
      testdens=paramIntegration(&SuffExp,param);
      if (setP->calcLoglik==1 && setP->iter>1) loglik+=getLogLikelihood(param);

   //report error E1 if E[W1],E[W2] is not on the tomography line
  if (fabs(caseP->W[0]-getW1FromW2(caseP->X, caseP->Y,caseP->W[1]))>0.01)
    Rprintf("E1 %d %5g %5g %5g %5g %5g %5g %5g %5g \n", i, caseP->X, caseP->Y, caseP->normcT, Wstar[i][0],Wstar[i][1],Wstar[i][2],Wstar[i][3],Wstar[i][4]);
  //report error E2 if Jensen's inequality doesn't hold
  if (Wstar[i][4]<pow(Wstar[i][1],2) || Wstar[i][2]<pow(Wstar[i][0],2))
     Rprintf("E2 %d %5g %5g %5g %5g %5g %5g %5g %5g\n", i, caseP->X, caseP->Y, caseP->normcT, Wstar[i][0],Wstar[i][1],Wstar[i][2],Wstar[i][3],Wstar[i][4]);
  //used for debugging if necessary
  if (verbose>=2 && i<20)
     Rprintf("%d %4g %4g %4g %4g %4g %4g %4g %4g %4g\n", i, caseP->X, caseP->Y, caseP->mu[0], param->setP->Sigma[0][1], caseP->normcT, caseP->W[0],caseP->W[1],Wstar[i][2],Wstar[i][3]);
    }
  }

    /* analytically compute E{W2_i|Y_i} given W1_i, mu and Sigma in x1 homeogeneous areas */
    for (i=n_samp; i<n_samp+x1_samp; i++) {
      temp0=params[i].caseP.Wstar[0];
      temp1=params[i].caseP.mu[1]+setP->Sigma[0][1]/setP->Sigma[0][0]*(temp0-params[i].caseP.mu[0]);
      Wstar[i][0]=temp0;
      Wstar[i][1]=temp1;
      Wstar[i][2]=temp0*temp0;
      Wstar[i][3]=temp0*temp1;
      Wstar[i][4]=temp1*temp1;
    }

  /*analytically compute E{W1_i|Y_i} given W2_i, mu and Sigma in x0 homeogeneous areas */
    for (i=n_samp+x1_samp; i<n_samp+x1_samp+x0_samp; i++) {
      temp1=params[i].caseP.Wstar[1];
      temp0=params[i].caseP.mu[0]+setP->Sigma[0][1]/setP->Sigma[1][1]*(temp1-params[i].caseP.mu[1]);
      Wstar[i][0]=temp0;
      Wstar[i][1]=temp1;
      Wstar[i][2]=temp0*temp0;
      Wstar[i][3]=temp0*temp1;
      Wstar[i][4]=temp1*temp1;
    }

    /* Use the values given by the survey data */
    for (i=n_samp+x1_samp+x0_samp; i<n_samp+x1_samp+x0_samp+s_samp; i++) {
      Wstar[i][0]=params[i].caseP.Wstar[0];
      Wstar[i][1]=params[i].caseP.Wstar[1];
      Wstar[i][2]=Wstar[i][0]*Wstar[i][0];
      Wstar[i][3]=Wstar[i][0]*Wstar[i][1];
      Wstar[i][4]=Wstar[i][1]*Wstar[i][1];
    }


  /*Calculate sufficient statistics */
  for (j=0; j<5; j++)
    suff[j]=0;

  /* compute sufficient statistics */
  for (i=0; i<t_samp; i++) {
    suff[0]+=Wstar[i][0];  /* sumE(W_i1|Y_i) */
    suff[1]+=Wstar[i][1];  /* sumE(W_i2|Y_i) */
    suff[2]+=Wstar[i][2];  /* sumE(W_i1^2|Y_i) */
    suff[3]+=Wstar[i][4];  /* sumE(W_i2^2|Y_i) */
    suff[4]+=Wstar[i][3];  /* sumE(W_i1^W_i2|Y_i) */
  }

  for(j=0; j<5; j++)
    suff[j]=suff[j]/t_samp;

  //if(verbose>=1) Rprintf("Log liklihood %15g\n",loglik);
  suff[5]=loglik;

FreeMatrix(Wstar,t_samp);

}

//Standard M-Step
//input: Suff
//output or mutated: pdTheta, params
void ecoMStep(double* Suff, double* pdTheta, Param* params) {

int i;
setParam* setP=params[0].setP;

  pdTheta[0]=Suff[0];  /*mu1*/
  pdTheta[1]=Suff[1];  /*mu2*/
  if (!setP->fixedRho) { //standard
    pdTheta[2]=Suff[2]-2*Suff[0]*pdTheta[0]+pdTheta[0]*pdTheta[0];  //sigma11
    pdTheta[3]=Suff[3]-2*Suff[1]*pdTheta[1]+pdTheta[1]*pdTheta[1];  //sigma22
    pdTheta[4]=Suff[4]-Suff[0]*pdTheta[1]-Suff[1]*pdTheta[0]+pdTheta[0]*pdTheta[1]; //sigma12
    pdTheta[4]=pdTheta[4]/sqrt(pdTheta[2]*pdTheta[3]); /*rho*/
  }
  else { //fixed rho
    double Imat[2][2];
    Imat[0][0]=Suff[2]-2*pdTheta[0]*Suff[0]+pdTheta[0]*pdTheta[0];  //I_11
    Imat[1][1]=Suff[3]-2*Suff[1]*pdTheta[1]+pdTheta[1]*pdTheta[1];  //I_22
    Imat[0][1]=Suff[4]-Suff[0]*pdTheta[1]-Suff[1]*pdTheta[0]+pdTheta[0]*pdTheta[1];  //I_12
    pdTheta[2]=(Imat[0][0]-pdTheta[4]*Imat[0][1]*pow(Imat[0][0]/Imat[1][1],0.5))/(1-pdTheta[4]*pdTheta[4]); //sigma11
    pdTheta[3]=(Imat[1][1]-pdTheta[4]*Imat[0][1]*pow(Imat[1][1]/Imat[0][0],0.5))/(1-pdTheta[4]*pdTheta[4]); //sigma22
    //sigma12 will be determined below by rho
  }

    //set Sigma
  setP->Sigma[0][0] = pdTheta[2];
  setP->Sigma[1][1] = pdTheta[3];
  setP->Sigma[0][1] = pdTheta[4]*sqrt(pdTheta[2]*pdTheta[3]);
  setP->Sigma[1][0] = setP->Sigma[0][1];

  //if(setP->verbose>=3) Rprintf("Sigma mstep: %5g %5g %5g %5g\n",setP->Sigma[0][0],setP->Sigma[0][1],setP->Sigma[1][0],setP->Sigma[1][1]);
  dinv2D((double*)(&(setP->Sigma[0][0])), 2, (double*)(&(setP->InvSigma[0][0])),"regular M-step");

  /* assign each data point the new mu (same for all points) */
  for(i=0;i<setP->t_samp;i++) {
    params[i].caseP.mu[0]=pdTheta[0];
    params[i].caseP.mu[1]=pdTheta[1];
  }
}


//M-Step under NCAR
void ecoMStepNCAR(double* Suff, double* pdTheta, Param* params) {

  setParam* setP=params[0].setP;
  //double Sigma[2][2]=setP->Sigma;
  //double[2][2] InvSigma=setP->InvSigma;
  //double[3][3] Sigma3=setP->Sigma3;   /* covariance matrix*/
  //double[3][3] InvSigma3=setP->Sigma3;   /* inverse covariance matrix*/
  int i,verbose,t_samp;
  verbose=setP->verbose;
  t_samp=setP->t_samp;

  //find mean of X
  double mu3=0; double mu3sq=0; double XW1=0; double XW2=0;
  char ebuffer[30];

  for(i=0;i<setP->t_samp;i++) {
    sprintf(ebuffer, "mstep X %i", i);
    double lx= logit(params[i].caseP.X,ebuffer);
    mu3 += lx;
    mu3sq += lx*lx;
    XW1 += params[i].caseP.Wstar[0]*lx;
    XW2 += params[i].caseP.Wstar[1]*lx;
  }
  mu3 = mu3/t_samp; mu3sq = mu3sq/t_samp;
  XW1 = XW1/t_samp; XW2 = XW2/t_samp;
  setP->mu3=mu3;

  pdTheta[0]=Suff[0];  /*mu1*/
  pdTheta[1]=Suff[1];  /*mu2*/

  //set variances and correlations
  pdTheta[2]=Suff[2]-2*Suff[0]*pdTheta[0]+pdTheta[0]*pdTheta[0]; //s11
  pdTheta[3]=Suff[3]-2*Suff[1]*pdTheta[1]+pdTheta[1]*pdTheta[1]; //s22
  pdTheta[4]=Suff[4]-Suff[0]*pdTheta[1]-Suff[1]*pdTheta[0]+pdTheta[0]*pdTheta[1]; //sigma12
  pdTheta[4]=pdTheta[4]/sqrt(pdTheta[2]*pdTheta[3]); //rho_12
  pdTheta[5]=(XW1 - mu3*Suff[0])/sqrt((Suff[2] - Suff[0]*Suff[0])*(setP->Sigma3[2][2] - mu3*mu3)); //rho_13
  pdTheta[6]=(XW2 - mu3*Suff[1])/sqrt((Suff[4] - Suff[1]*Suff[1])*(setP->Sigma3[2][2] - mu3*mu3)); //rho_23


  if (!setP->fixedRho) { //variable rho
    //variances
    setP->Sigma3[0][0] = pdTheta[2];
    setP->Sigma3[1][1] = pdTheta[3];
    setP->Sigma3[2][2] = mu3sq-mu3*mu3;

    //covariances
    setP->Sigma3[0][1] = pdTheta[4]*sqrt(pdTheta[2]*pdTheta[3]);
    setP->Sigma3[0][2] = pdTheta[5]*sqrt(pdTheta[2]*mu3);
    setP->Sigma3[1][2] = pdTheta[6]*sqrt(pdTheta[3]*mu3);

    //symmetry
    setP->Sigma3[1][0] = setP->Sigma3[0][1];
    setP->Sigma3[2][0] = setP->Sigma3[0][2];
    setP->Sigma3[2][1] = setP->Sigma3[1][2];

  }
  else { //fixed rho NEEDS WORK
    double Imat[2][2];
    Imat[0][0]=setP->Sigma3[0][0];
    Imat[1][1]=setP->Sigma3[1][1];
    Imat[0][1]=setP->Sigma3[0][1];
    pdTheta[2]=(Imat[0][0]-pdTheta[6]*Imat[0][1]*pow(Imat[0][0]/Imat[1][1],0.5))/(1-pdTheta[6]*pdTheta[6]); //sigma11
    pdTheta[3]=(Imat[1][1]-pdTheta[6]*Imat[0][1]*pow(Imat[1][1]/Imat[0][0],0.5))/(1-pdTheta[6]*pdTheta[6]); //sigma22
    pdTheta[4]=setP->Sigma3[0][1];
    pdTheta[5]=setP->Sigma3[0][2];
    setP->Sigma3[0][0]=pdTheta[2];
    setP->Sigma3[1][1]=pdTheta[3];
    setP->Sigma3[1][0] = pdTheta[6]*sqrt(pdTheta[2]*pdTheta[3]);
    setP->Sigma3[1][0] = setP->Sigma3[0][1];
  }

  initNCAR(params,pdTheta[0],pdTheta[1]);

}

//NCAR initialize
void initNCAR(Param* params,double mu0, double mu1) {
  setParam* setP=params[0].setP;
  /* OLD CODE
  //Set 2x2 Sigma
  setP->Sigma[0][0]= setP->Sigma3[0][0] - setP->Sigma3[0][2]*setP->Sigma3[0][2]/setP->Sigma3[2][2];
  setP->Sigma[0][1]= setP->Sigma3[0][1] - setP->Sigma3[0][2]*setP->Sigma3[1][2]/setP->Sigma3[2][2];
  setP->Sigma[1][0]= setP->Sigma[0][1];
  setP->Sigma[1][1]= setP->Sigma3[1][1] - setP->Sigma3[1][2]*setP->Sigma3[1][2]/setP->Sigma3[2][2];

  dinv2D((double*)(&(setP->Sigma[0][0])), 2, (double*)(&(setP->InvSigma[0][0])),"NCAR M-step S2");
  dinv2D((double*)(&(setP->Sigma3[0][0])), 3, (double*)(&(setP->InvSigma3[0][0])),"NCAR M-step S3");

  //assign each data point the new mu (different for each point)
  for(i=0;i<t_samp;i++) {
    params[i].caseP.mu[0]=pdTheta[0] + (setP->Sigma3[0][2]/setP->Sigma3[2][2])*(params[i].caseP.X-mu3);
    params[i].caseP.mu[1]=pdTheta[1] + (setP->Sigma3[1][2]/setP->Sigma3[2][2])*(params[i].caseP.X-mu3);
  }*/

  //NEW CODE, PLEASE CHECK
  setP->Sigma[0][0]= setP->Sigma3[0][0]*(1 - setP->Sigma3[0][2]*setP->Sigma3[0][2]);
  setP->Sigma[0][1]= setP->Sigma3[0][1]*(1 - setP->Sigma3[1][2]*setP->Sigma3[1][2])*pow((1 - setP->Sigma3[0][2])*(1 - setP->Sigma3[1][2]),0.5);
  setP->Sigma[1][0]= setP->Sigma[0][1];
  setP->Sigma[1][1]= setP->Sigma3[1][1]*(1 - setP->Sigma3[1][2]*setP->Sigma3[1][2]);

  dinv2D((double*)(&(setP->Sigma[0][0])), 2, (double*)(&(setP->InvSigma[0][0])),"NCAR M-step S2");
  dinv2D((double*)(&(setP->Sigma3[0][0])), 3, (double*)(&(setP->InvSigma3[0][0])),"NCAR M-step S3");

  //assign each data point the new mu (different for each point)
  int i;
  for(i=0;i<setP->t_samp;i++) {
    params[i].caseP.mu[0]=mu0 + setP->Sigma3[0][2]*pow(setP->Sigma3[0][0]/setP->Sigma3[2][2],0.5)*(logit(params[i].caseP.X,"initNCAR mu0")-setP->mu3);
    params[i].caseP.mu[1]=mu1 + setP->Sigma3[1][2]*pow(setP->Sigma3[1][1]/setP->Sigma3[2][2],0.5)*(logit(params[i].caseP.X,"initNCAR mu1")-setP->mu3);
  }
}

/*
 * Read in the data set and population params
 */
 void readData(Param* params, int n_dim, double* pdX, double* sur_W, double* x1_W1, double* x0_W2,
                int n_samp, int s_samp, int x1_samp, int x0_samp) {
     /* read the data set */
int itemp,i,j;
double dtemp;
setParam* setP=params[0].setP;

  /** Packing Y, X  **/
  itemp = 0;
  for (j = 0; j < n_dim; j++)
    for (i = 0; i < n_samp; i++) {
      params[i].caseP.data[j] = pdX[itemp++];
    }

  for (i = 0; i < n_samp; i++) {
    params[i].caseP.X=params[i].caseP.data[0];
    params[i].caseP.Y=params[i].caseP.data[1];
    //fix X edge cases
    params[i].caseP.X=(params[i].caseP.X >= 1) ? .9999 : ((params[i].caseP.X <= 0) ? 0.0001 : params[i].caseP.X);
  }

  if (setP->verbose>=2) {
    printf("Y X\n");
    for(i=0;i<10;i++)
      Rprintf("%5d%14g%14g\n",i,params[i].caseP.Y,params[i].caseP.X);
      }

  /*read homeogenous areas information */
    for (i=n_samp; i<n_samp+x1_samp; i++) {
      params[i].caseP.W[0]=(x1_W1[i] == 1) ? .9999 : ((x1_W1[i]==0) ? .0001 : x1_W1[i]);
      params[i].caseP.Wstar[0]=logit(params[i].caseP.W[0],"X1 read");
    }

    for (i=n_samp+x1_samp; i<n_samp+x1_samp+x0_samp; i++) {
      params[i].caseP.W[1]=(x0_W2[i] == 1) ? .9999 : ((x0_W2[i]==0) ? .0001 : x0_W2[i]);
      params[i].caseP.Wstar[1]=logit(params[i].caseP.W[1],"X0 read");
    }


  /*read the survey data */
     itemp=0;
    for (j=0; j<n_dim; j++)
      for (i=n_samp+x1_samp+x0_samp; i<n_samp+x1_samp+x0_samp+s_samp; i++) {
        dtemp=sur_W[itemp++];
        params[i].caseP.W[j]=(dtemp == 1) ? .9999 : ((dtemp==0) ? .0001 : dtemp);
        params[i].caseP.Wstar[j]=logit(params[i].caseP.W[j],"Survey read");
      }
 }

/*
 * Parameterizes the elements of theta
 * Input: pdTheta
 * Mutates: t_pdTheta
 */
void transformTheta(double* pdTheta, double* t_pdTheta, int len) {
  t_pdTheta[0]=pdTheta[0];
  t_pdTheta[1]=pdTheta[1];
  t_pdTheta[2]=log(pdTheta[2]);
  t_pdTheta[3]=log(pdTheta[3]);
  t_pdTheta[4]=.5*(log(1+pdTheta[4])-log(1-pdTheta[4]));
  if (len>5) {
    t_pdTheta[5]=.5*(log(1+pdTheta[5])-log(1-pdTheta[6]));
    t_pdTheta[6]=.5*(log(1+pdTheta[5])-log(1-pdTheta[6]));
  }
}

void untransformTheta(double* t_pdTheta,double* pdTheta, int len) {
  pdTheta[0]=t_pdTheta[0];
  pdTheta[1]=t_pdTheta[1];
  pdTheta[2]=exp(t_pdTheta[2]);
  pdTheta[3]=exp(t_pdTheta[3]);
  pdTheta[4]=(exp(2*t_pdTheta[4])-1)/(exp(2*t_pdTheta[4])+1);
  if (len>5) {
    pdTheta[5]=(exp(2*t_pdTheta[5])-1)/(exp(2*t_pdTheta[5])+1);
    pdTheta[6]=(exp(2*t_pdTheta[6])-1)/(exp(2*t_pdTheta[6])+1);
  }
}

/*
 * Input transformed theta, loglikelihood, iteration
 * Mutates: history_full
 */
void setHistory(double* t_pdTheta, double loglik, int iter,int len,double history_full[][8]) {
  int j;
  for(j=0;j<len;j++)
    history_full[iter][j]=t_pdTheta[j];
  history_full[iter][len]=loglik;
}

/*
 * Determines whether we have converged
 * Takes in the current and old (one step previous) array of theta values
 * maxerr is the maximum difference two corresponding values can have before the
 *  function returns false
 */
int closeEnough(double* pdTheta, double* pdTheta_old, int len, double maxerr) {
  int j;
  for(j = 0; j<len; j++)
    if (fabs(pdTheta[j]-pdTheta_old[j])>=maxerr) return 0;
  return 1;
}

int semDoneCheck(setParam* setP) {
  int j,len;
  len=setP->param_len - setP->fixedRho; //4 if there is a fixed rho
  for(j=0;j<len;j++)
    if(setP->semDone[j]==0) return 0;
  return 1;
}

void gridEStep(Param* params, int n_samp, int s_samp, int x1_samp, int x0_samp, double* suff, int verbose, double minW1, double maxW1) {

  int n_dim=2;
  int n_step=5000;    /* The default size of grid step */
  int ndraw=10000;
  int trapod=0;       /* 1 if use trapozodial ~= in numer. int.*/
  int *n_grid=intArray(n_samp);                /* grid size */
  double **W1g=doubleMatrix(n_samp, n_step);   /* grids for W1 */
  double **W2g=doubleMatrix(n_samp, n_step);   /* grids for W2 */
  double *vtemp=doubleArray(n_dim);
  int *mflag=intArray(n_step);
  double *prob_grid=doubleArray(n_step);
  double *prob_grid_cum=doubleArray(n_step);
  double **X=doubleMatrix(n_samp,n_dim);     /* Y and covariates */

  int itemp,i,j,k,t_samp;
  double dtemp,dtemp1,temp0,temp1;

  t_samp=n_samp+x1_samp+x0_samp+s_samp;

  double **W=doubleMatrix(t_samp,n_dim);     /* W1 and W2 matrix */
  double **Wstar=doubleMatrix(t_samp,5);     /* pseudo data(transformed*/

  for (i=0;i<t_samp;i++)
    for(j=0;j<n_dim;j++)
      X[i][j]=params[i].caseP.data[j];

  GridPrep(W1g, W2g, (double**) params[i].caseP.data, (double*)&maxW1, (double*)&minW1, n_grid, n_samp, n_step);

    for (i=0; i<n_step; i++) {
    mflag[i]=0;
  }


  //update W, Wstar given mu, Sigma in regular areas
  for (i=0;i<n_samp;i++){
    if ( params[i].caseP.Y!=0 && params[i].caseP.Y!=1 ) {
      // project BVN(mu, Sigma) on the inth tomo line
      dtemp=0;
      for (j=0;j<n_grid[i];j++){
        vtemp[0]=log(W1g[i][j])-log(1-W1g[i][j]);
        vtemp[1]=log(W2g[i][j])-log(1-W2g[i][j]);
        prob_grid[j]=dMVN(vtemp, params[i].caseP.mu, (double**)(params[i].setP->InvSigma), 2, 1) -
          log(W1g[i][j])-log(W2g[i][j])-log(1-W1g[i][j])-log(1-W2g[i][j]);
        prob_grid[j]=exp(prob_grid[j]);
        dtemp+=prob_grid[j];
        prob_grid_cum[j]=dtemp;
      }
      for (j=0;j<n_grid[i];j++){
        prob_grid_cum[j]/=dtemp; //standardize prob.grid
      }
      // MC numerical integration, compute E(W_i|Y_i, X_i, theta)
      //2 sample ndraw W_i on the ith tomo line
      //   use inverse CDF method to draw
      //   0-1 by 1/ndraw approx uniform distribution
      //3 compute Wsta_i from W_i
      j=0;
      itemp=1;

      for (k=0; k<ndraw; k++){
        j=findInterval(prob_grid_cum, n_grid[i],
		      (double)(1+k)/(ndraw+1), 1, 1, itemp, mflag);
        itemp=j-1;


        if ((W1g[i][j]==0) || (W1g[i][j]==1))
          Rprintf("W1g%5d%5d%14g", i, j, W1g[i][j]);
        if ((W2g[i][j]==0) || (W2g[i][j]==1))
          Rprintf("W2g%5d%5d%14g", i, j, W2g[i][j]);

        if (j==0 || trapod==0) {
          W[i][0]=W1g[i][j];
          W[i][1]=W2g[i][j];
        }
        else if (j>=1 && trapod==1) {
          if (prob_grid_cum[j]!=prob_grid_cum[(j-1)]) {
            dtemp1=((double)(1+k)/(ndraw+1)-prob_grid_cum[(j-1)])/(prob_grid_cum[j]-prob_grid_cum[(j-1)]);
            W[i][0]=dtemp1*(W1g[i][j]-W1g[i][(j-1)])+W1g[i][(j-1)];
            W[i][1]=dtemp1*(W2g[i][j]-W2g[i][(j-1)])+W2g[i][(j-1)];
          }
          else if (prob_grid_cum[j]==prob_grid_cum[(j-1)]) {
            W[i][0]=W1g[i][j];
            W[i][1]=W2g[i][j];
          }
        }
        temp0=log(W[i][0])-log(1-W[i][0]);
        temp1=log(W[i][1])-log(1-W[i][1]);
        Wstar[i][0]+=temp0;
        Wstar[i][1]+=temp1;
        Wstar[i][2]+=temp0*temp0;
        Wstar[i][3]+=temp0*temp1;
        Wstar[i][4]+=temp1*temp1;
      }
    }
  }

  // compute E_{W_i|Y_i} for n_samp
  for (i=0; i<n_samp; i++) {
    if ( X[i][1]!=0 && X[i][1]!=1 ) {
      Wstar[i][0]/=ndraw;  //E(W1i)
      Wstar[i][1]/=ndraw;  //E(W2i)
      Wstar[i][2]/=ndraw;  //E(W1i^2)
      Wstar[i][3]/=ndraw;  //E(W1iW2i)
      Wstar[i][4]/=ndraw;  //E(W2i^2)
    }
  } //for x0type, x1type and survey data, E-step is either the observed value or the analytical expectation

  /* compute sufficient statistics */
  for (j=0; j<5; j++)
    suff[j]=0;

  for (i=0; i<t_samp; i++) {
    suff[0]+=Wstar[i][0];  /* sumE(W_i1|Y_i) */
    suff[1]+=Wstar[i][1];  /* sumE(W_i2|Y_i) */
    suff[2]+=Wstar[i][2];  /* sumE(W_i1^2|Y_i) */
    suff[3]+=Wstar[i][4];  /* sumE(W_i2^2|Y_i) */
    suff[4]+=Wstar[i][3];  /* sumE(W_i1^W_i2|Y_i) */
  }


  for(j=0; j<5; j++)
    suff[j]=suff[j]/t_samp;

  Free(n_grid);Free(vtemp);Free(mflag);Free(prob_grid);Free(prob_grid_cum);
  FreeMatrix(W1g,n_samp);FreeMatrix(W2g,n_samp);FreeMatrix(X,n_samp);
  FreeMatrix(W,t_samp);FreeMatrix(Wstar,t_samp);

}
