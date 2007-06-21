#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
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
void ecoMStepCCAR(double* pdTheta, Param* params);
void MStepHypTest(Param* params, double* pdTheta);
void initTheta(double* pdTheta_in,Param* params, double* pdTheta);
void initNCAR(Param* params, double* pdTheta);
void setHistory(double* t_pdTheta, double loglik, int iter,setParam* setP,double history_full[][10]);
int closeEnough(double* pdTheta, double* pdTheta_old, int len, double maxerr);
int semDoneCheck(setParam* setP);
void gridEStep(Param* params, int n_samp, int s_samp, int x1_samp, int x0_samp, double* suff, int verbose, double minW1, double maxW1);
void transformTheta(double* pdTheta, double* t_pdTheta, int len, setParam* setP);
void untransformTheta(double* t_pdTheta,double* pdTheta, int len, setParam* setP);
void ncarFixedRhoTransform(double* pdTheta);
void ncarFixedRhoUnTransform(double* pdTheta);
void printColumnHeader(int main_loop, int iteration_max, setParam* setP, int finalTheta);

/**
 * Main function.
 * Important mutations (i.e., outputs): pdTheta, Suff, DMmatrix, history
 * See internal comments for details.
 */

void cEMeco(
	    /*data input */
	    double *pdX,         /* data (X, Y) */
	    double *pdTheta_in,  /* Theta^ t
				    CAR: mu1, mu2, var1, var2, rho
				    NCAR: mu1, mu2, var1, var2, p13,p13,p12*/
	    int *pin_samp,       /* sample size */

	    /* loop vairables */
	    int *iteration_max,          /* number of maximum iterations */
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
	    int *hypTest_L,   /* number of hypothesis constraints */
	    double *optTheta,  /*optimal theta obtained from previous EM result; if set, then we're doing SEM*/

	    /* storage */
      //Theta under CAR: mu1,mu2,s1,s2,p12
      //Theta under NCAR: mu_3, mu_1, mu_2, sig_3, sig_1, sig_2, r_13, r_23, r_12
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
  //int t_samp=n_samp+s_samp+x1_samp+x0_samp;  /* total sample size*/
  int t_samp=n_samp+s_samp;  /* total sample size, ignoring homog data*/
  int n_dim=2;        /* dimensions */

  setParam setP;
  //set options
  setP.ncar=bit(*flag,0);
  setP.fixedRho=bit(*flag,1);
  setP.sem=bit(*flag,2) & (optTheta[2]!=-1.1);
  setP.ccar=0; setP.ccar_nvar=0;

  //hard-coded hypothesis test
  //hypTest is the number of constraints.  hyptTest==0 when we're not checking a hypothesis
  setP.hypTest=(*hypTest_L);
  if (setP.hypTest>1) error("Unable to do hypothesis testing with more than one constraint");
  if (setP.hypTest==1) {
    setP.hypTestCoeff=doubleMatrix(setP.ncar ? 3 : 2,setP.hypTest);
    setP.hypTestCoeff[0][0]=1; setP.hypTestCoeff[1][0]=-1;
    if (setP.ncar) setP.hypTestCoeff[2][0]=0;
    setP.hypTestResult=0;
  }

  setP.verbose=*verbosiosity;
  if (setP.verbose>=1) Rprintf("OPTIONS::  Ncar: %s; Fixed Rho: %s; SEM: %s\n",setP.ncar==1 ? "Yes" : "No",
   setP.fixedRho==1 ? "Yes" : "No",setP.sem==1 ? "Second run" : (bit(*flag,2)==1 ? "First run" : "No"));
  setP.calcLoglik=*calcLoglik;
  setP.convergence=*convergence;
  setP.t_samp=t_samp; setP.n_samp=n_samp; setP.s_samp=s_samp; setP.x1_samp=x1_samp; setP.x0_samp=x0_samp;
  int param_len=setP.ccar ? setP.ccar_nvar : (setP.ncar ? 9 : 5);
  setP.param_len=param_len;
  setP.pdTheta=doubleArray(param_len);
  setP.suffstat_len=(setP.ncar ? 9 : 5);
  setP.SigmaK=doubleMatrix(param_len,param_len); //CCAR
  setP.InvSigmaK=doubleMatrix(param_len,param_len); //CCAR

  /* model parameters */
  //double **Sigma=doubleMatrix(n_dim,n_dim);/* inverse covariance matrix*/
  //double **InvSigma=doubleMatrix(n_dim,n_dim);/* inverse covariance matrix*/

  double *pdTheta_old=doubleArray(param_len);
  double *t_pdTheta=doubleArray(param_len); //transformed theta
  double *t_pdTheta_old=doubleArray(param_len);
  double Rmat_old[7][7];
  double Rmat[7][7];
  double history_full[*iteration_max+1][10];

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
      initTheta(pdTheta_in,params,pdTheta);
      transformTheta(pdTheta,t_pdTheta,param_len, &setP);
      setHistory(t_pdTheta,0,0,(setParam*)&setP,history_full);
      if (!setP.ncar) {
        for(i=0;i<t_samp;i++) {
          params[i].caseP.mu[0] = pdTheta[0];
          params[i].caseP.mu[1] = pdTheta[1];
        }
        setP.Sigma[0][0] = pdTheta[2];
        setP.Sigma[1][1] = pdTheta[3];
        setP.Sigma[0][1] = pdTheta[4]*sqrt(pdTheta[2]*pdTheta[3]);
        setP.Sigma[1][0] = setP.Sigma[0][1];
        dinv2D((double*)&setP.Sigma[0][0], 2, (double*)&setP.InvSigma[0][0], "Start of main loop");
      }
      else {
        if (setP.fixedRho) ncarFixedRhoTransform(pdTheta);
        initNCAR(params,pdTheta);
        if (setP.fixedRho) ncarFixedRhoUnTransform(pdTheta);
      }
      start=0;
    }
    for(i=0;i<param_len;i++) setP.pdTheta[i]=pdTheta[i];

    if (setP.verbose>=1) {
      if ((main_loop - 1) % 15 == 0) printColumnHeader(main_loop,*iteration_max,&setP,0);

      Rprintf("cycle %d/%d:",main_loop,*iteration_max);
      for(i=0;i<param_len;i++)
        if (setP.varParam[i]) {
          if (pdTheta[i]>=0) Rprintf("% 5.3f",pdTheta[i]);
          else Rprintf(" % 5.2f",pdTheta[i]);
        }
      if (setP.calcLoglik==1 && main_loop>2)
        Rprintf(" Prev LL: %5.2f",Suff[setP.suffstat_len]);
      Rprintf("\n");
    }
    //keep the old theta around for comaprison
    for(i=0;i<param_len;i++) pdTheta_old[i]=pdTheta[i];
    transformTheta(pdTheta_old,t_pdTheta_old,param_len,&setP);


    ecoEStep(params, Suff);
    if (!setP.ncar)
      ecoMStep(Suff,pdTheta,params);
    else
      ecoMStepNCAR(Suff,pdTheta,params);
    transformTheta(pdTheta,t_pdTheta,param_len,&setP);
    //char ch;
    //scanf(" %c", &ch );

    //if we're in the second run through of SEM
    if (setP.sem==1) {
      ecoSEM(optTheta, pdTheta, params, Rmat_old, Rmat);
    }
    else {
      setHistory(t_pdTheta,(main_loop<=1) ? 0 : Suff[setP.suffstat_len],main_loop,(setParam*)&setP,history_full);
    }


    if (setP.verbose>=2) {
      Rprintf("theta and suff\n");
      if (param_len>5) {
        Rprintf("%10g%10g%10g%10g%10g%10g%10g%10g%10g\n",pdTheta[0],pdTheta[1],pdTheta[2],pdTheta[3],pdTheta[4],pdTheta[5],pdTheta[6],pdTheta[7],pdTheta[8]);
      }
      else {
        Rprintf("%10g%10g%10g%10g%10g (%10g)\n",pdTheta[0],pdTheta[1],pdTheta[2],pdTheta[3],pdTheta[4],pdTheta[4]*sqrt(pdTheta[2]*pdTheta[3]));
      }
      Rprintf("%10g%10g%10g%10g%10g\n",Suff[0],Suff[1],Suff[2],Suff[3],Suff[4]);
      Rprintf("Sig: %10g%10g%10g\n",setP.Sigma[0][0],setP.Sigma[1][1],setP.Sigma[0][1]);
      if (setP.ncar) Rprintf("Sig3: %10g%10g%10g%10g\n",setP.Sigma3[0][0],setP.Sigma3[1][1],setP.Sigma3[2][2]);
      //char x;
      //R_ReadConsole("hit enter\n",(char*)&x,4,0);
    }
    main_loop++;
    R_FlushConsole();
    R_CheckUserInterrupt();
  }

  /***End main loop ***/
  //finish up: record results and loglik
  Param* param;
  Suff[setP.suffstat_len]=0.0;
  for(i=0;i<param_len;i++) setP.pdTheta[i]=pdTheta[i];
  for(i=0;i<t_samp;i++) {
     param=&(params[i]);
    if(i<n_samp) {
     for(j=0;j<2;j++) inSample[i*2+j]=param->caseP.W[j];
      //setBounds(param);
      //setNormConst(param);
    }
    Suff[setP.suffstat_len]+=getLogLikelihood(param);
  }

  if (setP.verbose>=1) {
    printColumnHeader(main_loop,*iteration_max,&setP,1);
    Rprintf("Final Theta:");
      for(i=0;i<param_len;i++) {
        if (pdTheta[i]>=0) Rprintf("% 5.3f",pdTheta[i]);
        else Rprintf(" % 5.2f",pdTheta[i]);
      }
      if (setP.calcLoglik==1 && main_loop>2) {
        Rprintf(" Final LL: %5.2f",Suff[setP.suffstat_len]);
        history_full[main_loop-1][param_len]=Suff[setP.suffstat_len];
      }
      Rprintf("\n");
    }

  //set the DM matrix (only matters for SEM)
  if (setP.sem==1) {
    int DMlen=0;
    for(i=0; i<param_len;i++)
      if(setP.varParam[i]) DMlen++;
    for(i=0;i<DMlen;i++)
      for(j=0;j<DMlen;j++)
        DMmatrix[i*DMlen+j]=Rmat[i][j];
  }

  *itersUsed=main_loop;
  for(i=0;i<(*itersUsed);i++) {
    for(j=0;j<(param_len+1);j++)
      history[i*(param_len+1)+j]=history_full[i][j];
  }


  /* write out the random seed */
  PutRNGstate();

  /* Freeing the memory */
  Free(pdTheta_old);
  //FreeMatrix(Rmat_old,5);
  //FreeMatrix(Rmat,5);
  }

/**
 * initializes Theta, varParam, and semDone
 * input: pdTheta_in,params
 * mutates: params.setP, pdTheta
 * CAR theta: mu_1, mu_2, sig_1, sig_2, rho_12
 * NCAR theta: mu_3, mu_1, mu_2, sig_3, sig_1, sig_2, r_13, r_23, r_12
 */
void initTheta(double* pdTheta_in,Param* params, double* pdTheta) {
  setParam* setP=params[0].setP;
  int param_len=setP->param_len;
  int i;
  if (!setP->ncar) {
    for(i=0;i<param_len;i++) {
      pdTheta[i]=pdTheta_in[i];
      setP->varParam[i]=1;
    }
    if (setP->fixedRho) setP->varParam[4]=0;
  }
  else {
    //constants
    double lx,mu3sq;
    pdTheta[0]=0; mu3sq=0;
    for(i=0;i<setP->t_samp;i++) {
      lx=logit(params[i].caseP.X,"initpdTheta0");
      pdTheta[0] += lx;
      mu3sq += lx*lx;
    }
    pdTheta[0] = pdTheta[0]/setP->t_samp;
    mu3sq = mu3sq/setP->t_samp;
    pdTheta[3] = mu3sq-pdTheta[0]*pdTheta[0]; //variance
    //fill from pdTheta_in
    pdTheta[1]=pdTheta_in[0];
    pdTheta[2]=pdTheta_in[1];
    pdTheta[4]=pdTheta_in[2];
    pdTheta[5]=pdTheta_in[3];
    pdTheta[6]=pdTheta_in[4];
    pdTheta[7]=pdTheta_in[5];
    pdTheta[8]=pdTheta_in[6];
    for(i=0;i<param_len;i++) setP->varParam[i]=1;
    setP->varParam[0]=0;setP->varParam[3]=0;
    //if (setP->fixedRho) setP->varParam[8]=0;
  }
  int varlen=0;
  for(i=0; i<param_len;i++)
    if(setP->varParam[i]) varlen++;
  for(i=0; i<varlen;i++)
      setP->semDone[i]=0;
}

/**
  * The E-step for parametric ecological inference
  * Takes in a Param array of length n_samp + t_samp + x0_samp + x1_samp
  * Suff should be an array with the same length as the number of params (+1)
  * On exit: suff holds the sufficient statistics and loglik as follows
  * CAR: (0) E[W1*] (1) E[W2*] (2) E[W1*^2] (3) E[W2*^2] (4) E[W1*W2*] (5) loglik
  * NCAR: (0) X, (1) W1, (2) W2, (3) X^2, (4) W1^2, (5) W2^2, (6) x*W1, (7) X*W2, (8) W1*W2, (9) loglik
 **/

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
  if (verbose>=3 && !setP->sem) Rprintf("E-step start\n");
  for (i = 0; i<n_samp; i++) {
    param = &(params[i]);
    caseP=&(param->caseP);
    if (caseP->Y>=.990 || caseP->Y<=.010) { //if Y is near the edge, then W1 and W2 are very constrained
      Wstar[i][0]=logit(caseP->Y,"Y maxmin W1");
      Wstar[i][1]=logit(caseP->Y,"Y maxmin W2");
      Wstar[i][2]=Wstar[i][0]*Wstar[i][0];
      Wstar[i][3]=Wstar[i][0]*Wstar[i][1];
      Wstar[i][4]=Wstar[i][1]*Wstar[i][1];
      caseP->Wstar[0]=Wstar[i][0];
      caseP->Wstar[1]=Wstar[i][1];
      caseP->W[0]=caseP->Y;
      caseP->W[1]=caseP->Y;
      if (setP->calcLoglik==1 && setP->iter>1) loglik+=getLogLikelihood(param);
      //Rprintf("Skipping %d, Y=%5g",i,caseP->Y);
    }
    else {
      setBounds(param); //I think you only have to do this once...check later
      /*if (verbose>=2 && setP->iter==12 && i==422) {
        Rprintf("Bounds: %5g %5g %5g %5g\n",caseP->Wbounds[0][0],caseP->Wbounds[0][1],caseP->Wbounds[1][0],caseP->Wbounds[1][1]);
        setP->weirdness=1;
      }
      else setP->weirdness=0;*/

      setNormConst(param);
      for (j=0;j<5;j++) {
        caseP->suff=j;
        Wstar[i][j]=paramIntegration(&SuffExp,param);
        if (j<2)
          caseP->Wstar[j]=Wstar[i][j];
      }
      caseP->suff=SS_W1;
      caseP->W[0]=paramIntegration(&SuffExp,param);
      caseP->suff=SS_W2;
      caseP->W[1]=paramIntegration(&SuffExp,param);
      caseP->suff=SS_Test;
      testdens=paramIntegration(&SuffExp,param);
      if (setP->calcLoglik==1 && setP->iter>1) loglik+=getLogLikelihood(param);

      //report error E1 if E[W1],E[W2] is not on the tomography line
      if (fabs(caseP->W[0]-getW1FromW2(caseP->X, caseP->Y,caseP->W[1]))>0.011) {
        Rprintf("E1 %d %5g %5g %5g %5g %5g %5g %5g %5g err:%5g\n", i, caseP->X, caseP->Y, caseP->mu[0], caseP->mu[1], caseP->normcT,Wstar[i][0],Wstar[i][1],Wstar[i][2],fabs(caseP->W[0]-getW1FromW2(caseP->X, caseP->Y,caseP->W[1])));
        char ch;
        scanf("Hit enter to continue %c\n", &ch );
      }
      //report error E2 if Jensen's inequality doesn't hold
      if (Wstar[i][4]<pow(Wstar[i][1],2) || Wstar[i][2]<pow(Wstar[i][0],2))
        Rprintf("E2 %d %5g %5g %5g %5g %5g %5g %5g %5g\n", i, caseP->X, caseP->Y, caseP->normcT, caseP->mu[1],Wstar[i][0],Wstar[i][1],Wstar[i][2],Wstar[i][4]);
      //used for debugging if necessary
      if (verbose>=2 && !setP->sem && ((i<10 && verbose>=3) || (caseP->mu[1] < -1.7 && caseP->mu[0] > 1.4)))
        Rprintf("%d %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n", i, caseP->X, caseP->Y, caseP->mu[0], caseP->mu[1], param->setP->Sigma[0][1], caseP->normcT, caseP->W[0],caseP->W[1],Wstar[i][2]);
    }
  }


  /* Use the values given by the survey data */
  //Calculate loglik also
  for (i=n_samp; i<n_samp+s_samp; i++) {
    param = &(params[i]);
    caseP=&(param->caseP);
    Wstar[i][0]=caseP->Wstar[0];
    Wstar[i][1]=caseP->Wstar[1];
    Wstar[i][2]=Wstar[i][0]*Wstar[i][0];
    Wstar[i][3]=Wstar[i][0]*Wstar[i][1];
    Wstar[i][4]=Wstar[i][1]*Wstar[i][1];
    if (setP->calcLoglik==1 && setP->iter>1) loglik+=getLogLikelihood(param);
  }

  /* analytically compute E{W2_i|Y_i} given W1_i, mu and Sigma in x1 homeogeneous areas */
  for (i=n_samp+s_samp; i<n_samp+s_samp+x1_samp; i++) {
    /*temp0=params[i].caseP.Wstar[0];
    temp1=params[i].caseP.mu[1]+setP->Sigma[0][1]/setP->Sigma[0][0]*(temp0-params[i].caseP.mu[0]);
    Wstar[i][0]=temp0;
    Wstar[i][1]=temp1;
    Wstar[i][2]=temp0*temp0;
    Wstar[i][3]=temp0*temp1;
    Wstar[i][4]=temp1*temp1;*/
  }

  /*analytically compute E{W1_i|Y_i} given W2_i, mu and Sigma in x0 homeogeneous areas */
  for (i=n_samp+s_samp+x1_samp; i<n_samp+s_samp+x1_samp+x0_samp; i++) {
    /*temp1=params[i].caseP.Wstar[1];
    temp0=params[i].caseP.mu[0]+setP->Sigma[0][1]/setP->Sigma[1][1]*(temp1-params[i].caseP.mu[1]);
    Wstar[i][0]=temp0;
    Wstar[i][1]=temp1;
    Wstar[i][2]=temp0*temp0;
    Wstar[i][3]=temp0*temp1;
    Wstar[i][4]=temp1*temp1;*/
  }


  /*Calculate sufficient statistics */
  for (j=0; j<setP->suffstat_len; j++)
    suff[j]=0;


  //CAR: (0) E[W1*] (1) E[W2*] (2) E[W1*^2] (3) E[W2*^2] (4) E[W1*W2*] (5) loglik
  //NCAR: (0) X, (1) W1, (2) W2, (3) X^2, (4) W1^2, (5) W2^2, (6) x*W1, (7) X*W2, (8) W1*W2, (9) loglik
  /* compute sufficient statistics */
  for (i=0; i<t_samp; i++) {
    if (!setP->ncar) {
      suff[0] += Wstar[i][0];  /* sumE(W_i1|Y_i) */
      suff[1] += Wstar[i][1];  /* sumE(W_i2|Y_i) */
      suff[2] += Wstar[i][2];  /* sumE(W_i1^2|Y_i) */
      suff[3] += Wstar[i][4];  /* sumE(W_i2^2|Y_i) */
      suff[4] += Wstar[i][3];  /* sumE(W_i1*W_i2|Y_i) */
    }
    else if (setP->ncar) {
      double lx= logit(params[i].caseP.X,"mstep X");
      suff[0] += lx;
      suff[1] += Wstar[i][0];
      suff[2] += Wstar[i][1];
      suff[3] += lx*lx;
      suff[4] += Wstar[i][2];
      suff[5] += Wstar[i][4];
      suff[6] += params[i].caseP.Wstar[0]*lx;
      suff[7] += params[i].caseP.Wstar[1]*lx;
      suff[8] += Wstar[i][3];
    }
  }

  for(j=0; j<setP->suffstat_len; j++)
    suff[j]=suff[j]/t_samp;
  //Rprintf("%5g suff0,2,4 %5g %5g %5g\n",setP->pdTheta[6],suff[0],suff[2],suff[4]);
  //if(verbose>=1) Rprintf("Log liklihood %15g\n",loglik);
  suff[setP->suffstat_len]=loglik;

  FreeMatrix(Wstar,t_samp);
}

/**
 * CAR M-Step
 * inputs: Suff (sufficient statistics)
 *    CAR Suff: E[W1], E[W2], E[W1^2], E[W2^2], E[W1W2]
 * mutated (i.e., output): pdTheta, params
 */
void ecoMStep(double* Suff, double* pdTheta, Param* params) {

  int i;
  setParam* setP=params[0].setP;

  pdTheta[0]=Suff[0];  /*mu1*/
  pdTheta[1]=Suff[1];  /*mu2*/

  if (setP->hypTest>0) {
    MStepHypTest(params,pdTheta);
  }

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


/**
 * M-Step under NCAR
 * Input: Suff (sufficient statistics)
 *    (0) X, (1) W1, (2) W2, (3) X^2, (4) W1^2, (5) W2^2, (6) x*W1, (7) X*W2, (8) W1*W2, (9) loglik
 * mutated (i.e., output): pdTheta, params
 */
void ecoMStepNCAR(double* Suff, double* pdTheta, Param* params) {

  setParam* setP=params[0].setP;
  //double Sigma[2][2]=setP->Sigma;
  //double[2][2] InvSigma=setP->InvSigma;
  //double[3][3] Sigma3=setP->Sigma3;   /* covariance matrix*/
  //double[3][3] InvSigma3=setP->Sigma3;   /* inverse covariance matrix*/
  int ii,i,j,verbose,t_samp;
  verbose=setP->verbose;
  t_samp=setP->t_samp;


  //set E[XW*]
  double XW1=Suff[6];
  double XW2=Suff[7];



  //for(i = 0;i<9; i++) Rprintf("%f5.2\n",pdTheta[i]);
  if (!setP->fixedRho) { //variable rho


    //pdTheta[0] is const
    pdTheta[1]=Suff[1];  /*mu1*/
    pdTheta[2]=Suff[2];  /*mu2*/

    //set variances and correlations
    //pdTheta[3] is const
    pdTheta[4]=Suff[4]-2*Suff[1]*pdTheta[1]+pdTheta[1]*pdTheta[1]; //s11
    pdTheta[5]=Suff[5]-2*Suff[2]*pdTheta[2]+pdTheta[2]*pdTheta[2]; //s22
    pdTheta[6]=(XW1 - pdTheta[0]*Suff[1])/sqrt((Suff[4] - Suff[1]*Suff[1])*pdTheta[3]); //rho_13
    pdTheta[7]=(XW2 - pdTheta[0]*Suff[2])/sqrt((Suff[5] - Suff[2]*Suff[2])*pdTheta[3]); //rho_23
    pdTheta[8]=Suff[8]-Suff[1]*pdTheta[2]-Suff[2]*pdTheta[1]+pdTheta[1]*pdTheta[2]; //sigma12
    pdTheta[8]=pdTheta[8]/sqrt(pdTheta[4]*pdTheta[5]); //rho_12


    //reference: (0) mu_3, (1) mu_1, (2) mu_2, (3) sig_3, (4) sig_1, (5) sig_2, (6) r_13, (7) r_23, (8) r_12
    //variances
    setP->Sigma3[0][0] = pdTheta[4];
    setP->Sigma3[1][1] = pdTheta[5];
    setP->Sigma3[2][2] = pdTheta[3];

    //covariances
    setP->Sigma3[0][1] = pdTheta[8]*sqrt(pdTheta[4]*pdTheta[5]);
    setP->Sigma3[0][2] = pdTheta[6]*sqrt(pdTheta[4]*pdTheta[3]);
    setP->Sigma3[1][2] = pdTheta[7]*sqrt(pdTheta[5]*pdTheta[3]);

    //symmetry
    setP->Sigma3[1][0] = setP->Sigma3[0][1];
    setP->Sigma3[2][0] = setP->Sigma3[0][2];
    setP->Sigma3[2][1] = setP->Sigma3[1][2];
              //if (verbose>=2) {
            //Rprintf("Sigma3: %5g %5g %5g %5g %5g\n",setP->Sigma3[0][0],setP->Sigma3[0][1],setP->Sigma3[1][1],setP->Sigma3[1][2],setP->Sigma3[2][2]);
          //}

  }
  else { //fixed rho
    //reference: (0) mu_3, (1) mu_1, (2) mu_2, (3) sig_3, (4) sig_1 | 3, (5) sig_2 | 3, (6) beta1, (7) beta2, (8) r_12 | 3

    ncarFixedRhoTransform(pdTheta); //need the fixed param (pdTheta[8]) to be the conditional correlation

    //CODE BLOCK D
    //compute beta based on previous sigma
    //beta is mu1,beta1,mu2,beta, which are pdTheta 1,2,6,7
    double **InvSigma=doubleMatrix(2,2);
    double **Zmat=doubleMatrix(4,2);
    double **Zmat_t=doubleMatrix(2,4);
    double **tmp41=doubleMatrix(4,1);
    double **tmp42=doubleMatrix(4,2);
    double **tmp44=doubleMatrix(4,4);
    double **tmp21=doubleMatrix(2,1);
    double **denom=doubleMatrix(4,4);
    double **numer=doubleMatrix(4,1);
    for (i=0;i<4;i++) {
      for(j=0;j<4;j++) {
        if (j<2) {
          if (i<2) InvSigma[i][j]=setP->InvSigma[i][j];
          Zmat[i][j]=0; Zmat_t[j][i]=0;
        }
        denom[i][j]=0;
      }
      numer[i][0]=0;
    }
    //Rprintf("InvSigma %5g %5g %5g\n",InvSigma[0][0],InvSigma[1][1],InvSigma[0][1]);
    for(ii=0;ii<setP->t_samp;ii++) {
        double lx=logit(params[ii].caseP.X,"NCAR beta");
        for(j=0;j<2;j++) {
          Zmat_t[j][j*2+1]=lx - pdTheta[0];
          Zmat_t[j][j*2]=1;
          Zmat[j*2+1][j]=lx - pdTheta[0];
          Zmat[j*2][j]=1;
        }
        matrixMul(Zmat,InvSigma,4,2,2,2,tmp42);
        matrixMul(tmp42,Zmat_t,4,2,2,4,tmp44);
        for (i=0;i<4;i++)
          for(j=0;j<4;j++)
            denom[i][j]+=tmp44[i][j];
        //for (i=0;i<2;i++) tmp21[i][0]=(params[ii].caseP.Wstar[i] - pdTheta[i+1]); //Wtilde ??
        for (i=0;i<2;i++) tmp21[i][0]=params[ii].caseP.Wstar[i]; //Wstar
        //matrixMul(Zmat,InvSigma,4,2,2,2,tmp42);  //no need to repeat calculation
        matrixMul(tmp42,tmp21,4,2,2,1,tmp41);
        for (i=0;i<4;i++) numer[i][0]+=tmp41[i][0];
    }
    dinv(denom,4,denom);
    matrixMul(denom,numer,4,4,4,1,numer);

    pdTheta[1]=numer[0][0]; //mu1
    pdTheta[6]=numer[1][0]; //beta1
    pdTheta[2]=numer[2][0]; //mu2
    pdTheta[7]=numer[3][0]; //beta2
    //pdTheta[8] is constant
    //Rprintf("Compare Suff1 %5g to pdT1 %5g \n",Suff[1],pdTheta[1]);
    //Rprintf("Compare Suff2 %5g to pdT2 %5g \n",Suff[2],pdTheta[2]);

    if (setP->hypTest>0) {
      MStepHypTest(params,pdTheta);
    }

    //CAR: (0) E[W1*] (1) E[W2*] (2) E[W1*^2] (3) E[W2*^2] (4) E[W1*W2*] (5) loglik
    //NCAR: (0) X, (1) W1, (2) W2, (3) X^2, (4) W1^2, (5) W2^2, (6) x*W1, (7) X*W2, (8) W1*W2, (9) loglik
    //0->1, 1->2, 2->4, 3->5, 4->8


    //CODE BLOCK C
    //Compute sigma conditional on beta
    //reference: (0) mu_3, (1) mu_1, (2) mu_2, (3) sig_3, (4) sig_1 | 3, (5) sig_2 | 3, (6) beta1, (7) beta2, (8) r_12 | 3
    double Smat[2][2]; //the S matrix (divided by n) in the paper
    double Tmat[2][2]; //the T matrix (divided by n) in the paper
    double S1=Suff[1]; //S_1 = Sufficient stat of W1* - beta1 * (sum of [(X_i - \mu3)]) ; second term goes to zero
    double S2=Suff[2]; //S_2 =  Sufficient stat of W2*

    Smat[0][0]=Suff[4] - 2*pdTheta[6]*(XW1 - pdTheta[0]*Suff[1]) + pdTheta[6]*pdTheta[6]*pdTheta[3];  //S_11
    Smat[1][1]=Suff[5] - 2*pdTheta[7]*(XW2 - pdTheta[0]*Suff[2]) + pdTheta[7]*pdTheta[7]*pdTheta[3];  //S_22
    Smat[0][1]=Suff[8] - pdTheta[6]*(XW2 - pdTheta[0]*Suff[2]) - pdTheta[7]*(XW1 - pdTheta[0]*Suff[1]) + pdTheta[6]*pdTheta[7]*pdTheta[3] ;  //S_12

    Tmat[0][0]=Smat[0][0] - S1*S1;
    Tmat[1][1]=Smat[1][1] - S2*S2;
    Tmat[0][1]=Smat[0][1] - S1*S2;

    pdTheta[4]=(Tmat[0][0]-pdTheta[8]*Tmat[0][1]*pow(Tmat[0][0]/Tmat[1][1],0.5))/(1-pdTheta[8]*pdTheta[8]); //sigma11 | 3
    pdTheta[5]=(Tmat[1][1]-pdTheta[8]*Tmat[0][1]*pow(Tmat[1][1]/Tmat[0][0],0.5))/(1-pdTheta[8]*pdTheta[8]); //sigma22 | 3

    //variances
    //CODE BLOCK B
    setP->Sigma3[0][0] = pdTheta[4] + pdTheta[6]*pdTheta[6]*pdTheta[3];
    setP->Sigma3[1][1] = pdTheta[5] + pdTheta[7]*pdTheta[7]*pdTheta[3];
    setP->Sigma3[2][2] = pdTheta[3];

    //covariances
    setP->Sigma3[0][1] = (pdTheta[8]*sqrt(pdTheta[4]*pdTheta[5]) + pdTheta[6]*pdTheta[7]*pdTheta[3])/
                          (sqrt((pdTheta[4] + pdTheta[6]*pdTheta[6]*pdTheta[3])*(pdTheta[5] + pdTheta[7]*pdTheta[7]*pdTheta[3])));//rho_12 unconditional
    setP->Sigma3[0][1] = setP->Sigma3[0][1]*sqrt(setP->Sigma3[0][0]*setP->Sigma3[1][1]); //sig_12
    setP->Sigma3[0][2] = pdTheta[6]*sqrt((pdTheta[3])/(pdTheta[4] + pdTheta[6]*pdTheta[6]*pdTheta[3]))*sqrt(setP->Sigma3[0][0]*setP->Sigma3[2][2]);
    setP->Sigma3[1][2] = pdTheta[7]*sqrt((pdTheta[3])/(pdTheta[5] + pdTheta[7]*pdTheta[7]*pdTheta[3]))*sqrt(setP->Sigma3[1][1]*setP->Sigma3[2][2]);

    //symmetry
    setP->Sigma3[1][0] = setP->Sigma3[0][1];
    setP->Sigma3[2][0] = setP->Sigma3[0][2];
    setP->Sigma3[2][1] = setP->Sigma3[1][2];
  }
  dinv2D((double*)(&(setP->Sigma3[0][0])), 3, (double*)(&(setP->InvSigma3[0][0])),"NCAR M-step S3");
  initNCAR(params,pdTheta);
  if (setP->fixedRho) ncarFixedRhoUnTransform(pdTheta);
}

/**
 * M-Step under CCAR
 * Input: params
 * mutated (i.e., output): pdTheta, params
 */
void ecoMStepCCAR(double* pdTheta, Param* params) {
  setParam* setP=params[0].setP;
  int k=setP->ccar_nvar;
  int ii,i,j,verbose,t_samp;
  verbose=setP->verbose;
  t_samp=setP->t_samp;
  double **InvSigma=doubleMatrix(2,2);
  double **Z_i=doubleMatrix(k,2);
  double **Z_i_t=doubleMatrix(2,k);
  double **tmpk1=doubleMatrix(k,1);
  double **tmpk2=doubleMatrix(k,2);
  double **tmpkk=doubleMatrix(k,k);
  double **tmp21=doubleMatrix(2,1);
  double **tmp21_b=doubleMatrix(2,1);
  double **tmp12=doubleMatrix(1,2);
  double **tmp22=doubleMatrix(2,2);
  double **denom=doubleMatrix(k,k);
  double **numer=doubleMatrix(k,1);
  //betas
  for (i=0;i<k;i++) {
    for(j=0;j<k;j++) {
      if (j<2) {
        if (i<2) InvSigma[i][j]=setP->InvSigma[i][j];
      }
      denom[i][j]=0;
    }
    numer[i][0]=0;
  }
  //Rprintf("InvSigma %5g %5g %5g\n",InvSigma[0][0],InvSigma[1][1],InvSigma[0][1]);
  for(ii=0;ii<setP->t_samp;ii++) {
    for (i=0;i<k;i++) {
      for(j=0;j<k;j++) {
        Z_i[i][j]=params[ii].caseP.Z_i[i][j];
        Z_i_t[i][j]=params[ii].caseP.Z_i[j][i];
      }
    }
      matrixMul(Z_i,InvSigma,k,2,2,2,tmpk2);
      matrixMul(tmpk2,Z_i_t,k,2,2,k,tmpkk);
      for (i=0;i<k;i++)
        for(j=0;j<k;j++)
          denom[i][j]+=tmpkk[i][j];
      for (i=0;i<2;i++) tmp21[i][0]=params[ii].caseP.Wstar[i]; //Wstar
      matrixMul(tmpk2,tmp21,k,2,2,1,tmpk1);
      for (i=0;i<k;i++) numer[i][0]+=tmpk1[i][0];
  }
  dinv(denom,k,denom);
  matrixMul(denom,numer,k,k,k,1,numer);
  for(i=0; i<k;i++) pdTheta[i]=numer[i][0]; //betas


  if (setP->hypTest>0) {
    MStepHypTest(params,pdTheta);
  }

  //conditional Sigma
  //start at 0
  for(i=0; i<2;i++)
    for(j=0; j<2;j++)
      setP->Sigma[i][j] = 0;


  for(ii=0;ii<setP->t_samp;ii++) {
    for (i=0;i<k;i++) {
      for(j=0;j<k;j++) {
        Z_i_t[i][j]=params[ii].caseP.Z_i[j][i];
      }
    }
    matrixMul(Z_i_t,numer,2,k,k,1,tmp21_b);
    for (i=0;i<2;i++) tmp21[i][0]=params[ii].caseP.Wstar[i]; //Wstar
    for (i=0;i<2;i++) tmp21[i][0] = tmp21[i][0] - tmp21_b[i][0]; //Wstar - Z_t*B
    for (i=0;i<2;i++) tmp12[0][i] = tmp21[i][0]; //invserse
    matrixMul(tmp21,tmp12,2,1,1,2,tmp22);
    for(i=0; i<2;i++)
      for(j=0; j<2;j++)
        setP->Sigma[i][j] += tmp22[i][j];
  }
  dinv2D((double*)(&(setP->Sigma[0][0])), 2, (double*)(&(setP->InvSigma[0][0])),"CCAR M-step S2");

  //variances
  //CODE BLOCK B
  setP->Sigma3[0][0] = pdTheta[4] + pdTheta[6]*pdTheta[6]*pdTheta[3];
  setP->Sigma3[1][1] = pdTheta[5] + pdTheta[7]*pdTheta[7]*pdTheta[3];
  setP->Sigma3[2][2] = pdTheta[3];

  //covariances
  setP->Sigma3[0][1] = (pdTheta[8]*sqrt(pdTheta[4]*pdTheta[5]) + pdTheta[6]*pdTheta[7]*pdTheta[3])/
                        (sqrt((pdTheta[4] + pdTheta[6]*pdTheta[6]*pdTheta[3])*(pdTheta[5] + pdTheta[7]*pdTheta[7]*pdTheta[3])));//rho_12 unconditional
  setP->Sigma3[0][1] = setP->Sigma3[0][1]*sqrt(setP->Sigma3[0][0]*setP->Sigma3[1][1]); //sig_12
  setP->Sigma3[0][2] = pdTheta[6]*sqrt((pdTheta[3])/(pdTheta[4] + pdTheta[6]*pdTheta[6]*pdTheta[3]))*sqrt(setP->Sigma3[0][0]*setP->Sigma3[2][2]);
  setP->Sigma3[1][2] = pdTheta[7]*sqrt((pdTheta[3])/(pdTheta[5] + pdTheta[7]*pdTheta[7]*pdTheta[3]))*sqrt(setP->Sigma3[1][1]*setP->Sigma3[2][2]);

  //symmetry
  setP->Sigma3[1][0] = setP->Sigma3[0][1];
  setP->Sigma3[2][0] = setP->Sigma3[0][2];
  setP->Sigma3[2][1] = setP->Sigma3[1][2];

  dinv2D((double*)(&(setP->Sigma3[0][0])), 3, (double*)(&(setP->InvSigma3[0][0])),"NCAR M-step S3");
  initNCAR(params,pdTheta);

}

/**
 * Exta M-Step for hypothesis testing
 * Input: params
 * Mutates pdTheta
 */
void MStepHypTest(Param* params, double* pdTheta) {
  setParam* setP=params[0].setP;
  double offset,denom;
  int dim,i,j,l,k;
  dim=setP->ncar ? 3 : 2;
  l=setP->hypTest;
  double** Sigma=doubleMatrix(dim,dim);
  double** temp_LbyD=doubleMatrix(l,dim);
  double** temp_DbyL=doubleMatrix(dim,l);
  double** temp_LbyL=doubleMatrix(l,l);

  for(i=0;i<dim;i++)
    for(j=0;j<dim;j++) {
      if (dim==3) {
        Sigma[i][j]=setP->Sigma3[i][j];
      }
      else {
        Sigma[i][j]=setP->Sigma[i][j];
      }
    }
  //transpose
  double** hypTestCoeffT=doubleMatrix(l,dim);
  for(i=0;i<dim;i++) hypTestCoeffT[0][i]=setP->hypTestCoeff[i][0];

  //numerator
  for(k=0;k<2;k++) temp_DbyL[k][0]=0;
  for(i=0;i<setP->t_samp;i++) {
    temp_DbyL[0][0]+=params[i].caseP.Wstar[0];
    temp_DbyL[1][0]+=params[i].caseP.Wstar[1];
  }
  matrixMul(hypTestCoeffT,temp_DbyL,l,dim,dim,l,temp_LbyL);
  temp_LbyL[0][0]=temp_LbyL[0][0]-(setP->t_samp*setP->hypTestResult);
  matrixMul(Sigma,setP->hypTestCoeff,dim,dim,dim,l,temp_DbyL);
  for(k=0;k<2;k++) temp_DbyL[k][0]*=temp_LbyL[0][0];

  //denominator
  //matrixMul(hypTestCoeffT,InvSigma,l,dim,dim,dim,temp_LbyD);
  matrixMul(hypTestCoeffT,Sigma,l,dim,dim,dim,temp_LbyD);
  matrixMul(temp_LbyD,setP->hypTestCoeff,l,dim,dim,l,temp_LbyL);
  denom=setP->t_samp*temp_LbyL[0][0];

  //offset theta
  for(k=0;k<2;k++) {
    offset=temp_DbyL[k][0]/denom;
    int kindex= (setP->ncar) ? (k+1) : k;
    pdTheta[kindex]=pdTheta[kindex]-offset;
  }

}


/**
 * NCAR initialize
 * note that for fixed rho, the input is the UNTRANSFORMED PARAMETERS
 * input: pdTheta
 * mutates: params
 */
void initNCAR(Param* params, double* pdTheta) {
  setParam* setP=params[0].setP;
  int i;
  if (!setP->fixedRho) { //variable rho
    //reference: (0) mu_3, (1) mu_1, (2) mu_2, (3) sig_3, (4) sig_1, (5) sig_2, (6) r_13, (7) r_23, (8) r_12

    setP->Sigma[0][0]= pdTheta[4]*(1 - pdTheta[6]*pdTheta[6]);
    setP->Sigma[1][1]= pdTheta[5]*(1 - pdTheta[7]*pdTheta[7]);
    setP->Sigma[0][1]= (pdTheta[8] - pdTheta[6]*pdTheta[7])/sqrt((1 - pdTheta[6]*pdTheta[6])*(1 - pdTheta[7]*pdTheta[7])); //correlation
    setP->Sigma[0][1]= setP->Sigma[0][1]*sqrt(setP->Sigma[0][0]*setP->Sigma[1][1]); //covar
    setP->Sigma[1][0]= setP->Sigma[0][1]; //symmetry
    dinv2D((double*)(&(setP->Sigma[0][0])), 2, (double*)(&(setP->InvSigma[0][0])),"NCAR M-step S2");

    //assign each data point the new mu (different for each point)
    for(i=0;i<setP->t_samp;i++) {
      params[i].caseP.mu[0]=pdTheta[1] + pdTheta[6]*sqrt(pdTheta[4]/pdTheta[3])*(logit(params[i].caseP.X,"initNCAR mu0")-pdTheta[0]);
      params[i].caseP.mu[1]=pdTheta[2] + pdTheta[7]*sqrt(pdTheta[5]/pdTheta[3])*(logit(params[i].caseP.X,"initNCAR mu1")-pdTheta[0]);
      if(setP->verbose>=2 && !setP->sem && (i<3 || i==422))
      //if(setP->verbose>=2  && i<3)
        Rprintf("mu primes for %d: %5g %5g (mu2: %5g p7: %5g p5: %5g X-T: %5g)\n",i,params[i].caseP.mu[0],params[i].caseP.mu[1],pdTheta[2],pdTheta[7],pdTheta[5],logit(params[i].caseP.X,"initNCAR mu0")-pdTheta[0]);
    }
  }
  else { //fixed rho
    //reference: (0) mu_3, (1) mu_1, (2) mu_2, (3) sig_3, (4) sig_1 | 3, (5) sig_2 | 3, (6) beta1, (7) beta2, (8) r_12 | 3
    //CODE BLOCK A
    setP->Sigma[0][0]= pdTheta[4];
    setP->Sigma[1][1]= pdTheta[5];
    setP->Sigma[0][1]= pdTheta[8]*sqrt(pdTheta[4]*pdTheta[5]); //covar
    setP->Sigma[1][0]= setP->Sigma[0][1]; //symmetry
    dinv2D((double*)(&(setP->Sigma[0][0])), 2, (double*)(&(setP->InvSigma[0][0])),"NCAR M-step S2");

    for(i=0;i<setP->t_samp;i++) {
      params[i].caseP.mu[0]=pdTheta[1] + pdTheta[6]*(logit(params[i].caseP.X,"initNCAR mu0")-pdTheta[0]);
      params[i].caseP.mu[1]=pdTheta[2] + pdTheta[7]*(logit(params[i].caseP.X,"initNCAR mu1")-pdTheta[0]);
      if(setP->verbose>=2 && !setP->sem && (i<3 || i==422))
      //if(setP->verbose>=2  && i<3)
        Rprintf("mu primes for %d: %5g %5g (mu2: %5g p7: %5g p5: %5g X-T: %5g)\n",i,params[i].caseP.mu[0],params[i].caseP.mu[1],pdTheta[2],pdTheta[7],pdTheta[5],logit(params[i].caseP.X,"initNCAR mu0")-pdTheta[0]);
    }

  }
}

/**
 * CCAR initialize
 * Note that fixed rho is currently unimplemented
 * input: pdTheta
 * mutates: params
 */
void initCCAR(Param* params, double* pdTheta) {
  setParam* setP=params[0].setP;
  int i;
  if (!setP->fixedRho) { //variable rho
    //reference: (0) mu_3, (1) mu_1, (2) mu_2, (3) sig_3, (4) sig_1, (5) sig_2, (6) r_13, (7) r_23, (8) r_12

    setP->Sigma[0][0]= pdTheta[4]*(1 - pdTheta[6]*pdTheta[6]);
    setP->Sigma[1][1]= pdTheta[5]*(1 - pdTheta[7]*pdTheta[7]);
    setP->Sigma[0][1]= (pdTheta[8] - pdTheta[6]*pdTheta[7])/sqrt((1 - pdTheta[6]*pdTheta[6])*(1 - pdTheta[7]*pdTheta[7])); //correlation
    setP->Sigma[0][1]= setP->Sigma[0][1]*sqrt(setP->Sigma[0][0]*setP->Sigma[1][1]); //covar
    setP->Sigma[1][0]= setP->Sigma[0][1]; //symmetry
    dinv2D((double*)(&(setP->Sigma[0][0])), 2, (double*)(&(setP->InvSigma[0][0])),"NCAR M-step S2");
    //assign each data point the new mu (different for each point)
    for(i=0;i<setP->t_samp;i++) {
      params[i].caseP.mu[0]=pdTheta[1] + pdTheta[6]*sqrt(pdTheta[4]/pdTheta[3])*(logit(params[i].caseP.X,"initNCAR mu0")-pdTheta[0]);
      params[i].caseP.mu[1]=pdTheta[2] + pdTheta[7]*sqrt(pdTheta[5]/pdTheta[3])*(logit(params[i].caseP.X,"initNCAR mu1")-pdTheta[0]);
      if(setP->verbose>=2 && !setP->sem && (i<3 || i==422))
      //if(setP->verbose>=2  && i<3)
        Rprintf("mu primes for %d: %5g %5g (mu2: %5g p7: %5g p5: %5g X-T: %5g)\n",i,params[i].caseP.mu[0],params[i].caseP.mu[1],pdTheta[2],pdTheta[7],pdTheta[5],logit(params[i].caseP.X,"initNCAR mu0")-pdTheta[0]);
    }
  }
  else { //fixed rho
  }
}

/**
 * input: optTheta,pdTheta,params,Rmat
 * mutate/output: matrices Rmat and Rmat_old (dimensions of param_len x param_len)
 * optTheta is optimal theta
 * pdTheta is current theta
 * Rmat_old contains the input Rmat
 */
 void ecoSEM(double* optTheta, double* pdTheta, Param* params, double Rmat_old[7][7], double Rmat[7][7]) {
  //assume we have optTheta, ie \hat{phi}
  //pdTheta is phi^{t+1}
  int i,j,verbose,len,param_len;
  setParam setP_sem=*(params[0].setP);
  param_len=setP_sem.param_len;
  double *SuffSem=doubleArray(setP_sem.suffstat_len+1); //sufficient stats
  double phiTI[param_len]; //phi^t_i
  double phiTp1I[param_len]; //phi^{t+1}_i
  double t_optTheta[param_len]; //transformed optimal
  double t_phiTI[param_len]; //transformed phi^t_i
  double t_phiTp1I[param_len]; //transformed phi^{t+1}_i
  Param* params_sem=(Param*) Calloc(params->setP->t_samp,Param);
  verbose=setP_sem.verbose;
  //determine length of R matrix
  len=0;
  for(j=0; j<param_len;j++)
    if(setP_sem.varParam[j]) len++;

  //first, save old Rmat
  for(i=0;i<len;i++)
    for(j=0;j<len;j++)
      Rmat_old[i][j]=Rmat[i][j];

  for(i=0;i<len;i++) {
    if (!setP_sem.semDone[i]) { //we're not done with this row
      //step 1: set phi^t_i
      if (verbose>=2) Rprintf("Theta(%d):",(i+1));
      int switch_index_ir=0; int switch_index_it;
      for(j=0;j<param_len;j++) {
        if (!setP_sem.varParam[j]) //const
          phiTI[j]=optTheta[j];
        else {
          if (i==switch_index_ir) {
            phiTI[j]=pdTheta[j]; //current value
            switch_index_it=j;
          }
          else phiTI[j]=optTheta[j]; //optimal value
          switch_index_ir++;
        }
        if (verbose>=2) Rprintf(" %5g ", phiTI[j]);
      }
      //if (setP_sem.fixedRho) {
      //  phiTI[len-1]=pdTheta[len-1];
      //  phiTp1I[len-1]=pdTheta[len-1];
      // if (verbose>=2) Rprintf(" %5g ", phiTI[len-1]);
      //}
      if (verbose>=2) Rprintf("\n");
      for(j=0;j<param_len;j++) phiTp1I[j]=phiTI[j]; //init next iteration

      //step 2: run an E-step and an M-step with phi^t_i
      //initialize params
      if (!setP_sem.ncar) {
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
        dinv2D((double*)(&(setP_sem.Sigma[0][0])), 2, (double*)(&(setP_sem.InvSigma[0][0])), "SEM: CAR init ");
      }
      else {
        for(j=0;j<setP_sem.t_samp;j++) {
          params_sem[j].setP=&setP_sem;
          params_sem[j].caseP=params[j].caseP;
        }
        setP_sem.Sigma3[0][0] = phiTI[4];
        setP_sem.Sigma3[1][1] = phiTI[5];
        setP_sem.Sigma3[2][2] = phiTI[3];

        //covariances
        setP_sem.Sigma3[0][1] = phiTI[8]*sqrt(phiTI[4]*phiTI[5]);
        setP_sem.Sigma3[0][2] = phiTI[6]*sqrt(phiTI[4]*phiTI[3]);
        setP_sem.Sigma3[1][2] = phiTI[7]*sqrt(phiTI[5]*phiTI[3]);

        //symmetry
        setP_sem.Sigma3[1][0] = setP_sem.Sigma3[0][1];
        setP_sem.Sigma3[2][0] = setP_sem.Sigma3[0][2];
        setP_sem.Sigma3[2][1] = setP_sem.Sigma3[1][2];
        if (verbose>=2) {
          Rprintf("Sigma3: %5g %5g %5g %5g %5g %5g; %5g %5g\n",setP_sem.Sigma3[0][0],setP_sem.Sigma3[0][1],setP_sem.Sigma3[1][1],setP_sem.Sigma3[0][2],setP_sem.Sigma3[1][2],setP_sem.Sigma3[2][2],*(&(setP_sem.Sigma3[0][0])+0),*(&(setP_sem.Sigma3[0][0])+8));
        }
        dinv2D((double*)(&(setP_sem.Sigma3[0][0])), 3, (double*)(&(setP_sem.InvSigma3[0][0])),"SEM: NCAR Sig3 init");
        if (verbose>=2) {
          Rprintf("Check 1");
        }
        if (setP_sem.fixedRho) ncarFixedRhoTransform(phiTI);
        initNCAR(params_sem,phiTI);
        if (setP_sem.fixedRho) ncarFixedRhoUnTransform(phiTI);
        if (verbose>=2) {
          Rprintf("Check 2");
        }
      }

      //if (verbose>=2) {
      //  Rprintf("Sigma: %5g %5g %5g %5g\n",setP_sem.Sigma[0][0],setP_sem.Sigma[0][1],setP_sem.Sigma[1][0],setP_sem.Sigma[1][1]);
      //}

      ecoEStep(params_sem, SuffSem);
      if (!params[0].setP->ncar)
        ecoMStep(SuffSem,phiTp1I,params_sem);
      else
        ecoMStepNCAR(SuffSem,phiTp1I,params_sem);

      //step 3: create new R matrix row
      transformTheta(phiTp1I,t_phiTp1I,setP_sem.param_len,&setP_sem);
      transformTheta(optTheta,t_optTheta,setP_sem.param_len,&setP_sem);
      transformTheta(phiTI,t_phiTI,setP_sem.param_len,&setP_sem);
      /*if (verbose>=2) {
        Rprintf("T+1:");
        for (j=0;j<param_len;j++) Rprintf(" %5g ", phiTp1I[j]);
        Rprintf("\nOpt:");
        for (j=0;j<param_len;j++) Rprintf(" %5g ", optTheta[j]);
        Rprintf("\n 2nd item: %5g %5g %5g %5g", t_phiTp1I[2], t_optTheta[2], t_phiTI[switch_index_it], t_optTheta[switch_index_it]);
      }*/
      int index_jr=0;
      for(j = 0; j<param_len; j++) {
        if (setP_sem.varParam[j]) {
          Rmat[i][index_jr]=(t_phiTp1I[j]-t_optTheta[j])/(t_phiTI[switch_index_it]-t_optTheta[switch_index_it]);
          index_jr++;
        }
      }

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
      Rprintf("\nR Matrix row %d (%s): ", (i+1), (params[0].setP->semDone[i]) ? "    Done" : "Not done");
      for(j=0;j<len;j++) {
        Rprintf(" %5.2f ",Rmat[i][j]);
      }
    }
    Rprintf("\n\n");
  }
  Free(SuffSem);
  Free(params_sem);
}



/**
 * Read in the data set and population params
 * inputs:
 *   ndim: number of dimensions
 *   pdX: non-survey, non-homogenous data (length n_samp)
 *   sur_W: survey data (length s_samp)
 *   x1_W1: homogenous data (X==1) (length x1_samp)
 *   x0_W2: homogenous data (X==0) (length x0_samp)
 * mutates: params
 */
 void readData(Param* params, int n_dim, double* pdX, double* sur_W, double* x1_W1, double* x0_W2,
                int n_samp, int s_samp, int x1_samp, int x0_samp) {
     /* read the data set */
  int itemp,i,j,surv_dim;
  double dtemp;
  setParam* setP=params[0].setP;

  /** Packing Y, X  **/
  itemp = 0;
  for (j = 0; j < n_dim; j++)
    for (i = 0; i < n_samp; i++) {
      params[i].caseP.data[j] = pdX[itemp++];
    }

  for (i = 0; i < n_samp; i++) {
    params[i].caseP.dataType=DPT_General;
    params[i].caseP.X=params[i].caseP.data[0];
    params[i].caseP.Y=params[i].caseP.data[1];
    //fix X edge cases
    params[i].caseP.X=(params[i].caseP.X >= 1) ? .9999 : ((params[i].caseP.X <= 0) ? 0.0001 : params[i].caseP.X);
    //fix Y edge cases
    params[i].caseP.Y=(params[i].caseP.Y >= 1) ? .9999 : ((params[i].caseP.Y <= 0) ? 0.0001 : params[i].caseP.Y);
  }

  /*read the survey data */
  itemp=0;
  surv_dim=n_dim + (setP->ncar ? 1 : 0); //if NCAR, the survey data will include X's
  for (j=0; j<surv_dim; j++) {
    for (i=n_samp; i<n_samp+s_samp; i++) {
      dtemp=sur_W[itemp++];
      params[i].caseP.dataType=DPT_Survey;
      if (j<n_dim) {
        params[i].caseP.W[j]=(dtemp == 1) ? .9999 : ((dtemp==0) ? .0001 : dtemp);
        params[i].caseP.Wstar[j]=logit(params[i].caseP.W[j],"Survey read");
      }
      else { //if given the X (NCAR), we set the X and contruct Y
        params[i].caseP.X=(dtemp == 1) ? .9999 : ((dtemp==0) ? .0001 : dtemp);
        params[i].caseP.Y=params[i].caseP.X*params[i].caseP.W[0]+(1-params[i].caseP.X);
      }
    }
  }

  /*read homeogenous areas information */
  for (i=n_samp+s_samp; i<n_samp+s_samp+x1_samp; i++) {
    /*params[i].caseP.dataType=DPT_Homog_X1;
    params[i].caseP.W[0]=(x1_W1[i] == 1) ? .9999 : ((x1_W1[i]==0) ? .0001 : x1_W1[i]);
    params[i].caseP.Wstar[0]=logit(params[i].caseP.W[0],"X1 read");*/
  }

  for (i=n_samp+s_samp+x1_samp; i<n_samp+s_samp+x1_samp+x0_samp; i++) {
    /*params[i].caseP.dataType=DPT_Homog_X0;
    params[i].caseP.W[1]=(x0_W2[i] == 1) ? .9999 : ((x0_W2[i]==0) ? .0001 : x0_W2[i]);
    params[i].caseP.Wstar[1]=logit(params[i].caseP.W[1],"X0 read");*/
  }

  /*Current version of program does not handle homogenous data*/
  if ((x1_samp+x0_samp)>0) {
    printf("WARNING: Homogenous data is ignored and not handled by the current version of eco.");
  }


  if (setP->verbose>=2) {
    Rprintf("Y X\n");
    for(i=0;i<5;i++) Rprintf("%5d%14g%14g\n",i,params[i].caseP.Y,params[i].caseP.X);
    if (s_samp>0) {
      Rprintf("SURVEY data\nY X\n");
      int s_max=fmin2(n_samp+x1_samp+x0_samp+s_samp,n_samp+x1_samp+x0_samp+5);
      for(i=n_samp+x1_samp+x0_samp; i<s_max; i++) Rprintf("%5d%14g%14g\n",i,params[i].caseP.Y,params[i].caseP.X);
    }
  }

}

/**
 * During the main E-M loop, this function prints the output column headers
 * finalTheta: 1 if this is for the final theta -- include static variables
 **/
void printColumnHeader(int main_loop, int iteration_max, setParam* setP, int finalTheta) {
  int i;
  int param_len;
  param_len = setP->param_len;

  char temp[50]; int hlen;
  if (!finalTheta) hlen=sprintf(temp, "cycle %d/%d:",main_loop,iteration_max); //Length of cycle text
  else hlen=sprintf(temp, "Final Theta:");
  for (i=0;i<hlen;i++) Rprintf(" ");

  if (param_len<=5) { //CAR
    Rprintf("  mu_1  mu_2 sig_1 sig_2");
    if (!setP->fixedRho || finalTheta) Rprintf("  r_12");
  } else { //NCAR
    if (finalTheta) {
      Rprintf("  mu_3  mu_1  mu_2 sig_3 sig_1 sig_2  r_13  r_23  r_12");
    }
    else {
      Rprintf("  mu_1  mu_2 sig_1 sig_2  r_13  r_23  r_12");
    }
  }
  Rprintf("\n");
}

/**
 * Parameterizes the elements of theta
 * Input: pdTheta
 * Mutates: t_pdTheta
 */
void transformTheta(double* pdTheta, double* t_pdTheta, int len, setParam* setP) {
  if (len<=5) {
    t_pdTheta[0]=pdTheta[0];
    t_pdTheta[1]=pdTheta[1];
    t_pdTheta[2]=log(pdTheta[2]);
    t_pdTheta[3]=log(pdTheta[3]);
    t_pdTheta[4]=.5*(log(1+pdTheta[4])-log(1-pdTheta[4]));
  }
  else {
    t_pdTheta[0]=pdTheta[0];
    t_pdTheta[1]=pdTheta[1];
    t_pdTheta[2]=pdTheta[2];
    t_pdTheta[3]=log(pdTheta[3]);
    t_pdTheta[4]=log(pdTheta[4]);
    t_pdTheta[5]=log(pdTheta[5]);
    t_pdTheta[6]=.5*(log(1+pdTheta[6])-log(1-pdTheta[6]));
    t_pdTheta[7]=.5*(log(1+pdTheta[7])-log(1-pdTheta[7]));
    t_pdTheta[8]=.5*(log(1+pdTheta[8])-log(1-pdTheta[8]));
  }
}

/**
 * Un-parameterizes the elements of theta
 * Input: t_pdTheta
 * Mutates: pdTheta
 */
void untransformTheta(double* t_pdTheta,double* pdTheta, int len, setParam* setP) {
  if (len<=5) {
    pdTheta[0]=t_pdTheta[0];
    pdTheta[1]=t_pdTheta[1];
    pdTheta[2]=exp(t_pdTheta[2]);
    pdTheta[3]=exp(t_pdTheta[3]);
    pdTheta[4]=(exp(2*t_pdTheta[4])-1)/(exp(2*t_pdTheta[4])+1);
  }
  else {
    pdTheta[0]=t_pdTheta[0];
    pdTheta[1]=t_pdTheta[1];
    pdTheta[2]=t_pdTheta[2];
    pdTheta[3]=exp(t_pdTheta[3]);
    pdTheta[4]=exp(t_pdTheta[4]);
    pdTheta[5]=exp(t_pdTheta[5]);
    if (!setP->fixedRho) {
      pdTheta[6]=(exp(2*t_pdTheta[6])-1)/(exp(2*t_pdTheta[6])+1);
      pdTheta[7]=(exp(2*t_pdTheta[7])-1)/(exp(2*t_pdTheta[7])+1);
    }
    else {
      pdTheta[6]=t_pdTheta[6];
      pdTheta[7]=t_pdTheta[7];
    }
    pdTheta[8]=(exp(2*t_pdTheta[8])-1)/(exp(2*t_pdTheta[8])+1);
  }
}

/**
 * untransforms theta under ncar -- fixed rho
 * input reference:  (0) mu_3, (1) mu_1, (2) mu_2, (3) sig_3, (4) sig_1 | 3, (5) sig_2 | 3, (6) beta1, (7) beta2, (8) r_12 | 3
 * output reference: (0) mu_3, (1) mu_1, (2) mu_2, (3) sig_3, (4) sig_1, (5) sig_2, (6) r_13, (7) r_23, (8) r_12
 * mutates: pdTheta
 **/
void ncarFixedRhoUnTransform(double* pdTheta) {
  double* tmp=doubleArray(9);
  int i;
  for (i=0;i<9;i++) tmp[i]=pdTheta[i];
  pdTheta[0]=tmp[0];
  pdTheta[1]=tmp[1];
  pdTheta[2]=tmp[2];
  pdTheta[3]=tmp[3];
  pdTheta[4]=tmp[4] + tmp[6]*tmp[6]*tmp[3];
  pdTheta[5]=tmp[5] + tmp[7]*tmp[7]*tmp[3];
  pdTheta[6]=(tmp[6]*sqrt(tmp[3]))/(sqrt(pdTheta[4]));
  pdTheta[7]=(tmp[7]*sqrt(tmp[3]))/(sqrt(pdTheta[5]));
  pdTheta[8]=(tmp[8]*sqrt(tmp[4]*tmp[5]) + tmp[6]*tmp[7]*tmp[3])/(sqrt(pdTheta[4]*pdTheta[5]));
  Free(tmp);
}

/**
 * transforms theta under ncar -- fixed rho
 * input reference:  (0) mu_3, (1) mu_1, (2) mu_2, (3) sig_3, (4) sig_1, (5) sig_2, (6) r_13, (7) r_23, (8) r_12
 * output reference: (0) mu_3, (1) mu_1, (2) mu_2, (3) sig_3, (4) sig_1 | 3, (5) sig_2 | 3, (6) beta1, (7) beta2, (8) r_12 | 3
 * mutates: pdTheta
 **/
void ncarFixedRhoTransform(double* pdTheta) {
  double* tmp=doubleArray(9);
  int i;
  for (i=0;i<9;i++) tmp[i]=pdTheta[i];
  pdTheta[0]=tmp[0];
  pdTheta[1]=tmp[1];
  pdTheta[2]=tmp[2];
  pdTheta[3]=tmp[3];
  pdTheta[4]=tmp[4] - tmp[6]*tmp[6]*tmp[4];
  pdTheta[5]=tmp[5] - tmp[7]*tmp[7]*tmp[5];
  pdTheta[6]=tmp[6]*sqrt(tmp[4]/tmp[3]);
  pdTheta[7]=tmp[7]*sqrt(tmp[5]/tmp[3]);
  pdTheta[8]=(tmp[8] - tmp[6]*tmp[7])/(sqrt((1 - tmp[6]*tmp[6])*(1 - tmp[7]*tmp[7])));
  Free(tmp);
}


/**
 * Input transformed theta, loglikelihood, iteration
 * Mutates: history_full
 **/
void setHistory(double* t_pdTheta, double loglik, int iter,setParam* setP,double history_full[][10]) {
  int len=setP->param_len;
  int j;
  for(j=0;j<len;j++)
    history_full[iter][j]=t_pdTheta[j];
  if (iter>0) history_full[iter-1][len]=loglik;
}

/**
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

/**
 * Is the SEM process completely done.
 **/
int semDoneCheck(setParam* setP) {
  int varlen=0; int j;
  for(j=0; j<setP->param_len;j++)
    if(setP->varParam[j]) varlen++;
  for(j=0;j<varlen;j++)
    if(setP->semDone[j]==0) return 0;
  return 1;
}

/**
 * Older function. No longer used.
 **/
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
