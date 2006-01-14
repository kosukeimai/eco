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
                int n_samp, int s_samp, int x1_samp, int x0_samp, int verbose);
void ecoEStep(Param* params, int n_samp, int s_samp, int x1_samp, int x0_samp, double* suff, int verbose);
void ecoMStep(double* Suff, double* pdTheta, double** Sigma, double** InvSigma, Param* params, int t_samp, int verbose);
void ecoMStepNCAR(double* Suff, double* pdTheta, double** Sigma, double** InvSigma, Param* params, int t_samp, int verbose);
int closeEnough(double* pdTheta, double* pdTheta_old, double maxerr);
void gridEStep(Param* params, int n_samp, int s_samp, int x1_samp, int x0_samp, double* suff, int verbose, double minW1, double maxW1);

void cEMeco(
	    /*data input */
	    double *pdX,         /* data (X, Y) */
	    double *pdTheta_in,  /* Theta^ t
				    mu1, mu2, var1, var2, rho */
	    int *pin_samp,       /* sample size */

	    /* loop vairables */
	    int *iteration_max,          /* number of maximum interations */
	    int *convergence,          /* abs value limit before stopping */

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

	    /* storage */
	    double *pdTheta,  /*EM result for Theta^(t+1) */
	    double *Suff      /*out put suffucient statistics (E(W_1i|Y_i),
				E(E_1i*W_1i|Y_i..) when  conveges */
	    ){

  int n_samp  = *pin_samp;    /* sample size */
  int s_samp  = *survey ? *sur_samp : 0;     /* sample size of survey data */
  int x1_samp = *x1 ? *sampx1 : 0;       /* sample size for X=1 */
  int x0_samp = *x0 ? *sampx0 : 0;       /* sample size for X=0 */
  int t_samp=n_samp+s_samp+x1_samp+x0_samp;  /* total sample size*/
  int n_dim=2;        /* dimensions */
  int verbose=1;         /* one to print cycles, two to print the data */

  /* model parameters */
  double **Sigma=doubleMatrix(n_dim,n_dim);/* inverse covariance matrix*/
  double **InvSigma=doubleMatrix(n_dim,n_dim);/* inverse covariance matrix*/

  double *pdTheta_old=doubleArray(5);

  /* misc variables */
  int i, main_loop, start;   /* used for various loops */
  double rho;

  /* get random seed */
  GetRNGstate();

  //assign param
  Param* params=(Param*) R_alloc(t_samp,sizeof(Param));
  readData(params, n_dim, pdX, sur_W, x1_W1, x0_W2, n_samp, s_samp, x1_samp, x0_samp, verbose);

/***Begin main loop ***/
main_loop=1;start=1;
while (main_loop<=*iteration_max && (start==1 || !closeEnough(pdTheta,pdTheta_old,*convergence))) {

  if (start) {
    for(i=0;i<5;i++) pdTheta[i]=pdTheta_in[i];
    start=0;
  }
  //keep the old theta around for comaprison
  for(i=0;i<5;i++) pdTheta_old[i]=pdTheta[i];


  rho=Sigma[1][0]/sqrt(Sigma[0][0]*Sigma[1][1]);
  Rprintf("cycle %d/%d: %5g %5g %5g %5g %5g rho: %5g\n",main_loop,*iteration_max,pdTheta[0],pdTheta[1],Sigma[0][0],Sigma[1][1],Sigma[1][0],rho);

  ecoEStep(params, n_samp, s_samp, x1_samp, x0_samp, Suff,verbose);
  if (1==1)
    ecoMStep(Suff,pdTheta,Sigma,InvSigma,params,t_samp,verbose);
  else
    ecoMStepNCAR(Suff,pdTheta,Sigma,InvSigma,params,t_samp,verbose);
  //char ch;
  //scanf(" %c", &ch );


  if (verbose>=2) {
    Rprintf("theta and suff\n");
    Rprintf("%10g%10g%10g%10g%10g (%10g)\n",pdTheta[0],pdTheta[1],pdTheta[2],pdTheta[3],pdTheta[4],pdTheta[4]*sqrt(pdTheta[2]*pdTheta[3]));
    Rprintf("%10g%10g%10g%10g%10g\n",Suff[0],Suff[1],Suff[2],Suff[3],Suff[4]);
  }
  main_loop++;
}

/***End main loop ***/

Rprintf("End loop PDT:%5g %5g %5g %5g %5g\n",pdTheta[0],pdTheta[1],pdTheta[2],pdTheta[3],pdTheta[4]);

/* write out the random seed */
PutRNGstate();

/* Freeing the memory */
Free(pdTheta_old);
FreeMatrix(Sigma,n_dim);
FreeMatrix(InvSigma,n_dim);
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

 */

void ecoEStep(Param* params, int n_samp, int s_samp, int x1_samp, int x0_samp, double* suff, int verbose) {

int t_samp,i,j,temp0,temp1, loglik;
double testw1, testw2,testll;
Param param;

loglik=0;

t_samp=n_samp+x1_samp+x0_samp+s_samp;

  double **Wstar=doubleMatrix(t_samp,5);     /* pseudo data(transformed)*/
testll=0;
  for (i = 0; i<n_samp; i++) {
    param = params[i];
    if (param.Y>=.9999 || param.Y<=.0001) { //if Y is near the edge, then W1 and W2 are very constrained
      Wstar[i][0]=param.Y;
      Wstar[i][1]=param.Y;
      Wstar[i][2]=Wstar[i][0]*Wstar[i][0];
      Wstar[i][3]=Wstar[i][0]*Wstar[i][1];
      Wstar[i][4]=Wstar[i][1]*Wstar[i][1];
    }
    else {
      setBounds((Param*)&param);
      setNormConst((Param*)&param);

      for (j=0;j<5;j++) {
        param.suff=j;
        Wstar[i][j]=paramIntegration(&SuffExp,(void *)&param);
      }
      param.suff=5;
      testw1=paramIntegration(&SuffExp,(void *)&param);
      param.suff=6;
      testw2=paramIntegration(&SuffExp,(void *)&param);
      if (loglik) {
        param.suff=7;
        testll+=paramIntegration(&SuffExp,(void *)&param);
      }

   //report error E1 if E[W1],E[W2] is not on the tomography line
  if (fabs(testw1-getW1FromW2(param.X, param.Y,testw2))>0.01)
    Rprintf("E1 %d %5g %5g %5g %5g %5g %5g %5g %5g \n", i, param.X, param.Y, param.normcT, Wstar[i][0],Wstar[i][1],Wstar[i][2],Wstar[i][3],Wstar[i][4]);
  //report error E2 if Jensen's inequality doesn't hold
  if (Wstar[i][4]<pow(Wstar[i][1],2) || Wstar[i][2]<pow(Wstar[i][0],2))
     Rprintf("E2 %d %5g %5g %5g %5g %5g %5g %5g %5g\n", i, param.X, param.Y, param.normcT, Wstar[i][0],Wstar[i][1],Wstar[i][2],Wstar[i][3],Wstar[i][4]);
  //used for debugging if necessary
  if (verbose>=2 && i<20)
     Rprintf("%d %4g %4g %4g %4g %4g %4g %4g %4g %4g %4g\n", i, param.X, param.Y, param.normcT, Wstar[i][0],Wstar[i][1],testw1,testw2,Wstar[i][2],Wstar[i][3],Wstar[i][4]);
    }
  }

    /* analytically compute E{W2_i|Y_i} given W1_i, mu and Sigma in x1 homeogeneous areas */
    for (i=n_samp; i<n_samp+x1_samp; i++) {
      temp0=log(params[i].W[0])-log(1-params[i].W[0]);
      temp1=params[i].mu[1]+params[i].Sigma[0][1]/params[i].Sigma[0][0]*(temp0-params[i].mu[0]);
      Wstar[i][0]=temp0;
      Wstar[i][1]=temp1;
      Wstar[i][2]=temp0*temp0;
      Wstar[i][3]=temp0*temp1;
      Wstar[i][4]=temp1*temp1;
    }

  /*analytically compute E{W1_i|Y_i} given W2_i, mu and Sigma in x0 homeogeneous areas */
    for (i=n_samp+x1_samp; i<n_samp+x1_samp+x0_samp; i++) {
      temp1=log(params[i].W[1])-log(1-params[i].W[1]);
      temp0=params[i].mu[0]+params[i].Sigma[0][1]/params[i].Sigma[1][1]*(temp1-params[i].mu[1]);
      Wstar[i][0]=temp0;
      Wstar[i][1]=temp1;
      Wstar[i][2]=temp0*temp0;
      Wstar[i][3]=temp0*temp1;
      Wstar[i][4]=temp1*temp1;
    }

    /* Use the values given by the survey data */
    for (i=n_samp+x1_samp+x0_samp; i<n_samp+x1_samp+x0_samp+s_samp; i++) {
      Wstar[i][0]=log(params[i].W[0])-log(1-params[i].W[0]);
      Wstar[i][1]=log(params[i].W[1])-log(1-params[i].W[1]);
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

    if (loglik)
      if(verbose>=1) Rprintf("Log liklihood %15g",testll);

FreeMatrix(Wstar,t_samp);

}

//Standard M-Step
void ecoMStep(double* Suff, double* pdTheta, double** Sigma, double** InvSigma, Param* params, int t_samp, int verbose) {

int i,j,k;

  pdTheta[0]=Suff[0];  /*mu1*/
  pdTheta[1]=Suff[1];  /*mu2*/
  pdTheta[2]=Suff[2]-2*Suff[0]*pdTheta[0]+pdTheta[0]*pdTheta[0];  /*sigma11*/
  pdTheta[3]=Suff[3]-2*Suff[1]*pdTheta[1]+pdTheta[1]*pdTheta[1];  /*sigma22*/
  pdTheta[4]=Suff[4]-Suff[0]*pdTheta[1]-Suff[1]*pdTheta[0]+pdTheta[0]*pdTheta[1]; /*sigma12*/
  pdTheta[4]=pdTheta[4]/sqrt(pdTheta[2]*pdTheta[3]); /*rho*/


    //set Sigma
  Sigma[0][0] = pdTheta[2];
  Sigma[1][1] = pdTheta[3];
  Sigma[0][1] = pdTheta[4]*sqrt(pdTheta[2]*pdTheta[3]);
  Sigma[1][0] = Sigma[0][1];

  dinv(Sigma, 2, InvSigma);

  /* assign each data set the new mu and sigma */
  for(i=0;i<t_samp;i++) {
    params[i].mu[0]=pdTheta[0];
    params[i].mu[1]=pdTheta[1];
    for(j=0; j<2;j++)
      for(k=0; k<2;k++) {
        params[i].Sigma[j][k] = Sigma[j][k];
        params[i].InvSigma[j][k] = InvSigma[j][k];
      }
  }

}


//M-Step under NCAR
//NEEDS WORK
void ecoMStepNCAR(double* Suff, double* pdTheta,  double** Sigma, double** InvSigma, Param* params, int t_samp, int verbose) {

  double **Sigma3=doubleMatrix(3,3);   /* covariance matrix*/
  int i,j,k;

  //find mean of X
  double mu3=0;
  double mu3sq=0;
  for(i=0;i<t_samp;i++) {
    mu3 += logit(params[i].X);
    mu3sq += pow(logit(params[i].X),2);
  }
  mu3 = mu3/t_samp; mu3sq = mu3sq/t_samp;


  pdTheta[0]=Suff[0];  /*mu1*/
  pdTheta[1]=Suff[1];  /*mu2*/

  //set Sigma3 (3x3 sigma) ..probably needs fixing by using prior mus
  Sigma3[0][0] = Suff[2]-2*Suff[0]*pdTheta[0]+pdTheta[0]*pdTheta[0];
  Sigma3[0][1] = Suff[4]-Suff[0]*pdTheta[1]-Suff[1]*pdTheta[0]+pdTheta[0]*pdTheta[1];
  Sigma3[0][2] = 0; //What to do here?
  Sigma3[1][0] = Sigma[0][1];
  Sigma3[1][1] = Suff[3]-2*Suff[1]*pdTheta[1]+pdTheta[1]*pdTheta[1];
  Sigma3[1][2] = 0;
  Sigma3[2][0] = 0;
  Sigma3[2][1] = 0;
  Sigma3[2][2] = mu3sq-mu3*mu3;

  Sigma[0][0]= Sigma3[0][0] - Sigma3[0][2]*Sigma3[0][2]/Sigma3[2][2];
  Sigma[0][1]= Sigma3[0][1] - Sigma3[0][2]*Sigma3[1][2]/Sigma3[2][2];
  Sigma[1][0]= Sigma3[1][0] - Sigma3[0][2]*Sigma3[1][2]/Sigma3[2][2];
  Sigma[1][1]= Sigma3[1][1] - Sigma3[1][2]*Sigma3[1][2]/Sigma3[2][2];

  dinv(Sigma, 2, InvSigma);

  pdTheta[2]=Sigma[0][0];
  pdTheta[3]=Sigma[1][1];
  pdTheta[4]=Sigma[1][0]/(sqrt(Sigma[0][0]*Sigma[1][1]));

  /* assign each data point the new mu and sigma */
  for(i=0;i<t_samp;i++) {
    params[i].mu[0]=pdTheta[0] + (Sigma3[0][2]/Sigma3[2][2])*(params[i].X-mu3);
    params[i].mu[1]=pdTheta[1] + (Sigma3[1][2]/Sigma3[2][2])*(params[i].X-mu3);
    for(j=0; j<2;j++)
      for(k=0; k<2;k++) {
        params[i].Sigma[j][k] = Sigma[j][k];
        params[i].InvSigma[j][k] = InvSigma[j][k];
      }
  }

  FreeMatrix(Sigma3,3);

}

/*
 * Read in the data set and population params
 */
 void readData(Param* params, int n_dim, double* pdX, double* sur_W, double* x1_W1, double* x0_W2,
                int n_samp, int s_samp, int x1_samp, int x0_samp, int verbose) {
     /* read the data set */
int itemp,i,j;
double dtemp;
  /** Packing Y, X  **/
  itemp = 0;
  for (j = 0; j < n_dim; j++)
    for (i = 0; i < n_samp; i++) {
      params[i].data[j] = pdX[itemp++];
    }

  for (i = 0; i < n_samp; i++) {
    params[i].X=params[i].data[0];
    params[i].Y=params[i].data[1];
    //fix X edge cases
    params[i].X=(params[i].X == 1) ? .9999 : ((params[i].X == 0) ? 0.0001 : params[i].X);
  }

  if (verbose>=2) {
    printf("Y X\n");
    for(i=0;i<n_samp;i++)
      Rprintf("%5d%14g%14g\n",i,params[i].Y,params[i].X);
      }

  /*read homeogenous areas information */
    for (i=n_samp; i<n_samp+x1_samp; i++)
      params[i].W[0]=(x1_W1[i] == 1) ? .9999 : ((x1_W1[i]==0) ? .0001 : x1_W1[i]);

    for (i=0; i<x0_samp; i++)
      params[i].W[1]=(x0_W2[i] == 1) ? .9999 : ((x0_W2[i]==0) ? .0001 : x0_W2[i]);


  /*read the survey data */
     itemp=0;
    for (j=0; j<n_dim; j++)
      for (i=n_samp+x1_samp+x0_samp; i<n_samp+x1_samp+x0_samp+s_samp; i++) {
        dtemp=sur_W[itemp++];
        params[i].W[j]=(dtemp == 1) ? .9999 : ((dtemp==0) ? .0001 : dtemp);
      }
 }

/*
 * Determines whether we have converged
 * Takes in the current and old (one step previous) array of theta values
 * maxerr is the maximum difference two corresponding values can have before the
 *  function returns false
 */
int closeEnough(double* pdTheta, double* pdTheta_old, double maxerr) {
  int j;
  for(j = 0; j<5; j++)
    if (fabs(pdTheta[j]-pdTheta_old[j])>=maxerr) return 0;
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
      X[i][j]=params[i].data[j];

  GridPrep(W1g, W2g, (double**) params[i].data, (double*)&maxW1, (double*)&minW1, n_grid, n_samp, n_step);

    for (i=0; i<n_step; i++) {
    mflag[i]=0;
  }


  //update W, Wstar given mu, Sigma in regular areas
  for (i=0;i<n_samp;i++){
    if ( params[i].Y!=0 && params[i].Y!=1 ) {
      // project BVN(mu, Sigma) on the inth tomo line
      dtemp=0;
      for (j=0;j<n_grid[i];j++){
        vtemp[0]=log(W1g[i][j])-log(1-W1g[i][j]);
        vtemp[1]=log(W2g[i][j])-log(1-W2g[i][j]);
        prob_grid[j]=dMVN(vtemp, params[i].mu, (double**)(params[i].InvSigma), 2, 1) -
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
