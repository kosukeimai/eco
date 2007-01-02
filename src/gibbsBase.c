
#include <string.h>
#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "bayes.h"
#include "sample.h"

/* Normal Parametric Model for 2x2 Tables */
void cBaseeco(
	      /*data input */
	      double *pdX,     /* data (X, Y) */
	      int *pin_samp,   /* sample size */

	      /*MCMC draws */
	      int *n_gen,      /* number of gibbs draws */
	      int *burn_in,    /* number of draws to be burned in */
	      int *pinth,      /* keep every nth draw */
	      int *verbose,    /* 1 for output monitoring */

	      /* prior specification*/
	      int *pinu0,      /* prior df parameter for InvWish */
	      double *pdtau0,  /* prior scale parameter for Sigma */
	      double *mu0,     /* prior mean for mu */
	      double *pdS0,    /* prior scale for Sigma */
	      double *mustart, /* starting values for mu */
	      double *Sigmastart, /* starting values for Sigma */

	      /* incorporating survey data */
	      int *survey,     /*1 if survey data available (set of W_1, W_2)
				 0 not*/
	      int *sur_samp,   /*sample size of survey data*/
	      double *sur_W,   /*set of known W_1, W_2 */ 
				  
	      /* incorporating homeogenous areas */
	      int *x1,         /* 1 if X=1 type areas available 
				  W_1 known, W_2 unknown */
	      int *sampx1,     /* number X=1 type areas */
	      double *x1_W1,   /* values of W_1 for X1 type areas */
	      int *x0,         /* 1 if X=0 type areas available 
				  W_2 known, W_1 unknown */
	      int *sampx0,     /* number X=0 type areas */
	      double *x0_W2,   /* values of W_2 for X0 type areas */

	      /* bounds of W1 */
	      double *minW1, double *maxW1,

	      /* flags */
	      int *parameter,  /* 1 if save population parameter */
	      int *Grid,       /* 1 if Grid algorithm is used; 0 for
				  Metropolis */

	      /* storage for Gibbs draws of mu/sigmat*/
	      double *pdSMu0, double *pdSMu1, 
	      double *pdSSig00, double *pdSSig01, double *pdSSig11,
           
	      /* storage for Gibbs draws of W*/
	      double *pdSW1, double *pdSW2
	      ){	   
  
  /* some integers */
  int n_samp = *pin_samp;    /* sample size */
  int s_samp = *sur_samp;    /* sample size of survey data */ 
  int x1_samp = *sampx1;     /* sample size for X=1 */
  int x0_samp = *sampx0;     /* sample size for X=0 */
  int t_samp = n_samp+s_samp+x1_samp+x0_samp;  /* total sample size */
  int nth = *pinth;  
  int n_dim = 2;             /* dimension */
  int n_step = 1000;         /* 1/The default size of grid step */  

  /* prior parameters */ 
  double tau0 = *pdtau0;                          /* prior scale */
  int nu0 = *pinu0;                               /* prior degrees of freedom */   
  double **S0 = doubleMatrix(n_dim, n_dim);       /* The prior S parameter for InvWish */

  /* data */
  double **X = doubleMatrix(n_samp, n_dim);       /* The Y and covariates */
  double **W = doubleMatrix(t_samp, n_dim);       /* The W1 and W2 matrix */
  double **Wstar = doubleMatrix(t_samp, n_dim);   /* logit tranformed W */       
  double **S_W = doubleMatrix(s_samp, n_dim);     /* The known W1 and W2 matrix*/
  double **S_Wstar = doubleMatrix(s_samp, n_dim); /* logit transformed S_W*/

  /* grids */
  double **W1g = doubleMatrix(n_samp, n_step);    /* grids for W1 */
  double **W2g = doubleMatrix(n_samp, n_step);    /* grids for W2 */
  int *n_grid = intArray(n_samp);                 /* grid size */

  /* model parameters */
  double *mu = doubleArray(n_dim);                /* The mean */
  double **Sigma = doubleMatrix(n_dim, n_dim);    /* The covariance matrix */
  double **InvSigma = doubleMatrix(n_dim, n_dim); /* The inverse covariance matrix */

  /* misc variables */
  int i, j, k, main_loop;   /* used for various loops */
  int itemp, itempS, itempC, itempA;
  int progress = 1, itempP = ftrunc((double) *n_gen/10);
  double dtemp, dtemp1;

  /* get random seed */
  GetRNGstate();
  
  /* read the priors */
  itemp=0;
  for(k=0;k<n_dim;k++)
    for(j=0;j<n_dim;j++) S0[j][k]=pdS0[itemp++];

  /* read the data */
  itemp = 0;
  for (j = 0; j < n_dim; j++) 
    for (i = 0; i < n_samp; i++) 
      X[i][j] = pdX[itemp++];

  /* Initialize W, Wstar for n_samp */
  for (i=0; i< n_samp; i++) {
    if (X[i][1]!=0 && X[i][1]!=1) {
      W[i][0]=runif(minW1[i], maxW1[i]);
      W[i][1]=(X[i][1]-X[i][0]*W[i][0])/(1-X[i][0]);
    }
    if (X[i][1]==0) 
      for (j=0; j<n_dim; j++) W[i][j]=0.0001;
    if (X[i][1]==1) 
      for (j=0; j<n_dim; j++) W[i][j]=0.9999;
    for (j=0; j<n_dim; j++)
      Wstar[i][j]=log(W[i][j])-log(1-W[i][j]);
  }

  /* read homeogenous areas information */
  if (*x1==1) 
    for (i=0; i<x1_samp; i++) {
      W[(n_samp+i)][0]=x1_W1[i];
      if (W[(n_samp+i)][0]==0) W[(n_samp+i)][0]=0.0001;
      if (W[(n_samp+i)][0]==1) W[(n_samp+i)][0]=0.9999;
      Wstar[(n_samp+i)][0]=log(W[(n_samp+i)][0])-log(1-W[(n_samp+i)][0]);
    }
  if (*x0==1) 
    for (i=0; i<x0_samp; i++) {
      W[(n_samp+x1_samp+i)][1]=x0_W2[i];
      if (W[(n_samp+x1_samp+i)][1]==0) W[(n_samp+x1_samp+i)][1]=0.0001;
      if (W[(n_samp+x1_samp+i)][1]==1) W[(n_samp+x1_samp+i)][1]=0.9999;
      Wstar[(n_samp+x1_samp+i)][1]=log(W[(n_samp+x1_samp+i)][1])-log(1-W[(n_samp+x1_samp+i)][1]);
    }

  /* read the survey data */
  if (*survey==1) {
    itemp = 0;
    for (j=0; j<n_dim; j++)
      for (i=0; i<s_samp; i++) {
	S_W[i][j]=sur_W[itemp++];
	if (S_W[i][j]==0) S_W[i][j]=0.0001;
	if (S_W[i][j]==1) S_W[i][j]=0.9999;
	S_Wstar[i][j]=log(S_W[i][j])-log(1-S_W[i][j]);
	W[(n_samp+x1_samp+x0_samp+i)][j]=S_W[i][j];
	Wstar[(n_samp+x1_samp+x0_samp+i)][j]=S_Wstar[i][j];
      }
  }

  /* counters */
  itempA=0; /* for alpha */
  itempS=0; /* for storage */
  itempC=0; /* control nth draw */

  /*** calculate grids ***/
  if (*Grid) 
    GridPrep(W1g, W2g, X, maxW1, minW1, n_grid, n_samp, n_step);
    
  /* starting vales of mu and Sigma */
  itemp = 0;
  for(j=0;j<n_dim;j++){
    mu[j] = mustart[j];
    for(k=0;k<n_dim;k++)
      Sigma[j][k]=Sigmastart[itemp++];
  }
  dinv(Sigma, n_dim, InvSigma);
  
  /*** Gibbs sampler! ***/
  if (*verbose)
    Rprintf("Starting Gibbs Sampler...\n");
  for(main_loop=0; main_loop<*n_gen; main_loop++){
    /** update W, Wstar given mu, Sigma in regular areas **/
    for (i=0;i<n_samp;i++){
      if ( X[i][1]!=0 && X[i][1]!=1 ) {
	if (*Grid)
	  rGrid(W[i], W1g[i], W2g[i], n_grid[i], mu, InvSigma, n_dim);
	else 
	  rMH(W[i], X[i], minW1[i], maxW1[i], mu, InvSigma, n_dim);
      } 
      /*3 compute Wsta_i from W_i*/
      Wstar[i][0]=log(W[i][0])-log(1-W[i][0]);
      Wstar[i][1]=log(W[i][1])-log(1-W[i][1]);
    }
    
    /* update W2 given W1, mu and Sigma in x1 homeogeneous areas */
    if (*x1==1)
      for (i=0; i<x1_samp; i++) {
	dtemp=mu[1]+Sigma[0][1]/Sigma[0][0]*(Wstar[n_samp+i][0]-mu[0]);
	dtemp1=Sigma[1][1]*(1-Sigma[0][1]*Sigma[0][1]/(Sigma[0][0]*Sigma[1][1]));
	dtemp1=sqrt(dtemp1);
	Wstar[n_samp+i][1]=rnorm(dtemp, dtemp1);
	W[n_samp+i][1]=exp(Wstar[n_samp+i][1])/(1+exp(Wstar[n_samp+i][1]));
      }
    
    /* update W1 given W2, mu and Sigma in x0 homeogeneous areas */
    if (*x0==1)
      for (i=0; i<x0_samp; i++) {
	dtemp=mu[0]+Sigma[0][1]/Sigma[1][1]*(Wstar[n_samp+x1_samp+i][1]-mu[1]);
	dtemp1=Sigma[0][0]*(1-Sigma[0][1]*Sigma[0][1]/(Sigma[0][0]*Sigma[1][1]));
	dtemp1=sqrt(dtemp1);
	Wstar[n_samp+x1_samp+i][0]=rnorm(dtemp, dtemp1);
	W[n_samp+x1_samp+i][0]=exp(Wstar[n_samp+x1_samp+i][0])/(1+exp(Wstar[n_samp+x1_samp+i][0]));
      }
    
    /* update mu, Sigma given wstar using effective sample of Wstar */
    NIWupdate(Wstar, mu, Sigma, InvSigma, mu0, tau0, nu0, S0, t_samp, n_dim);
    
    /*store Gibbs draw after burn-in and every nth draws */      
    if (main_loop>=*burn_in){
      itempC++;
      if (itempC==nth){
	pdSMu0[itempA]=mu[0];
	pdSMu1[itempA]=mu[1];
	pdSSig00[itempA]=Sigma[0][0];
	pdSSig01[itempA]=Sigma[0][1];
	pdSSig11[itempA]=Sigma[1][1];
	itempA++;
	for(i=0; i<(n_samp+x1_samp+x0_samp); i++){
	  pdSW1[itempS]=W[i][0];
	  pdSW2[itempS]=W[i][1];
	  itempS++;
	}
	itempC=0;
      }
    } 
    if (*verbose)
      if (itempP == main_loop) {
	Rprintf("%3d percent done.\n", progress*10);
	itempP+=ftrunc((double) *n_gen/10); progress++;
	R_FlushConsole();
      }
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */ 

  if(*verbose)
    Rprintf("100 percent done.\n");

  /** write out the random seed **/
  PutRNGstate();

  /* Freeing the memory */
  FreeMatrix(X, n_samp);
  FreeMatrix(W, t_samp);
  FreeMatrix(Wstar, t_samp);
  FreeMatrix(S_W, s_samp);
  FreeMatrix(S_Wstar, s_samp);
  FreeMatrix(S0, n_dim);
  FreeMatrix(W1g, n_samp);
  FreeMatrix(W2g, n_samp);
  free(n_grid);
  free(mu);
  FreeMatrix(Sigma,n_dim);
  FreeMatrix(InvSigma, n_dim);
  
} /* main */

