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
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "bayes.h"
#include "sample.h"

/* Normal Parametric Model for 2xC (with C > 2) Tables */
void cBase2C(
	     /*data input */
	     double *pdX,     /* X: matrix */
	     double *Y,       /* Y: vector */
	     double *pdWmin,  /* lower bounds */
	     double *pdWmax,  /* uppwer bounds */
	     int *pin_samp,   /* sample size */
	     int *pin_dim,    /* number of columns */
	     
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
	     
	     /* storage */
	     int *parameter,  /* 1 if save population parameter */
	     double *pdSmu, 
	     double *pdSSigma,
	     double *pdSW
	     ){	   
  
  /* some integers */
  int n_samp = *pin_samp;    /* sample size */
  int nth = *pinth;          /* keep every pth draw */
  int n_dim = *pin_dim;      /* dimension */

  /* prior parameters */ 
  double tau0 = *pdtau0;                          /* prior scale */
  int nu0 = *pinu0;                               /* prior degrees of freedom */   
  double **S0 = doubleMatrix(n_dim, n_dim);       /* The prior S parameter for InvWish */

  /* data */
  double **X = doubleMatrix(n_samp, n_dim);       /* X */
  double **W = doubleMatrix(n_samp, n_dim);       /* The W1 and W2 matrix */
  double **Wstar = doubleMatrix(n_samp, n_dim);   /* logit tranformed W */       

  /* The lower and upper bounds of Z = W*X/Y **/
  double **minZ = doubleMatrix(n_samp, n_dim);
  double **maxZ = doubleMatrix(n_samp, n_dim);    

  /* model parameters */
  double *mu = doubleArray(n_dim);                /* The mean */
  double **Sigma = doubleMatrix(n_dim, n_dim);    /* The covariance matrix */
  double **InvSigma = doubleMatrix(n_dim, n_dim); /* The inverse covariance matrix */

  /* misc variables */
  int i, j, k, main_loop;   /* used for various loops */
  int itemp;
  int itempM = 0; /* for mu */
  int itempS = 0; /* for Sigma */
  int itempW = 0; /* for W */
  int itempC = 0; /* control nth draw */
  int progress = 1, itempP = ftrunc((double) *n_gen/10);
  double dtemp, dtemp1;
  double *param = doubleArray(n_dim);   /* Dirichlet parameters */
  double *dvtemp = doubleArray(n_dim);

  /* get random seed */
  GetRNGstate();
  
  /* read X */
  itemp = 0;
  for (j = 0; j < n_dim; j++) 
    for (i = 0; i < n_samp; i++) 
      X[i][j] = pdX[itemp++];

  /* compute bounds on Z */
  itemp = 0;
  for (j = 0; j < n_dim; j++) 
    for (i = 0; i < n_samp; i++) 
      minZ[i][j] = fmax2(0, pdWmin[itemp++]*X[i][j]/Y[i]);
  itemp = 0;
  for (j = 0; j < n_dim; j++) 
    for (i = 0; i < n_samp; i++) 
      maxZ[i][j] = fmin2(1, pdWmax[itemp++]*X[i][j]/Y[i]);

  /* initial values for W */
  for (j = 0; j < n_dim; j++)
    param[j] = 1;
  for (i = 0; i < n_samp; i++) {
    k = 0; itemp = 1;
    while (itemp > 0) {
    rDirich(dvtemp, param, n_dim);
    itemp = 0;
    for (j = 0; j < n_dim; j++)
      if (dvtemp[j] > maxZ[i][j] || dvtemp[j] < minZ[i][j])
        itemp++;
    k++;
    if (k > 100000)
      error("gibbs sampler cannot start because bounds are too short.");
    }
    for (j = 0; j < n_dim; j++) {
      W[i][j] = dvtemp[j]*Y[i]/X[i][j];
      Wstar[i][j] = log(W[i][j])-log(1-W[i][j]);
    }
  }

  /* read the prior */
  itemp = 0;
  for(k = 0; k < n_dim; k++)
    for(j = 0; j < n_dim; j++) 
      S0[j][k] = pdS0[itemp++];

  /* initialize vales of mu and Sigma */
  for(j = 0; j < n_dim; j++){
    mu[j] = mu0[j];
    for(k = 0; k < n_dim; k++)
      Sigma[j][k] = S0[j][k];
  }
  dinv(Sigma, n_dim, InvSigma);
  
  /*** Gibbs sampler! ***/
  for(main_loop = 0; main_loop < *n_gen; main_loop++){
    /** update W, Wstar given mu, Sigma **/
    for (i = 0; i < n_samp; i++){
      rMHrc(W[i], X[i], Y[i], minZ[i], maxZ[i], mu, InvSigma, n_dim);
      for (j = 0; j < n_dim; j++) 
	Wstar[i][j] = log(W[i][j])-log(1-W[i][j]);
    }
    
    /* update mu, Sigma given wstar using effective sample of Wstar */
    NIWupdate(Wstar, mu, Sigma, InvSigma, mu0, tau0, nu0, S0, n_samp, n_dim);
    
    /*store Gibbs draw after burn-in and every nth draws */      
    if (main_loop>=*burn_in){
      itempC++;
      if (itempC==nth){
	for (j = 0; j < n_dim; j++) {
	  pdSmu[itempM++]=mu[j];
	  for (k = 0; k < n_dim; k++)
	    if (j <=k)
	      pdSSigma[itempS++]=Sigma[j][k];
	}
	for(i = 0; i < n_samp; i++)
	  for (j = 0; j < n_dim; j++)
	    pdSW[itempW++] = W[i][j];
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
  FreeMatrix(W, n_samp);
  FreeMatrix(Wstar, n_samp);
  FreeMatrix(minZ, n_samp);
  FreeMatrix(maxZ, n_samp);
  FreeMatrix(S0, n_dim);
  free(mu);
  FreeMatrix(Sigma, n_dim);
  FreeMatrix(InvSigma, n_dim);
  free(dvtemp);
} /* main */

