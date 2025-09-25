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

/* Normal Parametric Model for 2xC (with C > 2) Tables */
void cBase2C(
	     /*data input */
	     double *pdX,     /* X: matrix */
	     double *Y,       /* Y: vector */
	     double *pdWmin,  /* lower bounds */
	     double *pdWmax,  /* uppwer bounds */
	     int *pin_samp,   /* sample size */
	     int *pin_col,    /* number of columns */
	     
	     /*MCMC draws */
	     int *reject,     /* whether to use rejection sampling */
	     int *maxit,      /* max number of iterations for
				 rejection sampling */
	     int *n_gen,      /* number of gibbs draws */
	     int *burn_in,    /* number of draws to be burned in */
	     int *pinth,      /* keep every nth draw */
	     int *verbose,    /* 1 for output monitoring */
	     
	     /* prior specification*/
	     int *pinu0,      /* prior df parameter for InvWish */
	     double *pdtau0,  /* prior scale parameter for Sigma */
	     double *mu0,     /* prior mean for mu */
	     double *pdS0,    /* prior scale for Sigma */

	     /* starting values */
	     double *mu,
	     double *SigmaStart,
	     
	     /* storage */
	     int *parameter,  /* 1 if save population parameter */
	     double *pdSmu, 
	     double *pdSSigma,
	     double *pdSW
	     ){	   
  
  /* some integers */
  int n_samp = *pin_samp;    /* sample size */
  int nth = *pinth;          /* keep every pth draw */
  int n_col = *pin_col;      /* dimension */

  /* prior parameters */ 
  double tau0 = *pdtau0;                          /* prior scale for mu */
  int nu0 = *pinu0;                               /* prior degrees of freedom */   
  double **S0 = doubleMatrix(n_col, n_col);       /* prior scale for Sigma */

  /* data */
  double **X = doubleMatrix(n_samp, n_col);       /* X */
  double **W = doubleMatrix(n_samp, n_col);       /* The W matrix */
  double **Wstar = doubleMatrix(n_samp, n_col);   /* logit(W) */     

  /* The lower and upper bounds of U = W*X/Y **/
  double **minU = doubleMatrix(n_samp, n_col);
  double **maxU = doubleMatrix(n_samp, n_col);    

  /* model parameters */
  double **Sigma = doubleMatrix(n_col, n_col);    /* The covariance matrix */
  double **InvSigma = doubleMatrix(n_col, n_col); /* The inverse covariance matrix */

  /* misc variables */
  int i, j, k, main_loop;   /* used for various loops */
  int itemp;
  int itempM = 0; /* for mu */
  int itempS = 0; /* for Sigma */
  int itempW = 0; /* for W */
  int itempC = 0; /* control nth draw */
  int progress = 1, itempP = ftrunc((double) *n_gen/10);
  double dtemp, dtemp1;
  double *param = doubleArray(n_col);   /* Dirichlet parameters */
  double *dvtemp = doubleArray(n_col);

  /* get random seed */
  GetRNGstate();
  
  /* read X */
  itemp = 0;
  for (j = 0; j < n_col; j++) 
    for (i = 0; i < n_samp; i++) 
      X[i][j] = pdX[itemp++];

  /* read initial values of Sigma */
  itemp = 0;
  for (k = 0; k < n_col; k++) 
    for (j = 0; j < n_col; j++) 
      Sigma[j][k] = SigmaStart[itemp++];
  dinv(Sigma, n_col, InvSigma);

  /* compute bounds on U */
  itemp = 0;
  for (j = 0; j < n_col; j++) 
    for (i = 0; i < n_samp; i++) 
      minU[i][j] = fmax2(0, pdWmin[itemp++]*X[i][j]/Y[i]);
  itemp = 0;
  for (j = 0; j < n_col; j++) 
    for (i = 0; i < n_samp; i++) 
      maxU[i][j] = fmin2(1, pdWmax[itemp++]*X[i][j]/Y[i]);

  /* initial values for W */
  for (j = 0; j < n_col; j++)
    param[j] = 1;
  for (i = 0; i < n_samp; i++) {
    k = 0; itemp = 1;
    while (itemp > 0) { /* rejection sampling */
      rDirich(dvtemp, param, n_col);
      itemp = 0; k++;
      for (j = 0; j < n_col; j++)
	if (dvtemp[j] > maxU[i][j] || dvtemp[j] < minU[i][j])
	  itemp++;
      if (itemp == 0)
	for (j = 0; j < n_col; j++)
	  W[i][j] = dvtemp[j]*Y[i]/X[i][j];
      if (k > *maxit) { /* if rejection sampling fails, then use
			   midpoits of bounds sequentially */
	itemp = 0;
	dtemp = Y[i]; dtemp1 = 1;
	for (j = 0; j < n_col-1; j++) {
	  W[i][j] = 0.5*(fmax2(0,(X[i][j]/dtemp1+dtemp-1)*dtemp1/X[i][j])+
			 fmin2(1,dtemp*dtemp1/X[i][j]));
	  dtemp -= W[i][j]*X[i][j]/dtemp1;
	  dtemp1 -= X[i][j];
	}
	W[i][n_col-1] = dtemp;
      }
    }
    for (j = 0; j < n_col; j++) 
      Wstar[i][j] = log(W[i][j])-log(1-W[i][j]);
  }

  /* read the prior */
  itemp = 0;
  for(k = 0; k < n_col; k++)
    for(j = 0; j < n_col; j++) 
      S0[j][k] = pdS0[itemp++];

  /*** Gibbs sampler! ***/
  if (*verbose)
    Rprintf("Starting Gibbs sampler...\n");
  for(main_loop = 0; main_loop < *n_gen; main_loop++){
    /** update W, Wstar given mu, Sigma **/
    for (i = 0; i < n_samp; i++){
      rMH2c(W[i], X[i], Y[i], minU[i], maxU[i], mu, InvSigma, n_col,
	    *maxit, *reject);
      for (j = 0; j < n_col; j++) 
	Wstar[i][j] = log(W[i][j])-log(1-W[i][j]);
    }
    
    /* update mu, Sigma given wstar using effective sample of Wstar */
    NIWupdate(Wstar, mu, Sigma, InvSigma, mu0, tau0, nu0, S0, n_samp, n_col);
    
    /*store Gibbs draw after burn-in and every nth draws */      
    if (main_loop>=*burn_in){
      itempC++;
      if (itempC==nth){
	for (j = 0; j < n_col; j++) {
	  pdSmu[itempM++]=mu[j];
	  for (k = 0; k < n_col; k++)
	    if (j <=k)
	      pdSSigma[itempS++]=Sigma[j][k];
	}
	for(i = 0; i < n_samp; i++)
	  for (j = 0; j < n_col; j++)
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
  FreeMatrix(S0, n_col);
  FreeMatrix(X, n_samp);
  FreeMatrix(W, n_samp);
  FreeMatrix(Wstar, n_samp);
  FreeMatrix(minU, n_samp);
  FreeMatrix(maxU, n_samp);
  FreeMatrix(Sigma, n_col);
  FreeMatrix(InvSigma, n_col);
  free(dvtemp);
  free(param);
} /* main */

