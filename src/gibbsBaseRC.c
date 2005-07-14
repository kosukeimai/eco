/******************************************************************
  This file is a part of eco: R Package for Ecological Inference
  by Kosuke Imai and Ying Lu
  Copyright: GPL version 2 or later.
*******************************************************************/

#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "bayes.h"
#include "sample.h"

/* Normal Parametric Model for RxC (with R >= 2, C >= 2) Tables */
void cBaseRC(
	     /*data input */
	     double *pdX,     /* X */
	     double *pdY,     /* Y */
	     double *pdWmin,  /* lower bounds */
	     double *pdWmax,  /* uppwer bounds */
	     int *pin_samp,   /* sample size */
	     int *pin_col,    /* number of columns */
	     int *pin_row,    /* number of rows */

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
	     double *pdMu,
	     double *pdSigma,
	     
	     /* storage */
	     int *parameter,  /* 1 if save population parameter */
	     double *pdSmu, 
	     double *pdSSigma,
	     double *pdSW
	     ){	   
  
  /* some integers */
  int n_samp = *pin_samp;    /* sample size */
  int nth = *pinth;          /* keep every pth draw */
  int n_col = *pin_col;      /* number of columns */
  int n_dim = *pin_row-1;    /* number of rows - 1 */

  /* prior parameters */ 
  double tau0 = *pdtau0;                     /* prior scale */
  int nu0 = *pinu0;                          /* prior degrees of freedom */   
  double **S0 = doubleMatrix(n_col, n_col);  /* prior scale for InvWish */

  /* data */
  double **Y = doubleMatrix(n_samp, n_dim);               /* Y */
  double **X = doubleMatrix(n_samp, n_col);               /* X */
  double ***W = doubleMatrix3D(n_samp, n_dim, n_col);     /* W */
  double ***Wstar = doubleMatrix3D(n_col, n_samp, n_dim); /* logratio(W) */       
  double **Wsum = doubleMatrix(n_samp, n_col);            /* sum_{r=1}^{R-1} W_{irc} */
  double **SWstar = doubleMatrix(n_col, n_dim);

  /* The lower and upper bounds of U = W*X/Y **/
  double ***minU = doubleMatrix3D(n_samp, n_dim, n_col);
  double *maxU = doubleArray(n_col);

  /* model parameters */
  double **mu = doubleMatrix(n_col, n_dim);                 /* mean */
  double ***Sigma = doubleMatrix3D(n_col, n_dim, n_dim);    /* covariance */
  double ***InvSigma = doubleMatrix3D(n_col, n_dim, n_dim); /* inverse */

  /* misc variables */
  int i, j, k, l, main_loop;   /* used for various loops */
  int itemp, counter;
  int itempM = 0;           /* for mu */
  int itempS = 0;           /* for Sigma */
  int itempW = 0;           /* for W */
  int itempC = 0;           /* control nth draw */
  int progress = 1, itempP = ftrunc((double) *n_gen/10);
  double dtemp, dtemp1;
  double *param = doubleArray(n_col);   /* Dirichlet parameters */
  double *dvtemp = doubleArray(n_col);
  double *dvtemp1 = doubleArray(n_col);

  /* get random seed */
  GetRNGstate();
  
  /* read X */
  itemp = 0;
  for (k = 0; k < n_col; k++) 
    for (i = 0; i < n_samp; i++) 
      X[i][k] = pdX[itemp++];

  /* read Y */
  itemp = 0;
  for (j = 0; j < n_dim; j++) 
    for (i = 0; i < n_samp; i++) 
      Y[i][j] = pdY[itemp++];

  /* compute bounds on U */
  itemp = 0; 
  for (k = 0; k < n_col; k++) 
    for (j = 0; j < n_dim; j++) 
      for (i = 0; i < n_samp; i++) 
	minU[i][j][k] = fmax2(0, (X[i][k]+Y[i][j]-1)/Y[i][j]);

  /* initial values for mu and Sigma */
  itemp = 0;
  for (k = 0; k < n_col; k++)
    for (j = 0; j < n_dim; j++)
      mu[k][j] = pdMu[itemp++]; 
  itemp = 0;
  for (k = 0; k < n_col; k++)
    for (j = 0; j < n_dim; j++) 
      for (i = 0; i < n_dim; i++) 
	Sigma[k][j][i] = pdSigma[itemp++];
  for (k = 0; k < n_col; k++)
    dinv(Sigma[k], n_dim, InvSigma[k]);
  
  /* initial values for W */
  for (k = 0; k < n_col; k++)
    param[k] = 1.0;
  for (i = 0; i < n_samp; i++) {
    for (k = 0; k < n_col; k++)
      Wsum[i][k] = 0.0;
    for (j = 0; j < n_dim; j++) {
      counter = 0; itemp = 1; 
      while (itemp > 0) { /* first try rejection sampling */
	rDirich(dvtemp, param, n_col);
	itemp = 0;
	for (k = 0; k < n_col; k++) {
	  if (dvtemp[k] < minU[i][j][k] || 
	      dvtemp[k] > fmin2(1, X[i][k]*(1-Wsum[i][k])/Y[i][j]))
	    itemp++;
	}
	if (itemp < 1) 
	  for (k = 0; k < n_col; k++) {
	    W[i][j][k] = dvtemp[k]*Y[i][j]/X[i][k];
	    Wsum[i][k] += W[i][j][k];
	  }
	counter++;
	if (counter > *maxit && itemp > 0) { /* if rejection sampling fails, then
				   use midpoints of bounds */
	  itemp = 0;
	  dtemp = Y[i][j]; dtemp1 = 1;
	  for (k = 0; k < n_col-1; k++) {
	    W[i][j][k] = 0.25*(fmax2(0,(X[i][k]/dtemp1+dtemp-1)*dtemp1/X[i][k])+
			      fmin2(1-Wsum[i][k],dtemp*dtemp1/X[i][k]));
	    dtemp -= W[i][j][k]*X[i][k]/dtemp1;
	    dtemp1 -= X[i][k];
	    Wsum[i][k] += W[i][j][k];
	  }
	  W[i][j][n_col-1] = dtemp;
	  Wsum[i][n_col-1] += dtemp;
	}
	R_CheckUserInterrupt();
      }
      for (l = 0; l < n_dim; l++) 
	for (k = 0; k < n_col; k++)
	  Wstar[k][i][l] = log(W[i][l][k])-log(1-Wsum[i][k]);
    }
  }

  /* read the prior */
  itemp = 0;
  for(k = 0; k < n_dim; k++)
    for(j = 0; j < n_dim; j++) 
      S0[j][k] = pdS0[itemp++];

  /*** Gibbs sampler! ***/
  if (*verbose)
    Rprintf("Starting Gibbs sampler...\n");
  for(main_loop = 0; main_loop < *n_gen; main_loop++){
    /** update W, Wstar given mu, Sigma **/
    for (i = 0; i < n_samp; i++) {
      /* sampling W through Metropolis Step for each row */
      for (j = 0; j < n_dim; j++) {
	/* computing upper bounds for U */
	for (k = 0; k < n_col; k++) {
	  Wsum[i][k] -= W[i][j][k];
	  maxU[k] = fmin2(1, X[i][k]*(1-Wsum[i][k])/Y[i][j]);
	}
	/** MH step **/
	/* Sample a candidate draw of W from truncated Dirichlet */
	l = 0; itemp = 1;
	while (itemp > 0) {
	  rDirich(dvtemp, param, n_col);
	  itemp = 0;
	  for (k = 0; k < n_col; k++) 
	    if (dvtemp[k] > maxU[k] || dvtemp[k] < minU[i][j][k])
	      itemp++;
	  l++;
	  if (l > *maxit)
	    error("rejection algorithm failed because bounds are too tight.\n increase maxit or use gibbs sampler instead.");
	}
	/* get W and its log-ratio transformation */
	for (k = 0; k < n_col; k++) {
	  dvtemp[k] = dvtemp[k]*Y[i][j]/X[i][k];
	  dvtemp1[k] = Wsum[i][k]+dvtemp[k];
	}
	for (k = 0; k < n_col; k++) 
	  for (l = 0; l < n_dim; l++) 
	    if (l == j)
	      SWstar[k][l] = log(dvtemp[k])-log(1-dvtemp1[k]);
	    else
	      SWstar[k][l] = log(W[i][j][k])-log(1-dvtemp1[k]);
	/* computing acceptance ratio */
	dtemp = 0; dtemp1 = 0;
	for (k= 0; k < n_col; k++) {
	  dtemp += dMVN(SWstar[k], mu[k], InvSigma[k], n_dim, 1);
	  dtemp1 += dMVN(Wstar[k][i], mu[k], InvSigma[k], n_dim, 1);
	  dtemp -= log(dvtemp[k]);
	  dtemp1 -= log(W[i][j][k]);
	}
	if (unif_rand() < fmin2(1, exp(dtemp-dtemp1))) 
	  for (k = 0; k < n_col; k++)
	    W[i][j][k] = dvtemp[k]; 
	/* updating Wsum and Wstar with new draws */
	for (k = 0; k < n_col; k++) {
	  Wsum[i][k] += W[i][j][k];
	  for (l = 0; l < n_dim; l++) 
	    Wstar[k][i][l] = log(W[i][l][k])-log(1-Wsum[i][k]);
	}
      }
    }    
    
    /* update mu, Sigma given wstar using effective sample of Wstar */
    for (k = 0; k < n_col; k++)
      NIWupdate(Wstar[k], mu[k], Sigma[k], InvSigma[k], mu0, tau0,
		nu0, S0, n_samp, n_dim); 
    
    /*store Gibbs draw after burn-in and every nth draws */     
    if (main_loop >= *burn_in){
      itempC++;
      if (itempC==nth){
	for (k = 0; k < n_col; k++) 
	  for (j = 0; j < n_dim; j++) {
	    pdSmu[itempM++]=mu[k][j];
	    for (i = 0; i < n_dim; i++)
	      if (j <= i)
		pdSSigma[itempS++]=Sigma[k][j][i];
	  }
	for(i = 0; i < n_samp; i++)
	  for (k = 0; k < n_col; k++)
	    for (j = 0; j < n_dim; j++)
	      pdSW[itempW++] = W[i][j][k];
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
  if (*verbose)
    Rprintf("100 percent done.\n");

  /** write out the random seed **/
  PutRNGstate();

  /* Freeing the memory */
  FreeMatrix(S0, n_col);
  FreeMatrix(X, n_samp);
  FreeMatrix(Y, n_samp);
  Free3DMatrix(W, n_samp, n_dim);
  Free3DMatrix(Wstar, n_col, n_samp);
  FreeMatrix(Wsum, n_samp);
  Free3DMatrix(minU, n_samp, n_dim);
  FreeMatrix(mu, n_col);
  Free3DMatrix(Sigma, n_col, n_dim);
  Free3DMatrix(InvSigma, n_col, n_dim);
  free(param);
  free(dvtemp);
} /* main */

