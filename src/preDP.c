/******************************************************************
  This file is a part of eco: R Package for Fitting Bayesian Models 
  of Ecological Inference for 2x2 Tables
  by Kosuke Imai and Ying Lu
  Copyright: GPL version 2 or later.
*******************************************************************/

#include <stddef.h>
#include <string.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rinterface.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "bayes.h"
#include "sample.h"

/* Prediction for Nonparametric Model for 2x2 Tables */
void preDP(
	   double *pdmu, 
	   double *pdSigma,
	   int *pin_samp,
	   int *pin_draw,
	   int *pin_dim,
	   int *verbose,    /* 1 for output monitoring */
	   double *pdStore
	   ){	   
  
  /* some integers */
  int n_samp = *pin_samp;    /* sample size */
  int n_draw = *pin_draw;    /* sample size of survey data */ 
  int n_dim = *pin_dim;      /* dimension */

  double *mu = doubleArray(n_dim);                /* The mean */
  double *Wstar = doubleArray(n_dim);
  double **Sigma = doubleMatrix(n_dim, n_dim);    /* The covariance matrix */

  /* misc variables */
  int i, j, k, main_loop;   /* used for various loops */
  int itemp = 0;
  int itempM = 0;
  int itempS = 0;
  int progress = 1, itempP = ftrunc((double) n_draw/10);

  /* get random seed */
  GetRNGstate();
  
  for(main_loop=0; main_loop<n_draw; main_loop++){
    for(i=0; i<n_samp; i++) {
      for (j=0;j<n_dim;j++) {
	mu[j] = pdmu[itempM++];
	for (k=j;k<n_dim;k++) {
	  Sigma[j][k] = pdSigma[itempS++];
	  Sigma[k][j] = Sigma[j][k];
	}
      }
      rMVN(Wstar, mu, Sigma, n_dim);
      for (j=0; j<n_dim; j++)
	pdStore[itemp++] = exp(Wstar[j])/(1+exp(Wstar[j]));
    }
    if (*verbose)
      if (itempP == main_loop) {
        Rprintf("%3d percent done.\n", progress*10);
        itempP+=ftrunc((double) n_draw/10); progress++;
        R_FlushConsole();
      }
    R_CheckUserInterrupt();
  }
  
  if(*verbose)
    Rprintf("100 percent done.\n");

  /** write out the random seed **/
  PutRNGstate();

  /* Freeing the memory */
  free(mu);
  free(Wstar);
  FreeMatrix(Sigma,n_dim);
  
} /* main */

