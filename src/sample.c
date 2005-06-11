/******************************************************************
  This file is a part of ECO: R Package for Estimating Fitting Bayesian 
  Models of Ecological Inference for 2X2 tables
  by Ying Lu and Kosuke Imai
  Copyright: GPL version 2 or later.
*******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "sample.h"


/* Grid method samping from tomography line*/
void rGrid(
	  double *Sample,         /* W_i sampled from each tomography line */                 
	  double *W1gi,           /* The grid lines of W1[i] */
	  double *W2gi,           /* The grid lines of W2[i] */
	  int ni_grid,            /* number of grids for observation i*/
          double *mu,             /* mean vector for normal */ 
          double **InvSigma,         /* Inverse covariance matrix for normal */
          int n_dim)               /* dimension of parameters */
{
  int j;
  double dtemp;
  double *vtemp=doubleArray(n_dim);
  double *prob_grid=doubleArray(ni_grid);  /*density by grid */
  double *prob_grid_cum=doubleArray(ni_grid); /* cumulative density by grid */
    
  dtemp=0;
  for (j=0;j<ni_grid;j++){
    vtemp[0]=log(W1gi[j])-log(1-W1gi[j]);
    vtemp[1]=log(W2gi[j])-log(1-W2gi[j]);
    prob_grid[j]=dMVN(vtemp, mu, InvSigma, n_dim, 1) -
      log(W1gi[j])-log(W2gi[j])-log(1-W1gi[j])-log(1-W2gi[j]);
    prob_grid[j]=exp(prob_grid[j]);
    dtemp+=prob_grid[j];
    prob_grid_cum[j]=dtemp;
  }
  for (j=0;j<ni_grid;j++)
    prob_grid_cum[j]/=dtemp; /*standardize prob.grid */
  /*2 sample W_i on the ith tomo line */

  j=0;
  dtemp=unif_rand();
  while (dtemp > prob_grid_cum[j]) j++;
  Sample[0]=W1gi[j];
  Sample[1]=W2gi[j];

  free(vtemp);
  free(prob_grid);
  free(prob_grid_cum);

}

/* sample W via MH */
void rMH(
	  double *Sample,         /* W_i sampled from each tomography line */                 
          double *W,              /* the previous draws */
	  double *XY,               /*  X_i and Y_i */
	  double W1min,           /* The grid lines of W1[i] */
	  double W1max,           /* The grid lines of W2[i] */
          double *mu0,             /* mean vector for normal */ 
          double **InvSigma0,         /* Inverse covariance matrix for normal */
          int n_dim)               /* dimension of parameters */
{
  double dens1, dens2;
  double ratio;
  int j;
  
  /* draw Sample[0] (W_1) from unif(W1min, W1max) */
  Sample[0]=W1min+unif_rand()*(W1max-W1min);
  Sample[1]=XY[1]/(1-XY[0])-Sample[0]*XY[0]/(1-XY[0]);
  
  dens1 = dMVN(Sample, mu0, InvSigma0, n_dim, 1) -
    log(Sample[0])-log(Sample[1])-log(1-Sample[0])-log(1-Sample[1]);
  dens2 = dMVN(W, mu0, InvSigma0, n_dim, 1) -
    log(W[0])-log(W[1])-log(1-W[0])-log(1-W[1]);

  ratio=fmin2(1, exp(dens1-dens2));

  /* translate p(acceptance=1) into p(ratio< runif(1)) */
  if (ratio < unif_rand() ) {
    for (j=0; j<n_dim; j++)
      Sample[j]=W[j];
  }  
}

