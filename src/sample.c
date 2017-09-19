/******************************************************************
  This file is a part of eco: R Package for Fitting Bayesian Models 
  of Ecological Inference for 2x2 Tables
  by Kosuke Imai and Ying Lu
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


/* Grid method samping from tomography line*/
void rGrid(
	   double *Sample,         /* W_i sampled from each tomography line */                 
	   double *W1gi,           /* The grid lines of W1[i] */
	   double *W2gi,           /* The grid lines of W2[i] */
	   int ni_grid,            /* number of grids for observation i*/
	   double *mu,             /* mean vector for normal */ 
	   double **InvSigma,      /* Inverse covariance matrix for normal */
	   int n_dim)              /* dimension of parameters */
{
  int j;
  double dtemp;
  double *vtemp=doubleArray(n_dim);
  double *prob_grid=doubleArray(ni_grid);     /* density by grid */
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

  Free(vtemp);
  Free(prob_grid);
  Free(prob_grid_cum);

}

/* preparation for Grid */
void GridPrep(
	      double **W1g,  /* grids holder for W1 */
	      double **W2g,  /* grids holder for W2 */
	      double **X,    /* data: [X Y] */
	      double *maxW1, /* upper bound for W1 */
	      double *minW1, /* lower bound for W1 */
	      int *n_grid,   /* number of grids */
	      int  n_samp,   /* sample size */
	      int  n_step    /* step size */
)
{
  int i, j;
  double dtemp;
  double *resid = doubleArray(n_samp);

  for(i=0; i<n_samp; i++)
    for (j=0; j<n_step; j++){
      W1g[i][j]=0;
      W2g[i][j]=0;
    }
  for(i=0;i<n_samp;i++) {
    if (X[i][1]!=0 && X[i][1]!=1) {
      /* 1/n_step is the length of the grid */
      dtemp=(double)1/n_step;
      if ((maxW1[i]-minW1[i]) > (2*dtemp)) { 
	n_grid[i]=ftrunc((maxW1[i]-minW1[i])*n_step);
	resid[i]=(maxW1[i]-minW1[i])-n_grid[i]*dtemp;
	/*if (maxW1[i]-minW1[i]==1) resid[i]=dtemp/4; */
	j=0; 
	while (j<n_grid[i]) {
	  W1g[i][j]=minW1[i]+(j+1)*dtemp-(dtemp+resid[i])/2;
	  if ((W1g[i][j]-minW1[i])<resid[i]/2) W1g[i][j]+=resid[i]/2;
	  if ((maxW1[i]-W1g[i][j])<resid[i]/2) W1g[i][j]-=resid[i]/2;
	  W2g[i][j]=(X[i][1]-X[i][0]*W1g[i][j])/(1-X[i][0]);
	  j++;
	}
      }
      else {
	W1g[i][0]=minW1[i]+(maxW1[i]-minW1[i])/3;
	W2g[i][0]=(X[i][1]-X[i][0]*W1g[i][0])/(1-X[i][0]);
	W1g[i][1]=minW1[i]+2*(maxW1[i]-minW1[i])/3;
	W2g[i][1]=(X[i][1]-X[i][0]*W1g[i][1])/(1-X[i][0]);
	n_grid[i]=2;
      }
    }
  }

  Free(resid);
}

/* sample W via MH for 2x2 table */
void rMH(
	 double *W,              /* previous draws */
	 double *XY,             /* X_i and Y_i */
	 double W1min,           /* lower bound for W1 */
	 double W1max,           /* upper bound for W1 */
	 double *mu,            /* mean vector for normal */ 
	 double **InvSigma,     /* Inverse covariance matrix for normal */
	 int n_dim)              /* dimension of parameters */
{
  int j;
  double dens1, dens2, ratio;
  double *Sample = doubleArray(n_dim);
  double *vtemp = doubleArray(n_dim);
  double *vtemp1 = doubleArray(n_dim);
  
  /* sample W_1 from unif(W1min, W1max) */
  Sample[0] = runif(W1min, W1max);
  Sample[1] = XY[1]/(1-XY[0])-Sample[0]*XY[0]/(1-XY[0]);
  for (j = 0; j < n_dim; j++) {
    vtemp[j] = log(Sample[j])-log(1-Sample[j]);
    vtemp1[j] = log(W[j])-log(1-W[j]);
  }
  
  /* acceptance ratio */
  dens1 = dMVN(vtemp, mu, InvSigma, n_dim, 1) -
    log(Sample[0])-log(Sample[1])-log(1-Sample[0])-log(1-Sample[1]);
  dens2 = dMVN(vtemp1, mu, InvSigma, n_dim, 1) -
    log(W[0])-log(W[1])-log(1-W[0])-log(1-W[1]);
  ratio = fmin2(1, exp(dens1-dens2));
  
  /* accept */
  if (unif_rand() < ratio) 
    for (j=0; j<n_dim; j++) 
      W[j]=Sample[j];
  
  Free(Sample);
  Free(vtemp);
  Free(vtemp1);
}


/* sample W via MH for 2xC table */
void rMH2c(
	   double *W,              /* W */
	   double *X,              /* X_i */
	   double Y,               /* Y_i */
	   double *minU,           /* lower bound for U */
	   double *maxU,           /* upper bound for U */
	   double *mu,             /* mean vector for normal */ 
	   double **InvSigma,      /* Inverse covariance matrix for normal */
	   int n_dim,              /* dimension of parameters */
	   int maxit,              /* max number of iterations for
				      rejection sampling */
	   int reject)             /* if 1, use rejection sampling to
				      draw from the truncated Dirichlet
				      if 0, use Gibbs sampling
				   */  
{
  int iter = 100;   /* number of Gibbs iterations */
  int i, j, exceed;
  double dens1, dens2, ratio, dtemp;
  double *Sample = doubleArray(n_dim);
  double *param = doubleArray(n_dim);
  double *vtemp = doubleArray(n_dim);
  double *vtemp1 = doubleArray(n_dim);
  
  /* set parent Dirichlet parameter to 1 */
  for (j = 0; j < n_dim; j++)
    param[j] = 1.0;

  /* Sample a candidate draw of W from truncated Dirichlet */
  if (reject) { /* rejection sampling */
    i = 0; exceed = 1;
    while (exceed > 0) {
      rDirich(vtemp, param, n_dim);
      exceed = 0;
      for (j = 0; j < n_dim; j++) 
	if (vtemp[j] > maxU[j] || vtemp[j] < minU[j])
	  exceed++;
      i++;
      if (i > maxit)
	error("rMH2c: rejection algorithm failed because bounds are too tight.\n increase maxit or use gibbs sampler instead.");
    }
  }
  else { /* gibbs sampler */
    for (j = 0; j < n_dim; j++) 
      vtemp[j] = W[j]*X[j]/Y;
    for (i = 0; i < iter; i++) {
      dtemp = vtemp[n_dim-1];
      for (j = 0; j < n_dim-1; j++) {
	dtemp += vtemp[j];
	vtemp[j] = runif(fmax2(minU[j], dtemp-maxU[n_dim-1]), 
			 fmin2(maxU[j], dtemp-minU[n_dim-1]));
	dtemp -= vtemp[j];
      }
      vtemp[n_dim-1] = dtemp;
    }
  }
  /* calcualte W and its logit transformation */
  for (j = 0; j < n_dim; j++) {
    Sample[j] = vtemp[j]*Y/X[j];
    vtemp[j] = log(Sample[j])-log(1-Sample[j]);
    vtemp1[j] = log(W[j])-log(1-W[j]);
  }
  
  /* acceptance ratio */
  dens1 = dMVN(vtemp, mu, InvSigma, n_dim, 1);
  dens2 = dMVN(vtemp1, mu, InvSigma, n_dim, 1);
  for (j=0; j<n_dim; j++) {
    dens1 -= (log(Sample[j])+log(1-Sample[j]));
    dens2 -= (log(W[j])+log(1-W[j]));
  }
  ratio=fmin2(1, exp(dens1-dens2));
  
  /* accept */
  if (unif_rand() < ratio) 
    for (j = 0; j < n_dim; j++)
      W[j] = Sample[j];
  
  Free(Sample);
  Free(param);
  Free(vtemp);
  Free(vtemp1);
}


