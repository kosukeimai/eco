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

/** Normal-InvWishart updating 
    Y|mu, Sigma ~ N(mu, Sigma) 
       mu|Sigma ~ N(mu0, Sigma/tau0) 
          Sigma ~ InvWish(nu0, S0^{-1}) **/
void NIWupdate(
	       double **Y,         /* data */
	       double *mu,         /* mean */
	       double **Sigma,     /* variance */
	       double **InvSigma,  /* precision */
	       double *mu0,        /* prior mean */
	       double tau0,        /* prior scale */
	       int nu0,            /* prior df */
	       double **S0,        /* prior scale */
	       int n_samp,         /* sample size */
	       int n_dim)          /* dimension */
{
  int i,j,k;
  double *Ybar = doubleArray(n_dim);
  double *mun = doubleArray(n_dim);
  double **Sn = doubleMatrix(n_dim, n_dim);
  double **mtemp = doubleMatrix(n_dim, n_dim);

  /*read data */
  for (j=0; j<n_dim; j++) {
    Ybar[j] = 0;
    for (i=0; i<n_samp; i++)
      Ybar[j] += Y[i][j];
    Ybar[j] /= n_samp;
    for (k=0; k<n_dim; k++)
      Sn[j][k] = S0[j][k];
  }

  /* posterior updating*/

  for (j=0; j<n_dim; j++) 
    {
      mun[j] = (tau0*mu0[j]+n_samp*Ybar[j])/(tau0+n_samp);
      for (k=0; k<n_dim; k++) 
	{
	  Sn[j][k] += (tau0*n_samp)*(Ybar[j]-mu0[j])*(Ybar[k]-mu0[k])/(tau0+n_samp);
	  for (i=0; i<n_samp; i++)
	    Sn[j][k] += (Y[i][j]-Ybar[j])*(Y[i][k]-Ybar[k]);
	}
    }

  dinv(Sn, n_dim, mtemp);
  rWish(InvSigma, mtemp, nu0+n_samp, n_dim);
  dinv(InvSigma, n_dim, Sigma);
 
  for (j=0; j<n_dim; j++)
    for (k=0; k<n_dim; k++)
      mtemp[j][k] = Sigma[j][k]/(tau0+n_samp);

  rMVN(mu, mun, mtemp, n_dim);

  Free(Ybar);
  Free(mun);
  FreeMatrix(Sn, n_dim);
  FreeMatrix(mtemp, n_dim);
}
