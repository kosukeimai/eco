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
#include "sample.h"

/* Multivariate Normal density */
double dMVN(			
	double *Y,		/* The data */
	double *MEAN,		/* The parameters */
	double **SIG_INV,         /* inverse of the covariance matrix */	
	int dim,                /* dimension */
	int give_log){          /* 1 if log_scale 0 otherwise */
  
  int j,k;
  double value=0.0;

  for(j=0;j<dim;j++){
    for(k=0;k<j;k++)
      value+=2*(Y[k]-MEAN[k])*(Y[j]-MEAN[j])*SIG_INV[j][k];
    value+=(Y[j]-MEAN[j])*(Y[j]-MEAN[j])*SIG_INV[j][j];
  }

  value=-0.5*value-0.5*dim*log(2*M_PI)+0.5*ddet(SIG_INV, dim, 1);


  if(give_log)  
    return(value);
  else
    return(exp(value));

}

/* the density of Multivariate T-distribution */
double dMVT(
            double *Y,          /* The data */
            double *MEAN,       /* mean */
            double **SIG_INV,   /* inverse of scale matrix */
            int nu,             /* Degrees of freedom */
            int dim,            /* dimension */
            int give_log)       /* 1 if log_scale 0 otherwise */
{
  int j,k;
  double value=0;

  for(j=0;j<dim;j++){
    for(k=0;k<j;k++)
      value+=2*(Y[k]-MEAN[k])*(Y[j]-MEAN[j])*SIG_INV[j][k];
    value+=(Y[j]-MEAN[j])*(Y[j]-MEAN[j])*SIG_INV[j][j];
  }

  value=0.5*ddet(SIG_INV, dim,1) - 0.5*dim*(log((double)nu)+log(M_PI)) -
    0.5*((double)dim+nu)*log(1+value/(double)nu) +
    lgammafn(0.5*(double)(nu+dim)) - lgammafn(0.5*(double)nu);

  if(give_log)
    return(value);
  else
    return(exp(value));
}


/* Sample from the MVN dist */
void rMVN(                      
	  double *Sample,         /* Vector for the sample */
	  double *mean,           /* The vector of means */
	  double **Var,           /* The matrix Variance */
	  int size)               /* The dimension */
{
  int j,k;
  double **Model = doubleMatrix(size+1, size+1);
  double cond_mean;
    
  /* draw from mult. normal using SWP */
  for(j=1;j<=size;j++){       
    for(k=1;k<=size;k++)
      Model[j][k]=Var[j-1][k-1];
    Model[0][j]=mean[j-1];
    Model[j][0]=mean[j-1];
  }
  Model[0][0]=-1;
  Sample[0]=(double)norm_rand()*sqrt(Model[1][1])+Model[0][1];
  for(j=2;j<=size;j++){
    SWP(Model,j-1,size+1);
    cond_mean=Model[j][0];
    for(k=1;k<j;k++) cond_mean+=Sample[k-1]*Model[j][k];
    Sample[j-1]=(double)norm_rand()*sqrt(Model[j][j])+cond_mean;
  }
  
  FreeMatrix(Model,size+1);
}


/* Sample from a wish dist */
/* Odell, P. L. and Feiveson, A. H. ``A Numerical Procedure to Generate
   a Sample Covariance Matrix'' Journal of the American Statistical
   Association, Vol. 61, No. 313. (Mar., 1966), pp. 199-203. */

void rWish(                  
	   double **Sample,        /* The matrix with to hold the sample */
	   double **S,             /* The parameter */
	   int df,                 /* the degrees of freedom */
	   int size)               /* The dimension */
{
  int i,j,k;
  double *V = doubleArray(size);
  double **B = doubleMatrix(size, size);
  double **C = doubleMatrix(size, size);
  double **N = doubleMatrix(size, size);
  double **mtemp = doubleMatrix(size, size);
  
  for(i=0;i<size;i++) {
    V[i]=rchisq((double) df-i-1);
    B[i][i]=V[i];
    for(j=(i+1);j<size;j++)
      N[i][j]=norm_rand();
  }

  for(i=0;i<size;i++) {
    for(j=i;j<size;j++) {
      Sample[i][j]=0;
      Sample[j][i]=0;
      mtemp[i][j]=0;
      mtemp[j][i]=0;
      if(i==j) {
	if(i>0)
	  for(k=0;k<j;k++)
	    B[j][j]+=N[k][j]*N[k][j];
      }
      else { 
	B[i][j]=N[i][j]*sqrt(V[i]);
	if(i>0)
	  for(k=0;k<i;k++)
	    B[i][j]+=N[k][i]*N[k][j];
      }
      B[j][i]=B[i][j];
    }
  }
  
  dcholdc(S, size, C);
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      for(k=0;k<size;k++)
	mtemp[i][j]+=C[i][k]*B[k][j];
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      for(k=0;k<size;k++)
	Sample[i][j]+=mtemp[i][k]*C[j][k];

  free(V);
  FreeMatrix(B, size);
  FreeMatrix(C, size);
  FreeMatrix(N, size);
  FreeMatrix(mtemp, size);
}

/* Sample from a Dirichlet distribution */
void rDirich(
	     double *Sample, /* Vector for the sample */
	     double *theta,  /* parameters */
	     int size)       /* The dimension */
{
  int j;
  double dtemp=0;
  
  for (j=0; j<size; j++) {
    Sample[j] = rgamma(theta[j], 1.0);
    dtemp += Sample[j];
  }
  for (j=0 ; j<size; j++)
    Sample[j] /= dtemp;
}
