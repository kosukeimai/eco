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
#include "sample.h"
#include "bayes.h"
#include "fintegrate.h"

#include "macros.h"




/*integrate normalize constant term */
void test(double *W1, int n, void *param)
{
  int ii;
  double mu[2];
  double Sigma[2][2];
  double *W2;
  double X, Y;
  double dtemp;

  Param *pp=(Param *)param;
  mu[0]=pp->mu[0];
  Sigma[0][0]=pp->Sigma[0][0];
  X=pp->X;
  Y=pp->Y;

  for (ii=0; ii<n; ii++)
    {
      W1[ii]=1/sqrt(M_PI*Sigma[0][0])*exp(-0.5*(W1[ii]-mu[0])*(W1[ii]-mu[0])/Sigma[0][0]);
    }

}


void NormConst(double *W1, int n, void *param)
{
  int ii;
  double mu[2];
  double Sigma[2][2];
  double *W2;
  double X, Y, rho;
  double dtemp;

  W2 = (double *) malloc(n*sizeof(double));

  Param *pp=(Param *)param;
  mu[0]=pp->mu[0];
  mu[1]=pp->mu[1];
  Sigma[0][0]=pp->Sigma[0][0];
  Sigma[1][1]=pp->Sigma[1][1];
  Sigma[0][1]=pp->Sigma[0][1];
  Sigma[1][0]=pp->Sigma[1][0];
  rho=Sigma[0][1]/sqrt(Sigma[0][0]*Sigma[1][1]);
  X=pp->X;
  Y=pp->Y;
 
  for (ii=0; ii<n; ii++)
    {
      W2[ii]=Y/(1-X)-X/(1-X)/(1+exp(-W1[ii]));
      W2[ii]=log(W2[ii]/(1-W2[ii]));
    }

  dtemp=1/(2*M_PI*sqrt(Sigma[0][0]*Sigma[1][1]*(1-rho*rho)));

  for (ii=0; ii<n; ii++)
    {
      W1[ii]=exp(-1/(2*(1-rho*rho))*
		 ((W1[ii]-mu[0])*(W1[ii]-mu[0])/Sigma[0][0]+
		  (W2[ii]-mu[1])*(W2[ii]-mu[1])/Sigma[1][1]-
		  2*rho*(W1[ii]-mu[0])*(W2[ii]-mu[0])
		  /sqrt(Sigma[0][0]*Sigma[1][1])))/dtemp;
    }
 
  free(W2);
}


/*E(W1) */

void W1Exp(double *W1, int n, void *param)
{
  int ii;
  int dim=2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double *W2;
  double X, Y;
  double dtemp, rho;
  
  double *vtemp=doubleArray(dim);

  W2 = (double *) malloc(n*sizeof(double));
  
  Param *pp=(Param *)param;
  mu[0]=pp->mu[0];
  mu[1]=pp->mu[1];
  Sigma[0][0]=pp->Sigma[0][0];
  Sigma[1][1]=pp->Sigma[1][1];
  Sigma[0][1]=pp->Sigma[0][1];
  Sigma[1][0]=pp->Sigma[1][0];

  Rprintf("Sigma00 15%g\n", Sigma[0][0]);
  Rprintf("Sigma11 15%g\n", Sigma[1][1]);
  Rprintf("Sigma01 15%g\n", Sigma[0][1]);


  X=pp->X;
  Y=pp->Y;
  
  for (ii=0; ii<n; ii++)
    {
      W2[ii] = Y/(1-X)-X/(1-X)/(1+exp(-W1[ii]));
      W2[ii] = log(W2[ii]/(1-W2[ii]));
    }
  
  for (ii=0; ii<n; ii++)
    {
      vtemp[0] = W1[ii];
      vtemp[1] = W2[ii];
      W1[ii] = W1[ii]*dBVNtomo(vtemp, X, Y, mu, Sigma, 0);
    }
  free(W2);
  free(mu);
  FreeMatrix(Sigma,dim);
}


/**E(W1*W1)*/
void W1W1Exp(double *W1, int n, void *param)
{
  int ii;
  double mu[2];
  double Sigma[2][2];
  double *W2;
  double X, Y;
  double dtemp, rho;
  
  double vtemp[2];
  W2 = (double *) malloc(1000*sizeof(double));

  Param *pp=(Param *)param;
  mu[0]=pp->mu[0];
  mu[1]=pp->mu[1];
  Sigma[0][0]=pp->Sigma[0][0];
  Sigma[1][1]=pp->Sigma[1][1];
  Sigma[0][1]=pp->Sigma[0][1];
  Sigma[1][0]=pp->Sigma[1][0];
  rho=Sigma[0][1]/sqrt(Sigma[0][0]*Sigma[1][1]);
  X=pp->X;
  Y=pp->Y;
  
  for (ii=0; ii<n; ii++)
    {
      W2[ii] = Y/(1-X)-X/(1-X)/(1+exp(-W1[ii]));
      W2[ii] = log(W2[ii]/(1-W2[ii]));
    }
  
  for (ii=0; ii<n; ii++)
    {
      vtemp[0] = W1[ii];
      vtemp[1] = W2[ii];
      W1[ii] = W1[ii] * W1[ii] * dBVNtomo(vtemp, X,Y, mu, Sigma,0);
    }
  free(W2);
}

