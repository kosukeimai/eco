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
  double dtemp, inp;
  int imposs;

  W2 = (double *) malloc(n*sizeof(double));
  if (W2==NULL) Rprintf("Malloc Error");

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
  imposs=0;


  dtemp=1/(2*M_PI*sqrt(Sigma[0][0]*Sigma[1][1]*(1-rho*rho)));
  //dtemp=1/(sqrt(2*M_PI*Sigma[0][0]*Sigma[1][1]*(1-rho*rho)));

  for (ii=0; ii<n; ii++)
    {
      imposs=0; inp=W1[ii];
      W2[ii]=Y/(1-X)-X/(1-X)/(1+exp(-W1[ii]));
      if(W2[ii]>=1) imposs=1; //impossible pair of values
      else W2[ii]=log(W2[ii]/(1-W2[ii]));
    //}

  //Rprintf("n: %d W1[0] %5g W2[0] %5g Y %5g X %5g imposs %d\n", n, W1[0], W2[0],Y,X,imposs);

  //for (ii=0; ii<n; ii++)
    //{
      /*W1[ii]=exp(-1/(2*(1-rho*rho))*
		 ((W1[ii]-mu[0])*(W1[ii]-mu[0])/Sigma[0][0]+
		  (W2[ii]-mu[1])*(W2[ii]-mu[1])/Sigma[1][1]-
		  2*rho*(W1[ii]-mu[0])*(W2[ii]-mu[0])
		  /sqrt(Sigma[0][0]*Sigma[1][1])))/dtemp; */
		  if (imposs==1) W1[ii]=0;
      else W1[ii]=exp(-1/(2*(1-rho*rho))*
		 ((W1[ii]-mu[0])*(W1[ii]-mu[0])/Sigma[0][0]+
		  (W2[ii]-mu[1])*(W2[ii]-mu[1])/Sigma[1][1]-
		  2*rho*(W1[ii]-mu[0])*(W2[ii]-mu[1])
		  /sqrt(Sigma[0][0]*Sigma[1][1])))*dtemp;

		  //Rprintf("Normc... %d %d %5g -> %5g with %5g imposs %d\n", ii, n, inp, W1[ii],dtemp,imposs);

    }
  //Rprintf("W1[0] %15g\n", W1[0]);
  free(W2);
}


/*E(W1)
  assumes bivariate normal distribution
*/

void W1Exp(double *W1, int n, void *param)
{
  int ii,imposs;
  int dim=2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double *W2;
  double X, Y, normc;
  double dtemp, rho, inp;

  double *vtemp=doubleArray(dim);

  W2 = (double *) malloc(n*sizeof(double));
  if (W2==NULL) Rprintf("Malloc Error");

  Param *pp=(Param *)param;
  mu[0]=pp->mu[0];
  mu[1]=pp->mu[1];
  Sigma[0][0]=pp->Sigma[0][0];
  Sigma[1][1]=pp->Sigma[1][1];
  Sigma[0][1]=pp->Sigma[0][1];
  Sigma[1][0]=pp->Sigma[1][0];
  X=pp->X;
  Y=pp->Y;
  normc=pp->normc;

  //Rprintf("mu1 %15g\n", mu[0]);
  //Rprintf("mu2 %15g\n", mu[1]);
  //Rprintf("Sigma00 %15g\n", Sigma[0][0]);
  //Rprintf("Sigma11 %15g\n", Sigma[1][1]);
  //Rprintf("Sigma01 %15g\n", Sigma[0][1]);
  imposs=0;


  for (ii=0; ii<n; ii++)
  //for (ii=0; ii<1; ii++)
    {
      inp=W1[ii]; //just for debugging
      imposs=0;
      W2[ii] = Y/(1-X)-X/(1-X)/(1+exp(-W1[ii]));
      if(W2[ii]>=1) imposs=1; //impossible pair of values
      else W2[ii] = log(W2[ii]/(1-W2[ii]));
    //}


  //this is constant over all W1's, so just calc once
  //double normc=getNormConst((void*)pp);
  //Rprintf("NC %15g\n", normc);

  //for (ii=0; ii<n; ii++)
  //for (ii=0; ii<1; ii++)
    //{
      vtemp[0] = W1[ii];
      vtemp[1] = W2[ii];
      if (imposs==1) W1[ii]=0;
      else W1[ii] = W1[ii]*dBVNtomo(vtemp, pp, 0);
      //Rprintf("W1Exp... %d %d %5g -> %5g via %5g imposs %d\n", ii, n, inp, W1[ii],dBVNtomo(vtemp, pp, 0),imposs);
    }
    //Rprintf("W1[0] %15g\n", W1[0]);
      //char ch;
      //scanf(" %c", &ch );
  free(W2);
  //free(mu);
  //FreeMatrix(Sigma,dim);
}


//E(W1*W1)

void W1W1Exp(double *W1, int n, void *param)
{
  int ii;
  int dim=2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double *W2;
  double X, Y, normc;
  double dtemp, rho;
  int imposs;
  double vtemp[2];
  W2 = (double *) malloc(1000*sizeof(double));
  if (W2==NULL) Rprintf("Malloc Error");

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
  normc=pp->normc;
  imposs=0;

  //Rprintf("mu1 %15g mu2 %15g\n", mu[0],mu[1]);
  //Rprintf("X %15g Y %15g\n", X,Y);
  //Rprintf("Sigma00 %15g\n", Sigma[0][0]);
  //Rprintf("Sigma11 %15g\n", Sigma[1][1]);
  //Rprintf("Sigma01 %15g\n", Sigma[0][1]);

  for (ii=0; ii<n; ii++)
    {
      imposs=0;
      W2[ii] = Y/(1-X)-X/(1-X)/(1+exp(-W1[ii]));
      if(W2[ii]>=1) imposs=1;
      else W2[ii] = log(W2[ii]/(1-W2[ii]));
    }

  //Rprintf("%d W2[0] %5g\n", imposs,W2[0]);
  //this is constant over all W1's, so just calc once
  //double normc=getNormConst((void*)pp);
  //Rprintf("%d W2[0] %5g NC %5g\n", imposs,W2[0],normc);

  for (ii=0; ii<n; ii++)
    {
      //Rprintf("W1W1 ... %d %d %5g %5g\n", ii,n,W1[ii],W2[ii]);
      vtemp[0] = W1[ii];
      vtemp[1] = W2[ii];
      if (imposs==1) W1[ii]=0;
      else W1[ii] = W1[ii] * W1[ii] * dBVNtomo(vtemp, pp,0);
    }
      //wait so I can see output
      //take out
      //Rprintf("W1sq[0] %15g\n", W1[0]);
      //char ch;
      //scanf(" %c", &ch );
  free(W2);
  //free(mu);
  //FreeMatrix(Sigma,dim);
}

/**E(W2*W2)*/
void W2W2Exp(double *W2, int n, void *param)
{
  int ii;
  int dim=2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double *W1;
  double X, Y,normc;
  double dtemp, rho;

  double vtemp[2];
  W1 = (double *) malloc(1000*sizeof(double));
  if (W1==NULL) Rprintf("Malloc Error");

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
  normc=pp->normc;

  for (ii=0; ii<n; ii++)
    {
      W1[ii] = Y/X-(1-X)/X/(1+exp(-W2[ii]));
      if(W1[ii]>=1)
        W1[ii]=99; //an otherwise impossible value
      else W1[ii] = log(W1[ii]/(1-W1[ii]));
    }

  //this is constant over all W1's, so just calc once
  //double normc=getNormConst((void*)pp);

  for (ii=0; ii<n; ii++)
    {
      vtemp[0] = W2[ii];
      vtemp[1] = W1[ii];
      if (W1[ii]==99) W2[ii]=0;
      else W2[ii] = W2[ii] * W2[ii] * dBVNtomo(vtemp, pp,0);
    }
  free(W1);
  //free(mu);
  //FreeMatrix(Sigma,dim);
}

/* integrate normalizing constant */
double getNormConst(void* pp) {

    Param *param=(Param *)pp;
    int inf=2;
    double bound=0.0;
    double epsabs=0.0000001, epsrel=0.0000001;
    double result=9999, anserr=9999;
    int limit=100;
    int last, neval, ier;
    int lenw=4*limit;
    int *iwork =   (int *) R_alloc(limit, sizeof(int));
    double *work=(double *)R_alloc(lenw, sizeof(double));
  Rdqagi(&NormConst, pp, &bound, &inf, &epsabs, &epsrel, &result,
	   &anserr, &neval, &ier, &limit, &lenw, &last, iwork, work);

    //Rprintf("Normc %15g\n", result);
    return result;

}
