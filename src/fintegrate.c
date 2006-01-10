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
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "sample.h"
#include "bayes.h"
#include "macros.h"
#include "fintegrate.h"
//#include  <gsl/gsl_integration.h>

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

//Bivariate normal distribution, with W2* dependent on W1* (which we are integrating over)
//see: http://mathworld.wolfram.com/BivariateNormalDistribution.html
void NormConstW1(double *W1, int n, void *param)
{
  int ii;
  double mu[2];
  double Sigma[2][2];
  double *W2;
  double X, Y, rho;
  double dtemp, inp;
  int imposs;

  W2 = Calloc(n,double);
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

  for (ii=0; ii<n; ii++)
    {
      imposs=0; inp=W1[ii];
      W2[ii]=getW2starFromW1star(X,Y,W1[ii],&imposs);

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
      //char ch;
      //scanf(" %c", &ch );
    }
  //Rprintf("W1[0] %15g\n", W1[0]);
  Free(W2);
}

//Bivariate normal distribution, with W1* dependent on W2* (which we are integrating over)
//see: http://mathworld.wolfram.com/BivariateNormalDistribution.html
void NormConstW2(double *W2, int n, void *param)
{
  int ii;
  double mu[2];
  double Sigma[2][2];
  double *W1;
  double X, Y, rho;
  double dtemp, inp;
  int imposs;

  W1 = Calloc(n,double);
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
  imposs=0;


  dtemp=1/(2*M_PI*sqrt(Sigma[0][0]*Sigma[1][1]*(1-rho*rho)));

  for (ii=0; ii<n; ii++)
    {
      imposs=0; inp=W2[ii];
      W1[ii]=getW1starFromW2star(X,Y,W2[ii],&imposs);

  //Rprintf("n: %d W2[0] %5g W1[0] %5g Y %5g X %5g imposs %d\n", n, W2[0], W1[0],Y,X,imposs);

  //for (ii=0; ii<n; ii++)
    //{
      /*W2[ii]=exp(-1/(2*(1-rho*rho))*
		 ((W2[ii]-mu[0])*(W2[ii]-mu[0])/Sigma[0][0]+
		  (W1[ii]-mu[1])*(W1[ii]-mu[1])/Sigma[1][1]-
		  2*rho*(W2[ii]-mu[0])*(W1[ii]-mu[0])
		  /sqrt(Sigma[0][0]*Sigma[1][1])))/dtemp; */
		  if (imposs==1) W2[ii]=0;
      else W2[ii]=exp(-1/(2*(1-rho*rho))*
		 ((W1[ii]-mu[0])*(W1[ii]-mu[0])/Sigma[0][0]+
		  (W2[ii]-mu[1])*(W2[ii]-mu[1])/Sigma[1][1]-
		  2*rho*(W1[ii]-mu[0])*(W2[ii]-mu[1])
		  /sqrt(Sigma[0][0]*Sigma[1][1])))*dtemp;

		  //Rprintf("Normc... %d %d %5g -> %5g with %5g imposs %d\n", ii, n, inp, W2[ii],dtemp,imposs);
      //char ch;
      //scanf(" %c", &ch );
    }
  //Rprintf("W2[0] %15g\n", W2[0]);
  free(W1);
}

//Bivariate normal distribution, with parameterization
//see: http://mathworld.wolfram.com/BivariateNormalDistribution.html
//see for param: http://www.math.uconn.edu/~binns/reviewII210.pdf
void NormConstT(double *t, int n, void *param)
{
  int ii;
  int dim=2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double *W1,*W1p,*W2,*W2p;
  double X, Y, rho;
  double dtemp, inp, pfact;
  int imposs;

  W1 = Calloc(n,double); if (W1==NULL) Rprintf("Malloc Error");
  W1p = Calloc(n,double); if (W1p==NULL) Rprintf("Malloc Error");
  W2 = Calloc(n,double); if (W2==NULL) Rprintf("Malloc Error");
  W2p = Calloc(n,double); if (W2p==NULL) Rprintf("Malloc Error");

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

  for (ii=0; ii<n; ii++)
    {
      imposs=0; inp=t[ii];
      W1[ii]=getW1starFromT(t[ii],pp,&imposs);
      if (!imposs) W2[ii]=getW2starFromT(t[ii],pp,&imposs);
		  if (imposs==1) t[ii]=0;
      else {
          W1p[ii]=getW1starPrimeFromT(t[ii],pp);
          W2p[ii]=getW2starPrimeFromT(t[ii],pp);
          pfact=sqrt(W1p[ii]*W1p[ii]+W2p[ii]*W2p[ii]);
          t[ii]=exp(-1/(2*(1-rho*rho))*
       ((W1[ii]-mu[0])*(W1[ii]-mu[0])/Sigma[0][0]+
        (W2[ii]-mu[1])*(W2[ii]-mu[1])/Sigma[1][1]-
        2*rho*(W1[ii]-mu[0])*(W2[ii]-mu[1])
        /sqrt(Sigma[0][0]*Sigma[1][1])))*dtemp*pfact;

        //Rprintf("Normc... %d %d %5g -> %5g with %5g imposs %d\n", ii, n, inp, W1[ii],dtemp,imposs);
        //char ch;
        //scanf(" %c", &ch );
      }
    }
  Free(W1);
  Free(W2);
  free(mu);
  FreeMatrix(Sigma,dim);
}


void SuffExp(double *t, int n, void *param)
{
  int ii,imposs,suff;
  int dim=2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double *W1,*W1p,*W2,*W2p;
  double dtemp, rho, inp,density,pfact,normc;

  double *vtemp=doubleArray(dim);

  W1 = Calloc(n,double); if (W1==NULL) Rprintf("Malloc Error");
  W1p = Calloc(n,double); if (W1p==NULL) Rprintf("Malloc Error");
  W2 = Calloc(n,double); if (W2==NULL) Rprintf("Malloc Error");
  W2p = Calloc(n,double); if (W2p==NULL) Rprintf("Malloc Error");

  Param *pp=(Param *)param;
  mu[0]=pp->mu[0];
  mu[1]=pp->mu[1];
  Sigma[0][0]=pp->Sigma[0][0];
  Sigma[1][1]=pp->Sigma[1][1];
  Sigma[0][1]=pp->Sigma[0][1];
  Sigma[1][0]=pp->Sigma[1][0];
  normc=pp->normcT;
  suff=pp->suff;
  imposs=0;


  for (ii=0; ii<n; ii++)
    {
     imposs=0; inp=t[ii];
      W1[ii]=getW1starFromT(t[ii],pp,&imposs);
      if (!imposs) W2[ii]=getW2starFromT(t[ii],pp,&imposs);
		  if (imposs==1) t[ii]=0;
      else {
          W1p[ii]=getW1starPrimeFromT(t[ii],pp);
          W2p[ii]=getW2starPrimeFromT(t[ii],pp);
          pfact=sqrt(W1p[ii]*W1p[ii]+W2p[ii]*W2p[ii]);
          vtemp[0] = W1[ii];
          vtemp[1] = W2[ii];
          density=dBVNtomo(vtemp, pp, 0,normc);
          t[ii] = density*pfact;
          if (suff==0) t[ii]=W1[ii]*t[ii];
          else if (suff==1) t[ii]=W2[ii]*t[ii];
          else if (suff==2) t[ii]=W1[ii]*W1[ii]*t[ii];
          else if (suff==3) t[ii]=W1[ii]*W2[ii]*t[ii];
          else if (suff==4) t[ii]=W2[ii]*W2[ii]*t[ii];
          else if (suff!=-1) Rprintf("Error Suff= %d",suff);
        }
    }
  free(W2);free(mu);free(vtemp);
  FreeMatrix(Sigma,dim);
}


void W1Exp(double *W1, int n, void *param)
{
  int ii,imposs;
  int dim=2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double *W2;
  double X, Y, normc;
  double dtemp, rho, inp,density;

  double *vtemp=doubleArray(dim);

  W2 = (double *) Calloc(n,double);
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
  normc=pp->normcW1;

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
      //W2[ii] = Y/(1-X)-X/(1-X)/(1+exp(-W1[ii]));
      //if(W2[ii]>=1) imposs=1; //impossible pair of values
      //else W2[ii] = log(W2[ii]/(1-W2[ii]));
      W2[ii]=getW2starFromW1star(X,Y,W1[ii],&imposs);
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
      else {
        density=dBVNtomo(vtemp, pp, 0,normc);
        W1[ii] = W1[ii]*density;
      }
      //Rprintf("W1Exp... %d %d %5g -> %5g via %5g imposs %d\n", ii, n, inp, W1[ii],density,imposs);
    }
    //Rprintf("W1[0] %15g\n", W1[0]);
      //char ch;
      //scanf(" %c", &ch );
  free(W2);
  free(mu);
    FreeMatrix(Sigma,dim);
}

/*E(W1)
  assumes bivariate normal distribution
  integral over W2
*/

void W1ExpW2(double *W2, int n, void *param)
{
  int ii,imposs;
  int dim=2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double *W1;
  double X, Y, normc;
  double dtemp, rho, inp, density;

  double *vtemp=doubleArray(dim);

  W1 = (double *) Calloc(n,double);
  if (W1==NULL) Rprintf("Malloc Error");

  Param *pp=(Param *)param;
  mu[0]=pp->mu[0];
  mu[1]=pp->mu[1];
  Sigma[0][0]=pp->Sigma[0][0];
  Sigma[1][1]=pp->Sigma[1][1];
  Sigma[0][1]=pp->Sigma[0][1];
  Sigma[1][0]=pp->Sigma[1][0];
  X=pp->X;
  Y=pp->Y;
  normc=pp->normcW2;

  //Rprintf("mu1 %15g\n", mu[0]);
  //Rprintf("mu2 %15g\n", mu[1]);
  //Rprintf("Sigma00 %15g\n", Sigma[0][0]);
  //Rprintf("Sigma11 %15g\n", Sigma[1][1]);
  //Rprintf("Sigma01 %15g\n", Sigma[0][1]);
  imposs=0;


  for (ii=0; ii<n; ii++)
  //for (ii=0; ii<1; ii++)
    {
      inp=W2[ii]; //just for debugging
      imposs=0;
      //W1[ii] = Y/(1-X)-X/(1-X)/(1+exp(-W2[ii]));
      //if(W1[ii]>=1) imposs=1; //impossible pair of values
      //else W1[ii] = log(W1[ii]/(1-W1[ii]));
      W1[ii]=getW1starFromW2star(X,Y,W2[ii],&imposs);
    //}


  //this is constant over all W2's, so just calc once
  //double normc=getNormConst((void*)pp);
  //Rprintf("NC %15g\n", normc);

  //for (ii=0; ii<n; ii++)
  //for (ii=0; ii<1; ii++)
    //{
      vtemp[0] = W1[ii];
      vtemp[1] = W2[ii];
      if (imposs==1) W2[ii]=0;
      else {
        density=dBVNtomo(vtemp, pp, 0,normc);
        W2[ii] = W1[ii]*density;
      }
      //Rprintf("W2Exp... %d %d %5g -> %5g via %5g imposs %d\n", ii, n, inp, W2[ii],density,imposs);
    }
    //Rprintf("W2[0] %15g\n", W2[0]);
      //char ch;
      //scanf(" %c", &ch );
  free(W1);
  free(mu);
  FreeMatrix(Sigma,dim);
}


/*E(W2)
  assumes bivariate normal distribution
  integral over W2
*/

void W2Exp(double *W2, int n, void *param)
{
  int ii,imposs;
  int dim=2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double *W1;
  double X, Y, normc;
  double dtemp, rho, inp, density;

  double *vtemp=doubleArray(dim);

  W1 = (double *) Calloc(n,double);
  if (W1==NULL) Rprintf("Malloc Error");

  Param *pp=(Param *)param;
  mu[0]=pp->mu[0];
  mu[1]=pp->mu[1];
  Sigma[0][0]=pp->Sigma[0][0];
  Sigma[1][1]=pp->Sigma[1][1];
  Sigma[0][1]=pp->Sigma[0][1];
  Sigma[1][0]=pp->Sigma[1][0];
  X=pp->X;
  Y=pp->Y;
  normc=pp->normcW2;

  //Rprintf("mu1 %15g\n", mu[0]);
  //Rprintf("mu2 %15g\n", mu[1]);
  //Rprintf("Sigma00 %15g\n", Sigma[0][0]);
  //Rprintf("Sigma11 %15g\n", Sigma[1][1]);
  //Rprintf("Sigma01 %15g\n", Sigma[0][1]);
  imposs=0;


  for (ii=0; ii<n; ii++)
  //for (ii=0; ii<1; ii++)
    {
      inp=W2[ii]; //just for debugging
      imposs=0;
      //W1[ii] = Y/(1-X)-X/(1-X)/(1+exp(-W2[ii]));
      //if(W1[ii]>=1) imposs=1; //impossible pair of values
      //else W1[ii] = log(W1[ii]/(1-W1[ii]));
      W1[ii]=getW1starFromW2star(X,Y,W2[ii],&imposs);
    //}


  //this is constant over all W2's, so just calc once
  //double normc=getNormConst((void*)pp);
  //Rprintf("NC %15g\n", normc);

  //for (ii=0; ii<n; ii++)
  //for (ii=0; ii<1; ii++)
    //{
      vtemp[0] = W1[ii];
      vtemp[1] = W2[ii];
      if (imposs==1) W2[ii]=0;
      else {
        density=dBVNtomo(vtemp, pp, 0,normc);
        W2[ii] = W2[ii]*density;
      }
      //Rprintf("W2Exp... %d %d %5g -> %5g via %5g imposs %d\n", ii, n, inp, W2[ii],density,imposs);
    }
    //Rprintf("W2[0] %15g\n", W2[0]);
      //char ch;
      //scanf(" %c", &ch );
  free(W1);
  free(mu);
  FreeMatrix(Sigma,dim);
}


/*E(W2)
  assumes bivariate normal distribution
  integral over W1
*/

void W2ExpW1(double *W1, int n, void *param)
{
  int ii,imposs;
  int dim=2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double *W2;
  double X, Y, normc;
  double dtemp, rho, inp,density;

  double *vtemp=doubleArray(dim);

  W2 = (double *) Calloc(n,double);
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
  normc=pp->normcW1;

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
      //W2[ii] = Y/(1-X)-X/(1-X)/(1+exp(-W1[ii]));
      //if(W2[ii]>=1) imposs=1; //impossible pair of values
      //else W2[ii] = log(W2[ii]/(1-W2[ii]));
      W2[ii]=getW2starFromW1star(X,Y,W1[ii],&imposs);
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
      else {
        density=dBVNtomo(vtemp, pp, 0,normc);
        W1[ii] = W2[ii]*density;
      }
      //Rprintf("W1Exp... %d %d %5g -> %5g via %5g imposs %d\n", ii, n, inp, W1[ii],density,imposs);
    }
    //Rprintf("W1[0] %15g\n", W1[0]);
      //char ch;
      //scanf(" %c", &ch );
  free(W2);
  free(mu);
  FreeMatrix(Sigma,dim);
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
  double dtemp, rho, inp, density;
  int imposs;
  double vtemp[2];
  W2 = (double *) Calloc(n,double);
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
  normc=pp->normcW1;
  imposs=0;

  //Rprintf("mu1 %15g mu2 %15g\n", mu[0],mu[1]);
  //Rprintf("X %15g Y %15g\n", X,Y);
  //Rprintf("Sigma00 %15g\n", Sigma[0][0]);
  //Rprintf("Sigma11 %15g\n", Sigma[1][1]);
  //Rprintf("Sigma01 %15g\n", Sigma[0][1]);

  for (ii=0; ii<n; ii++)
    {
      imposs=0; inp=W1[ii];
      //W2[ii] = Y/(1-X)-X/(1-X)/(1+exp(-W1[ii]));
      //if(W2[ii]>=1) imposs=1;
      //else W2[ii] = log(W2[ii]/(1-W2[ii]));
      W2[ii]=getW2starFromW1star(X,Y,W1[ii],&imposs);
      //Rprintf(" %5g %5g\n", inp,W1[ii]);
    //}

  //Rprintf("%d W2[0] %5g\n", imposs,W2[0]);
  //this is constant over all W1's, so just calc once
  //double normc=getNormConst((void*)pp);
  //Rprintf("%d W2[0] %5g NC %5g\n", imposs,W2[0],normc);

  //for (ii=0; ii<n; ii++)
    //{
      vtemp[0] = W1[ii];
      vtemp[1] = W2[ii];
      if (imposs==1) W1[ii]=0;
      else {
        density=dBVNtomo(vtemp, pp,0,normc);
        //if (rho==0) W1[ii] = W1[ii] * W1[ii] * density;
        //else W1[ii] = density;
        W1[ii] = W1[ii] * W1[ii] * density;
      }
      //Rprintf("W1W1 ... %d %d %5g -> %5g via %5g imp:%d\n", ii,n,inp,W1[ii],density,imposs);
    }
      //wait so I can see output
      //take out
      //char ch;
      //scanf(" %c", &ch );
  free(W2);
  free(mu);
  FreeMatrix(Sigma,dim);
}

/*E(W1^2)
  assumes bivariate normal distribution
  integral over W2
*/

void W1W1ExpW2(double *W2, int n, void *param)
{
  int ii,imposs;
  int dim=2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double *W1;
  double X, Y, normc;
  double dtemp, rho, inp, density;

  double *vtemp=doubleArray(dim);

  W1 = (double *) Calloc(n,double);
  if (W1==NULL) Rprintf("Malloc Error");

  Param *pp=(Param *)param;
  mu[0]=pp->mu[0];
  mu[1]=pp->mu[1];
  Sigma[0][0]=pp->Sigma[0][0];
  Sigma[1][1]=pp->Sigma[1][1];
  Sigma[0][1]=pp->Sigma[0][1];
  Sigma[1][0]=pp->Sigma[1][0];
  X=pp->X;
  Y=pp->Y;
  normc=pp->normcW2;

  //Rprintf("mu1 %15g\n", mu[0]);
  //Rprintf("mu2 %15g\n", mu[1]);
  //Rprintf("Sigma00 %15g\n", Sigma[0][0]);
  //Rprintf("Sigma11 %15g\n", Sigma[1][1]);
  //Rprintf("Sigma01 %15g\n", Sigma[0][1]);
  imposs=0;


  for (ii=0; ii<n; ii++)
  //for (ii=0; ii<1; ii++)
    {
      inp=W2[ii]; //just for debugging
      imposs=0;
      //W1[ii] = Y/(1-X)-X/(1-X)/(1+exp(-W2[ii]));
      //if(W1[ii]>=1) imposs=1; //impossible pair of values
      //else W1[ii] = log(W1[ii]/(1-W1[ii]));
      W1[ii]=getW1starFromW2star(X,Y,W2[ii],&imposs);
    //}


  //this is constant over all W2's, so just calc once
  //double normc=getNormConst((void*)pp);
  //Rprintf("NC %15g\n", normc);

  //for (ii=0; ii<n; ii++)
  //for (ii=0; ii<1; ii++)
    //{
      vtemp[0] = W1[ii];
      vtemp[1] = W2[ii];
      if (imposs==1) W2[ii]=0;
      else {
        density=dBVNtomo(vtemp, pp, 0,normc);
        W2[ii] = W1[ii]*W1[ii]*density;
      }
      //Rprintf("W2Exp... %d %d %5g -> %5g via %5g imposs %d\n", ii, n, inp, W2[ii],density,imposs);
    }
    //Rprintf("W2[0] %15g\n", W2[0]);
      //char ch;
      //scanf(" %c", &ch );
  free(W1);
  free(mu);
  FreeMatrix(Sigma,dim);
}


/**E(W2*W2)*/
void W2W2Exp(double *W2, int n, void *param)
{
  int ii, imposs;
  int dim=2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double *W1;
  double X, Y,normc,inp;
  double dtemp, rho,density;

  double vtemp[2];
  W1 = (double *) Calloc(n,double);
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
  normc=pp->normcW2;
  imposs=0;

  for (ii=0; ii<n; ii++)
    {
      //W1[ii] = Y/X-(1-X)/X/(1+exp(-W2[ii]));
      //if(W1[ii]>=1)
      //  W1[ii]=99; //an otherwise impossible value
      //else W1[ii] = log(W1[ii]/(1-W1[ii]));
      imposs=0;inp=W2[ii];
      W1[ii]=getW1starFromW2star(X,Y,W2[ii],&imposs);
    //}

  //this is constant over all W1's, so just calc once
  //double normc=getNormConst((void*)pp);

  //for (ii=0; ii<n; ii++)
    //{
      vtemp[0] = W1[ii];
      vtemp[1] = W2[ii];
      if (imposs==1) W2[ii]=0;
      else {
        density=dBVNtomo(vtemp, pp,0,normc);
        //if (rho==0) W2[ii] = W2[ii] * W2[ii] * density;
        //else W2[ii] = density;
        W2[ii] = W2[ii] * W2[ii] * density;
      }
      //if (rho!=0) Rprintf("W2W2 ... %d %d %5g -> %5g via %5g %5g imp:%d\n", ii,n,inp,W2[ii],(imposs==1) ? 0 : density,rho,imposs);
    }
      //if (rho!=0) {char ch;
      //scanf(" %c", &ch );}
  free(W1);
  free(mu);
  FreeMatrix(Sigma,dim);
}


//E(W2*W2)
//integrate over W1
void W2W2ExpW1(double *W1, int n, void *param)
{
  int ii;
  int dim=2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double *W2;
  double X, Y, normc;
  double dtemp, rho, inp, density;
  int imposs;
  double vtemp[2];
  W2 = (double *) Calloc(n,double);
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
  normc=pp->normcW1;
  imposs=0;

  //Rprintf("mu1 %15g mu2 %15g\n", mu[0],mu[1]);
  //Rprintf("X %15g Y %15g\n", X,Y);
  //Rprintf("Sigma00 %15g\n", Sigma[0][0]);
  //Rprintf("Sigma11 %15g\n", Sigma[1][1]);
  //Rprintf("Sigma01 %15g\n", Sigma[0][1]);

  for (ii=0; ii<n; ii++)
    {
      imposs=0; inp=W1[ii];
      //W2[ii] = Y/(1-X)-X/(1-X)/(1+exp(-W1[ii]));
      //if(W2[ii]>=1) imposs=1;
      //else W2[ii] = log(W2[ii]/(1-W2[ii]));
      W2[ii]=getW2starFromW1star(X,Y,W1[ii],&imposs);
      //Rprintf(" %5g %5g\n", inp,W1[ii]);
    //}

  //Rprintf("%d W2[0] %5g\n", imposs,W2[0]);
  //this is constant over all W1's, so just calc once
  //double normc=getNormConst((void*)pp);
  //Rprintf("%d W2[0] %5g NC %5g\n", imposs,W2[0],normc);

  //for (ii=0; ii<n; ii++)
    //{
      vtemp[0] = W1[ii];
      vtemp[1] = W2[ii];
      if (imposs==1) W1[ii]=0;
      else {
        density=dBVNtomo(vtemp, pp,0,normc);
        //if (rho==0) W1[ii] = W1[ii] * W1[ii] * density;
        //else W1[ii] = density;
        W1[ii] = W2[ii] * W2[ii] * density;
      }
      //Rprintf("W1W1 ... %d %d %5g -> %5g via %5g imp:%d\n", ii,n,inp,W1[ii],density,imposs);
    }
      //wait so I can see output
      //take out
      //char ch;
      //scanf(" %c", &ch );
  free(W2);
  free(mu);
  FreeMatrix(Sigma,dim);
}


//E(W1*W2)
void W1W2Exp(double *W1, int n, void *param)
{
  int ii;
  int dim=2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double *W2;
  double X, Y, normc;
  double dtemp, rho, inp, density;
  int imposs;
  double vtemp[2];
  W2 = (double *) Calloc(n,double);
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
  normc=pp->normcW1;
  imposs=0;

  //Rprintf("mu1 %15g mu2 %15g\n", mu[0],mu[1]);
  //Rprintf("X %15g Y %15g\n", X,Y);
  //Rprintf("Sigma00 %15g\n", Sigma[0][0]);
  //Rprintf("Sigma11 %15g\n", Sigma[1][1]);
  //Rprintf("Sigma01 %15g\n", Sigma[0][1]);

  for (ii=0; ii<n; ii++)
    {
      imposs=0; inp=W1[ii];
      //W2[ii] = Y/(1-X)-X/(1-X)/(1+exp(-W1[ii]));
      //if(W2[ii]>=1) imposs=1;
      //else W2[ii] = log(W2[ii]/(1-W2[ii]));
      W2[ii]=getW2starFromW1star(X,Y,W1[ii],&imposs);
      //Rprintf(" %5g %5g\n", inp,W1[ii]);
    //}

  //Rprintf("%d W2[0] %5g\n", imposs,W2[0]);
  //this is constant over all W1's, so just calc once
  //double normc=getNormConst((void*)pp);
  //Rprintf("%d W2[0] %5g NC %5g\n", imposs,W2[0],normc);

  //for (ii=0; ii<n; ii++)
    //{
      vtemp[0] = W1[ii];
      vtemp[1] = W2[ii];
      if (imposs==1) W1[ii]=0;
      else {
        density=dBVNtomo(vtemp, pp,0,normc);
        W1[ii] = W1[ii] * W2[ii] * density;
        //W1[ii] = density;
      }
      //Rprintf("W1W1 ... %d %d %5g -> %5g via %5g imp:%d\n", ii,n,inp,W1[ii],density,imposs);
    }
      //wait so I can see output
      //take out
      //char ch;
      //scanf(" %c", &ch );
  free(W2);
  free(mu);
  FreeMatrix(Sigma,dim);
}


//E(W2*W1)
//for checking purposes
void W2W1Exp(double *W2, int n, void *param)
{
  int ii;
  int dim=2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double *W1;
  double X, Y, normc;
  double dtemp, rho, inp, density;
  int imposs;
  double vtemp[2];
  W1 = (double *) Calloc(n,double);
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
  normc=pp->normcW2;
  imposs=0;

  //Rprintf("mu1 %15g mu2 %15g\n", mu[0],mu[1]);
  //Rprintf("X %15g Y %15g\n", X,Y);
  //Rprintf("Sigma00 %15g\n", Sigma[0][0]);
  //Rprintf("Sigma11 %15g\n", Sigma[1][1]);
  //Rprintf("Sigma01 %15g\n", Sigma[0][1]);

  for (ii=0; ii<n; ii++)
    {
      imposs=0; inp=W2[ii];
      //W1[ii] = Y/(1-X)-X/(1-X)/(1+exp(-W2[ii]));
      //if(W1[ii]>=1) imposs=1;
      //else W1[ii] = log(W1[ii]/(1-W1[ii]));
      W1[ii]=getW1starFromW2star(X,Y,W2[ii],&imposs);
      //Rprintf(" %5g %5g\n", inp,W2[ii]);
    //}

  //Rprintf("%d W1[0] %5g\n", imposs,W1[0]);
  //this is constant over all W2's, so just calc once
  //double normc=getNormConst((void*)pp);
  //Rprintf("%d W1[0] %5g NC %5g\n", imposs,W1[0],normc);

  //for (ii=0; ii<n; ii++)
    //{
      vtemp[0] = W1[ii];
      vtemp[1] = W2[ii];
      if (imposs==1) W2[ii]=0;
      else {
        density=dBVNtomo(vtemp, pp,0,normc);
        W2[ii] = W2[ii] * W1[ii] * density;
        //W2[ii] = density;
      }
      //Rprintf("W2W2 ... %d %d %5g -> %5g via %5g imp:%d\n", ii,n,inp,W2[ii],density,imposs);
    }
      //wait so I can see output
      //take out
      //char ch;
      //scanf(" %c", &ch );
  free(W1);
  free(mu);
  FreeMatrix(Sigma,dim);
}


/* integrate normalizing constant */
void setNormConst(Param* param) {
    /*int inf;
    double ub,lb;
    lb=param->W1_lb;
    ub=param->W1_ub;
    inf=param->W1_inf;
    param->normcW1=numIntegration(&NormConstW1,(void*)param,inf,lb,ub);

    lb=param->W2_lb;
    ub=param->W2_ub;
    inf=param->W2_inf;
    param->normcW2=numIntegration(&NormConstW2,(void*)param,inf,lb,ub);*/

    param->normcT=paramIntegration(&NormConstT,(void*)param);

}

//Finds W2star, given the equation
//Y=XW1 + (1-X)W2 and the Wistar=logit(Wi)
//imposs is set to 1 if the equation cannot be satisfied
double getW2starFromW1star(double X, double Y, double W1, int* imposs) {
      if (W1>30) W1=1; //prevent overflow or underflow
      else W1=1/(1+exp(-1*W1));
      double W2=Y/(1-X)-X*W1/(1-X);
      if(W2>=1 || W2<=0) *imposs=1; //impossible pair of values
      else W2=log(W2/(1-W2));
      return W2;
}

double getW1starFromW2star(double X, double Y, double W2, int* imposs) {
      if (W2>30) W2=1; //prevent overflow or underflow
      else W2=1/(1+exp(-1*W2));
      double W1=(Y-(1-X)*W2)/X;
      //Rprintf(" %5g %5g %5g %5g\n", X,Y,W2, W1);
      if(W1>=1 || W1<=0) *imposs=1; //impossible pair of values
      else W1=log(W1/(1-W1));
      return W1;
}

//W1star(t)
//W1(t)=(W1_ub - W1_lb)*t + W1_lb
double getW1starFromT(double t, Param* param, int* imposs) {
    double W1=(param->W1_ub - param->W1_lb)*t + param->W1_lb;
    if (W1==1 || W1==0) *imposs=1;
    else W1=log(W1/(1-W1));
    return W1;
}
//W2star(t)
//W2(t)=(W2_lb - W2_ub)*t + W2_lb
double getW2starFromT(double t, Param* param, int* imposs) {
    double W2=(param->W2_lb - param->W2_ub)*t + param->W2_ub;
    if (W2==1 || W2==0) *imposs=1;
    else W2=log(W2/(1-W2));
    return W2;
}
//W1star'(t)
double getW1starPrimeFromT(double t, Param* param) {
    double c=(param->W1_ub - param->W1_lb);
    double W1=c*t + param->W1_lb;
    W1=((1-W1)/W1)*((c/(1-W1))-c*W1);
    return W1;
}
//W2star'(t)
double getW2starPrimeFromT(double t, Param* param) {
    double c=(param->W2_lb - param->W2_ub);
    double W2=c*t + param->W2_ub;
    W2=((1-W2)/W2)*((c/(1-W2))-c*W2);
    return W2;
}

//parameterized integration: bounds always from 0,1
double paramIntegration(integr_fn f, void *ex) {
    double bound;
  double epsabs=0.000000001, epsrel=0.000000001;
  double result=9999, anserr=9999;
  int limit=100;
  int last, neval, ier;
  int lenw=5*limit;
  int *iwork=(int *) Calloc(limit, int);
  double *work=(double *)Calloc(lenw, double);
  double lb=0; double ub=1;
    Rdqags(f, ex, &lb, &ub, &epsabs, &epsrel, &result,
      &anserr, &neval, &ier, &limit, &lenw, &last, iwork, work);

  Free(iwork);
  Free(work);
  if (ier==0) return result;
  else {
    Param* p = (Param*) ex;
    Rprintf("Integration error %d: X %5g Y %5g [%5g,%5g] -> %5g +- %5g\n",ier,p->X,p->Y,p->W1_lb,p->W1_ub,result,anserr);
    return result;
  }

}


//numerical integration with R's function
//inf: 0->(lb,ub), -1->(-inf,ub), 1->(lb,inf), 2->(-inf,inf)
double numIntegration(integr_fn f, void *ex, int inf, double lb, double ub) {

  double bound;
  double epsabs=0.000000001, epsrel=0.000000001;
  double result=9999, anserr=9999;
  int limit=100;
  int last, neval, ier;
  int lenw=5*limit;
  int *iwork=(int *) Calloc(limit, int);
  double *work=(double *)Calloc(lenw, double);
  if (inf==0) {
    Rdqags(f, ex, &lb, &ub, &epsabs, &epsrel, &result,
      &anserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
  }
  else {
    bound=(inf == 1) ? lb : ub;
    Rdqagi(f, ex, &bound, &inf, &epsabs, &epsrel, &result,
      &anserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
  }
  Free(iwork);
  Free(work);
  if (ier==0) return result;
  else {
    Param* p = (Param*) ex;
    Rprintf("Integration error %d: X %5g Y %5g [%5g,%5g]_%d -> %5g +- %5g\n",ier,p->X,p->Y,lb,ub,inf,result,anserr);
    return result;
  }
}

/*
//numerical integration with GNU's function
double numIntegration2(gsl_fn f, void *ex) {

    double result, error;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    gsl_function F;
    F.function = f;
    F.params = ex;

  gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,w, &result, &error);

  return result;
}
*/

void setBounds(Param* param) {
  double X,Y,w1_lb,w1_ub,w2_lb,w2_ub;
  int w1_inf,w2_inf;
  double tol0=0.0001;
  double tol1=0.9999;
  X=param->X;
  Y=param->Y;

  //find bounds for W1
  w1_ub=(Y-(1-X)*0)/X; //W2=0
  if (w1_ub>tol1) w1_ub=1;
  w1_lb=(Y-(1-X)*1)/X; //W2=1
  if (w1_lb<tol0) w1_lb=0;

  //find bounds for W2
  w2_ub=Y/(1-X)-X*0/(1-X); //W1=0
  if (w2_ub>tol1) w2_ub=1;
  w2_lb=Y/(1-X)-X*1/(1-X); //W1=1
  if (w2_lb<tol0) w2_lb=0;

//Rprintf("GB0: %5g %5g %5g %5g %5g %5g \n",X,Y,w1_lb, w1_ub, w2_lb, w2_ub);
/*
  if (w1_lb==0 && w1_ub==1) w1_inf=2;
  else if (w1_lb==0) w1_inf=-1;
  else if (w1_ub==1) w1_inf=1;
  else w1_inf=0;
  w1_lb=log(w1_lb/(1-w1_lb));
  w1_ub=log(w1_ub/(1-w1_ub));

  if (w2_lb==0 && w2_ub==1) w2_inf=2;
  else if (w2_lb==0) w2_inf=-1;
  else if (w2_ub==1) w2_inf=1;
  else w2_inf=0;
  w2_lb=log(w2_lb/(1-w2_lb));
  w2_ub=log(w2_ub/(1-w2_ub));

  //Rprintf("GB: %5g %5g %5g %5g %d %5g %5g %d\n",X,Y,w1_lb, w1_ub, w1_inf, w2_lb, w2_ub, w2_inf);
*/
  param->W1_lb=w1_lb;
  param->W1_ub=w1_ub;
  param->W2_lb=w2_lb;
  param->W2_ub=w2_ub;
  param->W1_inf=w1_inf;
  param->W2_inf=w2_inf;

}
