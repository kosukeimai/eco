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
  free(W2);
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

  W1 = (double *) malloc(n*sizeof(double));
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


/*
//for GSL integration
double NormConst2(double W1, void *param)
{
  double mu[2];
  double Sigma[2][2];
  double W2;
  double X, Y, rho;
  double dtemp, inp;
  int imposs;

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
  imposs=0;
  W2=getW2starFromW1star(X,Y,W1,&imposs);

  if (imposs==1) W1=0;
  else W1=exp(-1/(2*(1-rho*rho))*
		((W1-mu[0])*(W1-mu[0])/Sigma[0][0]+
		 (W2-mu[1])*(W2-mu[1])/Sigma[1][1]-
		  2*rho*(W1-mu[0])*(W2-mu[1])
		  /sqrt(Sigma[0][0]*Sigma[1][1])))*dtemp;
  return W1;
}
*/

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
  double dtemp, rho, inp,density;

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
  //free(mu);
  //FreeMatrix(Sigma,dim);
}



/*E(W2)
  assumes bivariate normal distribution
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

  W1 = (double *) malloc(n*sizeof(double));
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
  double dtemp, rho, inp, density;
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
  //free(mu);
  //FreeMatrix(Sigma,dim);
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
  //free(mu);
  //FreeMatrix(Sigma,dim);
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
  //free(mu);
  //FreeMatrix(Sigma,dim);
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
  //free(mu);
  //FreeMatrix(Sigma,dim);
}


/* integrate normalizing constant */
void setNormConst(Param* param) {
    int inf;
    double ub,lb;
    lb=param->W1_lb;
    ub=param->W1_ub;
    inf=param->W1_inf;
    param->normcW1=numIntegration(&NormConstW1,(void*)param,inf,lb,ub);

    lb=param->W2_lb;
    ub=param->W2_ub;
    inf=param->W2_inf;
    param->normcW2=numIntegration(&NormConstW2,(void*)param,inf,lb,ub);

}

//Finds W2star, given the equation
//Y=XW1 + (1-X)W2 and the Wistar=logit(Wi)
//imposs is set to 1 if the equation cannot be satisfied
double getW2starFromW1star(double X, double Y, double W1, int* imposs) {
      if (W1>30) W1=1;
      else W1=1/(1+exp(-1*W1));
      double W2=Y/(1-X)-X*W1/(1-X);
      if(W2>=1 || W2<=0) *imposs=1; //impossible pair of values
      else W2=log(W2/(1-W2));
      return W2;
}

double getW1starFromW2star(double X, double Y, double W2, int* imposs) {
      if (W2>30) W2=1;
      else W2=1/(1+exp(-1*W2));
      double W1=(Y-(1-X)*W2)/X;
      //Rprintf(" %5g %5g %5g %5g\n", X,Y,W2, W1);
      if(W1>=1 || W1<=0) *imposs=1; //impossible pair of values
      else W1=log(W1/(1-W1));
      return W1;
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
  int *iwork=(int *) R_alloc(limit, sizeof(int));
  double *work=(double *)R_alloc(lenw, sizeof(double));
  if (inf==0) {
    Rdqags(f, ex, &lb, &ub, &epsabs, &epsrel, &result,
      &anserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
  }
  else {
    bound=(inf == 1) ? lb : ub;
    Rdqagi(f, ex, &bound, &inf, &epsabs, &epsrel, &result,
      &anserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
  }
  if (ier==0) return result;
  else {
    Rprintf("Integration error: %d\n",ier);
    return 9999;
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

  param->W1_lb=w1_lb;
  param->W1_ub=w1_ub;
  param->W2_lb=w2_lb;
  param->W2_ub=w2_ub;
  param->W1_inf=w1_inf;
  param->W2_inf=w2_inf;

}
