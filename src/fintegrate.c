/******************************************************************
  This file is a part of eco: R Package for Estimating Fitting
  Bayesian Models of Ecological Inference for 2X2 tables
  by Kosuke Imai, Ying Lu, and Aaron Strauss
  Copyright: GPL version 2 or later.
*******************************************************************/

#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/PrtUtil.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "sample.h"
#include "bayes.h"
#include "macros.h"
#include "fintegrate.h"
//#include  <gsl/gsl_integration.h>

/**
 * Bivariate normal distribution, with parameterization
 * see: http://mathworld.wolfram.com/BivariateNormalDistribution.html
 * see for param: http://www.math.uconn.edu/~binns/reviewII210.pdf
 */
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

  W1 = doubleArray(n);
  W1p = doubleArray(n);
  W2 = doubleArray(n);
  W2p = doubleArray(n);

  Param *pp=(Param *)param;
  mu[0]= pp->caseP.mu[0];
  mu[1]= pp->caseP.mu[1];
  Sigma[0][0]=pp->setP->Sigma[0][0];
  Sigma[1][1]=pp->setP->Sigma[1][1];
  Sigma[0][1]=pp->setP->Sigma[0][1];
  Sigma[1][0]=pp->setP->Sigma[1][0];
  rho=Sigma[0][1]/sqrt(Sigma[0][0]*Sigma[1][1]);
  //Rprintf("TESTING: %4g %4g %4g %4g", pp->caseP.mu[0], pp->caseP.mu[1], pp->setP->Sigma[0][0],pp->setP->Sigma[0][1]);
  X=pp->caseP.X;
  Y=pp->caseP.Y;
  imposs=0;

  dtemp=1/(2*M_PI*sqrt(Sigma[0][0]*Sigma[1][1]*(1-rho*rho)));

  for (ii=0; ii<n; ii++) {
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
  //if (pp->setP->weirdness)
   //   Rprintf("Normc... %d %d %5g -> %5g %5g => %5g with %5g imposs %d\n", ii, n, inp, W1[ii], W2[ii],t[ii],pfact,imposs);
      //char ch;
      //scanf(" %c", &ch );
    }
  }
  Free(W1);
  Free(W1p);
  Free(W2);
  Free(W2p);
  Free(mu);
  FreeMatrix(Sigma,dim);
}

/**
 * Integrand for computing sufficient statistic
 * Which statistic to estimate depends on param->suff (see macros.h)
 */
void SuffExp(double *t, int n, void *param)
{
  int ii,imposs,i,j;
  sufficient_stat suff;
  Param *pp=(Param *)param;
  int dim = (pp->setP->ncar==1) ? 3 : 2;
  double *mu=doubleArray(dim);
  double **Sigma=doubleMatrix(dim,dim);
  double **InvSigma=doubleMatrix(dim,dim);/* inverse covariance matrix*/
  //double Sigma[dim][dim];
  //double InvSigma[dim][dim];
  double *W1,*W1p,*W2,*W2p,*vtemp;
  double inp,density,pfact,normc;

  vtemp=doubleArray(dim);
  W1 = doubleArray(n);
  W1p = doubleArray(n);
  W2 = doubleArray(n);
  W2p = doubleArray(n);
  mu[0]= pp->caseP.mu[0];
  mu[1]= pp->caseP.mu[1];
  for(i=0;i<dim;i++) {
    for(j=0;j<dim;j++) {
      if (dim==3) {
        Sigma[i][j]=pp->setP->Sigma3[i][j];
        InvSigma[i][j]=pp->setP->InvSigma3[i][j];
      }
      else {
        Sigma[i][j]=pp->setP->Sigma[i][j];
        InvSigma[i][j]=pp->setP->InvSigma[i][j];
      }
    }
  }
  normc=pp->caseP.normcT;
  suff=pp->caseP.suff;
  imposs=0;

  for (ii=0; ii<n; ii++) {
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
      if (suff==SS_W1star) t[ii]=W1[ii]*t[ii];
      else if (suff==SS_W2star) t[ii]=W2[ii]*t[ii];
      else if (suff==SS_W1star2) t[ii]=W1[ii]*W1[ii]*t[ii];
      else if (suff==SS_W1W2star) t[ii]=W1[ii]*W2[ii]*t[ii];
      else if (suff==SS_W2star2) t[ii]=W2[ii]*W2[ii]*t[ii];
      else if (suff==SS_W1) t[ii]=invLogit(W1[ii])*t[ii];
      else if (suff==SS_W2) t[ii]=invLogit(W2[ii])*t[ii];
      else if (suff==SS_Loglik) {
        if (dim == 3) {
          //if(pp->setP->verbose>=2 && dim==3) Rprintf("InvSigma loglik: %5g %5g %5g %5g %5g %5g\n",InvSigma[0][0],InvSigma[0][1],InvSigma[1][0],InvSigma[1][1],InvSigma[1][2],InvSigma[2][2]);
          vtemp[2]=logit(pp->caseP.X,"log-likelihood");
          mu[0]=pp->setP->pdTheta[1];
          mu[1]=pp->setP->pdTheta[2];
          mu[2]=pp->setP->pdTheta[0];
        }
        t[ii]=dMVN(vtemp,mu,InvSigma,dim,0)*pfact;
        //t[ii]=dMVN3(vtemp,mu,(double*)(&(InvSigma[0][0])),dim,0)*pfact;
      }
      else if (suff!=SS_Test) Rprintf("Error Suff= %d",suff);
    }
  }
  Free(W1);Free(W1p);Free(W2);Free(W2p);Free(mu);Free(vtemp);
  FreeMatrix(Sigma,dim); FreeMatrix(InvSigma,dim);
}


/**
 * Returns the log likelihood of a particular case (i.e, record, datapoint)
 */
double getLogLikelihood(Param* param) {
  if (param->caseP.dataType==DPT_General  && !(param->caseP.Y>=.990 || param->caseP.Y<=.010)) {
    //non-survey data: do integration to find likelihood
    param->caseP.suff=SS_Loglik;
    return log(paramIntegration(&SuffExp,(void*)param));


  } else if (param->caseP.dataType==DPT_Homog_X1 || param->caseP.dataType==DPT_Homog_X0) {
      //Homogenenous data: just do normal likelihood on one dimension
      double lik,sigma2,val,mu;
      val = (param->caseP.dataType==DPT_Homog_X1) ? param->caseP.Wstar[0] : param->caseP.Wstar[1];
      if (!param->setP->ncar) {
        mu = (param->caseP.dataType==DPT_Homog_X1) ? param->setP->pdTheta[0] : param->setP->pdTheta[1];
        sigma2 = (param->caseP.dataType==DPT_Homog_X1) ? param->setP->pdTheta[2] : param->setP->pdTheta[3];
      } else {
        mu = (param->caseP.dataType==DPT_Homog_X1) ? param->setP->pdTheta[1] : param->setP->pdTheta[2];
        sigma2 = (param->caseP.dataType==DPT_Homog_X1) ? param->setP->pdTheta[4] : param->setP->pdTheta[5];
      }
      lik=(1/(sqrt(2*M_PI*sigma2)))*exp(-(.5/sigma2)*(val - mu)*(val - mu));
      //return log(lik);
      return 0; //fix later

  } else if (param->caseP.dataType==DPT_Survey || (param->caseP.Y>=.990 || param->caseP.Y<=.010)) {
    //Survey data (or v tight bounds): multi-variate normal
    int dim=param->setP->ncar ? 3 : 2;
    double *mu=doubleArray(dim);
    double *vtemp=doubleArray(dim);
    double **InvSig=doubleMatrix(dim,dim);/* inverse covariance matrix*/
    int i,j;
    for(i=0;i<dim;i++) {
      for(j=0;j<dim;j++) {
        if (dim==3) {
          InvSig[i][j]=param->setP->InvSigma3[i][j];
        }
        else {
          InvSig[i][j]=param->setP->InvSigma[i][j];
        }
      }
    }
    double loglik;
    vtemp[0] = param->caseP.Wstar[0];
    vtemp[1] = param->caseP.Wstar[1];
    mu[0]= param->caseP.mu[0];
    mu[1]= param->caseP.mu[1];
    if (param->setP->ncar) {
      vtemp[2]=logit(param->caseP.X,"log-likelihood survey");
      mu[0]=param->setP->pdTheta[1];
      mu[1]=param->setP->pdTheta[2];
      mu[2]=param->setP->pdTheta[0];
      loglik=dMVN(vtemp,mu,InvSig,dim,1);
    }
    else {
      loglik=dMVN(vtemp,mu,InvSig,dim,1);
    }
    Free(mu); Free(vtemp); FreeMatrix(InvSig,dim);
    return loglik;
  }
  else { //Unknown type
    Rprintf("Error; unkown type: %d\n",param->caseP.dataType);
    return 0;
  }
}

/**
 **********
 * Line integral helper function
 **********
 */

/**
 * Returns W2star from W1star, given the following equalities
 * Y=XW1 + (1-X)W2 and the Wi-star=logit(Wi)
 * mutation: imposs is set to 1 if the equation cannot be satisfied
 */
double getW2starFromW1star(double X, double Y, double W1star, int* imposs) {
  double W1;
  if (W1star>30) W1=1; //prevent overflow or underflow
  else W1=1/(1+exp(-1*W1star));
  double W2=Y/(1-X)-X*W1/(1-X);

  if(W2>=1 || W2<=0) *imposs=1; //impossible pair of values
  else W2=log(W2/(1-W2));
  return W2;
}

/**
 * Returns W1star from W2star, given the following equalities
 * Y=XW1 + (1-X)W2 and the Wi-star=logit(Wi)
 * mutation: imposs is set to 1 if the equation cannot be satisfied
 */
double getW1starFromW2star(double X, double Y, double W2star, int* imposs) {
  double W2;
  if (W2star>30) W2=1; //prevent overflow or underflow
  else W2=1/(1+exp(-1*W2star));
  double W1=(Y-(1-X)*W2)/X;

  if(W1>=1 || W1<=0) *imposs=1; //impossible pair of values
  else W1=log(W1/(1-W1));
  return W1;
}

/**
 * Returns W1 from W2, X, and Y given
 * Y=XW1 + (1-X)W2
 */
double getW1FromW2(double X, double Y, double W2) {
  return (Y-(1-X)*W2)/X;
}


 /**
 * W1star(t)
 * W1(t)=(W1_ub - W1_lb)*t + W1_lb
 * mutates impossible to true if W1 is non-finite at t
 */
double getW1starFromT(double t, Param* param, int* imposs) {
  double W1=(param->caseP.Wbounds[0][1] - param->caseP.Wbounds[0][0])*t + param->caseP.Wbounds[0][0];
  if (W1==1 || W1==0) *imposs=1;
  else W1=log(W1/(1-W1));
  return W1;
}

/**
 * W2star(t)
 * W2(t)=(W2_lb - W2_ub)*t + W2_lb
 */
double getW2starFromT(double t, Param* param, int* imposs) {
  double W2=(param->caseP.Wbounds[1][0] - param->caseP.Wbounds[1][1])*t + param->caseP.Wbounds[1][1];
  if (W2==1 || W2==0) *imposs=1;
  else W2=log(W2/(1-W2));
  return W2;
}

/**
 * W1star'(t)
 * see paper for derivation: W1*(t) = (1/W1)*((w1_ub - w1_lb)/(1-W1)
 */
double getW1starPrimeFromT(double t, Param* param) {
  double m=(param->caseP.Wbounds[0][1] - param->caseP.Wbounds[0][0]);
  double W1=m*t + param->caseP.Wbounds[0][0];
  W1=(1/W1)*(m/(1-W1));
  return W1;
}

/**
 * W2star'(t)
 * see paper for derivation: W2*(t) = (1/W2)*((w2_lb - w2_ub)/(1-W2)
 */
double getW2starPrimeFromT(double t, Param* param) {
  double m=(param->caseP.Wbounds[1][0] - param->caseP.Wbounds[1][1]);
  double W2=m*t + param->caseP.Wbounds[1][1];
  W2=(1/W2)*(m/(1-W2));
  return W2;
}

/**
 * parameterized line integration
 * lower bound is t=0, upper bound is t=1
 */
double paramIntegration(integr_fn f, void *ex) {
  double epsabs=pow(10,-11), epsrel=pow(10,-11);
  double result=9999, anserr=9999;
  int limit=100;
  int last, neval, ier;
  int lenw=5*limit;
  int *iwork=(int *) Calloc(limit, int);
  double *work=(double *)Calloc(lenw, double);
  double lb=0.00001; double ub=.99999;
  Rdqags(f, ex, &lb, &ub, &epsabs, &epsrel, &result,
    &anserr, &neval, &ier, &limit, &lenw, &last, iwork, work);

  Free(iwork);
  Free(work);
  if (ier==0) return result;
  else {
    Param* p = (Param*) ex;
    Rprintf("Integration error %d: Sf %d X %5g Y %5g [%5g,%5g] -> %5g +- %5g\n",ier,p->caseP.suff,p->caseP.X,p->caseP.Y,p->caseP.Wbounds[0][0],p->caseP.Wbounds[0][1],result,anserr);
    char ch;
    scanf("Hit enter to continue %c", &ch );
    return result;
  }

}

/**
 * integrate normalizing constant and set it in param
 */
void setNormConst(Param* param) {
  param->caseP.normcT=paramIntegration(&NormConstT,(void*)param);
}


/**
 * Set the bounds on W1 and W2 in their parameter
 */
void setBounds(Param* param) {
  double X,Y,w1_lb,w1_ub,w2_lb,w2_ub;
  //int w1_inf,w2_inf;
  double tol0=0.0001;
  double tol1=0.9999;
  X=param->caseP.X;
  Y=param->caseP.Y;

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


  /*
  if (w1_lb==0 && w1_ub==1) w1_inf=2;
  else if (w1_lb==0) w1_inf=-1;
  else if (w1_ub==1) w1_inf=1;
  else w1_inf=0;
  //w1_lb=log(w1_lb/(1-w1_lb));
  //w1_ub=log(w1_ub/(1-w1_ub));

  if (w2_lb==0 && w2_ub==1) w2_inf=2;
  else if (w2_lb==0) w2_inf=-1;
  else if (w2_ub==1) w2_inf=1;
  else w2_inf=0;
  //w2_lb=log(w2_lb/(1-w2_lb));
  //w2_ub=log(w2_ub/(1-w2_ub));
  */
  param->caseP.Wbounds[0][0]=w1_lb;
  param->caseP.Wbounds[0][1]=w1_ub;
  param->caseP.Wbounds[1][0]=w2_lb;
  param->caseP.Wbounds[1][1]=w2_ub;
  //param->W1_inf=w1_inf;
  //param->W2_inf=w2_inf;

}
