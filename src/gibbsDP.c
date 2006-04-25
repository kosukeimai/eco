/******************************************************************
  This file is a part of eco: R Package for Fitting Bayesian Models 
  of Ecological Inference for 2x2 Tables
  by Kosuke Imai and Ying Lu
  Copyright: GPL version 2 or later.
*******************************************************************/

#include <string.h>
#include <stddef.h>
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

void cDPeco(
	    /*data input */
	    double *pdX,     /* data (X, Y) */
	    int *pin_samp,   /* sample size */

	    /*MCMC draws */
	    int *n_gen,      /* number of gibbs draws */ 
	    int *burn_in,    /* number of draws to be burned in */
	    int *pinth,      /* keep every nth draw */
	    int *verbose,    /* 1 for output monitoring */

	    /* prior specification*/
	    int *pinu0,      /* prior df parameter for InvWish */
	    double *pdtau0,  /* prior scale parameter for Sigma under G0*/ 
	    double *mu0,     /* prior mean for mu under G0 */
	    double *pdS0,    /* prior scale for Sigma */

	    /* DP prior specification */
	    double *alpha0,  /* precision parameter, can be fixed or updated*/
	    int *pinUpdate,  /* 1 if alpha gets updated */
	    double *pda0, double *pdb0, /* prior for alpha if alpha updated*/  

	    /*incorporating survey data */
	    int *survey,     /* 1 if survey data available (set of W_1, W_2) */
	                     /* 0 otherwise*/
	    int *sur_samp,   /* sample size of survey data*/
	    double *sur_W,   /* set of known W_1, W_2 */

	    /*incorporating homeogenous areas */
	    int *x1,         /* 1 if X=1 type areas available 
				W_1 known, W_2 unknown */
	    int *sampx1,     /* number X=1 type areas */
	    double *x1_W1,   /* values of W_1 for X1 type areas */

	    int *x0,         /* 1 if X=0 type areas available 
				W_2 known, W_1 unknown */
	    int *sampx0,     /* number X=0 type areas */
	    double *x0_W2,   /* values of W_2 for X0 type areas */

	    /* bounds of W1 */
	    double *minW1, double *maxW1,

	    /* storage */
	    int *parameter,  /* 1 if save population parameter */
	    int *Grid,       /* 1 if Grid algorithm used; \
				0 if Metropolis algorithm used*/

	    /* storage for Gibbs draws of mu/sigmat*/
	    double *pdSMu0, double *pdSMu1, 
	    double *pdSSig00, double *pdSSig01, double *pdSSig11,           
	    /* storage for Gibbs draws of W*/
	    double *pdSW1, double *pdSW2,
	    /* storage for Gibbs draws of alpha */
	    double *pdSa,
	    /* storage for nstar at each Gibbs draw*/
	    int *pdSn
 	    ){	   
  /*some integers */
  int n_samp = *pin_samp;    /* sample size */
  int s_samp = *sur_samp;    /* sample size of survey data */
  int x1_samp = *sampx1;     /* sample size for X=1 */
  int x0_samp = *sampx0;     /* sample size for X=0 */
  int t_samp = n_samp+x1_samp+x0_samp+s_samp; /* total sample size */
  int nth = *pinth;          /* keep every nth draw */ 
  int n_dim = 2;             /* dimension */
  int n_step=1000;           /* The default size of grid step */  

  /*prior parameters */
  double tau0 = *pdtau0;     /* prior scale */ 
  int nu0 = *pinu0;          /* prior degree of freedom*/ 
  double **S0 = doubleMatrix(n_dim,n_dim);/*The prior S parameter for InvWish*/
  double alpha = *alpha0;      /* precision parameter*/
  double a0 = *pda0, b0 = *pdb0; /* hyperprior for alpha */ 
  
  /* data */
  double **X = doubleMatrix(n_samp,n_dim);     /* The Y and covariates */
  double **W = doubleMatrix(t_samp,n_dim);     /* The W1 and W2 matrix */
  double **Wstar = doubleMatrix(t_samp,n_dim); /* The pseudo data  */
  double **S_W = doubleMatrix(s_samp,n_dim);    /* The known W1 and W2 matrix*/
  double **S_Wstar = doubleMatrix(s_samp,n_dim); /* The logit transformed S_W*/

  /* grids */
  double **W1g = doubleMatrix(n_samp, n_step); /* grids for W1 */
  double **W2g = doubleMatrix(n_samp, n_step); /* grids for W2 */
  int *n_grid = intArray(n_samp);              /* grids size */

  /* Model parameters */
  /* Dirichlet variables */
  double **mu = doubleMatrix(t_samp,n_dim);                /* mean matrix  */
  double ***Sigma = doubleMatrix3D(t_samp,n_dim,n_dim);    /*covarince matrix*/
  double ***InvSigma = doubleMatrix3D(t_samp,n_dim,n_dim); /* inv of Sigma*/
  
  int nstar;		           /* # clusters with distict theta values */
  int *C = intArray(t_samp);       /* vector of cluster membership */
  double *q = doubleArray(t_samp); /* Weights of posterior of Dirichlet */
  double *qq = doubleArray(t_samp); /* cumulative weight vector of q */
  double **S_bvt = doubleMatrix(n_dim,n_dim); /* S paramter for BVT in q0 */

  /* variables defined in remixing step: cycle through all clusters */
  double **Wstarmix = doubleMatrix(t_samp,n_dim);  /*data matrix used */ 
  double *mu_mix = doubleArray(n_dim);             /*updated MEAN parameter */
  double **Sigma_mix = doubleMatrix(n_dim,n_dim);  /*updated VAR parameter */
  double **InvSigma_mix = doubleMatrix(n_dim,n_dim); /* Inv of Sigma_mix */
  int nj;                            /* record # of obs in each cluster */
  int *sortC = intArray(t_samp);     /* record (sorted)original obs id */
  int *indexC = intArray(t_samp);   /* record  original obs id */
  int *label = intArray(t_samp);    /* store index values */

  /* misc variables */
  int i, j, k, l, main_loop;   /* used for various loops */
  int itemp;
  int itempA=0; /* counter for alpha */
  int itempS=0; /* counter for storage */
  int itempC=0; /* counter to control nth draw */
  int progress = 1, itempP = ftrunc((double) *n_gen/10);
  double dtemp, dtemp1;
  double *vtemp = doubleArray(n_dim);
  double **mtemp = doubleMatrix(n_dim,n_dim); 
  double **mtemp1 = doubleMatrix(n_dim,n_dim); 
  double **onedata = doubleMatrix(1, n_dim);

  /* get random seed */
  GetRNGstate();

  /* read priors under G0*/
  itemp=0;
  for(k=0;k<n_dim;k++)
    for(j=0;j<n_dim;j++) S0[j][k]=pdS0[itemp++];

  /* read the data set */
  itemp = 0;
  for (j = 0; j < n_dim; j++) 
    for (i = 0; i < n_samp; i++) X[i][j] = pdX[itemp++];

  /*Intialize W, Wsatr for n_samp */
  for (i=0; i< n_samp; i++) {
    if (X[i][1]!=0 && X[i][1]!=1) {
      W[i][0]=runif(minW1[i], maxW1[i]);
      W[i][1]=(X[i][1]-X[i][0]*W[i][0])/(1-X[i][0]);
    }
    if (X[i][1]==0) 
      for (j=0; j<n_dim; j++) W[i][j]=0.0001;
    if (X[i][1]==1) 
      for (j=0; j<n_dim; j++) W[i][j]=0.9999;
    for (j=0; j<n_dim; j++)
      Wstar[i][j]=log(W[i][j])-log(1-W[i][j]);
  }

  /*read homeogenous areas information */
  if (*x1==1)
    for (i=0; i<x1_samp; i++) {
      W[(n_samp+i)][0]=x1_W1[i];
      if (W[(n_samp+i)][0]==0) W[(n_samp+i)][0]=0.0001;
      if (W[(n_samp+i)][0]==1) W[(n_samp+i)][0]=0.9999;
      Wstar[(n_samp+i)][0]=log(W[(n_samp+i)][0])-log(1-W[(n_samp+i)][0]);
    }

  if (*x0==1)
    for (i=0; i<x0_samp; i++) {
      W[(n_samp+x1_samp+i)][1]=x0_W2[i];
      if (W[(n_samp+x1_samp+i)][1]==0) W[(n_samp+x1_samp+i)][1]=0.0001;
      if (W[(n_samp+x1_samp+i)][1]==1) W[(n_samp+x1_samp+i)][1]=0.9999;
      Wstar[(n_samp+x1_samp+i)][1]=log(W[(n_samp+x1_samp+i)][1])-log(1-W[(n_samp+x1_samp+i)][1]);
    }

  /*read the survey data */
  if (*survey==1) {
    itemp = 0;
    for (j=0; j<n_dim; j++)
      for (i=0; i<s_samp; i++) {
        S_W[i][j]=sur_W[itemp++];
        if (S_W[i][j]==0) S_W[i][j]=0.0001;
        if (S_W[i][j]==1) S_W[i][j]=0.9999;
        S_Wstar[i][j]=log(S_W[i][j])-log(1-S_W[i][j]);
	W[n_samp+x1_samp+x0_samp+i][j]=S_W[i][j];
	Wstar[n_samp+x1_samp+x0_samp+i][j]=S_Wstar[i][j];
      }
  }

  /* Calcualte grids */
  if (*Grid)
    GridPrep(W1g,W2g, X, maxW1, minW1, n_grid, n_samp, n_step);

  /* parmeters for Bivaraite t-distribution-unchanged in MCMC */
  for (j=0;j<n_dim;j++)
    for(k=0;k<n_dim;k++)
      mtemp[j][k]=S0[j][k]*(1+tau0)/(tau0*(nu0-n_dim+1));
  dinv(mtemp, n_dim, S_bvt);

  /**draw initial values of mu_i, Sigma_i under G0  for all effective sample**/
  /*1. Sigma_i under InvWish(nu0, S0^-1) with E(Sigma)=S0/(nu0-3)*/
  /*   InvSigma_i under Wish(nu0, S0^-1 */
  /*2. mu_i|Sigma_i under N(mu0, Sigma_i/tau0) */
  dinv(S0, n_dim, mtemp);
  for(i=0;i<t_samp;i++){
    /*draw from wish(nu0, S0^-1) */
    rWish(InvSigma[i], mtemp, nu0, n_dim);
    dinv(InvSigma[i], n_dim, Sigma[i]);
    for (j=0;j<n_dim;j++)
      for(k=0;k<n_dim;k++) mtemp1[j][k]=Sigma[i][j][k]/tau0;
    rMVN(mu[i], mu0, mtemp1, n_dim);
  }

  /* initialize the cluster membership */
  nstar=t_samp;  /* the # of disticnt values */
  for(i=0;i<t_samp;i++)
    C[i]=i; /*cluster is from 0...n_samp-1 */
  
  if (*verbose)
    Rprintf("Starting Gibbs Sampler...\n");
  for(main_loop=0; main_loop<*n_gen; main_loop++){
    /**update W, Wstar given mu, Sigma only for the unknown W/Wstar**/
    for (i=0;i<n_samp;i++){
      if (X[i][1]!=0 && X[i][1]!=1) {
	if (*Grid) 
	  rGrid(W[i], W1g[i], W2g[i], n_grid[i], mu[i], InvSigma[i], n_dim);
	else
	  rMH(W[i], X[i], minW1[i], maxW1[i],  mu[i], InvSigma[i], n_dim);
      }
      /*3 compute Wsta_i from W_i*/
      Wstar[i][0]=log(W[i][0])-log(1-W[i][0]);
      Wstar[i][1]=log(W[i][1])-log(1-W[i][1]);
    }
  
  if (*x1==1)
    for (i=0; i<x1_samp; i++) {
      dtemp=mu[n_samp+i][1]+Sigma[n_samp+i][0][1]/Sigma[n_samp+i][0][0]*(Wstar[n_samp+i][0]-mu[n_samp+i][0]);
      dtemp1=Sigma[n_samp+i][1][1]*(1-Sigma[n_samp+i][0][1]*Sigma[n_samp+i][0][1]/(Sigma[n_samp+i][0][0]*Sigma[n_samp+i][1][1]));
      Wstar[n_samp+i][1]=norm_rand()*sqrt(dtemp1)+dtemp;
      W[n_samp+i][1]=exp(Wstar[n_samp+i][1])/(1+exp(Wstar[n_samp+i][1]));
    }

  /*update W1 given W2, mu_ord and Sigma_ord in x0 homeogeneous areas */
  if (*x0==1)
    for (i=0; i<x0_samp; i++) {
      dtemp=mu[n_samp+x1_samp+i][0]+Sigma[n_samp+x1_samp+i][0][1]/Sigma[n_samp+x1_samp+i][1][1]*(Wstar[n_samp+x1_samp+i][1]-mu[n_samp+x1_samp+i][1]);
      dtemp1=Sigma[n_samp+x1_samp+i][0][0]*(1-Sigma[n_samp+x1_samp+i][0][1]*Sigma[n_samp+x1_samp+i][0][1]/(Sigma[n_samp+x1_samp+i][0][0]*Sigma[n_samp+x1_samp+i][1][1]));
      Wstar[n_samp+i][0]=norm_rand()*sqrt(dtemp1)+dtemp;
      W[n_samp+x1_samp+i][0]=exp(Wstar[n_samp+x1_samp+i][0])/(1+exp(Wstar[n_samp+x1_samp+i][0]));
    }

  /**updating mu, Sigma given Wstar uisng effective sample size W_star**/
  for (i=0; i<t_samp; i++){
    /* generate weight vector q */
    dtemp=0;
    for (j=0; j<t_samp; j++){
      if (j!=i)
	q[j]=dMVN(Wstar[i], mu[j], InvSigma[j], n_dim, 0);
      else
	q[j]=alpha*dMVT(Wstar[i], mu0, S_bvt, nu0-n_dim+1, 2, 0);
      dtemp+=q[j]; 
      qq[j]=dtemp; /*compute qq, the cumlative of q*/    
    }
    /*standardize q and qq */
    for (j=0; j<t_samp; j++) qq[j]/=dtemp;
    
    /** draw the configuration parameter **/
    j=0; dtemp=unif_rand();
    while (dtemp > qq[j]) j++;
    /** Dirichlet update Sigma_i, mu_i|Sigma_i **/
    /* j=i: posterior update given Wstar[i] */
    if (j==i){
      onedata[0][0] = Wstar[i][0];
      onedata[0][1] = Wstar[i][1];
      NIWupdate(onedata, mu[i], Sigma[i], InvSigma[i], mu0, tau0,nu0, S0, 1, n_dim);
      C[i]=nstar;
      nstar++;
    }
    /* j=i': replace with i' obs */
    else {
      /*1. mu_i=mu_j, Sigma_i=Sigma_j*/
      /*2. update C[i]=C[j] */
      for(k=0;k<n_dim;k++) {
	mu[i][k]=mu[j][k];
	for(l=0;l<n_dim;l++) {
	  Sigma[i][k][l]=Sigma[j][k][l];
	  InvSigma[i][k][l]=InvSigma[j][k][l];
	}
      }
      C[i]=C[j];
    }
    sortC[i]=C[i];
  } /* end of i loop*/
  
  /** remixing step using effective sample of Wstar**/
  for(i=0;i<t_samp;i++) indexC[i]=i;
  R_qsort_int_I(sortC, indexC, 1, t_samp);
  
  nstar=0; i=0;
  while (i<t_samp){
    j=sortC[i]; /*saves the first element in a block of same values */
    nj=0;       /* counter for a block of same values */
    
    /* get data for remixing */
    while ((sortC[i]==j) && (i<t_samp)) {
      label[nj]=indexC[i];
      for (k=0; k<n_dim; k++)
	Wstarmix[nj][k]=Wstar[label[nj]][k];
      nj++; i++;
    } /* i records the current position in IndexC */
      /* nj records the # of obs in Psimix */
    
    /** posterior update for mu_mix, Sigma_mix based on Psimix **/
    NIWupdate(Wstarmix, mu_mix,Sigma_mix, InvSigma_mix, mu0, tau0, nu0, S0, nj, n_dim);     
    
    /**update mu, Simgat with mu_mix, Sigmat_mix via label**/
    for (j=0;j<nj;j++){
      C[label[j]]=nstar;  /* updating C vector with no gap */
      for (k=0; k<n_dim; k++){
	mu[label[j]][k]=mu_mix[k];
	for (l=0;l<n_dim;l++){
	  Sigma[label[j]][k][l]=Sigma_mix[k][l];
	  InvSigma[label[j]][k][l]=InvSigma_mix[k][l];
	}
      }
    }
    nstar++; /*finish update one distinct value*/
  } /* nstar is the number of distinct values */
  
  /** updating alpha **/
  if(*pinUpdate) {
    dtemp=b0-log(rbeta(alpha+1, (double) t_samp));
    dtemp1=(double)(a0+nstar-1)/(t_samp*dtemp);
    if(unif_rand() < dtemp1)
      alpha=rgamma(a0+nstar, 1/dtemp);
    else 
      alpha=rgamma(a0+nstar-1, 1/dtemp);
  }
  
  /*store Gibbs draws after burn_in */
   R_CheckUserInterrupt();
  if (main_loop>=*burn_in) {
     itempC++;
    if (itempC==nth){
      if(*pinUpdate) {
	pdSa[itempA]=alpha;
     }
	pdSn[itempA]=nstar;     
      itempA++;
      
      for(i=0; i<(n_samp+x1_samp+x0_samp); i++) {
	pdSMu0[itempS]=mu[i][0];
	pdSMu1[itempS]=mu[i][1];
	pdSSig00[itempS]=Sigma[i][0][0];
	pdSSig01[itempS]=Sigma[i][0][1];
	pdSSig11[itempS]=Sigma[i][1][1];
	pdSW1[itempS]=W[i][0];
	pdSW2[itempS]=W[i][1];
	itempS++;
      }
      itempC=0; 
    }
  }

  if (*verbose)
    if (itempP == main_loop) {
      Rprintf("%3d percent done.\n", progress*10);
      itempP+=ftrunc((double) *n_gen/10); progress++;
       R_FlushConsole();
    }
  } /*end of MCMC for DP*/
  
  if (*verbose)
    Rprintf("100 percent done.\n");
  
  /** write out the random seed **/
   PutRNGstate();
  
  /* Freeing the memory */
  FreeMatrix(S0, n_dim);
  FreeMatrix(X, n_samp);
  FreeMatrix(W, t_samp);
  FreeMatrix(Wstar, t_samp);
  FreeMatrix(S_W, s_samp);
  FreeMatrix(S_Wstar, s_samp);
  FreeMatrix(W1g, n_samp);
  FreeMatrix(W2g, n_samp);
  free(n_grid);
  FreeMatrix(mu, t_samp);
  Free3DMatrix(Sigma, t_samp,n_dim);
  Free3DMatrix(InvSigma, t_samp, n_dim);
  free(C);
  free(q);
  free(qq);
  FreeMatrix(S_bvt, n_dim);
  FreeMatrix(Wstarmix, t_samp);
  free(mu_mix);
  FreeMatrix(Sigma_mix, n_dim);
  FreeMatrix(InvSigma_mix, n_dim);
  free(sortC);
  free(indexC);
  free(label);
  free(vtemp);
  FreeMatrix(mtemp, n_dim);
  FreeMatrix(mtemp1, n_dim);
  free(onedata);
} /* main */


