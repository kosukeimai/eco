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

void cEMeco(
	    /*data input */
	    double *pdX,         /* data (X, Y) */
	    double *pdTheta_in,  /* Theta^ t 
				    mu1, mu2, var1, var2, rho */
	    int *pin_samp,       /* sample size */

	    /*MCMC draws?? */
	    int *n_gen,          /* number of gibbs draws */

	    /*incorporating survey data */
	    int *survey,         /*1 if survey data available(W_1, W_2)
				   0 not*/
	    int *sur_samp,       /*sample size of survey data*/
	    double *sur_W,       /*set of known W_1, W_2 */ 
	    
	    /*incorporating homeogenous areas */
	    int *x1,       /* 1 if X=1 type areas available W_1 known, 
			      W_2 unknown */
	    int *sampx1,   /* number X=1 type areas */
	    double *x1_W1, /* values of W_1 for X1 type areas */
	    
	    int *x0,       /* 1 if X=0 type areas available W_2 known, 
			      W_1 unknown */
	    int *sampx0,   /* number X=0 type areas */
	    double *x0_W2, /* values of W_2 for X0 type areas */

	    /* bounds of W1 */
	    double *minW1, double *maxW1,

	    /* flags */
	    int *Grid,    /*1 if intergration is approximated by grids
			    0 R function */
	    
	    /* storage */
	    double *pdTheta,  /*EM result for Theta^(t+1) */
	    double *Suff      /*out put suffucient statistics (E(W_1i|Y_i), 
				E(E_1i*W_1i|Y_i..) when  conveges */
	    ){	   
  
  int n_samp = *pin_samp;    /* sample size */
  int s_samp= *sur_samp;     /* sample size of survey data */ 
  int x1_samp=*sampx1;       /* sample size for X=1 */
  int x0_samp=*sampx0;       /* sample size for X=0 */
  int t_samp=n_samp+s_samp+x1_samp+x0_samp;  /* total sample size*/
  int n_dim=2;        /* dimensions */
  int trapod=0;       /* 1 if use trapozodial ~= in numer. int.*/
  int n_step=5000;    /* The default size of grid step */  
  int data=0;         /* one to print the data */
  int ndraw=*n_gen;   /* number of draws */

  /* data */

  double **X=doubleMatrix(n_samp,n_dim);     /* Y and covariates */
  double **W=doubleMatrix(t_samp,n_dim);     /* W1 and W2 matrix */  
  double **Wstar=doubleMatrix(t_samp,5);     /* pseudo data(transformed*/
  double **S_W=doubleMatrix(s_samp, n_dim);  /* known W1 and W2 matrix*/
  double **S_Wstar=doubleMatrix(s_samp, n_dim);   /* transformed S_W*/

  /* grids */
  double **W1g=doubleMatrix(n_samp, n_step);   /* grids for W1 */
  double **W2g=doubleMatrix(n_samp, n_step);   /* grids for W2 */
  int *n_grid=intArray(n_samp);                /* grid size */

  /* numerical intergration */
  double *prob_grid=doubleArray(n_step);
  double *prob_grid_cum=doubleArray(n_step);

  /* model parameters */
  double *mu=doubleArray(n_dim);        /* mean vector*/
  double **Sigma=doubleMatrix(n_dim,n_dim);   /* covariance matrix*/
  double **InvSigma=doubleMatrix(n_dim,n_dim);/* inverse covariance matrix*/

  /* misc variables */
  int i, j, k, l, main_loop;   /* used for various loops */
  int itemp, itemp0;
  double dtemp, dtemp1, temp0, temp1;
  double *vtemp=doubleArray(n_dim);
  int *mflag=intArray(n_step);

  /* get random seed */
  GetRNGstate();

  /* read the data set */
  /** Packing Y, X  **/
  itemp = 0;
  for (j = 0; j < n_dim; j++) 
    for (i = 0; i < n_samp; i++) {
      X[i][j] = pdX[itemp++];
    }
  if (data) { 
    printf("Y X\n");
    for(i=0;i<n_samp;i++)
      Rprintf("%5d%14g%14g\n",i,X[i][1],X[i][0]);
      }

  /** read mu_org, Simga_ord **/
  mu[0] = pdTheta_in[0];
  mu[1] = pdTheta_in[1];
  Sigma[0][0] = pdTheta_in[2];
  Sigma[1][1] = pdTheta_in[3];
  Sigma[0][1] = pdTheta_in[4]*sqrt(pdTheta_in[2]*pdTheta_in[3]);
  Sigma[1][0] = Sigma[0][1];
  dinv(Sigma, n_dim, InvSigma);

  /* initialize W, Wstar for t_samp*/
  for (i=0; i< t_samp; i++) 
    for (j=0; j<5; j++) {
      if (j<2)
	W[i][j]=0;
      Wstar[i][j]=0;
    }

  for (i=0; i<n_samp; i++) {
    if (X[i][1]==0) {
      W[i][0]=0.0001;
      W[i][1]=0.0001;
      temp0=log(W[i][0])-log(1-W[i][0]);
      temp1=log(W[i][1])-log(1-W[i][1]);
      Wstar[i][0]=temp0;
      Wstar[i][1]=temp1;
      Wstar[i][2]=temp0*temp0;
      Wstar[i][3]=temp0*temp1;
      Wstar[i][4]=temp1*temp1;
    }
    else if (X[i][1]==1) {
      W[i][0]=0.9999; 
      W[i][1]=0.9999;
      temp0=log(W[i][0])-log(1-W[i][0]);
      temp1=log(W[i][1])-log(1-W[i][1]);
      Wstar[i][0]=temp0;
      Wstar[i][1]=temp1;
      Wstar[i][2]=temp0*temp0;
      Wstar[i][3]=temp0*temp1;
      Wstar[i][4]=temp1*temp1;
    }
  } 
  
  /*read homeogenous areas information */
  if (*x1) 
    for (i=0; i<x1_samp; i++) {
      W[(n_samp+i)][0]=x1_W1[i];
      if (W[(n_samp+i)][0]==0) W[(n_samp+i)][0]=0.0001;
      if (W[(n_samp+i)][0]==1) W[(n_samp+i)][0]=0.9999;
    }
  
  if (*x0) 
    for (i=0; i<x0_samp; i++) {
      W[(n_samp+x1_samp+i)][1]=x0_W2[i];
      if (W[(n_samp+x1_samp+i)][1]==0) W[(n_samp+x1_samp+i)][1]=0.0001;
      if (W[(n_samp+x1_samp+i)][1]==1) W[(n_samp+x1_samp+i)][1]=0.9999;
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
	W[(n_samp+x1_samp+x0_samp+i)][j]=S_W[i][j];
	Wstar[(n_samp+x1_samp+x0_samp+i)][j]=S_Wstar[i][j];
      }
    for (i=0; i<s_samp; i++) {
      Wstar[(n_samp+x1_samp+x0_samp+i)][2]=S_Wstar[i][0]*S_Wstar[i][0];
      Wstar[(n_samp+x1_samp+x0_samp+i)][3]=S_Wstar[i][0]*S_Wstar[i][1];
      Wstar[(n_samp+x1_samp+x0_samp+i)][4]=S_Wstar[i][1]*S_Wstar[i][1];
    }
  }    
  /*    if (data) { 
      Rprintf("survey W1 W2 W1* W2*\n");
      for(i=0;i<(n_samp+s_samp+x1_samp+x0_samp);i++)
	Rprintf("%5d%14g%14g%14g%14g\n",i,W[i][0],W[i][1],Wstar[i][0], Wstar[i][1]);
    }
  */
  
  /*** calculate grids ***/
  if (*Grid)
    GridPrep(W1g, W2g, X, maxW1, minW1, n_grid, n_samp, n_step);

  for (i=0; i<n_step; i++) {
    mflag[i]=0;
  }


  /**update W, Wstar given mu, Sigma in regular areas**/
  for (i=0;i<n_samp;i++){
    if ( X[i][1]!=0 && X[i][1]!=1 ) {
      /*1 project BVN(mu, Sigma) on the inth tomo line */
      dtemp=0;
      for (j=0;j<n_grid[i];j++){
	vtemp[0]=log(W1g[i][j])-log(1-W1g[i][j]);
	vtemp[1]=log(W2g[i][j])-log(1-W2g[i][j]);
	prob_grid[j]=dMVN(vtemp, mu, InvSigma, 2, 1) -
	  log(W1g[i][j])-log(W2g[i][j])-log(1-W1g[i][j])-log(1-W2g[i][j]);
	prob_grid[j]=exp(prob_grid[j]);
	dtemp+=prob_grid[j];
	prob_grid_cum[j]=dtemp;
      }
      for (j=0;j<n_grid[i];j++){
	prob_grid_cum[j]/=dtemp; /*standardize prob.grid */ 
      }
      /*2 sample W_i on the ith tomo line */
      /*3 compute Wsta_i from W_i*/
      j=0;
      itemp=1;
      
      for (k=0; k<ndraw; k++){
	j=findInterval(prob_grid_cum, n_grid[i],
		       (double)(1+k)/(ndraw+1), 1, 1, itemp, mflag); 
	itemp=j;
      if ((W1g[i][j]==0) || (W1g[i][j]==1)) 
	Rprintf("W1g%5d%5d%14g", i, j, W1g[i][j]);
      if ((W2g[i][j]==0) || (W2g[i][j]==1)) 
	Rprintf("W2g%5d%5d%14g", i, j, W2g[i][j]);

	if (j==0 || trapod==0) {
	  W[i][0]=W1g[i][j];
	  W[i][1]=W2g[i][j];
	}
	else if (j>=1 && trapod==1) {
	   if (prob_grid_cum[j]!=prob_grid_cum[(j-1)]) {
	  dtemp1=((double)(1+k)/(ndraw+1)-prob_grid_cum[(j-1)])/(prob_grid_cum[j]-prob_grid_cum[(j-1)]);
	    W[i][0]=dtemp1*(W1g[i][j]-W1g[i][(j-1)])+W1g[i][(j-1)];
	    W[i][1]=dtemp1*(W2g[i][j]-W2g[i][(j-1)])+W2g[i][(j-1)];
	   }
          else if (prob_grid_cum[j]==prob_grid_cum[(j-1)]) {
	  W[i][0]=W1g[i][j];
	  W[i][1]=W2g[i][j];
	  }
	}
	temp0=log(W[i][0])-log(1-W[i][0]);
	temp1=log(W[i][1])-log(1-W[i][1]);
	Wstar[i][0]+=temp0;
	Wstar[i][1]+=temp1;
	Wstar[i][2]+=temp0*temp0;
	Wstar[i][3]+=temp0*temp1;
	Wstar[i][4]+=temp1*temp1;
      }
    }
  }
  
  /* compute E_{W_i|Y_i} for n_samp*/
  for (i=0; i<n_samp; i++) {
    if ( X[i][1]!=0 && X[i][1]!=1 ) {  
      Wstar[i][0]/=ndraw;  /*E(W1i) */
      Wstar[i][1]/=ndraw;  /*E(W2i) */
      Wstar[i][2]/=ndraw;  /*E(W1i^2) */
      Wstar[i][3]/=ndraw;  /*E(W1iW2i) */
      Wstar[i][4]/=ndraw;  /*E(W2i^2) */
    }
  } /*for x0type, x1type and survey data, E-step is either the observed value or the analytical expectation*/
  
  
    /* analytically compute E{W2_i|Y_i} given W1_i, mu and Sigma in x1 homeogeneous areas */
  if (*x1==1)
    for (i=0; i<x1_samp; i++) {
      dtemp=mu[1]+Sigma[0][1]/Sigma[0][0]*(Wstar[n_samp+i][0]-mu[0]);
      temp0=log(W[n_samp+i][0])-log(1-W[n_samp+i][0]);
      temp1=dtemp;
      Wstar[n_samp+i][0]=temp0;
      Wstar[n_samp+i][1]=temp1;
      Wstar[n_samp+i][2]=temp0*temp0;
      Wstar[n_samp+i][3]=temp0*temp1;
      Wstar[n_samp+i][4]=temp1*temp1;
    }
    
  /*analytically compute E{W1_i|Y_i} given W2_i, mu and Sigma in x0 homeogeneous areas */
  if (*x0==1)
    for (i=0; i<x0_samp; i++) {
      dtemp=mu[0]+Sigma[0][1]/Sigma[1][1]*(Wstar[n_samp+x1_samp+i][1]-mu[1]);
      temp0=dtemp;
      temp1=log(W[n_samp+x1_samp+i][1])-log(1-W[n_samp+x1_samp+i][1]);
      Wstar[n_samp+x1_samp+i][0]=temp0;
      Wstar[n_samp+x1_samp+i][1]=temp1;
      Wstar[n_samp+x1_samp+i][2]=temp0*temp0;
      Wstar[n_samp+x1_samp+i][3]=temp0*temp1;
      Wstar[n_samp+x1_samp+i][4]=temp1*temp1;
    }

    if (data) { 
      Rprintf("survey W1 W2 W1* W2*\n");
      for(i=0;i<(n_samp+s_samp+x1_samp+x0_samp);i++)
	Rprintf("%5d%14g%14g%14g%14g\n",i,W[i][0],W[i][1],Wstar[i][0], Wstar[i][1]);
    }
       
  /*M-step: same procedure as normal model */
  for (j=0; j<5; j++) 
    Suff[j]=0;
  
  /* compute sufficient statistics */
  for (i=0; i<t_samp; i++) {
    Suff[0]+=Wstar[i][0];  /* sumE(W_i1|Y_i) */
    Suff[1]+=Wstar[i][1];  /* sumE(W_i2|Y_i) */
    Suff[2]+=Wstar[i][2];  /* sumE(W_i1^2|Y_i) */
    Suff[3]+=Wstar[i][4];  /* sumE(W_i2^2|Y_i) */
    Suff[4]+=Wstar[i][3];  /* sumE(W_i1^W_i2|Y_i) */
  }
  
  pdTheta[0]=Suff[0]/t_samp;  /*mu1*/
  pdTheta[1]=Suff[1]/t_samp;  /*mu2*/
  pdTheta[2]=(Suff[2]-2*Suff[0]*pdTheta[0]+t_samp*pdTheta[0]*pdTheta[0])/t_samp;  /*sigma11*/
  pdTheta[3]=(Suff[3]-2*Suff[1]*pdTheta[1]+t_samp*pdTheta[1]*pdTheta[1])/t_samp;  /*sigma22*/
  pdTheta[4]=(Suff[4]-Suff[0]*pdTheta[1]-Suff[1]*pdTheta[0]+t_samp*pdTheta[0]*pdTheta[1])/t_samp; /*sigma12*/
  pdTheta[4]=pdTheta[4]/sqrt(pdTheta[2]*pdTheta[3]); /*rho*/

  if (data) { 
      Rprintf("theta and suff\n");

	Rprintf("%14g%14g%14g%14g\n",pdTheta[0],pdTheta[1],pdTheta[2],pdTheta[3]);
	Rprintf("%14g%14g%14g%14g\n",Suff[0],Suff[1],Suff[2],Suff[3]);

    }  
  /* write out the random seed */
  PutRNGstate();

  /* Freeing the memory */
  FreeMatrix(X, n_samp);
  FreeMatrix(W, t_samp);
  FreeMatrix(Wstar, t_samp);
  FreeMatrix(S_W, s_samp);
  FreeMatrix(S_Wstar, s_samp);
  free(n_grid);
  FreeMatrix(W1g, n_samp);
  FreeMatrix(W2g, n_samp);
  free(prob_grid);
  free(prob_grid_cum);
  free(mu);
  FreeMatrix(Sigma,n_dim);
  FreeMatrix(InvSigma, n_dim);
  free(vtemp);
  free(mflag);
} /* main */

