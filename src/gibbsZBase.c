#include <string.h>
#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "sample.h"

void cBaseecoZ(
	      /*data input */
	      double *pdX,     /* data (X, Y) */
	      double *pdZ,     /* covariates Z */
	      int *pinZp,      /* dimension of Z
				  if =1, =gibbsBase
			             =2 and Z=X, gibbsXBase
			             >2 or Z!=X, regression*/
	      int *pin_samp,   /* sample size */
	      /*MCMC draws */
	      int *n_gen,      /* number of gibbs draws */
	      int *burn_in,    /* number of draws to be burned in */
	      int *pinth,        /* keep every nth draw */
	      int *verbose,    /* 1 for output monitoring */

	      /* prior specification for imputation, (beta, Sigma)~N-InvWish*/
	      /* prior for Sigma~InvWish(nu, S)*/
	      int *pinu0,       /* prior df parameter for InvWish */
	      double *pdS0,     /* prior scale for Sigma */
	      /* prior for beta~N(b0, Sigma*A0^-1) */
	      double *pdbeta0,      /* prior mean for beta*/
	      double *pdA0,   /* prior PRECISION=1/SCALE parameter for beta*/

	      /* staring values */
	      double *betastart, double *Sigmastart,

	      /*incorporating survey data */
	      int *survey,      /*1 if survey data available (set of W_1, W_2)
				  0 not*/
	      int *sur_samp,    /*sample size of survey data*/
	      double *sur_W,    /*set of known W_1, W_2 */ 
	      
	      /*incorporating homeogenous areas */
	      int *x1,       /* 1 if X=1 type areas available W_1 known, W_2 unknown */
	      int *sampx1,   /* number X=1 type areas */
	      double *x1_W1, /* values of W_1 for X1 type areas */

	      int *x0,       /* 1 if X=0 type areas available W_2 known, W_1 unknown */
	      int *sampx0,   /* number X=0 type areas */
	      double *x0_W2, /* values of W_2 for X0 type areas */

	      /* bounds of W1 */
	      double *minW1, double *maxW1,

	      /* storage */
	      int *parameter,/* 1 if save population parameter */
	      int *Grid,

	      /* storage for Gibbs draws of beta and Sigam, packed */
	      double *pdSBeta, double *pdSSigma,
	      /* storage for Gibbs draws of W*/
	      double *pdSW1, double *pdSW2
	      ){	   
  
  int n_samp = *pin_samp; /* sample size */
  int nth = *pinth;  
  int s_samp = *sur_samp; /* sample size of survey data */ 
  int x1_samp = *sampx1;
  int x0_samp = *sampx0;
  int t_samp = n_samp+s_samp+x1_samp+x0_samp;  /* total sample size */ 
  int n_dim = 2;          /* The dimension of the ecological table */
  int n_cov = *pinZp;     /* The dimension of the covariates */
  int n_step = 1000;
  
  /* priors */
  double *beta0 = doubleArray(n_cov); /* prior mean of beta */
  double **S0 = doubleMatrix(n_dim, n_dim); /* prior scale for Sigma */
  double **A0 = doubleMatrix(n_cov, n_cov); /* prior precision for beta */ 
  int nu0 = *pinu0;                         /* prior df for Sigma */   

  /* data */
  double **X = doubleMatrix(n_samp, n_dim); /* The Y and X */

  /*The known W1 and W2 matrix*/
  double **S_W = doubleMatrix(s_samp, n_dim); 
  double **S_Wstar=doubleMatrix(s_samp, n_dim); 

  /* pseudo data Wstar */
  double **W = doubleMatrix(t_samp, n_dim);
  double **Wstar = doubleMatrix(t_samp, n_dim);
  double *Wstar_bar = doubleArray(n_dim);

  /* The covariates and W */ 
  double **Z = doubleMatrix(t_samp*n_dim+n_cov, n_cov+1);
  /* Z*cholesky factor of covaraince matrix*/ 
  double **Zstar = doubleMatrix(t_samp*n_dim+n_cov, n_cov+1);

  /* grids */
  double **W1g = doubleMatrix(n_samp, n_step); /* grids for W1 */
  double **W2g = doubleMatrix(n_samp, n_step); /* grids for W2 */
  int *n_grid = intArray(n_samp);              /* grid size */

  /* paramters for Wstar under Normal baseline model */
  double *beta = doubleArray(n_cov); /* vector of regression coefficients */
  double **mu = doubleMatrix(t_samp, n_dim); 
  double **Sigma = doubleMatrix(n_dim, n_dim);
  double **InvSigma = doubleMatrix(n_dim, n_dim);

  /*posterior parameters for beta and Sigma*/
  double *mbeta = doubleArray(n_cov);         /* posterior mean of beta*/
  double **Vbeta = doubleMatrix(n_cov,n_cov); /* posterior varaince of beta */ 

  /* matrices used for sweep */
  /* quantities used in sweep */
  double **SS = doubleMatrix(n_cov+1, n_cov+1); /* the sum of square matrix */
  double *epsilon = doubleArray(t_samp*n_dim);  /* The error term */
  double **R = doubleMatrix(n_dim, n_dim);      /* ee' */
 
  /* misc variables */
  int i, j, k, l, main_loop;   /* used for various loops */
  int itemp;
  int itempA=0; /* counter for alpha */
  int itempB=0; 
  int itempC=0; /* counter to control nth draw */
  int itempS=0; /* counter for storage */

  int progress = 1, itempP = ftrunc((double) *n_gen/10);
  double dtemp, dtemp1;
  double *vtemp = doubleArray(n_dim);
  double **mtemp = doubleMatrix(n_dim, n_dim);
  double **mtemp1 = doubleMatrix(n_dim, n_dim);
  double **mtemp2 = doubleMatrix(n_cov, n_cov);

  /* get random seed */
  GetRNGstate();

  /**read prior information*/
  itemp=0;
  for (k=0; k<n_cov; k++) {
    beta0[k]=pdbeta0[k];
    for (j=0; j<n_cov; j++) 
      A0[j][k]=pdA0[itemp++];
  }
  itemp=0;
  for (k=0; k<n_dim; k++)
    for (j=0; j<n_dim; j++) 
      S0[j][k]=pdS0[itemp++];
    
  /* read the data set */
  /** Packing Y, X  **/
  itemp = 0;
  for (j = 0; j < n_dim; j++) 
    for (i = 0; i < n_samp; i++) 
      X[i][j] = pdX[itemp++];
  
  /**read Z **/
  for (i=0; i<t_samp*n_dim+n_cov; i++)
    for(j=0; j<=n_cov;j++) {
      Z[i][j]=0;
      Zstar[i][j]=0;
    }
  itemp = 0;
  for (k=0; k<n_cov; k++)
    for (j=0; j<n_dim; j++)
      for (i=0; i<t_samp; i++)
	Z[j*t_samp+i][k]=pdZ[itemp++];
  
  /* add prior information to Z*/
  dcholdc(A0, n_cov, mtemp2);  /*Cholesky decomopsition*/

  for (j=0; j<n_cov; j++) {
    Zstar[t_samp*n_dim+j][n_cov]=beta0[j];
    for (k=0; k<n_cov; k++){
      Zstar[t_samp*n_dim+j][n_cov]+=mtemp2[j][k]*beta0[j];
      Zstar[t_samp*n_dim+j][k]=mtemp2[j][k];
    }
  }      
    
  /* initialize W, Wstar for n_samp*/
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
	W[(n_samp+x1_samp+x0_samp+i)][j]=S_W[i][j];
	Wstar[(n_samp+x1_samp+x0_samp+i)][j]=S_Wstar[i][j];
	Z[(i+n_samp+x1_samp+x0_samp)*n_dim+j][n_cov]=Wstar[(i+n_samp+x1_samp+x0_samp)][j];
      }
  }

  /* calculate grids */
  if (*Grid)
    GridPrep(W1g, W2g, X, maxW1, minW1, n_grid, n_samp, n_step);

  /* starting vales of mu and Sigma */
  itemp = 0;
  for(j=0;j<n_cov;j++)
    beta[j] = betastart[j];
  for(j=0;j<n_dim;j++)
    for(k=0;k<n_dim;k++)
      Sigma[j][k]=Sigmastart[itemp++];
  dinv(Sigma, n_dim, InvSigma);

  /***Gibbs for  normal prior ***/
  for(main_loop=0; main_loop<*n_gen; main_loop++){
    for (i=0; i<t_samp; i++)
      for (j=0; j<n_dim; j++) 
	mu[i][j]=0;
    /**update W, Wstar given mu, Sigma in regular areas**/
    for (i=0;i<t_samp;i++)
      for (j=0; j<n_dim; j++)
	for (k=0; k<n_cov; k++) 
	  mu[i][j]+=Z[i*n_dim+j][k]*beta[k];
    
    for (i=0; i<n_samp; i++) {
      if ( X[i][1]!=0 && X[i][1]!=1 ) {
	/*1 project BVN(mu, Sigma) on the inth tomo line */
	/*2 sample W_i on the ith tomo line */
	if (*Grid)
	  rGrid(W[i], W1g[i], W2g[i], n_grid[i], mu[i], InvSigma, n_dim);
	else
	  rMH(W[i], X[i], minW1[i], maxW1[i], mu[i], InvSigma, n_dim);
      } 
      /*3 compute Wsta_i from W_i*/
      Wstar[i][0]=log(W[i][0])-log(1-W[i][0]);
      Wstar[i][1]=log(W[i][1])-log(1-W[i][1]);
      Z[i*n_dim][n_cov]=Wstar[i][0];
      Z[i*n_dim+1][n_cov]=Wstar[i][1];
    }
    
    /*update W2 given W1, mu and Sigma in x1 homeogeneous areas */
    if (*x1==1)
      for (i=0; i<x1_samp; i++) {
	dtemp=mu[n_samp+i][1]+Sigma[0][1]/Sigma[0][0]*(Wstar[n_samp+i][0]-mu[n_samp+i][0]);
	dtemp1=Sigma[1][1]*(1-Sigma[0][1]*Sigma[0][1]/(Sigma[0][0]*Sigma[1][1]));
	Rprintf("\n%14g%14g\n", dtemp, dtemp1);
	dtemp1=sqrt(dtemp1);
	Wstar[n_samp+i][1]=rnorm(dtemp, dtemp1);
	W[n_samp+i][1]=exp(Wstar[n_samp+i][1])/(1+exp(Wstar[n_samp+i][1]));
	Z[(i+n_samp)*n_dim][n_cov]=Wstar[(i+n_samp)][0];
	Z[(i+n_samp)*n_dim+1][n_cov]=Wstar[(i+n_samp)][1];
      }

    /*update W1 given W2, mu and Sigma in x0 homeogeneous areas */
    if (*x0==1)
      for (i=0; i<x0_samp; i++) {
	dtemp=mu[n_samp+x1_samp+i][0]+Sigma[0][1]/Sigma[1][1]*(Wstar[n_samp+x1_samp+i][1]-mu[n_samp+x1_samp+i][1]);
	dtemp1=Sigma[0][0]*(1-Sigma[0][1]*Sigma[0][1]/(Sigma[0][0]*Sigma[1][1]));
	dtemp1=sqrt(dtemp1);
	Wstar[n_samp+x1_samp+i][0]=rnorm(dtemp, dtemp1);
	W[n_samp+x1_samp+i][0]=exp(Wstar[n_samp+x1_samp+i][0])/(1+exp(Wstar[n_samp+x1_samp+i][0]));
	Z[(i+n_samp+x1_samp)*n_dim][n_cov]=Wstar[(i+n_samp+x1_samp)][0];
	Z[(i+n_samp+x1_samp)*n_dim+1][n_cov]=Wstar[(i+n_samp+x1_samp)][1];
      }

    dcholdc(InvSigma, n_dim, mtemp);
    for (i=0; i<t_samp*n_dim; i++)
      for (j=0; j<=n_cov; j++)
	Zstar[i][j]=0;
    for (i=0; i<t_samp; i++)
      for(j=0; j<n_dim; j++)
	for(k=0; k<n_dim; k++)
	  for(l=0; l<=n_cov; l++)
	    Zstar[i*n_dim+k][l]+=mtemp[k][j]*Z[i*n_dim+j][l];

    /*construct SS matrix for SWEEP */
    for (j=0; j<=n_cov; j++)
      for (k=0; k<=n_cov; k++)
	SS[j][k]=0;
    for(i=0; i<(t_samp*n_dim); i++)
      for(k=0; k<=n_cov; k++)
	for(l=0; l<=n_cov; l++)
	  SS[k][l]+=Zstar[i][k]*Zstar[i][l];
    for(j=0; j<n_cov; j++)
      for(k=0; k<=n_cov; k++)
	for(l=0; l<=n_cov; l++)
	  SS[k][l]+=Zstar[n_samp*n_dim+j][k]*Zstar[n_samp*n_dim+j][l];

    /*SWEEP to get posterior mean anf variance for beta */
    for (j=0; j<n_cov; j++) 
      SWP(SS,j,n_cov+1);

    /*draw beta given Sigma and W */
    for (j=0; j<n_cov; j++) {
      mbeta[j]=SS[j][n_cov]; 
      for (k=0; k<n_cov; k++)
	Vbeta[j][k]=-SS[j][k];
    }
    rMVN(beta, mbeta, Vbeta, n_cov);

    /*draw Sigmar give beta and Wstar */
    for(i=0; i<t_samp; i++)
      for(j=0; j<n_dim; j++) {
	epsilon[i*n_dim+j]=Z[i*n_dim+j][n_cov];
	for (k=0;k<n_cov; k++)
	  epsilon[i*n_dim+j]-=Z[i*n_dim+j][k]*beta[k];
      }
    for(j=0; j<n_dim; j++)
      for(k=0; k<n_dim; k++) 
	R[j][k]=0;
    for (i=0; i<t_samp; i++)
      for(j=0; j<n_dim; j++)
	for(k=0; k<n_dim; k++) 
	  R[j][k]+=epsilon[i*n_dim+j]*epsilon[i*n_dim+k];
    for(j=0; j<n_dim; j++)
      for (k=0; k<n_dim; k++)
	mtemp[j][k]=S0[j][k]+R[j][k];
    dinv(mtemp, n_dim, mtemp1);
    rWish(InvSigma, mtemp1, nu0+t_samp, n_dim);
    dinv(InvSigma, n_dim, Sigma);
    
    /*store Gibbs draw after burn-in and every nth draws */      
    R_CheckUserInterrupt();
    if (main_loop>=*burn_in){
      itempC++;
      if (itempC==nth){
	for (j=0; j<n_cov; j++)
	  pdSBeta[itempA++]=beta[j];
	for (j=0; j<n_dim; j++)
	  for (k=j; k<n_dim; k++)
	    pdSSigma[itempB++]=Sigma[j][k];
	for(i=0; i<(n_samp+x1_samp+x0_samp); i++){
	  pdSW1[itempS]=W[i][0];
	  pdSW2[itempS]=W[i][1];
	  itempS++;
	}
	itempC=0;
      }
    } /*end of stroage *burn_in*/
    if (*verbose)
      if (itempP == main_loop) {
	Rprintf("%3d percent done.\n", progress*10);
	itempP+=ftrunc((double) *n_gen/10); progress++;
      R_FlushConsole();
      }
  } /*end of MCMC for normal */ 
  
  /** write out the random seed **/
  PutRNGstate();

  /* Freeing the memory */
  FreeMatrix(X, n_samp);
  FreeMatrix(W, t_samp);
  FreeMatrix(Wstar, t_samp);
  FreeMatrix(S_W, s_samp);
  FreeMatrix(S_Wstar, s_samp);
  Free(minW1);
  Free(maxW1);
  free(n_grid);
  FreeMatrix(S0, n_dim);
  FreeMatrix(W1g, n_samp);
  FreeMatrix(W2g, n_samp);
  FreeMatrix(mu,t_samp);
  FreeMatrix(Sigma,n_dim);
  FreeMatrix(InvSigma, n_dim);
  FreeMatrix(Z, t_samp*n_dim+n_cov);
  FreeMatrix(Zstar, t_samp*n_dim+n_cov);
  Free(Wstar_bar);
  Free(vtemp);
  FreeMatrix(mtemp, n_dim);
  FreeMatrix(mtemp1, n_dim);
  FreeMatrix(mtemp2, n_cov);
  Free(beta);
  Free(beta0);
  FreeMatrix(A0, n_cov);
  FreeMatrix(SS, n_cov+1);
  Free(mbeta);
  FreeMatrix(Vbeta, n_cov);
  Free(epsilon);
  FreeMatrix(R, n_dim);
  
} /* main */

