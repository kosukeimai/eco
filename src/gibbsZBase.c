#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"

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

	      /*incorporating survey data */
	      int *survey,      /*1 if survey data available (set of W_1, W_2) */
	                       /*0 not*/
	      int *sur_samp,     /*sample size of survey data*/
	      double *sur_W,    /*set of known W_1, W_2 */ 
				  
	      /*incorporating homeogenous areas */
	      int *x1,       /* 1 if X=1 type areas available W_1 known, W_2 unknown */
	      int *sampx1,  /* number X=1 type areas */
	      double *x1_W1, /* values of W_1 for X1 type areas */

	      int *x0,       /* 1 if X=0 type areas available W_2 known, W_1 unknown */
	      int *sampx0,  /* number X=0 type areas */
	      double *x0_W2, /* values of W_2 for X0 type areas */

	      /* storage */
	      int *pred,       /* 1 if draw posterior prediction */
	      int *parameter,   /* 1 if save population parameter */

	      /* storage for Gibbs draws of beta and Sigam, packed */
	      double *pdSBeta, 
	      double *pdSSigma,
	      /* storage for Gibbs draws of W*/
	      double *pdSW1, double *pdSW2,
	      /* storage for posterior predictions of W */
	      double *pdSWt1, double *pdSWt2
	      ){	   
  
  int n_samp = *pin_samp;    /* sample size */
  int nth=*pinth;  

  int n_dim=2;           /* The dimension of the ecological table */
  int n_cov=*pinZp;      /* The dimension of the covariates */

  double *beta;        /* vector of regression coefficients */
  double *beta0;          /* prior mean of beta */
  double **A0;             /* prior precision matrix for beta */ 
  int nu0 = *pinu0;    /* prior scale parameter for Sigma */   
  double **S0;         /* The prior S parameter for Sigma */

  double **X;	    	 /* The Y and X */
  double **Z;	    	 /* The covariates and W matrix for sweep , X can be in Z*/
  double **Zstar;	 /* Z*cholesky factor of covaraince matrix*/

  int s_samp= *sur_samp;   /* sample size of survey data */ 
  double **S_W;            /*The known W1 and W2 matrix*/
  double **S_Wstar;        /*The inverse logit transformation of S_W*/

  int x1_samp=*sampx1;
  int x0_samp=*sampx0;

  int t_samp;              /* total effective sample size =n_samp+s_samp;*/

  /*bounds condition variables */
  double **W;            /* The W1 and W2 matrix */
  double *minW1, *maxW1; /* The lower and upper bounds of W_1i */
  int n_step=1000;    /* 1/The default size of grid step */  
  int *n_grid;           /* The number of grids for sampling on tomoline */
  double **W1g, **W2g;   /* The grids taken for W1 and W2 on tomoline */
  double *prob_grid;     /* The projected density on tomoline */
  double *prob_grid_cum; /* The projected cumulative density on tomoline */
  double *resid;         /* The centralizing vector for grids */

  /* pseudo data Wstar */
  double **Wstar;        /* The pseudo data  */
  double *Wstar_bar;

  /* paramters for Wstar under Normal baseline model */
  double **Sigma_ord;   /* The posterior estimation of covariance matrix of Wstar*/
  double **InvSigma_ord;
  double **mu_ord;        /* The posterior estimation of mean of Wstar=Z*beta*/

  /*posterior parameters for beta and Sigma*/
  double *mbeta;           /* The posterior mean of beta*/
  double **Vbeta;          /* The posterior varaince of beta */ 
  /* unused: double **Sn;               The posterior S parameter for InvWish */

  /* matrices used for sweep */
  double **SS;           /* the sum of square matrix used for sweep operator */
  double *epsilon;       /* The error term of the regression */
  double **R;            /* The sum of square matrix of ee' */

  /* misc variables */
  int i, j, k, t, l, main_loop;   /* used for various loops */
  int itemp, itempS, itempC, itempA, itempB;
  int progress = 1, itempP = ftrunc((double) *n_gen/10);
  double dtemp, dtemp1;
  double *vtemp;
  double **mtemp, **mtemp2, **mtemp3;

  /* get random seed */
  GetRNGstate();

  /* defining vectors and matricies */
  /* data */
  X=doubleMatrix(n_samp,n_dim);
  W=doubleMatrix((n_samp+s_samp+x0_samp+x1_samp),n_dim);
  Wstar=doubleMatrix((n_samp+s_samp+x0_samp+x1_samp),n_dim);

  Z=doubleMatrix((n_samp+s_samp+x0_samp+x1_samp)*n_dim+n_cov, n_cov+1);
  Zstar=doubleMatrix((n_samp+s_samp+x0_samp+x1_samp)*n_dim+n_cov, n_cov+1);

  S_W=doubleMatrix(s_samp, n_dim);
  S_Wstar=doubleMatrix(s_samp, n_dim);


  /* bounds */
  minW1=doubleArray(n_samp);
  maxW1=doubleArray(n_samp);
  n_grid=intArray(n_samp);
  resid=doubleArray(n_samp);

  /*bounds condition */
  W1g=doubleMatrix(n_samp, n_step);
  W2g=doubleMatrix(n_samp, n_step);
  prob_grid=doubleArray(n_step);
  prob_grid_cum=doubleArray(n_step);


  /*priors of beta/sigma*/
  beta0=doubleArray(n_cov);
  S0=doubleMatrix(n_dim,n_dim);
  A0=doubleMatrix(n_cov,n_cov);

  /*posteriors of beta*/
  mbeta=doubleArray(n_cov);
  Vbeta=doubleMatrix(n_cov,n_cov);

  /*Normal model */
  beta=doubleArray(n_cov);

 /* mu for Wstar*/
  mu_ord=doubleMatrix((n_samp+s_samp+x1_samp+x0_samp), n_dim); 
  Sigma_ord=doubleMatrix(n_dim,n_dim);
  InvSigma_ord=doubleMatrix(n_dim,n_dim);
  Wstar_bar=doubleArray(n_dim);

  vtemp=doubleArray(n_dim);
  mtemp=doubleMatrix(n_dim,n_dim);
  mtemp2=doubleMatrix(n_dim,n_dim);
  mtemp3=doubleMatrix(n_cov,n_cov);

  /* quantities used in sweep */
  SS=doubleMatrix(n_cov+1, n_cov+1);
  epsilon=doubleArray((n_samp+s_samp+x1_samp+x0_samp)*n_dim);
  R=doubleMatrix(n_dim, n_dim);
 

  /**read prior information*/
  itemp=0;
  for (k=0; k<n_cov; k++) {
    beta0[k]=pdbeta0[k];
    beta[k]=pdbeta0[k];
    for (j=0; j<n_cov; j++) 
      A0[j][k]=pdA0[itemp++];
  }


  itemp=0;
  for (k=0; k<n_dim; k++)
    for (j=0; j<n_dim; j++) 
      S0[j][k]=pdS0[itemp++];
    

  t_samp=n_samp+s_samp+x1_samp+x0_samp;  

  /* read the data set */
  /** Packing Y, X  **/
  itemp = 0;
  for (j = 0; j < n_dim; j++) 
    for (i = 0; i < n_samp; i++) {
      X[i][j] = pdX[itemp++];
    }

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
	Z[i*n_dim+j][k]=pdZ[itemp++];


  /* add prior information to Z*/

  dcholdc(A0, n_cov, mtemp3);  /*Cholesky decomopsition*/

  /*  for (j=0; j<n_dim; j++)
    for (k=0; k<n_dim; k++)
      Rprintf("%5d%5d%14g%14g\n", j, k, A0[j][k], mtemp3[j][k]);
  */


  for (j=0; j<n_cov; j++) {
    Zstar[t_samp*n_dim+j][n_cov]=beta0[j];
    for (k=0; k<n_cov; k++){
      Zstar[t_samp*n_dim+j][n_cov]+=mtemp3[j][k]*beta0[j];
      Zstar[t_samp*n_dim+j][k]=mtemp3[j][k];
    }
  }      
    

  /* initialize W, Wstar for n_samp*/
  for (j=0; j<n_dim; j++)
    for (i=0; i< n_samp; i++) {
      W[i][j]=0;
      if (X[i][1]==0) W[i][j]=0.0001;
      else if (X[i][1]==1) W[i][j]=0.9999;
    }

  for (j=0; j<n_dim; j++)
    Wstar_bar[j]=0;

  for (i=0; i< n_samp; i++) {
    for (j=0; j<n_dim; j++)
      Wstar[i][j]=0;
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

  itempA=0; /* counter for alpha */
  itempS=0; /* counter for storage */
  itempB=0; 
  itempC=0; /* counter to control nth draw */

  /*initialize W and Wstar */
  
  /*initialize W1g and W2g */
  for(i=0; i<n_samp; i++)
    for (j=0; j<n_step; j++){
      W1g[i][j]=0;
      W2g[i][j]=0;
    }

  /*** calculate bounds and grids ***/
  for(i=0;i<n_samp;i++) {
    if (X[i][1]!=0 && X[i][1]!=1) {
      /* min and max for W1 */ 
      minW1[i]=fmax2(0.0, (X[i][0]+X[i][1]-1)/X[i][0]);
      maxW1[i]=fmin2(1.0, X[i][1]/X[i][0]);
      /* number of grid points */
      /* note: 1/n_step is the length of the grid */
      dtemp=(double)1/n_step;
      if ((maxW1[i]-minW1[i]) > (2*dtemp)) { 
	n_grid[i]=ftrunc((maxW1[i]-minW1[i])*n_step);
	resid[i]=(maxW1[i]-minW1[i])-n_grid[i]*dtemp;
	/*if (maxW1[i]-minW1[i]==1) resid[i]=dtemp/4;*/
	j=0; 
	while (j<n_grid[i]) {
	  W1g[i][j]=minW1[i]+(j+1)*dtemp-(dtemp+resid[i])/2;
	  if ((W1g[i][j]-minW1[i])<resid[i]/2) W1g[i][j]+=resid[i]/2;
	  if ((maxW1[i]-W1g[i][j])<resid[i]/2) W1g[i][j]-=resid[i]/2;
	  W2g[i][j]=(X[i][1]-X[i][0]*W1g[i][j])/(1-X[i][0]);
	  /*if (i<20) printf("\n%5d%5d%14g%14g", i, j, W1g[i][j], W2g[i][j]);*/
	  j++;
	}
      }
      else {
	W1g[i][0]=minW1[i]+(maxW1[i]-minW1[i])/3;
	W2g[i][0]=(X[i][1]-X[i][0]*W1g[i][0])/(1-X[i][0]);
	W1g[i][1]=minW1[i]+2*(maxW1[i]-minW1[i])/3;
	W2g[i][1]=(X[i][1]-X[i][0]*W1g[i][1])/(1-X[i][0]);
	n_grid[i]=2;
	
      }
    }
  }

  /* initialize vales of mu_ord and Sigma_ord */
  for(j=0;j<n_dim;j++)
    for(k=0;k<n_dim;k++)
      Sigma_ord[j][k]=S0[j][k];
  
  dinv(Sigma_ord, n_dim, InvSigma_ord);
  /*    for(j=0; j<n_dim; j++)
      for(k=j; k<n_dim; k++)
      Rprintf("%14g", InvSigma_ord[j][k]);*/

  /***Gibbs for  normal prior ***/
  for(main_loop=0; main_loop<*n_gen; main_loop++){

    for (i=0; i<t_samp; i++)
      for (j=0; j<n_dim; j++) 
	mu_ord[i][j]=0;

    /**update W, Wstar given mu, Sigma in regular areas**/
    for (i=0;i<t_samp;i++)
      for (j=0; j<n_dim; j++)
	for (k=0; k<n_cov; k++) 
	  mu_ord[i][j]+=Z[i*n_dim+j][k]*beta[k];

    for (i=0; i<n_samp; i++) {
      if ( X[i][1]!=0 && X[i][1]!=1 ) {
	/*1 project BVN(mu_ord, Sigma_ord) on the inth tomo line */
	dtemp=0;
	for (j=0;j<n_grid[i];j++){
	    vtemp[0]=log(W1g[i][j])-log(1-W1g[i][j]);
	    vtemp[1]=log(W2g[i][j])-log(1-W2g[i][j]);
	    prob_grid[j]=dMVN(vtemp, mu_ord[i], InvSigma_ord, 2, 1) -
	      log(W1g[i][j])-log(W2g[i][j])-log(1-W1g[i][j])-log(1-W2g[i][j]);
	  prob_grid[j]=exp(prob_grid[j]);
	  dtemp+=prob_grid[j];
	  prob_grid_cum[j]=dtemp;
	}
	for (j=0;j<n_grid[i];j++)
	  prob_grid_cum[j]/=dtemp; /*standardize prob.grid */ 
	/*    Rprintf("ok9\n");*/
	/*2 sample W_i on the ith tomo line */
	/*3 compute Wsta_i from W_i*/
	j=0;
	dtemp=unif_rand();
	while (dtemp > prob_grid_cum[j]) j++;
	W[i][0]=W1g[i][j];
	W[i][1]=W2g[i][j];
      } /* end of *1 */

	Wstar[i][0]=log(W[i][0])-log(1-W[i][0]);
	Wstar[i][1]=log(W[i][1])-log(1-W[i][1]);
	Z[i*n_dim][n_cov]=Wstar[i][0];
	Z[i*n_dim+1][n_cov]=Wstar[i][1];
    }

    /*update W2 given W1, mu_ord and Sigma_ord in x1 homeogeneous areas */
    /*printf("W2 draws\n");*/
    if (*x1==1)
      for (i=0; i<x1_samp; i++) {
	dtemp=mu_ord[n_samp+i][1]+Sigma_ord[0][1]/Sigma_ord[0][0]*(Wstar[n_samp+i][0]-mu_ord[n_samp+i][0]);
	dtemp1=Sigma_ord[1][1]*(1-Sigma_ord[0][1]*Sigma_ord[0][1]/(Sigma_ord[0][0]*Sigma_ord[1][1]));
	printf("\n%14g%14g\n", dtemp, dtemp1);
	dtemp1=sqrt(dtemp1);
	Wstar[n_samp+i][1]=rnorm(dtemp, dtemp1);
	W[n_samp+i][1]=exp(Wstar[n_samp+i][1])/(1+exp(Wstar[n_samp+i][1]));
	Z[(i+n_samp)*n_dim][n_cov]=Wstar[(i+n_samp)][0];
	Z[(i+n_samp)*n_dim+1][n_cov]=Wstar[(i+n_samp)][1];
      }

    /*update W1 given W2, mu_ord and Sigma_ord in x0 homeogeneous areas */
    /*printf("W1 draws\n");*/
    if (*x0==1)
      for (i=0; i<x0_samp; i++) {
	dtemp=mu_ord[n_samp+x1_samp+i][0]+Sigma_ord[0][1]/Sigma_ord[1][1]*(Wstar[n_samp+x1_samp+i][1]-mu_ord[n_samp+x1_samp+i][1]);
	dtemp1=Sigma_ord[0][0]*(1-Sigma_ord[0][1]*Sigma_ord[0][1]/(Sigma_ord[0][0]*Sigma_ord[1][1]));
	dtemp1=sqrt(dtemp1);
	Wstar[n_samp+x1_samp+i][0]=rnorm(dtemp, dtemp1);
	W[n_samp+x1_samp+i][0]=exp(Wstar[n_samp+x1_samp+i][0])/(1+exp(Wstar[n_samp+x1_samp+i][0]));
	Z[(i+n_samp+x1_samp)*n_dim][n_cov]=Wstar[(i+n_samp+x1_samp)][0];
	Z[(i+n_samp+x1_samp)*n_dim+1][n_cov]=Wstar[(i+n_samp+x1_samp)][1];
      }
  




    /***regression block, given Wstar, beta/sigma priors*****/
    /*ss=0;
    for (j=0; j<n_dim; j++)
      for (k<0; k<n_dim; k++) 
	mtemp[j][k]=0;
    for (i=0; i<n_dim; i++)
      for (j=0; j<n_dim; j++)
	for (k=0; k<n_dim; k++) 
	  mtemp[j][k]+=S0[j][i]*InvSigma_ord[i][k];
    for (j=0; j<n_dim; j++)
    ss+=mtemp[j][j];*/

    /* mutliply Z and W by the Inverse of hte Cholesky factor */
    /*    for(j=0; j<n_dim; j++)
      for(k=j; k<n_dim; k++)
      Rprintf("%14g", InvSigma_ord[j][k]);*/

    dcholdc(InvSigma_ord, n_dim, mtemp);

    for (i=0; i<t_samp*n_dim; i++)
      for (j=0; j<=n_cov; j++)
	Zstar[i][j]=0;
    for (i=0; i<t_samp; i++)
      for(j=0; j<n_dim; j++)
	for(k=0; k<n_dim; k++)
	  for(l=0; l<=n_cov; l++)
	    Zstar[i*n_dim+k][l]+=mtemp[k][j]*Z[i*n_dim+j][l];


    /*    for (i=0; i<t_samp*n_dim+n_cov; i++){
      Rprintf("\n%5d\n",i);
      for(j=0; j<=n_cov;j++)
	Rprintf("%14g", Z[i][j]);
      Rprintf("\n%5d\n",i);
      for(j=0; j<=n_cov;j++)
	Rprintf("%14g", Zstar[i][j]);
    }
  for (j=0; j<n_dim; j++)
    for (k=0; k<n_dim; k++)
      Rprintf("%5d%5d%14g%14g\n", j, k, InvSigma_ord[j][k], mtemp[j][k]);
    */

    /*construct SS matrix for SWEEP */
    for (j=0; j<=n_cov; j++)
      for (k=0; k<=n_cov; k++)
	SS[j][k]=0;
    /*    for(i=0; i<t_samp; i++)
      for(j=0; j<n_dim; j++)
	for(k=0; k<=n_cov; k++)
	  for(l=0; l<=n_cov; l++)
	  SS[k][l]+=Zstar[i*n_dim+j][k]*Zstar[i*n_dim+j][l]; */

    for(i=0; i<(t_samp*n_dim); i++)
      for(k=0; k<=n_cov; k++)
	for(l=0; l<=n_cov; l++)
	  SS[k][l]+=Zstar[i][k]*Zstar[i][l];

    for(j=0; j<n_cov; j++)
      for(k=0; k<=n_cov; k++)
	for(l=0; l<=n_cov; l++)
	  SS[k][l]+=Zstar[n_samp*n_dim+j][k]*Zstar[n_samp*n_dim+j][l];


    /*SWEEP to get posterior mean anf variance for beta */
    for (j=0; j<=n_cov; j++) 
      SWP(SS,j,n_cov+1);
    Rprintf("ok1\n");
    /*draw alpha2 given Sigma and W */
    /*ss+=SS[n_cov][n_cov];
      alpha2=ss/(double)rchisq((double)(n_samp+nu0)*n_dim);*/

    /*draw beta given Sigma and W */
    for (j=0; j<n_cov; j++) {
      beta[j]=SS[j][n_cov]; /*mbeta[j]=SS[j][n_cov];*/
      for (k=0; k<n_cov; k++)
	Vbeta[j][k]=-SS[j][k];
    }
    /*
    Rprintf("beta\n");
        for (j=0; j<n_cov; j++)
      Rprintf("%14g", mbeta[j]);
    Rprintf("vbeta\n");
    for (k=0; k<n_cov; k++) {
    Rprintf("\n");
      for (j=0; j<n_cov; j++)
      Rprintf("%14g", Vbeta[k][j]);
      }*/

    /*    rMVN(beta, mbeta, Vbeta, n_cov); */

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
    dinv(mtemp, n_dim, mtemp2);
    rWish(InvSigma_ord, mtemp2, nu0+t_samp, n_dim);
    dinv(InvSigma_ord, n_dim, Sigma_ord);


    
    /*store Gibbs draw after burn-in and every nth draws */      
    R_CheckUserInterrupt();
    if (main_loop>=*burn_in){
      itempC++;
      if (itempC==nth){
	for (j=0; j<n_cov; j++)
	  pdSBeta[itempA++]=beta[j];
	for (j=0; j<n_dim; j++)
	  for (k=j; k<n_dim; k++)
	    pdSSigma[itempB++]=Sigma_ord[j][k];
	for(i=0; i<(n_samp+x1_samp+x0_samp); i++){
	  pdSW1[itempS]=W[i][0];
	  pdSW2[itempS]=W[i][1];
	  /*Wstar prediction */
	  if (*pred) {
	    rMVN(vtemp, mu_ord[i], Sigma_ord, n_dim);
	      pdSWt1[itempS]=exp(vtemp[0])/(exp(vtemp[0])+1);
	      pdSWt2[itempS]=exp(vtemp[1])/(exp(vtemp[1])+1);
	  }
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
  free(minW1);
  free(maxW1);
  free(n_grid);
  free(resid);
  FreeMatrix(S0, n_dim);
  FreeMatrix(W1g, n_samp);
  FreeMatrix(W2g, n_samp);
  free(prob_grid);
  free(prob_grid_cum);
  FreeMatrix(mu_ord,t_samp);
  FreeMatrix(Sigma_ord,n_dim);
  FreeMatrix(InvSigma_ord, n_dim);
  FreeMatrix(Z, t_samp*n_dim+n_cov);
  FreeMatrix(Zstar, t_samp*n_dim+n_cov);
  free(Wstar_bar);
  free(vtemp);
  FreeMatrix(mtemp, n_dim);
  FreeMatrix(mtemp2, n_dim);
  FreeMatrix(mtemp3, n_cov);
  free(beta);
  free(beta0);
  FreeMatrix(A0, n_cov);
  FreeMatrix(SS, n_cov+1);
  free(mbeta);
  FreeMatrix(Vbeta, n_cov);
  free(epsilon);
  FreeMatrix(R, n_dim);
  
} /* main */

