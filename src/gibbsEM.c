#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"

void cEMeco(
	      /*data input */
	      double *pdX,     /* data (X, Y) */
              double *pdTheta_in,  /* Theta^ t */

	      int *pin_samp,   /* sample size */
	      /*MCMC draws */
	      int *n_gen,      /* number of gibbs draws */

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
	      double *pdTheta  /*EM result for Theta^(t+1) */
	      ){	   
  
  int n_samp = *pin_samp;    /* sample size */

  int data=0;            /* one to print the data */
  int keep=1;            /* keeps every #num draw */ 
  int n_cov=2;           /* The number of covariates */

  double **X;	    	 /* The Y and covariates */
  
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

  /* ordinary model variables */
  double **Sigma_ord;   /* The posterior covariance matrix of psi (oridinary)*/
  double **InvSigma_ord;
  double *mu_ord;        /* The posterior mean of psi (ordinary)*/

  double **Wstar;        /* The pseudo data  */
  double *Wstar_bar;


  /* misc variables */
  int i, j, k, l, main_loop;   /* used for various loops */
  int itemp, itempS, itempC, itempA;
  double dtemp, dtemp1, temp0, temp1;
  double *vtemp;
  double **mtemp;

  /* get random seed */
  GetRNGstate();

  /* defining vectors and matricies */
  /* data */
  X=doubleMatrix(n_samp,n_cov);
  W=doubleMatrix((n_samp+s_samp+x0_samp+x1_samp),n_cov);
  Wstar=doubleMatrix((n_samp+s_samp+x0_samp+x1_samp),n_cov);

  S_W=doubleMatrix(s_samp, n_cov);
  S_Wstar=doubleMatrix(s_samp, n_cov);


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

  /*ordinary model */
  mu_ord=doubleArray(n_cov);
  Sigma_ord=doubleMatrix(n_cov,n_cov);
  InvSigma_ord=doubleMatrix(n_cov,n_cov);

  Wstar_bar=doubleArray(n_cov);

  vtemp=doubleArray(n_cov);
  mtemp=doubleMatrix(n_cov,n_cov);

  t_samp=n_samp+s_samp+x1_samp+x0_samp;  



  /* read the data set */
  /** Packing Y, X  **/
  itemp = 0;
  for (j = 0; j < n_cov; j++) 
    for (i = 0; i < n_samp; i++) {
      X[i][j] = pdX[itemp++];
    }
  if(data==1) { 
    printf("Y X\n");
    for(i=0;i<n_samp;i++)
      printf("%5d%14g%14g\n",i,X[i][1],X[i][0]);
    fflush(stdout);
  }

  /** assigm mu_org, Simga_ord **/
  mu_ord[0] = pdTheta_in[0];
  mu_ord[1] = pdTheta_in[1];
  Sigma_ord[0][0] = pdTheta_in[2];
  Sigma_ord[0][1] = pdTheta_in[3];
  Sigma_ord[1][0] = pdTheta_in[3];
  Sigma_ord[1][1] = pdTheta_in[4];

  dinv(Sigma_ord, n_cov, InvSigma_ord);


  /* initialize W, Wstar for n_samp*/
  for (j=0; j<n_cov; j++)
    for (i=0; i< n_samp; i++) {
      W[i][j]=0;
      Wstar[i][j]=0;
      if (X[i][1]==0) W[i][j]=0.000001;
      else if (X[i][1]==1) W[i][j]=0.999999;

    }


  /*read homeogenous areas information */
  if (*x1==1) 
    for (i=0; i<x1_samp; i++) {
      W[(n_samp+i)][0]=x1_W1[i];
      if (W[(n_samp+i)][0]==0) W[(n_samp+i)][0]=0.000001;
      if (W[(n_samp+i)][0]==1) W[(n_samp+i)][0]=0.999999;
    }

  if (*x0==1) 
    for (i=0; i<x0_samp; i++) {
      W[(n_samp+x1_samp+i)][1]=x0_W2[i];
      if (W[(n_samp+x1_samp+i)][1]==0) W[(n_samp+x1_samp+i)][1]=0.000001;
      if (W[(n_samp+x1_samp+i)][1]==1) W[(n_samp+x1_samp+i)][1]=0.999999;
    }


  /*read the survey data */

  if (*survey==1) {
    itemp = 0;
    for (j=0; j<n_cov; j++)
      for (i=0; i<s_samp; i++) {
	S_W[i][j]=sur_W[itemp++];
	if (S_W[i][j]==0) S_W[i][j]=0.000001;
	if (S_W[i][j]==1) S_W[i][j]=0.999999;
	S_Wstar[i][j]=log(S_W[i][j])-log(1-S_W[i][j]);
	W[(n_samp+x1_samp+x0_samp+i)][j]=S_W[i][j];
	Wstar[(n_samp+x1_samp+x0_samp+i)][j]=S_Wstar[i][j];
      }


    if(data==1) { 
      printf("survey W1 W2 W1* W2*\n");
      for(i=0;i<(n_samp+s_samp+x1_samp+x0_samp);i++)
	printf("%5d%14g%14g%14g%14g\n",i,W[i][0],W[i][1],Wstar[i][0], Wstar[i][1]);
      fflush(stdout);
    }
  }

  itempA=0; /* counter for alpha */
  itempS=0; /* counter for storage */
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
      /*    if (i<n_samp){
	    printf("grids\n");
	    printf("minW1 maxW1 resid\n");
	    printf("%5d%14g%14g%14g\n", i, minW1[i], maxW1[i], resid[i]);
	    for (j=0;j<n_grid[i];j++){
	    if (j<5 | j>(n_grid[i]-5))
	    printf("%5d%5d%14g%14g\n", i, j, W1g[i][j], W2g[i][j]);
	    }
	    }*/
    }
  }
    
  
  /***Gibbs for  normal prior ***/
  for(main_loop=0; main_loop<*n_gen; main_loop++){
    /**update W, Wstar given mu, Sigma in regular areas**/
    for (i=0;i<n_samp;i++){
      if ( X[i][1]!=0 && X[i][1]!=1 ) {
	/*1 project BVN(mu_ord, Sigma_ord) on the inth tomo line */
	dtemp=0;
	for (j=0;j<n_grid[i];j++){

	    vtemp[0]=log(W1g[i][j])-log(1-W1g[i][j]);
	    vtemp[1]=log(W2g[i][j])-log(1-W2g[i][j]);
	    prob_grid[j]=dMVN(vtemp, mu_ord, InvSigma_ord, 2, 1) -
	      log(W1g[i][j])-log(W2g[i][j])-log(1-W1g[i][j])-log(1-W2g[i][j]);

	  prob_grid[j]=exp(prob_grid[j]);
	  dtemp+=prob_grid[j];
	  prob_grid_cum[j]=dtemp;
	}
	for (j=0;j<n_grid[i];j++)
	  prob_grid_cum[j]/=dtemp; /*standardize prob.grid */ 
	
	/*2 sample W_i on the ith tomo line */
	/*3 compute Wsta_i from W_i*/
	j=0;
	dtemp=unif_rand();
	while (dtemp > prob_grid_cum[j]) j++;
	W[i][0]=W1g[i][j];
	W[i][1]=W2g[i][j];
      } /* end of *1 */
        temp0=log(W[i][0])-log(1-W[i][0]);
        temp1=log(W[i][1])-log(1-W[i][1]);
    	Wstar[i][0]+=temp0;
	Wstar[i][1]+=temp1;
    }

    
    /*update W2 given W1, mu_ord and Sigma_ord in x1 homeogeneous areas */
    /*printf("W2 draws\n");*/
    if (*x1==1)
      for (i=0; i<x1_samp; i++) {
	dtemp=mu_ord[1]+Sigma_ord[0][1]/Sigma_ord[0][0]*(Wstar[n_samp+i][0]-mu_ord[0]);
	dtemp1=Sigma_ord[1][1]*(1-Sigma_ord[0][1]*Sigma_ord[0][1]/(Sigma_ord[0][0]*Sigma_ord[1][1]));
	dtemp1=sqrt(dtemp1);
        temp0=log(W[i][0])-log(1-W[i][0]);
        temp1=rnorm(dtemp, dtemp1);
    	Wstar[n_samp+i][0]+=temp0;
	Wstar[n_samp+i][1]+=temp1;
       }
    
    /*update W1 given W2, mu_ord and Sigma_ord in x0 homeogeneous areas */
    /*printf("W1 draws\n");*/
    if (*x0==1)
      for (i=0; i<x0_samp; i++) {
	dtemp=mu_ord[0]+Sigma_ord[0][1]/Sigma_ord[1][1]*(Wstar[n_samp+x1_samp+i][1]-mu_ord[1]);
	dtemp1=Sigma_ord[0][0]*(1-Sigma_ord[0][1]*Sigma_ord[0][1]/(Sigma_ord[0][0]*Sigma_ord[1][1]));
	/* printf("\n%14g%14g\n", dtemp, dtemp1);*/
	dtemp1=sqrt(dtemp1);
        temp0=rnorm(dtemp, dtemp1);
        temp1=log(W[i][1])-log(1-W[i][1]);
        Wstar[n_samp+x1_samp+i][0]+=temp0;
        Wstar[n_samp+x1_samp+i][1]+=temp1;
 	/*printf("\n%5d%14g%14g\n", i, Wstar[n_samp+x1_samp+i][0], W[n_samp+x1_samp+i][0]);*/ 
      }
    
  }

  /* compute E_{W_i|Y_i} */
  for (i=0; i<(n_samp+x1_samp+x0_samp); i++)
    {
        Wstar[i][0]/=*n_gen;
        Wstar[i][1]/=*n_gen;
    } /*for survey data, E-step is the same value*/
   
  /*M-step: same procedure as normal model */
  for (j=0; j<5; j++) 
    {
      pdTheta[j]=0;
    }  

  for (i=0; i<t_samp; i++)
    {
      pdTheta[0]+=Wstar[i][0]/t_samp;
      pdTheta[1]+=Wstar[i][1]/t_samp;
    }   

  for(i=0; i<t_samp; i++)
    {
      pdTheta[2]+=(Wstar[i][0]-pdTheta[0])*(Wstar[i][0]-pdTheta[0])/t_samp;
      pdTheta[3]+=(Wstar[i][0]-pdTheta[0])*(Wstar[i][1]-pdTheta[1])/t_samp;
      pdTheta[4]+=(Wstar[i][1]-pdTheta[1])*(Wstar[i][1]-pdTheta[1])/t_samp;
    }


  for (j=0; j<5; j++)
    {
      printf("\n%14g", pdTheta[j]);
    }

  /** write out the random seed **/
  PutRNGstate();

  /* Freeing the memory */
  FreeMatrix(X, n_samp);
  FreeMatrix(W, t_samp);
  FreeMatrix(Wstar, t_samp);
  free(minW1);
  free(maxW1);
  free(n_grid);
  free(resid);
  FreeMatrix(W1g, n_samp);
  FreeMatrix(W2g, n_samp);
  free(prob_grid);
  free(prob_grid_cum);
  free(mu_ord);
  FreeMatrix(Sigma_ord,n_cov);
  FreeMatrix(InvSigma_ord, n_cov);
  free(Wstar_bar);
  free(vtemp);
  FreeMatrix(mtemp, n_cov);
  
} /* main */

