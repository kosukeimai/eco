#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"

void cDPecoX(
	    /*data input */
	    double *pdX,     /* data (X, Y) */
	    int *pin_samp,   /* sample size */
	    /*MCMC draws */
	    int *n_gen,      /* number of gibbs draws */ 
	    int *burn_in,    /* number of draws to be burned in */
	    int *pinth,        /* keep every nth draw */
	    int *verbose,    /* 1 for output monitoring */
	    /* prior specification*/
	    int *pinu0,      /* prior df parameter for InvWish */
	    double *pdtau0,  /* prior scale parameter for Sigma under G0*/ 
	    double *mu0,     /* prior mean for mu under G0 (3x1) */
	    double *pdS0,    /* prior scale for Sigma (3x3) */

	    double *alpha0,   /* precision parameter, can be fixed or updated*/
	    int *pinUpdate,  /* 1 if alpha gets updated */
	    double *pda0, double *pdb0, /* prior for alpha if alpha updated*/  

	    /*incorporating survey data */
	    int *survey,      /*1 if survey data available (set of W_1, W_2) */
	                      /*0 otherwise*/
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

	    /*  unused:	    int *link,       one Logit transformation 
				two Probit transformation 
	                        three cloglog transformation */

	    /* storage for Gibbs draws of mu/sigmat*/
	    double *pdSMu0, double *pdSMu1, double *pdSMu2, 
	    double *pdSSig00, double *pdSSig01, double *pdSSig02, 
	    double *pdSSig11, double *pdSSig12, double *pdSSig22,          
	    /* storage for Gibbs draws of W*/
	    double *pdSW1, double *pdSW2,
	    /* storage for posterior predictions of W */
	    double *pdSWt1, double *pdSWt2,
	    /* storage for Gibbs draws of alpha */
	    double *pdSa,
	    /* storage for nstar at each Gibbs draw*/
	    int *pdSn
	    /* unused: storage for posterior predictions of Y 
	    double *pdY,*/
 	    ){	   
  
  int n_samp = *pin_samp;    /* sample size */
  int nu0 = *pinu0;          /* prior parameters */ 
  double tau0 = *pdtau0;   
  double a0=*pda0, b0=*pdb0;  
  int nth=*pinth;  
  int n_dim=2;           /* The number of covariates */

  double **X;	    	 /* The Y and covariates */
  double **S0;           /* The prior S parameter for InvWish */

  double alpha=*alpha0;  /*precision parameter*/
  int s_samp= *sur_samp;   /* sample size of survey data */
  double **S_W;            /*The known W1 and W2 matrix*/
  double **S_Wstar;        /*The inverse logit transformation of S_W*/

  int x1_samp=*sampx1;
  int x0_samp=*sampx0;

  int t_samp; /*effective sample size of W when survey data available*/

  /*bounds condition variables */
  double **W;            /* The W1 and W2 matrix */
  double *minW1, *maxW1; /* The lower and upper bounds of W_1i */
  int n_step=1000;    /* 1/The default size of grid step */  
  int *n_grid;           /* The number of grids for sampling on tomoline */
  double **W1g, **W2g;   /* The grids taken for W1 and W2 on tomoline */
  double *prob_grid;     /* The projected density on tomoline */
  double *prob_grid_cum; /* The projected cumulative density on tomoline */
  double *resid;         /* The centralizing vector for grids */

  /* dirichlet variables */
  double **Wstar;        /* The pseudo data  */
                         /*The unidentified parameters */
  double **Sn;           /* The posterior S parameter for InvWish */
  double *mun;           /* The posterior mean of mu under G0*/
  double ***Sigma;      /* The covarince matrix of psi (nsamp,n_dim+1,n_dim+1)*/
  double ***InvSigma;   /* The inverse of Sigma*/
  double **mu;           /* The mean of psi (nsamp,n_dim+1)  */

  /* conditional distribution of W1 and W2 on X, to draw W1 and W2 marginally */
  double **Sigma_w;      /* The covarince matrix of Wstar|X (n_dim,n_dim)*/
  double **InvSigma_w;   /* The inverse of Sigma_w*/
  double *mu_w;           /* The mean of Wstar|X (n_dim)  */

  int *C;    	         /* vector of cluster membership */
  int nstar;		 /* # of clusters with distict theta values */
  double *q;	         /* Weights used in drawing posterior from Dirichlet */
  double *qq;            /* cumulative dist of weight vector q */
  double **S_tvt;        /* The matrix paramter for TVT in q0 part */

  /* variables defined in remixing step */
  double **Snj;           /* The posterior S parameter for InvWish */
  double *munj;           /* The posterior mean of mu under G0*/
  int nj;                 /* track the number of obs in each cluster */
  int *sortC;             /* track original obs id */
  int *indexC;            /* track original obs id */
  double **Wstarmix;      /* extracted data matrix used in remix step */ 

  int *label;             /* store index values */
  double *mu_mix;         /* store mu update from remixing step */
  double **Sigma_mix;    /* store Sigma update from remixing step */
  double **InvSigma_mix; /* store InvSigma update from remixing step */
  double *Wstar_bar;      /* mean of Psi_nj at remixing step */

  /* misc variables */
  int i, j, k, l, main_loop;   /* used for various loops */
  int itemp, itempS, itempC, itempA;
  int progress = 1, itempP = ftrunc((double) *n_gen/10);
  double dtemp, dtemp1, dtemp2;
  double *vtemp, *vtemp1;
  double **mtemp, **mtemp1;

  /* get random seed */
  GetRNGstate();

  /* defining vectors and matricies */
  /* data */
  X=doubleMatrix(n_samp,n_dim);
  W=doubleMatrix((n_samp+x1_samp+x0_samp+s_samp),n_dim);
  Wstar=doubleMatrix((n_samp+x1_samp+x0_samp+s_samp),(n_dim+1));

  S_W=doubleMatrix(s_samp,n_dim+1);
  S_Wstar=doubleMatrix(s_samp,n_dim+1);


  /* bounds */
  minW1=doubleArray(n_samp);
  maxW1=doubleArray(n_samp);
  n_grid=intArray(n_samp);
  resid=doubleArray(n_samp);

  /*priors*/
  S0=doubleMatrix((n_dim+1),(n_dim+1));

  /*posteriors*/
  Sn=doubleMatrix((n_dim+1),(n_dim+1));
  mun=doubleArray(n_dim+1); 

  /*bounds condition */
  W1g=doubleMatrix(n_samp, n_step);
  W2g=doubleMatrix(n_samp, n_step);
  prob_grid=doubleArray(n_step);
  prob_grid_cum=doubleArray(n_step);

  /*Dirichlet variables*/
  Sigma=doubleMatrix3D((n_samp+x1_samp+x0_samp+s_samp),(n_dim+1),(n_dim+1));
  InvSigma=doubleMatrix3D((n_samp+x1_samp+x0_samp+s_samp),(n_dim+1),(n_dim+1));
  mu=doubleMatrix((n_samp+x1_samp+x0_samp+s_samp),(n_dim+1));

  /*conditional distribution parameter */
  Sigma_w=doubleMatrix(n_dim,n_dim);
  InvSigma_w=doubleMatrix(n_dim,n_dim);
  mu_w=doubleArray(n_dim);

  C=intArray((n_samp+x1_samp+x0_samp+s_samp));
  q=doubleArray((n_samp+x1_samp+x0_samp+s_samp));
  qq=doubleArray((n_samp+x1_samp+x0_samp+s_samp));
  S_tvt=doubleMatrix((n_dim+1),(n_dim+1));

  sortC=intArray((n_samp+x1_samp+x0_samp+s_samp));
  indexC=intArray((n_samp+x1_samp+x0_samp+s_samp));
  Wstarmix=doubleMatrix((n_samp+x1_samp+x0_samp+s_samp),(n_dim+1));
  Snj=doubleMatrix((n_dim+1),(n_dim+1));
  munj=doubleArray(n_dim+1);
 
  Wstar_bar=doubleArray(n_dim+1);

  mu_mix=doubleArray(n_dim+1);
  Sigma_mix=doubleMatrix((n_dim+1),(n_dim+1));
  InvSigma_mix=doubleMatrix((n_dim+1),(n_dim+1));
  label=intArray((n_samp+x1_samp+x0_samp+s_samp));

  vtemp=doubleArray(n_dim);
  vtemp1=doubleArray(n_dim+1);
  mtemp=doubleMatrix((n_dim+1),(n_dim+1));
  mtemp1=doubleMatrix((n_dim+1),(n_dim+1));


  t_samp=n_samp+x1_samp+x0_samp+s_samp;
  Rprintf("ok0\n");
  /* read the data set */
  /** Packing Y, X  **/
  itemp = 0;
  for (j = 0; j < n_dim; j++) 
    for (i = 0; i < n_samp; i++) X[i][j] = pdX[itemp++];


  /* priors under G0*/
  itemp=0;
  for(k=0;k<=n_dim;k++)
    for(j=0;j<=n_dim;j++) {
    S0[j][k]=pdS0[itemp++];
    Rprintf("%14g\n", S0[j][k]);
    }

  for (i=0; i<n_samp; i++) {
    for (j=0; j<n_dim; j++)
      {
	W[i][j]=0;
	Wstar[i][j]=0;
	if (X[i][1]==0) W[i][j]=0.0001;
	else if (X[i][1]==1) W[i][j]=0.9999;
      }
    if (X[i][0]==0) 
      { Wstar[i][2]=log(0.0001/0.9999); }
    else {
      if (X[i][0]==1) 
	{ Wstar[i][2]=log(0.9999/0.0001); }
      else 
	{  Wstar[i][2]=log(X[i][0]/(1-X[i][0])); }
    }
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

  Rprintf("ok1\n");
  /*read the survey data */

  if (*survey==1) {
    itemp = 0;
    for (j=0; j<=n_dim; j++)
      for (i=0; i<s_samp; i++) {
        S_W[i][j]=sur_W[itemp++];
        if (S_W[i][j]==0) S_W[i][j]=0.0001;
        if (S_W[i][j]==1) S_W[i][j]=0.9999;
        S_Wstar[i][j]=log(S_W[i][j])-log(1-S_W[i][j]);
        if (j<n_dim) W[n_samp+x1_samp+x0_samp+i][j]=S_W[i][j];
	Wstar[n_samp+x1_samp+x0_samp+i][j]=S_Wstar[i][j];
      }

  }

  for (i=0; i<t_samp; i++)
    Rprintf("%5d%14g%14g%14g\n", i, Wstar[i][0], Wstar[i][1], Wstar[i][2]);


  itempA=0; /* counter for alpha */
  itempS=0; /* counter for storage */
  itempC=0; /* counter to control nth draw */


  /*initialize W1g and W2g */
  for(i=0; i<n_samp; i++)
    for (j=0; j<n_step; j++){
      W1g[i][j]=0;
      W2g[i][j]=0;
    }

  /*** calculate bounds and grids ***/
  dtemp=(double)1/n_step;
  for(i=0;i<n_samp;i++) {
    if (X[i][1]!=0 && X[i][1]!=1) {
      /* min and max for W1 */
      minW1[i]=fmax2(0.0, (X[i][0]+X[i][1]-1)/X[i][0]);
      maxW1[i]=fmin2(1.0, X[i][1]/X[i][0]);
      /* number of grid points */
      /* note: 1/n_step is the length of the grid */
      if ((maxW1[i]-minW1[i]) > (2*dtemp)) {
        n_grid[i]=ftrunc((maxW1[i]-minW1[i])*n_step);
        resid[i]=(maxW1[i]-minW1[i])-n_grid[i]*dtemp;
        j=0;
        while (j<n_grid[i]) {
          W1g[i][j]=minW1[i]+(j+1)*dtemp-(dtemp+resid[i])/2;
          if ((W1g[i][j]-minW1[i])<resid[i]/2) W1g[i][j]+=resid[i]/2;
	  if ((maxW1[i]-W1g[i][j])<resid[i]/2) W1g[i][j]-=resid[i]/2;
	  W2g[i][j]=(X[i][1]-X[i][0]*W1g[i][j])/(1-X[i][0]);
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

  /* parmeters for Trivaraite t-distribution-unchanged in MCMC */
  for (j=0;j<=n_dim;j++)
    for(k=0;k<=n_dim;k++)
      mtemp[j][k]=tau0*(nu0-1)*S0[j][k]/(1+tau0);
  dinv(mtemp, (n_dim+1), S_tvt);

  Rprintf("S_tvt");
  Rprintf("%14g%14g%14g%14g%14g%14g\n", S_tvt[0][0], S_tvt[0][1], S_tvt[0][2], S_tvt[1][1], S_tvt[1][2], S_tvt[2][2]);

  /**draw initial values of mu_i, Sigma_i under G0  for all effective sample**/
  /*1. Sigma_i under InvWish(nu0, S0^-1) with E(Sigma)=S0/(nu0-3)*/
  /*   InvSigma_i under Wish(nu0, S0^-1 */
  /*2. mu_i|Sigma_i under N(mu0, Sigma_i/tau0) */
  dinv(S0, (n_dim+1), mtemp);
  Rprintf("S0inverse\n");
  Rprintf("%14g%14g%14g%14g%14g%14g\n", mtemp[0][0], mtemp[0][1], mtemp[0][2], mtemp[1][1], mtemp[1][2], mtemp[2][2]);
  for(i=0;i<t_samp;i++){
    /*draw from wish(nu0, S0^-1) */
    rWish(InvSigma[i], mtemp, nu0, (n_dim+1));
    dinv(InvSigma[i], (n_dim+1), Sigma[i]);
    for (j=0;j<=n_dim;j++)
      for(k=0;k<=n_dim;k++) mtemp1[j][k]=Sigma[i][j][k]/tau0;
    rMVN(mu[i], mu0, mtemp1, (n_dim+1));
  }

  
  Rprintf("ok3a\n");
  /* initialize the cluster membership */
  nstar=t_samp;  /* the # of disticnt values */
  for(i=0;i<t_samp;i++)
    C[i]=i; /*cluster is from 0...n_samp-1 */

  for (i=0; i<60; i++)
    Rprintf("%5d%14g%14g%14g%14g%14g%14g\n",i, mu[i][0], mu[i][1], mu[i][2], Sigma[i][0][0], Sigma[i][1][1], Sigma[i][2][2]);

  for(main_loop=0; main_loop<*n_gen; main_loop++){
    /**update W, Wstar given mu, Sigma only for the unknown W/Wstar**/
    for (i=0;i<t_samp;i++){
      for (j=0; j<n_dim; j++) {
        mu_w[j]=mu[i][j]+Sigma[i][n_dim][j]/Sigma[i][n_dim][n_dim]*(Wstar[i][n_dim]-mu[i][n_dim]);
      }
      for (j=0; j<n_dim; j++)
        for (k=0; k<n_dim; k++) {
          Sigma_w[j][k]=Sigma[i][j][k]-Sigma[i][n_dim][j]/Sigma[i][n_dim][n_dim]*Sigma[i][n_dim][k];
	}
      dinv(Sigma_w, n_dim, InvSigma_w);

      if (X[i][1]!=0 && X[i][1]!=1 && i<n_samp) {

        /*1 project BVN(mu_i, Sigma_i) on the inth tomo line */
        dtemp=0;
        for (j=0;j<n_grid[i];j++){
          /*  if (*link==1){*/
	  vtemp[0]=log(W1g[i][j])-log(1-W1g[i][j]);
	  vtemp[1]=log(W2g[i][j])-log(1-W2g[i][j]);
	  prob_grid[j]=dMVN(vtemp, mu_w, InvSigma_w, n_dim, 1) -
	    log(W1g[i][j])-log(W2g[i][j])-log(1-W1g[i][j])-log(1-W2g[i][j]);
	 
          /*{
          else if (*link==2){
            vtemp[0]=qnorm(W1g[i][j], 0, 1, 1, 0);
            vtemp[1]=qnorm(W2g[i][j], 0, 1, 1, 0);
            prob_grid[j]=dMVN(vtemp, mu[i], InvSigma[i], 2, 1) -
              dnorm(vtemp[0], 0, 1, 1)-dnorm(vtemp[1], 0, 1, 1);
          }
          else if (*link==3) {
            vtemp[0]=-log(-log(W1g[i][j]));
            vtemp[1]=-log(-log(W2g[i][j]));
            prob_grid[j]=dMVN(vtemp, mu[i], InvSigma[i], 2, 1) -
              log(W1g[i][j])-log(W2g[i][j])-log(-log(W1g[i][j]))-log(-log(W2g[i][j]));
              }*/
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
        while ((dtemp > prob_grid_cum[j]) && (j<(n_grid[i]-1))) j++;
        W[i][0]=W1g[i][j];
        W[i][1]=W2g[i][j];
      }
      /*      if (*link==1) {*/
      Wstar[i][0]=log(W[i][0])-log(1-W[i][0]);
      Wstar[i][1]=log(W[i][1])-log(1-W[i][1]);
      /* }
      else if (*link==2) {
        Wstar[i][0]=qnorm(W[i][0],0 ,1, 1, 0);
        Wstar[i][1]=qnorm(W[i][1],0 ,1, 1, 0);
      }
      else if (*link==3) {
        Wstar[i][0]=-log(-log(W[i][0]));
        Wstar[i][1]=-log(-log(W[i][1]));
        }*/


      if (*x1==1 && i>=n_samp && i<(n_samp+x1_samp)) {
    Rprintf("ok3c\n");
  Rprintf("%5d\n", i);
	dtemp=mu_w[1]+Sigma_w[0][1]/Sigma_w[0][0]*(Wstar[i][0]-mu_w[0]);
	dtemp1=Sigma_w[1][1]*(1-Sigma_w[0][1]*Sigma_w[0][1]/(Sigma_w[0][0]*Sigma_w[1][1]));
	/*dtemp1=sqrt(dtemp1);
	  Wstar[n_samp+i][1]=rnorm(dtemp, dtemp1);*/
	Wstar[i][1]=norm_rand()*sqrt(dtemp1)+dtemp;
	W[i][1]=exp(Wstar[i][1])/(1+exp(Wstar[i][1]));
      }

      /*update W1 given W2, mu_ord and Sigma_ord in x0 homeogeneous areas */
      /*printf("W1 draws\n");*/
      if (*x0==1  && i>=(n_samp+x1_samp) && i<(n_samp+x1_samp+x0_samp)) {
  Rprintf("ok3d\n");
  Rprintf("%5d\n", i);
        dtemp=mu_w[0]+Sigma_w[0][1]/Sigma_w[1][1]*(Wstar[i][1]-mu_w[1]);
        dtemp1=Sigma_w[0][0]*(1-Sigma_w[0][1]*Sigma_w[0][1]/(Sigma_w[0][0]*Sigma_w[1][1]));
        Wstar[i][0]=norm_rand()*sqrt(dtemp1)+dtemp;
        W[i][0]=exp(Wstar[i][0])/(1+exp(Wstar[i][0]));
      }
    }


  /**updating mu, Sigma given Wstar uisng effective sample size t_samp**/
  for (i=0; i<t_samp; i++){

       /* generate weight vector q */
       dtemp=0;
       for (j=0; j<t_samp; j++){
	 if (j!=i)
	   q[j]=dMVN(Wstar[i], mu[j], InvSigma[j], (n_dim+1), 0);
	 else
	   q[j]=alpha*dMVT(Wstar[i], mu0, S_tvt, (nu0-1), (n_dim+1), 0);
	 dtemp+=q[j];
	 qq[j]=dtemp;    /*compute qq, the cumlative of q*/
       }
       /*standardize q and qq */
       for (j=0; j<t_samp; j++) {
	 qq[j]/=dtemp;
       }

       /** draw the configuration parameter **/
       /* j=i means to draw from posterior baseline */
       /* j=i' means to replace with i' obs */
       j=0; dtemp=unif_rand();
       while (dtemp > qq[j]) j++;
       /** Dirichlet update Sigma_i, mu_i|Sigma_i **/
       if (j==i){
	 /*1. draw Sigma_i from Invwish(nun, Sn^-1) with E(Sigma_i)=Sn/(nun-3)*/
	 /*   draw InvSigma_i from Wish(nun, Sn^-1) */
	 /*2. draw mu_i from N(mun, Sigma_i) with Sigma_i/=taun */
	 /*3. also update c[i]=nstar*/
	 for (k=0; k<=n_dim; k++)
	   for (l=0; l<=n_dim; l++)
	     Sn[k][l]=S0[k][l]+tau0*(Wstar[i][k]-mu0[k])*(Wstar[i][l]-mu0[l])/(tau0+1);
	 dinv(Sn, (n_dim+1), mtemp);
	 rWish(InvSigma[i], mtemp, (nu0+1), (n_dim+1));
	 dinv(InvSigma[i], (n_dim+1), Sigma[i]);
	 for (k=0;k<=n_dim;k++){
	   mun[k]=(tau0*mu0[k]+Wstar[i][k])/(tau0+1);
	   for (l=0;l<=n_dim;l++)  mtemp[k][l]=Sigma[i][k][l]/(tau0+1);
	 }
	 rMVN(mu[i], mun, mtemp, (n_dim+1));
	 C[i]=nstar;
	 nstar++;
       }
       else {
	 /*1. mu_i=mu_j, Sigma_i=Sigma_j*/
	 /*2. update C[i]=C[j] */
	 for(k=0;k<=n_dim;k++) {
	   mu[i][k]=mu[j][k];
	   for(l=0;l<=n_dim;l++) {
	     Sigma[i][k][l]=Sigma[j][k][l];
	     InvSigma[i][k][l]=InvSigma[j][k][l];
	   }
	 }
	 C[i]=C[j];
       }
       sortC[i]=C[i];
  } /* end of i loop*/
  /** remixing step using effective sample**/
  for(i=0;i<t_samp;i++)
    indexC[i]=i;
  R_qsort_int_I(sortC, indexC, 1, t_samp);

  nstar=0;
  i=0;
  while (i<t_samp){
    /*initialize the vector and matrix */
    for(k=0; k<=n_dim; k++) {
      Wstar_bar[k]=0;
      for(l=0;l<=n_dim;l++) Snj[k][l]=S0[k][l];
    }

    j=sortC[i]; /*saves the first element in a block of same values */
    nj=0; /* counter for a block of same values */

    /* get data for remixing */
    while ((sortC[i]==j) && (i<t_samp)) {
      label[nj]=indexC[i];
      for (k=0; k<=n_dim; k++) {
	Wstarmix[nj][k]=Wstar[label[nj]][k];
	Wstar_bar[k]+=Wstarmix[nj][k];
      }
      nj++;
      i++;
    } /* i records the current position in IndexC */
    /* nj records the # of obs in Psimix */

    /** posterior update for mu_mix, Sigma_mix based on Psimix **/
    for (k=0; k<=n_dim; k++)
      Wstar_bar[k]/=nj;
    /* compute Snj */
    /*1. first two terms in Snj */
    /* Snj=S0 */
    for (j=0; j<nj; j++)
      for (k=0; k<=n_dim; k++)
	for (l=0; l<=n_dim; l++)
	  Snj[k][l]+=(Wstarmix[j][k]-Wstar_bar[k])*(Wstarmix[j][l]-Wstar_bar[l]);

    /*2. plus third term in Snj*/
    for (k=0;k<=n_dim;k++)
      for (l=0;l<=n_dim;l++)
	Snj[k][l]+=tau0*nj*(Wstar_bar[k]-mu0[k])*(Wstar_bar[l]-mu0[l])/(tau0+nj);


    /*darw unscaled InvSigma_mix ~Wish(Snj^-1) */
    dinv(Snj, (n_dim+1), mtemp);
    rWish(InvSigma_mix, mtemp, (nu0+nj), (n_dim+1));
    dinv(InvSigma_mix, (n_dim+1), Sigma_mix);


    /*2. draw mu_mix from N(munj, Sigma_mix) with Sigma_mix/=(tau0+nj) */
    for (k=0; k<=n_dim; k++){
      munj[k]=(tau0*mu0[k]+nj*Wstar_bar[k])/(tau0+nj);
      for (l=0; l<=n_dim; l++)  mtemp[k][l]=Sigma_mix[k][l]/(tau0+nj);
    }
    rMVN(mu_mix, munj, mtemp, (n_dim+1));

    /**update mu, Simgat with mu_mix, Sigmat_mix via label**/
    for (j=0;j<nj;j++){
      C[label[j]]=nstar;  /*updating C vector with no gap */
      for (k=0; k<=n_dim; k++){
	mu[label[j]][k]=mu_mix[k];
	for (l=0;l<=n_dim;l++){
	  Sigma[label[j]][k][l]=Sigma_mix[k][l];
	  InvSigma[label[j]][k][l]=InvSigma_mix[k][l];
	}
      }
    }
    nstar++; /*finish update one distinct value*/
  } /* nstar is the number of distinct values */
  Rprintf("ok5\n");
  /** updating alpha **/
  if(*pinUpdate) {
    dtemp1=(double)(alpha+1);
    dtemp2=(double)t_samp;
    dtemp=b0-log(rbeta(dtemp1, dtemp2));

    dtemp1=(double)(a0+nstar-1)/(t_samp*dtemp);
    if(unif_rand() < dtemp1) {
      dtemp2=(double)(a0+nstar);
      alpha=rgamma(dtemp2, 1/dtemp);
    }
    else {
      dtemp2=(double)(a0+nstar-1);
      alpha=rgamma(dtemp2, 1/dtemp);
    }
  }

  /*store Gibbs draws after burn_in */
  R_CheckUserInterrupt();
  if (main_loop>=*burn_in) {
    itempC++;
    if (itempC==nth){
      if(*pinUpdate) {
	pdSa[itempA]=alpha;
	pdSn[itempA]=nstar;
	itempA++;
      }

      for(i=0; i<(n_samp+x1_samp+x0_samp); i++) {
	pdSMu0[itempS]=mu[i][0];
	pdSMu1[itempS]=mu[i][1];
	pdSMu2[itempS]=mu[i][2];
	pdSSig00[itempS]=Sigma[i][0][0];
	pdSSig01[itempS]=Sigma[i][0][1];
	pdSSig02[itempS]=Sigma[i][0][2];
	pdSSig11[itempS]=Sigma[i][1][1];
	pdSSig12[itempS]=Sigma[i][1][2];
	pdSSig22[itempS]=Sigma[i][2][2];
	pdSW1[itempS]=W[i][0];
	pdSW2[itempS]=W[i][1];
	/* Wstar prediction */
	if (*pred) {
	  rMVN(vtemp1, mu[i], Sigma[i], n_dim+1);
	  /*  if (*link==1){ */
	  pdSWt1[itempS]=exp(vtemp1[0])/(exp(vtemp1[0])+1);
	  pdSWt2[itempS]=exp(vtemp1[1])/(exp(vtemp1[1])+1);
	  /* }
            else if (*link==2){
              pdSWt1[itempS]=pnorm(vtemp[0], 0, 1, 1, 0);
              pdSWt2[itempS]=pnorm(vtemp[1], 0, 1, 1, 0);
            }
           else if (*link==3){
              pdSWt1[itempS]=exp(-exp(-vtemp[0]));
              pdSWt2[itempS]=exp(-exp(-vtemp[1]));
              }*/
	}
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
  FreeMatrix(S0, n_dim+1);
  FreeMatrix(Sn, n_dim+1);
  free(mun);
  FreeMatrix(W1g, n_samp);
  FreeMatrix(W2g, n_samp);
  free(prob_grid);
  free(prob_grid_cum);
  Free3DMatrix(Sigma, t_samp,n_dim+1);
  Free3DMatrix(InvSigma, t_samp, n_dim+1);
  FreeMatrix(mu, t_samp);
  free(mu_w);
  FreeMatrix(Sigma_w, n_dim);
  FreeMatrix(InvSigma_w, n_dim);
  free(C);
  free(q);
  free(qq);
  FreeMatrix(S_tvt, n_dim+1);
  free(sortC);
  free(indexC);
  FreeMatrix(Wstarmix, t_samp);
  FreeMatrix(Snj, n_dim+1);
  free(munj);
  free(Wstar_bar);
  free(mu_mix);
  FreeMatrix(Sigma_mix, n_dim+1);
  FreeMatrix(InvSigma_mix, n_dim+1);
  free(label);
  free(vtemp);
  free(vtemp1);
  FreeMatrix(mtemp, n_dim+1);
  FreeMatrix(mtemp1, n_dim+1);
  FreeMatrix(S_W, s_samp);
  FreeMatrix(S_Wstar, s_samp);
} /* main */


