#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"

void cDPeco(int *n_gen,      /* number of gibbs draws */ 
	    int *link,       /* one Logit transformation 
				two Probit transformation 
	                        three cloglog transformation */
	    double *pdX,     /* data (X, Y) */
	    int *pin_samp,   /* sample size */
	    int *pinu0,      /* prior df parameter for InvWish */
	    double *pdtau0,  /* prior scale parameter for Sigma under G0*/ 
	    double *mu0,     /* prior mean for mu under G0 */
	    double *pdS0,    /* prior scale for Sigma */
	    int *pinUpdate,  /* 1 if alpha gets updated */
	    int *burn_in,    /* number of draws to be burned in */
	    int *pinth,        /* keep every nth draw */
	    int *pred,       /* 1 if draw posterior prediction */
	    double *pda0, double *pdb0, /* prior for alpha */  
	    /* storage for Gibbs draws of mu/sigmat*/
	    double *pdSMu0, double *pdSMu1, 
	    double *pdSSig00, double *pdSSig01, double *pdSSig11,           
	    /* storage for Gibbs draws of W*/
	    double *pdSW1, double *pdSW2,
	    /* storage for posterior predictions of W */
	    double *pdSWt1, double *pdSWt2,
	    /* storage for Gibbs draws of alpha */
	    double *pdSa,
	    /* storage for nstar at each Gibbs draw*/
	    int *pdSn
 	    ){	   
  
  int n_samp = *pin_samp;    /* sample size */
  int nu0 = *pinu0;          /* prior parameters */ 
  double tau0 = *pdtau0;   
  double a0=*pda0, b0=*pdb0;  
  int nth=*pinth;  

  int data=1;            /* one to print the data */
  int keep=1;            /* keeps every #num draw */ 
  int n_cov=2;           /* The number of covariates */

  double **X;	    	 /* The Y and covariates */
  double **S0;           /* The prior S parameter for InvWish */

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
  double ***Sigma;      /* The covarince matrix of psi (nsamp,cov,cov)*/
  double ***InvSigma;   /* The inverse of Sigma*/
  double **mu;           /* The mean of psi (nsamp,cov)  */

  int *C;    	         /* vector of cluster membership */
  int nstar;		 /* # of clusters with distict theta values */
  double *q;	         /* Weights used in drawing posterior from Dirichlet */
  double *qq;            /* cumulative dist of weight vector q */
  double alpha=1;        /* the precision parameter for Dirichlet */
  double **S_bvt;        /* The matrix paramter for BVT in q0 part */

  /* variables defined in remixing step */
  double **Snj;           /* The posterior S parameter for InvWish */
  double *munj;           /* The posterior mean of mu under G0*/
  int nj;                 /* track the number of obs in each cluster */
  int *sortC;             /* track original obs id */
  int *indexC;            /* track original obs id */
  double **Wstarmix;      /* extracted data matrix used in remix step */ 

  int which_obs;          /* pass index values */
  int *label;             /* store index values */
  double *mu_mix;         /* store mu update from remixing step */
  double **Sigma_mix;    /* store Sigma update from remixing step */
  double **InvSigma_mix; /* store InvSigma update from remixing step */
  double *Wstar_bar;      /* mean of Psi_nj at remixing step */

  /* misc variables */
  int i, j, k, l, main_loop;   /* used for various loops */
  int itemp, itempS, itempC, itempA;
  double dtemp;
  double *vtemp;
  double **mtemp, **mtemp1;

  /* get random seed */
  GetRNGstate();

  /* defining vectors and matricies */
  /* data */
  X=doubleMatrix(n_samp,n_cov);
  W=doubleMatrix(n_samp,n_cov);
  Wstar=doubleMatrix(n_samp,n_cov);

  /* bounds */
  minW1=doubleArray(n_samp);
  maxW1=doubleArray(n_samp);
  n_grid=intArray(n_samp);
  resid=doubleArray(n_samp);

  /*priors*/
  S0=doubleMatrix(n_cov,n_cov);

  /*posteriors*/
  Sn=doubleMatrix(n_cov,n_cov);
  mun=doubleArray(n_cov); 

  /*bounds condition */
  W1g=doubleMatrix(n_samp, n_step);
  W2g=doubleMatrix(n_samp, n_step);
  prob_grid=doubleArray(n_step);
  prob_grid_cum=doubleArray(n_step);

  /*Dirichlet variables*/
  Sigma=doubleMatrix3D(n_samp,n_cov,n_cov);
  InvSigma=doubleMatrix3D(n_samp,n_cov,n_cov);
  mu=doubleMatrix(n_samp,n_cov);

  C=intArray(n_samp);
  q=doubleArray(n_samp);
  qq=doubleArray(n_samp);
  S_bvt=doubleMatrix(n_cov,n_cov);

  sortC=intArray(n_samp);
  indexC=intArray(n_samp);
  Wstarmix=doubleMatrix(n_samp,n_cov);
  Snj=doubleMatrix(n_cov,n_cov);
  munj=doubleArray(n_cov);
 
  Wstar_bar=doubleArray(n_cov);

  mu_mix=doubleArray(n_cov);
  Sigma_mix=doubleMatrix(n_cov,n_cov);
  InvSigma_mix=doubleMatrix(n_cov,n_cov);
  label=intArray(n_samp);

  vtemp=doubleArray(n_cov);
  mtemp=doubleMatrix(n_cov,n_cov);
  mtemp1=doubleMatrix(n_cov,n_cov);

  /* priors under G0*/
  itemp=0;
  for(k=0;k<n_cov;k++)
    for(j=0;j<n_cov;j++) S0[j][k]=pdS0[itemp++];
  
  /* read the data set */
  /** Packing Y, X  **/
  itemp = 0;
  for (j = 0; j < n_cov; j++) 
    for (i = 0; i < n_samp; i++) X[i][j] = pdX[itemp++];
  if(data==1) { 
    printf("Y X\n");
    for(i=0;i<n_samp;i++)
      printf("%5d%14g%14g\n",i,X[i][1],X[i][0]);
    fflush(stdout);
  }

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
  for(i=0;i<n_samp;i++) {
    /* min and max for W1 */ 
    minW1[i]=fmax2(0.0, (X[i][0]+X[i][1]-1)/X[i][0]);
    maxW1[i]=fmin2(1.0, X[i][1]/X[i][0]);
    /* number of grid points */
    /* note: 1/n_step is the length of the grid */
    if ((maxW1[i]-minW1[i]) > (double)(2/n_step)) { 
      n_grid[i]=ftrunc((maxW1[i]-minW1[i])*n_step);
      resid[i]=(maxW1[i]-minW1[i])-(double)n_grid[i]/n_step;
      j=0; 
      while (j<n_grid[i]) {
	W1g[i][j]=minW1[i]+(double)(j+1)/n_step-((double)1/n_step+resid[i])/2;
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

  /* parmeters for Bivaraite t-distribution-unchanged in MCMC */
  for (j=0;j<n_cov;j++)
    for(k=0;k<n_cov;k++)
      mtemp[j][k]=tau0*(nu0-1)*S0[j][k]/(1+tau0);
  dinv(mtemp, n_cov, S_bvt); 
  
  /**draw initial values of mu_i, Sigma_i under G0 **/
  /*1. Sigma_i under InvWish(nu0, S0^-1) with E(Sigma)=S0/(nu0-3)*/
  /*   InvSigma_i under Wish(nu0, S0^-1 */
  /*2. mu_i|Sigma_i under N(mu0, Sigma_i/tau0) */
  dinv(S0, n_cov, mtemp);
  for(i=0;i<n_samp;i++){
    /*draw from wish(nu0, S0^-1) */
    rWish(InvSigma[i], mtemp, nu0, n_cov); 
    dinv(InvSigma[i], n_cov, Sigma[i]);
    for (j=0;j<n_cov;j++)
      for(k=0;k<n_cov;k++) mtemp1[j][k]=Sigma[i][j][k]/tau0;
    rMVN(mu[i], mu0, mtemp1, n_cov);
  }
  
  /* initialize the cluster membership */
  nstar=n_samp;  /* the # of disticnt values */
  for(i=0;i<n_samp;i++) 
    C[i]=i; /*cluster is from 0...n_samp-1 */
  
  for(main_loop=0; main_loop<*n_gen; main_loop++){
    
    /**update W, Wstar given mu, Sigma**/
    for (i=0;i<n_samp;i++){
      /*1 project BVN(mu_i, Sigma_i) on the inth tomo line */
      dtemp=0;
      for (j=0;j<n_grid[i];j++){
	if (*link==1){
	  vtemp[0]=log(W1g[i][j])-log(1-W1g[i][j]);
	  vtemp[1]=log(W2g[i][j])-log(1-W2g[i][j]);
	  prob_grid[j]=dMVN(vtemp, mu[i], InvSigma[i], 2, 1) -
	    log(W1g[i][j])-log(W2g[i][j])-log(1-W1g[i][j])-log(1-W2g[i][j]);
	}
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
	}
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
      if (*link==1) {
	Wstar[i][0]=log(W1g[i][j])-log(1-W1g[i][j]);
	Wstar[i][1]=log(W2g[i][j])-log(1-W2g[i][j]);
      }
      else if (*link==2) {
	Wstar[i][0]=qnorm(W1g[i][j],0 ,1, 1, 0);
	Wstar[i][1]=qnorm(W2g[i][j],0 ,1, 1, 0);
      }
      else if (*link==3) {
	Wstar[i][0]=-log(-log(W1g[i][j]));
	Wstar[i][1]=-log(-log(W2g[i][j]));
      }
    }
    
    /**updating mu, Sigma given Wstar**/
    for (i=0; i<n_samp; i++){
      /* generate weight vector q */
      dtemp=0;
      for (j=0; j<n_samp; j++){
	if (j!=i)
	  q[j]=dMVN(Wstar[i], mu[j], InvSigma[j], 2, 0);
	else 
	  q[j]=alpha*dMVT(Wstar[i], mu0, S_bvt, nu0-1, 2, 0);
	dtemp+=q[j];
	qq[j]=dtemp;	/*compute qq, the cumlative of q*/
      }
      /*standardize q and qq */
      for (j=0; j<n_samp; j++) {
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
	for (k=0; k<n_cov; k++)
	  for (l=0; l<n_cov; l++)
	    Sn[k][l]=S0[k][l]+tau0*(Wstar[i][k]-mu0[k])*(Wstar[i][l]-mu0[l])/(tau0+1); 
	dinv(Sn, n_cov, mtemp); 
	rWish(InvSigma[i], mtemp, nu0+1, n_cov);
	dinv(InvSigma[i], n_cov, Sigma[i]);
	for (k=0;k<n_cov;k++){
	  mun[k]=(tau0*mu0[k]+Wstar[i][k])/(tau0+1); 
	  for (l=0;l<n_cov;l++)  mtemp[k][l]=Sigma[i][k][l]/(tau0+1); 
	}
	rMVN(mu[i], mun, mtemp, n_cov);
	C[i]=nstar;
	nstar++;
      }
      else { 
	/*1. mu_i=mu_j, Sigma_i=Sigma_j*/
	/*2. update C[i]=C[j] */
	for(k=0;k<n_cov;k++) {       
	  mu[i][k]=mu[j][k];
	  for(l=0;l<n_cov;l++) {       
	    Sigma[i][k][l]=Sigma[j][k][l];
	    InvSigma[i][k][l]=InvSigma[j][k][l];
	  }
	}
	C[i]=C[j];
      }
      sortC[i]=C[i];
    } /* end of i loop*/
    
      /** remixing step **/
    for(i=0;i<n_samp;i++)
      indexC[i]=i;
    R_qsort_int_I(sortC, indexC, 1, n_samp);
    
    nstar=0;
    i=0;
    /*    printf("\n");
	  printf("C sortC indexC mu0 mu1 Sigmat\n");
	  for (i=0; i<n_samp; i++) {
	  printf("\n%5d%5d%5d", C[i], sortC[i], indexC[i]);
	  printf("%14g%14g%14g", mu[i][0], mu[i][0], Sigmat[i][0][0]);
	  }*/
    
    while (i<n_samp){
      /*initialize the vector and matrix */
      for(k=0; k<n_cov; k++) {
	Wstar_bar[k]=0;
	for(l=0;l<n_cov;l++) Snj[k][l]=S0[k][l];
      }
      
      j=sortC[i]; /*saves the first element in a block of same values */
      nj=0; /* counter for a block of same values */
      
      /* get data for remixing */
      while ((sortC[i]==j) && (i<n_samp)) {
	label[nj]=indexC[i];
	for (k=0; k<n_cov; k++) {
	  Wstarmix[nj][k]=Wstar[label[nj]][k];
	  Wstar_bar[k]+=Wstarmix[nj][k];
	}
	nj++;
	i++;
      } /* i records the current position in IndexC */
        /* nj records the # of obs in Psimix */
      
	/** posterior update for mu_mix, Sigma_mix based on Psimix **/
      for (k=0; k<n_cov; k++)
	Wstar_bar[k]/=nj;
      /* compute Snj */
      /*1. first two terms in Snj */
      /* Snj=S0 */
      for (j=0; j<nj; j++) 
	for (k=0; k<n_cov; k++)
	  for (l=0; l<n_cov; l++)
	    Snj[k][l]+=(Wstarmix[j][k]-Wstar_bar[k])*(Wstarmix[j][l]-Wstar_bar[l]);
      
      /*2. plus third term in Snj*/
      for (k=0;k<n_cov;k++)
	for (l=0;l<n_cov;l++) 
	  Snj[k][l]+=tau0*nj*(Wstar_bar[k]-mu0[k])*(Wstar_bar[l]-mu0[l])/(tau0+nj); 
      
      /*darw unscaled InvSigma_mix ~Wish(Snj^-1) */
      dinv(Snj, n_cov, mtemp);
      rWish(InvSigma_mix, mtemp, nu0+nj, n_cov); 
      dinv(InvSigma_mix, n_cov, Sigma_mix);
      
      /*2. draw mu_mix from N(munj, Sigma_mix) with Sigma_mix/=(tau0+nj) */ 
      for (k=0; k<n_cov; k++){
	munj[k]=(tau0*mu0[k]+nj*Wstar_bar[k])/(tau0+nj); 
	for (l=0; l<n_cov; l++)  mtemp[k][l]=Sigma_mix[k][l]/(tau0+nj);
      }
      rMVN(mu_mix, munj, mtemp, n_cov);
      
      /**update mu, Simgat with mu_mix, Sigmat_mix via label**/
      for (j=0;j<nj;j++){
	C[label[j]]=nstar;  /*updating C vector with no gap */
	for (k=0; k<n_cov; k++){
	  mu[label[j]][k]=mu_mix[k];
	  for (l=0;l<n_cov;l++){
	    Sigma[label[j]][k][l]=Sigma_mix[k][l];
	    InvSigma[label[j]][k][l]=InvSigma_mix[k][l];
	  }
	}
      }
      nstar++; /*finish update one distinct value*/
    } /* nstar is the number of distinct values */
    
      /** updating alpha **/
    if(*pinUpdate) {
      dtemp=b0-log(rbeta((double)alpha+1, (double)n_samp)); 
      if(unif_rand() < (double)(a0+nstar-1)/(n_samp*dtemp))
	alpha=rgamma((double)a0+nstar, 1/dtemp);
      else
	alpha=rgamma((double)a0+nstar-1, 1/dtemp);
      printf("%5d%14g%5d\n", main_loop, alpha, nstar); 
    }
    
    /*store Gibbs draws after burn_in */
    if (main_loop>=*burn_in) {
      itempC++;
      if (itempC==nth){
	if(*pinUpdate) { 
	  pdSa[itempA]=alpha;
	  pdSn[itempA]=nstar;
	  itempA++;
	}
	for(i=0; i<n_samp; i++) {
	  pdSMu0[itempS]=mu[i][0];
	  pdSMu1[itempS]=mu[i][1];
	  pdSSig00[itempS]=Sigma[i][0][0];
	  pdSSig01[itempS]=Sigma[i][0][1];
	  pdSSig11[itempS]=Sigma[i][1][1];
	  pdSW1[itempS]=W[i][0];
	  pdSW2[itempS]=W[i][1];
	  /* Wstar prediction */
	  if (*pred) {
	    rMVN(vtemp, mu[i], Sigma[i], n_cov);  
	    if (*link==1){
	      pdSWt1[itempS]=exp(vtemp[0])/(exp(vtemp[0])+1);
	      pdSWt2[itempS]=exp(vtemp[1])/(exp(vtemp[1])+1);
	    }
	    else if (*link==2){
	      pdSWt1[itempS]=pnorm(vtemp[0], 0, 1, 1, 0);
	      pdSWt2[itempS]=pnorm(vtemp[1], 0, 1, 1, 0);
	    }
	    else if (*link==3){
	      pdSWt1[itempS]=exp(-exp(-vtemp[0]));
	      pdSWt2[itempS]=exp(-exp(-vtemp[1]));
	    }
	  }
	  itempS++;
	}
	itempC=0;
      }
    }
  } /*end of MCMC for DP*/
  
  /** write out the random seed **/
  PutRNGstate();

  /* Freeing the memory */
  FreeMatrix(X, n_samp);
  FreeMatrix(W, n_samp);
  FreeMatrix(Wstar, n_samp);
  free(minW1);
  free(maxW1);
  free(n_grid);
  free(resid);
  FreeMatrix(S0, n_cov);
  FreeMatrix(Sn, n_cov);
  free(mun);
  FreeMatrix(W1g, n_samp);
  FreeMatrix(W2g, n_samp);
  free(prob_grid);
  free(prob_grid_cum);
  Free3DMatrix(Sigma, n_samp,n_cov);
  Free3DMatrix(InvSigma, n_samp, n_cov);
  FreeMatrix(mu, n_samp);
  free(C);
  free(q);
  free(qq);
  FreeMatrix(S_bvt, n_cov);
  free(sortC);
  free(indexC);
  FreeMatrix(Wstarmix, n_samp);
  FreeMatrix(Snj, n_cov);
  free(munj);
  free(Wstar_bar);
  free(mu_mix);
  FreeMatrix(Sigma_mix, n_cov);
  FreeMatrix(InvSigma_mix, n_cov);
  free(label);
  free(vtemp);
  FreeMatrix(mtemp, n_cov);
  FreeMatrix(mtemp1, n_cov);

  
} /* main */

