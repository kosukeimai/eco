#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"

void cEMeco(
	      /*data input */
	      double *pdX,     /* data (X, Y) */
              double *pdTheta_old,  /* Theta^ t */

	      int *pin_samp,   /* sample size */
	      int *n_gen,      /* number of gibbs draws */
	      int *by_draw,    /* increase of draws at each iteration */
	      int *max_draw,   /* max # of draws */
	      double *converge, /*convergence criterion */
	      int *max_iter, /*maximum # of iteration before convergence*/
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
	      double *pdTheta,  /*EM result for Theta^(t+1) */
	      double *Ioc /*Ioc matrix given converged theta values */
	      ){	   
  
  int n_samp = *pin_samp;    /* sample size */

  int data=0;            /* one to print the data */
  int cflag;            /* one if satisfy the convergence criterion */
  int n_cov=2;           /* The number of covariates */
  int trapod=0;          /* one if use trapozodial approx */
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
  int n_step=4000;    /* 1/The default size of grid step */  
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


  /* misc variables */
  int i, j, k, l, main_loop;   /* used for various loops */
  int ndraw, itemp, itemp0;
  double dtemp, dtemp1, temp0, temp1;
  double w1, w2, w11, w12, w22;
  double d1, d2, d11, d12, d22, rho, rho2, sig1, sig2;
  double *vtemp, *utemp, *ttemp;
  double **mtemp;
  int *mflag;

  /* get random seed */
  GetRNGstate();

  ndraw=*n_gen;

  /* defining vectors and matricies */
  /* data */
  X=doubleMatrix(n_samp,n_cov);
  W=doubleMatrix((n_samp+s_samp+x0_samp+x1_samp),n_cov);
  Wstar=doubleMatrix((n_samp+s_samp+x0_samp+x1_samp),5);

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

  vtemp=doubleArray(n_cov);
  utemp=doubleArray(ndraw);
  ttemp=doubleArray(5);
  mtemp=doubleMatrix(n_cov,n_cov);

  t_samp=n_samp+s_samp+x1_samp+x0_samp;  

  mflag=intArray(n_step);
  for (i=0; i<n_step; i++) {
    mflag[i]=0;
  }

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
  mu_ord[0] = pdTheta_old[0];
  mu_ord[1] = pdTheta_old[1];
  Sigma_ord[0][0] = pdTheta_old[2]*pdTheta_old[2];
  Sigma_ord[1][1] = pdTheta_old[3]*pdTheta_old[3];
  Sigma_ord[0][1] = pdTheta_old[4]*pdTheta_old[2]*pdTheta_old[3];
  Sigma_ord[1][0] = Sigma_ord[0][1];

  dinv(Sigma_ord, n_cov, InvSigma_ord);


  /* initialize W, Wstar for t_samp*/
  for (i=0; i< t_samp; i++) {
    W[i][0]=0;
    W[i][1]=0;
    Wstar[i][0]=0;
    Wstar[i][1]=0;
    Wstar[i][2]=0;
    Wstar[i][3]=0;
    Wstar[i][4]=0;
  }

  for (i=0; i<n_samp; i++)
    {
      if (X[i][1]==0) 
	{
	  W[i][0]=0.000001;
	  W[i][1]=0.000001;
	  temp0=log(W[i][0])-log(1-W[i][0]);
	  temp1=log(W[i][1])-log(1-W[i][1]);
	  Wstar[i][0]=temp0;
	  Wstar[i][1]=temp1;
	  Wstar[i][2]=temp0*temp0;
	  Wstar[i][3]=temp0*temp1;
	  Wstar[i][4]=temp1*temp1;
	}
      if (X[i][1]==1) 
	{
	  W[i][0]=0.999999; 
	  W[i][1]=0.999999;
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
    for (i=0; i<s_samp; i++) {
	Wstar[(n_samp+x1_samp+x0_samp+i)][2]=S_Wstar[i][0]*S_Wstar[i][0];
	Wstar[(n_samp+x1_samp+x0_samp+i)][3]=S_Wstar[i][0]*S_Wstar[i][1];
	Wstar[(n_samp+x1_samp+x0_samp+i)][4]=S_Wstar[i][1]*S_Wstar[i][1];
    }

    if(data==1) { 
      printf("survey W1 W2 W1* W2*\n");
      for(i=0;i<(n_samp+s_samp+x1_samp+x0_samp);i++)
	printf("%5d%14g%14g%14g%14g\n",i,W[i][0],W[i][1],Wstar[i][0], Wstar[i][1]);
      fflush(stdout);
    }
  }


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
    


     /* draw n_gen's from uniform 
	for (i=0; i<ndraw; i++) {
	utemp[i]=unif_rand();
	}
	R_rsort(utemp, ndraw); */

  cflag=0;
  main_loop=0;

  /* EM interation until cflag==1 or exceed max_iter*/
  while ((cflag==0) && (main_loop<*max_iter)) 
    {
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
	  for (j=0;j<n_grid[i];j++){
	    prob_grid_cum[j]/=dtemp; /*standardize prob.grid */ 
	  }
	  /*2 sample W_i on the ith tomo line */
	  /*3 compute Wsta_i from W_i*/
	  
	  j=0;
	  itemp=1;
	  
	  for (k=0; k<ndraw; k++){
	    /* dtemp=utemp[k]; */
	    j=findInterval(prob_grid_cum, n_grid[i],
			   (double)(1+k)/(ndraw+1), 1, 1, itemp, mflag); 
	    itemp=j;
	    if (j==0 || trapod==0) {
	      W[i][0]=W1g[i][j];
	      W[i][1]=W2g[i][j];
	    }
	    else if (j>=1 && trapod==1) {
	      dtemp1=((double)(1+k)/(ndraw+1)-prob_grid_cum[(j-1)])/(prob_grid_cum[j]-prob_grid_cum[j-1]);
	      W[i][0]=dtemp1*(W1g[i][j]-W1g[i][(j-1)])+W1g[i][(j-1)];
	      W[i][1]=dtemp1*(W2g[i][j]-W2g[i][(j-1)])+W2g[i][(j-1)];
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
      for (i=0; i<n_samp; i++)
	{
	  if ( X[i][1]!=0 && X[i][1]!=1 ) {  
	    Wstar[i][0]/=ndraw;  /*E(W1i) */
	    Wstar[i][1]/=ndraw;  /*E(W2i) */
	    Wstar[i][2]/=ndraw;  /*E(W1i^2) */
	    Wstar[i][3]/=ndraw;  /*E(W1iW2i) */
	    Wstar[i][4]/=ndraw;  /*E(W2i^2) */
	  }
	} /*for x0type, x1type and survey data, E-step is either the observed value or the analytical expectation*/
      
      
      /*analytically compute E{W2_i|Y_i} given W1_i, mu_ord and Sigma_ord in x1 homeogeneous areas */
      if (*x1==1)
	for (i=0; i<x1_samp; i++) {
	  dtemp=mu_ord[1]+Sigma_ord[0][1]/Sigma_ord[0][0]*(Wstar[n_samp+i][0]-mu_ord[0]);
	  temp0=log(W[n_samp+i][0])-log(1-W[n_samp+i][0]);
	  temp1=dtemp;
	  Wstar[n_samp+i][0]=temp0;
	  Wstar[n_samp+i][1]=temp1;
	  Wstar[n_samp+i][2]=temp0*temp0;
	  Wstar[n_samp+i][3]=temp0*temp1;
	  Wstar[n_samp+i][4]=temp1*temp1;
	}
      
      /*analytically compute E{W1_i|Y_i} given W2_i, mu_ord and Sigma_ord in x0 homeogeneous areas */
      
      if (*x0==1)
	for (i=0; i<x0_samp; i++) {
	  dtemp=mu_ord[0]+Sigma_ord[0][1]/Sigma_ord[1][1]*(Wstar[n_samp+x1_samp+i][1]-mu_ord[1]);
	  temp0=dtemp;
	  temp1=log(W[n_samp+x1_samp+i][1])-log(1-W[n_samp+x1_samp+i][1]);
	  Wstar[n_samp+x1_samp+i][0]=temp0;
	  Wstar[n_samp+x1_samp+i][1]=temp1;
	  Wstar[n_samp+x1_samp+i][2]=temp0*temp0;
	  Wstar[n_samp+x1_samp+i][3]=temp0*temp1;
	  Wstar[n_samp+x1_samp+i][4]=temp1*temp1;
	}
      
      
      
      
      
      /*M-step: same procedure as normal model */
      for (j=0; j<5; j++) 
	{
	  ttemp[j]=0;
	}
      w1=0;
      w2=0;
      w11=0;
      w12=0;
      w22=0;
      
      for (i=0; i<t_samp; i++)
	{
	  w1+=Wstar[i][0];   /*sum E(W_1i) */
	  w2+=Wstar[i][1];   /*sum E(W_2i) */
	  w11+=Wstar[i][0]*Wstar[i][0];   /*sum E(W_1i*W_1i) */
	  w12+=Wstar[i][0]*Wstar[i][1];   /*sum E(W_1i*W_2i) */
	  w22+=Wstar[i][1]*Wstar[i][1];   /*sum E(W_2i*W_2i) */
	}
      
      ttemp[0]=(double)w1/t_samp;  /*mu1*/
      ttemp[1]=(double)w2/t_samp;  /*mu2*/
      ttemp[2]=(double)sqrt((w11-2*ttemp[0]*w1+t_samp*ttemp[0]*ttemp[0])/t_samp); /*sigma11^0.5 */ 
      ttemp[3]=(double)sqrt((w22-2*ttemp[1]*w2+t_samp*ttemp[1]*ttemp[1])/t_samp); /*sigma22^0.5 */
      ttemp[4]=(double)((w12-ttemp[1]*w1-ttemp[0]*w2+t_samp*ttemp[0]*ttemp[1])/t_samp)/(ttemp[2]*ttemp[3]); /*pho */
   
      
      cflag=1;
      printf("\n%5d", main_loop);
      for (j=0; j<5; j++) {
	dtemp=ttemp[j]-pdTheta_old[j];
	if ((dtemp > *converge) || (-dtemp > *converge)) { cflag=0; }
	pdTheta_old[j]=ttemp[j];
	printf("%15g", ttemp[j]);
      }

      mu_ord[0] = pdTheta_old[0];
      mu_ord[1] = pdTheta_old[1];
      Sigma_ord[0][0] = pdTheta_old[2]*pdTheta_old[2];
      Sigma_ord[1][1] = pdTheta_old[3]*pdTheta_old[3];
      Sigma_ord[0][1] = pdTheta_old[4]*pdTheta_old[2]*pdTheta_old[3];
      Sigma_ord[1][0] = Sigma_ord[0][1];
      
      dinv(Sigma_ord, n_cov, InvSigma_ord);

      if (ndraw<(*max_draw-*by_draw)) ndraw+=*by_draw;
  
      main_loop++;
    }


  for (j=0; j<5; j++) {
    pdTheta[j]=ttemp[j];
  }

  d1=0;
  d2=0;
  d11=0;
  d12=0;
  d22=0;
  
  if (cflag==1) {
    printf("converged!\n"); 
    d1=(w1-t_samp*pdTheta[0])/pdTheta[2]; 
    d2=(w2-t_samp*pdTheta[1])/pdTheta[3];
    d11=(w11-2*pdTheta[0]*w1+t_samp*pdTheta[0]*pdTheta[0])/pdTheta[2]*pdTheta[2];
    d22=(w22-2*pdTheta[1]*w2+t_samp*pdTheta[1]*pdTheta[1])/pdTheta[3]*pdTheta[3];
    d12=(w12-pdTheta[0]*w2-pdTheta[1]*w1+t_samp*pdTheta[0]*pdTheta[1])/(pdTheta[2]*pdTheta[3]);

    sig1=pdTheta[2];
    sig2=pdTheta[3];
    rho=pdTheta[4];
    rho2=1-pdTheta[4]*pdTheta[4];

    Ioc[0]=(double)t_samp/(rho2*sig1*sig1);
    Ioc[1]=(double)-t_samp*rho/(rho2*sig1*sig2);
    Ioc[2]=(double)1/(rho2*sig1*sig1)*(2*d1-rho*d2);
    Ioc[3]=(double)-rho/rho2/(sig1*sig2)*d2;
    Ioc[4]=(double)-1/(rho2*rho2*sig1)*(2*rho*d1-(1+rho*rho)*d2);
    Ioc[5]=(double)t_samp/(rho2*sig2*sig2);
    Ioc[6]=(double)-rho/rho2/(sig1*sig2)*d1;
    Ioc[7]=(double)1/(rho2*sig2*sig2)*(2*d2-rho*d1);
    Ioc[8]=(double)-1/(rho2*rho2*sig2)*(2*rho*d2-(1+rho*rho)*d1);
    Ioc[9]=(double)-t_samp/(sig1*sig1)+1/(rho2*sig1*sig1)*(3*d11-2*rho*d12);
    Ioc[10]=(double)-rho/(rho2*sig1*sig2)*d12;
    Ioc[11]=(double)-1/(rho2*rho2*sig1)*(2*rho*d11-(1+rho*rho)*d12);
    Ioc[12]=(double)-t_samp/(sig2*sig2)+1/(rho2*sig2*sig2)*(3*d22-2*rho*d12);
    Ioc[13]=(double)-1/(rho2*rho2*sig2)*(2*rho*d22-(1+rho*rho)*d12);
    Ioc[14]=(double)-t_samp*(1+rho*rho)/(rho2*rho2)+1/(rho2*rho2*rho2)*((3*rho*rho+1)*(d11+d22)-2*rho*(3+rho*rho)*d12);
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
  free(vtemp);
  FreeMatrix(mtemp, n_cov);
  
} /* main */

