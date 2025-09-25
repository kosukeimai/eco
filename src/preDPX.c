#include <string.h>
#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "bayes.h"
#include "sample.h"

/* Conditional Prediction for Nonparametric Model for 2x2 Tables */
void preDPX(
	   double *pdmu, 
	   double *pdSigma,
	   double *X,
	   int *pin_samp,
	   int *pin_draw,
	   int *pin_dim,
	   int *verbose,    /* 1 for output monitoring */
	   double *pdStore
	   ){	   
  
  /* some integers */
  int n_samp = *pin_samp;    /* sample size */
  int n_draw = *pin_draw;    /* sample size of survey data */ 
  int n_dim = *pin_dim;      /* dimension */

  double *mu = doubleArray(n_dim);                /* The mean */
  double *Wstar = doubleArray(n_dim);
  double **Sigma = doubleMatrix(n_dim, n_dim);    /* The covariance matrix */

  /* misc variables */
  int i, j, main_loop;   /* used for various loops */
  int itemp = 0;
  int itempM = 0;
  int itempS = 0;
  int progress = 1, itempP = ftrunc((double) n_draw/10);

  /* get random seed */
  GetRNGstate();
  
  for(main_loop=0; main_loop<n_draw; main_loop++){
    for(i=0; i<n_samp; i++) {
      mu[0] = pdmu[itempM]+pdSigma[itempS+2]/pdSigma[itempS+5]*(X[i]-pdmu[itempM+2]);
      mu[1] = pdmu[itempM+1]+pdSigma[itempS+4]/pdSigma[itempS+5]*(X[i]-pdmu[itempM+2]);
      Sigma[0][0] = pdSigma[itempS]-pdSigma[itempS+2]*pdSigma[itempS+2]/pdSigma[itempS+5];
      Sigma[1][1] = pdSigma[itempS+3]-pdSigma[itempS+4]*pdSigma[itempS+4]/pdSigma[itempS+5];
      Sigma[0][1] = pdSigma[itempS+1]-pdSigma[itempS+2]*pdSigma[itempS+4]/pdSigma[itempS+5];
      Sigma[1][0] = Sigma[0][1];
      rMVN(Wstar, mu, Sigma, n_dim);
      for (j=0; j<n_dim; j++)
	pdStore[itemp++] = exp(Wstar[j])/(1+exp(Wstar[j]));
      itempS += 6;
      itempM += 3;
    }
    if (*verbose)
      if (itempP == main_loop) {
        Rprintf("%3d percent done.\n", progress*10);
        itempP+=ftrunc((double) n_draw/10); progress++;
        R_FlushConsole();
      }
    R_CheckUserInterrupt();
  }
  
  if(*verbose)
    Rprintf("100 percent done.\n");

  /** write out the random seed **/
  PutRNGstate();

  /* Freeing the memory */
  free(mu);
  free(Wstar);
  FreeMatrix(Sigma,n_dim);
  
} /* main */

