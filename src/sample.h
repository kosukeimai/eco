/******************************************************************
  This file is a part of eco: R Package for Fitting Bayesian Models 
  of Ecological Inference for 2x2 Tables
  by Kosuke Imai and Ying Lu
  Copyright: GPL version 2 or later.
*******************************************************************/

void rGrid(double *Sample, double *W1gi, double *W2gi, int ni_grid, 
	   double *mu, double **InvSigma, int n_dim); 
void GridPrep(double **W1g, double **W2g, double **X, double *maxW1,
	      double *minW1, int *n_grid, int n_samp, int n_step);
void rMH(double *W, double *XY, double W1min, double W1max, 
	 double *mu, double **InvSigma, int n_dim);
void rMH2c(double *W, double *X, double Y, double *minU, 
	   double *maxU, double *mu, double **InvSigma, int n_dim, 
	   int maxit, int reject);
