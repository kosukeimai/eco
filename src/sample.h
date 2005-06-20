/******************************************************************
  This file is a part of eco: R Package for Estimating Fitting 
  Bayesian Models of Ecological Inference for 2X2 tables
  by Ying Lu and Kosuke Imai
  Copyright: GPL version 2 or later.
*******************************************************************/

void rGrid(double *Sample, double *W1gi, double *W2gi, int ni_grid, 
	   double *mu, double **InvSigma, int n_dim); 
void GridPrep(double **W1g, double **W2g, double **X, double *maxW1,
	      double *minW1, int *n_grid, int n_samp, int n_step);
void rMH(double *W, double *XY, double W1min, double W1max, 
	 double *mu, double **InvSigma, int n_dim);
void rMHrc(double *W, double *X, double Y, double *minZ, 
	   double *maxZ, double *mu, double **InvSigma, int n_dim);
