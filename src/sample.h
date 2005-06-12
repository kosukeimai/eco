/******************************************************************
  This file is a part of eco: R Package for Estimating Fitting 
  Bayesian Models of Ecological Inference for 2X2 tables
  by Ying Lu and Kosuke Imai
  Copyright: GPL version 2 or later.
*******************************************************************/

void rGrid(double *Sample, double *W1gi, double *W2gi, int ni_grid, 
	   double *mu, double **InvSigma, int n_dim); 
void rMH(double *Sample, double *W, double *XY, double W1min, 
	 double W1max, double *mu, double **InvSigma, int n_dim);
void rMHrc(double *Sample, double *W, double *XY, double *Zmin, 
	   double *Zmax, double *mu, double **InvSigma, int n_dim);
