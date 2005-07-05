/******************************************************************
  This file is a part of eco: R Package for Fitting Bayesian Models 
  of Ecological Inference for 2x2 Tables
  by Kosuke Imai and Ying Lu
  Copyright: GPL version 2 or later.
*******************************************************************/

void NIWupdate(double **Y, double *mu, double **Sigma, double **InvSigma,
	       double *mu0, double tau0, int nu0, double **S0, 
	       int n_samp, int n_dim); 
