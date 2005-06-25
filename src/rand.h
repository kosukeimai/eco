/******************************************************************
  This file is a part of eco: R Package for Estimating Fitting 
  Bayesian Models of Ecological Inference for 2X2 tables
  by Ying Lu and Kosuke Imai
  Copyright: GPL version 2 or later.
*******************************************************************/

double dMVN(double *Y, double *MEAN, double **SIGMA, int dim, int give_log);
double dMVT(double *Y, double *MEAN, double **SIG_INV, int nu, int dim, int give_log);
void rMVN(double *Sample, double *mean, double **inv_Var, int size);
void rWish(double **Sample, double **S, int df, int size);
void rDirich(double *Sample, double *theta, int size);

