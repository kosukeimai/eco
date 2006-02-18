/******************************************************************
  This file is a part of eco: R Package for Fitting Bayesian Models
  of Ecological Inference for 2x2 Tables
  by Kosuke Imai and Ying Lu
  Copyright: GPL version 2 or later.
*******************************************************************/

double dMVN(double *Y, double *MEAN, double **SIG_INV, int dim, int give_log);
double dMVT(double *Y, double *MEAN, double **SIG_INV, int nu, int dim, int give_log);
void rMVN(double *Sample, double *mean, double **inv_Var, int size);
void rWish(double **Sample, double **S, int df, int size);
void rDirich(double *Sample, double *theta, int size);
double dBVNtomo(double *Wstar, void* pp, int give_log, double normc);
double invLogit(double x);
double logit(double x,char* emsg);
int bit(int t, int n);


