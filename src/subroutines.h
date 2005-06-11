/******************************************************************
  This file is a part of eco: R Package for Estimating Fitting 
  Bayesian Models of Ecological Inference for 2X2 tables
  by Ying Lu and Kosuke Imai
  Copyright: GPL version 2 or later.
*******************************************************************/

void SWP( double **X, int k, int size);
void dinv(double **X, int size, double **X_inv);
void dcholdc(double **X, int size, double **L);
double ddet(double **X, int size, int give_log);
