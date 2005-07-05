/******************************************************************
  This file is a part of eco: R Package for Fitting Bayesian Models 
  of Ecological Inference for 2x2 Tables
  by Kosuke Imai and Ying Lu
  Copyright: GPL version 2 or later.
*******************************************************************/

void SWP( double **X, int k, int size);
void dinv(double **X, int size, double **X_inv);
void dcholdc(double **X, int size, double **L);
double ddet(double **X, int size, int give_log);
