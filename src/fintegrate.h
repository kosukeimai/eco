/******************************************************************
  This file is a part of eco: R Package for Fitting Bayesian Models
  of Ecological Inference for 2x2 Tables
  by Kosuke Imai and Ying Lu
  Copyright: GPL version 2 or later.
*******************************************************************/
void test(double *W1, int n, void *param);
void NormConst(double *W1, int n, void *param);
void W1Exp(double *W1, int n, void *param);
void W1W1Exp(double *W1, int n, void *param);
void W2W2Exp(double *W2, int n, void *param);
double getNormConst(void* pp);
