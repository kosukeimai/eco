/******************************************************************
  This file is a part of eco: R Package for Fitting Bayesian Models
  of Ecological Inference for 2x2 Tables
  by Kosuke Imai and Ying Lu
  Copyright: GPL version 2 or later.
*******************************************************************/
void test(double *W1, int n, void *param);
void NormConstW1(double *W1, int n, void *param);
void NormConstW2(double *W1, int n, void *param);
void NormConstT(double *t, int n, void *param);
void SuffExp(double *t, int n, void *param);
void W1Exp(double *W1, int n, void *param);
void W1ExpW2(double *W2, int n, void *param);
void W2Exp(double *W2, int n, void *param);
void W2ExpW1(double *W1, int n, void *param);
void W1W1Exp(double *W1, int n, void *param);
void W1W1ExpW2(double *W2, int n, void *param);
void W2W2Exp(double *W2, int n, void *param);
void W2W2ExpW1(double *W1, int n, void *param);
void W1W2Exp(double *W1, int n, void *param);
void W2W1Exp(double *W2, int n, void *param);
void setNormConst(Param* param);
double getW2starFromW1star(double X, double Y, double W1, int* imposs);
double getW1starFromW2star(double X, double Y, double W2, int* imposs);
double getW1starFromT(double t, Param* param, int* imposs);
double getW2starFromT(double t, Param* param, int* imposs);
double getW1starPrimeFromT(double t, Param* param);
double getW2starPrimeFromT(double t, Param* param);
double numIntegration(integr_fn f, void *ex, int inf, double lb, double ub);
double paramIntegration(integr_fn f, void *ex);
//double numIntegration2(gsl_fn f, void *ex);
void setBounds(Param* param);

