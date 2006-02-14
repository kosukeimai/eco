/******************************************************************
  This file is a part of eco: R Package for Fitting Bayesian Models
  of Ecological Inference for 2x2 Tables
  by Kosuke Imai and Ying Lu
  Copyright: GPL version 2 or later.
*******************************************************************/
#include <R_ext/Applic.h>

void NormConstT(double *t, int n, void *param);
void SuffExp(double *t, int n, void *param);
double getLogLikelihood(Param* param) ;
void setNormConst(Param* param);
double getW2starFromW1star(double X, double Y, double W1, int* imposs);
double getW1starFromW2star(double X, double Y, double W2, int* imposs);
double getW1FromW2(double X, double Y, double W2);
double getW1starFromT(double t, Param* param, int* imposs);
double getW2starFromT(double t, Param* param, int* imposs);
double getW1starPrimeFromT(double t, Param* param);
double getW2starPrimeFromT(double t, Param* param);
double paramIntegration(integr_fn f, void *ex);
void setNormConst(Param* param);
void setBounds(Param* param);

