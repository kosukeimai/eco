/******************************************************************
  This file is a part of eco: R Package for Fitting Bayesian Models
  of Ecological Inference for 2x2 Tables
  by Kosuke Imai and Ying Lu
  Copyright: GPL version 2 or later.
*******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <R_ext/Lapack.h>
#include "vector.h"
#include "rand.h"
#include "subroutines.h"


/*
 * Computes the dot product of two vectors
 */
double dotProduct(double* a, double* b, int size) {
  int i; double ans=0;
  for (i=0; i<size; i++) {
    ans+=a[i]*b[i];
  }
  return ans;
}

/*
 * Multiply two matrices (A,B) with dims r1,c1,r2,c2
 * mutates C to return the answer
 */
void matrixMul(double** A, double** B, int r1, int c1, int r2, int c2, double** C) {
  int i,j,k;
  double tmp[r1][c2];
  if (c1!=r2) error("Matrix multiplication: %d != %d", c2, r1);
  else {
    for (i=0; i<r1; i++)
      for (j=0; j<c2; j++) {
        double entry=0;
        for(k=0;k<r2;k++) entry += A[i][k]*B[k][j];
        tmp[i][j]=entry;
      }
    for (i=0; i<r1; i++)
      for (j=0; j<c2; j++) {
        C[i][j]=tmp[i][j];
      }
  }
}

/*  The Sweep operator */
void SWP(
	 double **X,             /* The Matrix to work on */
	 int k,                  /* The row to sweep */
	 int size)               /* The dim. of X */
{
  int i,j;

  if (X[k][k] < 10e-20)
    error("SWP: singular matrix.\n");
  else
    X[k][k]=-1/X[k][k];
  for(i=0;i<size;i++)
    if(i!=k){
      X[i][k]=-X[i][k]*X[k][k];
      X[k][i]=X[i][k];
    }
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      if(i!=k && j!=k)
	X[i][j]=X[i][j]+X[i][k]*X[k][j]/X[k][k];

}


/* inverting a matrix */
void dinv(double **X,
	  int	size,
	  double **X_inv)
{
  int i,j, k, errorM;
  double *pdInv = doubleArray(size*size);

  for (i = 0, j = 0; j < size; j++)
    for (k = 0; k <= j; k++)
      pdInv[i++] = X[k][j];
  F77_CALL(dpptrf)("U", &size, pdInv, &errorM FCONE);
  if (!errorM) {
    F77_CALL(dpptri)("U", &size, pdInv, &errorM FCONE);
    if (errorM) {
      if (errorM>0) {
        Rprintf("The matrix being inverted is singular. Error code %d\n", errorM);
      } else {
        Rprintf("The matrix being inverted contained an illegal value. Error code %d.\n", errorM);
      }
      error("Exiting from dinv().\n");
    }
  }
  else {
    if (errorM>0) {
      Rprintf("The matrix being inverted was not positive definite. Error code %d\n", errorM);
    } else {
      Rprintf("The matrix being inverted contained an illegal value. Error code %d.\n", errorM);
    }
    error("Exiting from dinv().\n");
  }
  for (i = 0, j = 0; j < size; j++) {
    for (k = 0; k <= j; k++) {
      X_inv[j][k] = pdInv[i];
      X_inv[k][j] = pdInv[i++];
    }
  }

  Free(pdInv);
}

/* inverting a matrix, first tyring positive definite trick, and then symmetric
 * Uses special syntax since we don't know dimensions of array
 * Prevents memory errors for matrices created with double[][]
 */
void dinv2D(double* X,
	  int	size,
	  double* X_inv,char* emsg)
{
  int i,j, k, errorM, skip;
  double *pdInv = doubleArray(size*size);
  skip=0;

  for (i = 0, j = 0; j < size; j++)
    for (k = 0; k <= j; k++)
      //pdInv[i++] = X[k][j];
      pdInv[i++] = *(X+k*size+j);

//Rprintf("test: %5g %5g %d",pdInv[0],pdInv[(size == 3) ? 5 : 2],i);
  F77_CALL(dpptrf)("U", &size, pdInv, &errorM FCONE);
  if (!errorM) {
    F77_CALL(dpptri)("U", &size, pdInv, &errorM FCONE);
    if (errorM) {
      Rprintf(emsg);
    if (errorM>0) {
      Rprintf(": The matrix being inverted is singular. Error code %d\n", errorM);
    } else {
      Rprintf(": The matrix being inverted contained an illegal value. Error code %d.\n", errorM);
    }
      error("Exiting from dinv2D().\n");
    }
  }
  else {
    Rprintf(emsg);
    if (errorM>0) {
      /* The matrix is not positive definite.
       * This error does occur with proper data, when the likelihood curve is flat,
       * usually with the combination of NCAR and SEM.  At one point we tried
       * inverting the matrix via an alternative method that does not rely on
       * positive definiteness (see below), but that just led to further errors.
       * Instead, the program halts as gracefully as possible.
       */
      //Inverting the matrix anyway:
      //Rprintf(": Warning, the matrix being inverted was not positive definite on minor order %d.\n", errorM);
      //dinv2D_sym(X,size,X_inv,emsg);
      //skip=1;
      Rprintf(": Error, the matrix being inverted was not positive definite on minor order %d.\n", errorM);
      error("The program cannot continue; try a different model or including supplemental data.\n");
    } else {
      Rprintf(": The matrix being inverted contained an illegal value. Error code %d.\n", errorM);
      error("Exiting from dinv2D().\n");
    }
  }

  if (skip==0) {
  for (i = 0, j = 0; j < size; j++) {
    for (k = 0; k <= j; k++) {
      *(X_inv+size*j+k) = pdInv[i];
      *(X_inv+size*k+j) = pdInv[i++];
    }
  }
  }

  Free(pdInv);
}


/* inverting a matrix, assumes symmtretric, but not pos def
 * Uses special syntax since we don't know dimensions of array
 * Prevents memory errors for matrices created with double[][]
*/
void dinv2D_sym(double* X,
	  int	size,
	  double* X_inv,char* emsg)
{
  int i,j, k, errorM, size2;
  size2=size*size;
  double *pdInv = doubleArray(size2);
  double *B= doubleArray(size2);
  int *factor_out = intArray(size);

  //init pdInv and B.  B is identity
  for (i = 0, j = 0; j < size; j++)
    for (k = 0; k < size; k++) {
      if (j==k) B[i]=1;
      else B[i]=0;
      pdInv[i]=*(X+k*size+j);
      i++;
    }

  //for (i = 0, j = 0; j < size; j++)
  //  for (k = 0; k <= j; k++) {
  //    pdInv[i++] = *(X+k*size+j);
  //  }

  double *work0 = doubleArray(size2);
  int test=-1;
  F77_CALL(dsysv)("U", &size, &size, pdInv, &size, factor_out, B, &size, work0, &test, &errorM FCONE);
  int lwork=(int)work0[0];
  Free(work0);

  //Rprintf("work size %d\n",lwork);
  double *work = doubleArray(lwork);
  //Rprintf("In A: %5g %5g %5g %5g\n",pdInv[0],pdInv[1],pdInv[2],pdInv[3]);
  //Rprintf("In B: %5g %5g %5g %5g\n",B[0],B[1],B[2],B[3]);
  F77_CALL(dsysv)("U", &size, &size, pdInv, &size, factor_out, B, &size, work, &lwork, &errorM FCONE);
  Free(work);
  //Rprintf("Out1: %5g %5g %5g %5g %d\n",B[0],B[1],B[2],B[3],errorM);

  if (errorM) {
    Rprintf(emsg);
    if (errorM>0) {
      Rprintf(": The matrix being inverted is singular. Error code %d\n", errorM);
    } else {
      Rprintf(": The matrix being inverted contained an illegal value. Error code %d.\n", errorM);
    }
    error("Exiting from dinv2D_sym() (dsytrf).\n");
  }

  for (i = 0, j = 0; j < size; j++) {
    for (k = 0; k < size; k++) {
      *(X_inv+size*j+k) = B[i++];
    }
  }

  free(factor_out);
  Free(B);
  Free(pdInv);
}


/* Cholesky decomposition */
/* returns lower triangular matrix */
void dcholdc(double **X, int size, double **L)
{
  int i, j, k, errorM;
  double *pdTemp = doubleArray(size*size);

  for (j = 0, i = 0; j < size; j++)
    for (k = 0; k <= j; k++)
      pdTemp[i++] = X[k][j];
  F77_CALL(dpptrf)("U", &size, pdTemp, &errorM FCONE);
  if (errorM) {
    if (errorM>0) {
      Rprintf("The matrix being inverted was not positive definite. Error code %d\n", errorM);
    } else {
      Rprintf("The matrix being inverted contained an illegal value. Error code %d.\n", errorM);
    }
    error("Exiting from dcholdc().\n");
  }
  for (j = 0, i = 0; j < size; j++) {
    for (k = 0; k < size; k++) {
      if(j<k)
	L[j][k] = 0.0;
      else
	L[j][k] = pdTemp[i++];
    }
  }

  Free(pdTemp);
}

/* calculate the determinant of the positive definite symmetric matrix
   using the Cholesky decomposition  */
double ddet(double **X, int size, int give_log)
{
  int i;
  double logdet=0.0;
  double **pdTemp = doubleMatrix(size, size);

  dcholdc(X, size, pdTemp);
  for(i = 0; i < size; i++)
    logdet += log(pdTemp[i][i]);

  FreeMatrix(pdTemp, size);
  if(give_log)
    return(2.0*logdet);
  else
    return(exp(2.0*logdet));
}


/* calculate the determinant of the positive definite symmetric matrix
   using the Cholesky decomposition  use with double[][]*/
double ddet2D(double** X, int size, int give_log)
{
  int i;
  double logdet=0.0;
  double **pdTemp = doubleMatrix(size, size);

  dcholdc2D((double*)(&X[0][0]), size, (double*)(&pdTemp[0][0]));
  for(i = 0; i < size; i++)
    logdet += log(pdTemp[i][i]);

  FreeMatrix(pdTemp, size);
  if(give_log)
    return(2.0*logdet);
  else
    return(exp(2.0*logdet));
}

/*double ddet2Db(double* X, int size, int give_log)
{
  int i;
  double logdet=0.0;
  double **pdTemp = doubleMatrix(size, size);

  dcholdc2D(X, size, (double*)(&pdTemp[0][0]));
  for(i = 0; i < size; i++)
    logdet += log(pdTemp[i][i]);

  FreeMatrix(pdTemp, size);
  if(give_log)
    return(2.0*logdet);
  else
    return(exp(2.0*logdet));
}*/

/* Cholesky decomposition */
/* returns lower triangular matrix; use with double[][] */
void dcholdc2D(double *X, int size, double *L)
{
  int i, j, k, errorM;
  double *pdTemp = doubleArray(size*size);

  for (j = 0, i = 0; j < size; j++)
    for (k = 0; k <= j; k++)
      pdTemp[i++] = *(X+size*k+j); //pdTemp[i++] = X[k][j];
  F77_CALL(dpptrf)("U", &size, pdTemp, &errorM FCONE);
  if (errorM) {
    if (errorM>0) {
      Rprintf("The matrix being inverted was not positive definite. Error code %d\n", errorM);
    } else {
      Rprintf("The matrix being inverted contained an illegal value. Error code %d.\n", errorM);
    }
    error("Exiting from dcholdc2D().\n");
  }
  for (j = 0, i = 0; j < size; j++) {
    for (k = 0; k < size; k++) {
      if(j<k)
        *(L+size*j+k)=0.0; //L[j][k] = 0.0;
      else
        *(L+size*j+k)=pdTemp[i++]; //L[j][k] = pdTemp[i++];
    }
  }

  Free(pdTemp);
}

