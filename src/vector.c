/******************************************************************
  This file is a part of eco: R Package for Fitting Bayesian Models
  of Ecological Inference for 2x2 Tables
  by Kosuke Imai and Ying Lu
  Copyright: GPL version 2 or later.
*******************************************************************/

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <R_ext/Utils.h>
#include <R_ext/PrtUtil.h>
#include <R.h>

int* intArray(int num) {
  int *iArray = (int *)malloc(num * sizeof(int));
  if (iArray)
    return iArray;
  else {
    error("Out of memory error in intArray\n");
    return NULL;
  }
}

int** intMatrix(int row, int col) {
  int i;
  int **iMatrix = (int **)malloc(row * sizeof(int *));
  if (iMatrix) {
    for (i = 0; i < row; i++) {
      iMatrix[i] = (int *)malloc(col *  sizeof(int));
      if (!iMatrix[i])
	error("Out of memory error in intMatrix\n");
    }
    return iMatrix;
  }
  else {
    error("Out of memory error in intMatrix\n");
    return NULL;
  }
}

double* doubleArray(int num) {
  double *dArray = (double *)malloc(num * sizeof(double));
  //double *dArray = (double *)Calloc(num,double);
  if (dArray)
    return dArray;
  else {
    error("Out of memory error in doubleArray\n");
    return NULL;
  }
}

double** doubleMatrix(int row, int col) {
  int i;
  double **dMatrix = (double **)malloc((size_t)(row * sizeof(double *)));
  //double **dMatrix = Calloc(row,double*);
  if (dMatrix) {
    for (i = 0; i < row; i++) {
      /* dMatrix[i] = Calloc(col,double); */
      dMatrix[i] = (double *)malloc(col * sizeof(double));
      if (!dMatrix[i]) {
        error("Out of memory error in doubleMatrix\n");
        return NULL;
      }
    }
    return dMatrix;
  }
  else {
    error("Out of memory error in doubleMatrix\n");
    return NULL;
  }
}

double*** doubleMatrix3D(int x, int y, int z) {
  int i;
  double ***dM3 = (double ***)malloc(x * sizeof(double **));
  if (dM3) {
    for (i = 0; i < x; i++)
      dM3[i] = doubleMatrix(y, z);
    return dM3;
  }
  else {
    error("Out of memory error in doubleMatrix3D\n");
    return NULL;
  }
}

long* longArray(int num) {
  long *lArray = (long *)malloc(num * sizeof(long));
  if (lArray)
    return lArray;
  else {
    error("Out of memory error in longArray\n");
    return NULL;
  }
}

void FreeMatrix(double **Matrix, int row) {
  int i;
  for (i = 0; i < row; i++)
    free(Matrix[i]);
  free(Matrix);
}

void FreeintMatrix(int **Matrix, int row) {
  int i;
  for (i = 0; i < row; i++)
    free(Matrix[i]);
  free(Matrix);
}

void Free3DMatrix(double ***Matrix, int index, int row) {
  int i;
  for (i = 0; i < index; i++)
    FreeMatrix(Matrix[i], row);
  free(Matrix);
}




