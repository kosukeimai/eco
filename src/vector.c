#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

int* intArray(int num)
{
  int *iArray = (int *)malloc(num * sizeof(int));
  if (iArray)
    {
      return iArray;
    }
  else
    {
      fprintf(stderr, "Out of memory error in intArray\n");
      exit(-1);
    }
}

int** intMatrix(int row, int col)
{
  int i;
  int **iMatrix = (int **)malloc(row * sizeof(int *));
  if (iMatrix)
	{
	  for (i = 0; i < row; i++)
	    {
	      iMatrix[i] = (int *)malloc(col *  sizeof(int));
	      if (!iMatrix[i])
		{
		  fprintf(stderr, "Out of memory error in intMatrix\n");
		  exit(-1);
		}
	    }
	  return iMatrix;
	}
  fprintf(stderr, "Out of memory error in intMatrix\n");
  exit(-1);
}

double* doubleArray(int num)
{
  double *dArray = (double *)malloc(num * sizeof(double));
  if (dArray)
    {
      return dArray;
    }
  else
    {
      fprintf(stderr, "Out of memory error in doubleArray\n");
      exit(-1);
    }
}

double** doubleMatrix(int row, int col)
{
  int i;
  double **dMatrix = (double **)malloc((size_t)(row * sizeof(double *)));
  if (dMatrix)
    {
      for (i = 0; i < row; i++)
	{
	  dMatrix[i] = (double *)malloc((size_t)(col * sizeof(double)));
	  if (!dMatrix[i])
	    {
	      fprintf(stderr, "Out of memory error in doubleMatrix\n");
	      exit(-1);
	    }
	}
      return dMatrix;
    }
  fprintf(stderr, "Out of memory error in doubleMatrix\n");
  exit(-1);
}

double*** doubleMatrix3D(int x, int y, int z) {
  int i;
  double ***dM3 = (double ***)malloc(x * sizeof(double **));
  if (dM3) {
    for (i = 0; i < x; i++) {
      dM3[i] = doubleMatrix(y, z);
    }
    return dM3;
  }
  fprintf(stderr, "Out of memory error in doubleMatrix3D\n");
  exit(-1);
}

long* longArray(int num)
{
  long *lArray = (long *)malloc(num * sizeof(long));
  if (lArray)
    {
      return lArray;
    }
  else
    {
      fprintf(stderr, "Out of memory error in longArray\n");
      exit(-1);
    }
}

void FreeMatrix(double **Matrix, int row)
{
  int i;
  for (i = 0; i < row; i++)
    {
      free(Matrix[i]);
    }
  free(Matrix);
}

void FreeintMatrix(int **Matrix, int row)
{
  int i;
  for (i = 0; i < row; i++)
    {
      free(Matrix[i]);
    }
  free(Matrix);
}

void Free3DMatrix(double ***Matrix, int index, int row)
{
  int i;
  for (i = 0; i < index; i++)
    {
      FreeMatrix(Matrix[i], row);
    }
  free(Matrix);
}
		
	
		
			
