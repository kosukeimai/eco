#include <stdlib.h>
#include <assert.h>

int *intArray(int num);
int **intMatrix(int row, int col);

double *doubleArray(int num);
double **doubleMatrix(int row, int col);
double ***doubleMatrix3D(int x, int y, int z);

long *longArray(int num);

void FreeMatrix(double **Matrix, int row);
void FreeintMatrix(int **Matrix, int row);
void Free3DMatrix(double ***Matrix, int index, int row);
