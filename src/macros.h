# ifndef MACROS_H
# define MACROS_H


/****************/
/** structrues **/
/****************/
/* parameters and observed data */
struct Param{
  double mu[2];
  double Sigma[2][2];
  double X;
  double Y;
  double normc; //normalized const on tomog line
};

typedef struct Param Param;

/***************************/
/** typedef functions     **/
/***************************/

typedef void integr_fn(double *x, int n, void *ex);

# endif
