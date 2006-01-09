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
  double normcW1; //normalized const on tomog line (integrating over W1)
  double normcW2; //normalized const on tomog line (integrating over W2)
  double W1_lb; //lower and upper bounds for W1star and W2star
  double W1_ub;
  double W2_lb;
  double W2_ub;
  int W1_inf; //inf: 0->(lb,ub), -1->(-inf,ub), 1->(lb,inf), 2->(-inf,inf)
  int W2_inf;
};

typedef struct Param Param;

/***************************/
/** typedef functions     **/
/***************************/

typedef void integr_fn(double *x, int n, void *ex);
typedef double gsl_fn(double x, void *ex);

# endif
