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
  double normcT; //normalized const on tomog line (integrating with parameterization)
  double W1_lb; //lower and upper bounds for W1 and W2 (no starred)
  double W1_ub;
  double W2_lb;
  double W2_ub;
  int W1_inf; //inf: 0->(lb,ub), -1->(-inf,ub), 1->(lb,inf), 2->(-inf,inf)
  int W2_inf;
  int suff; //the sufficient stat we're calculating: 1->W1, 2->W2,3->W1^2,4->W1W2,5->W2^2, -1 ->test case
};

typedef struct Param Param;

/***************************/
/** typedef functions     **/
/***************************/

typedef void integr_fn(double *x, int n, void *ex);
typedef double gsl_fn(double x, void *ex);

# endif
