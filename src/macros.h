# ifndef MACROS_H
# define MACROS_H


/****************/
/** structrues **/
/****************/
/* parameters and observed data */
struct Param{
  double mu[2];
  double Sigma[2][2];
  double InvSigma[2][2];
  double data[2]; //collect the data
  double X; //X,Y here for ease of use
  double Y;
  double normcT; //normalized const on tomog line (integrating with parameterization)
  double W[2]; //if W is known, also handy place to store E[W1] when we calculate it each step
  double W1_lb; //lower and upper bounds for W1 and W2 (not starred)
  double W1_ub;
  double W2_lb;
  double W2_ub;
  int W1_inf; //inf: 0->(lb,ub), -1->(-inf,ub), 1->(lb,inf), 2->(-inf,inf)
  int W2_inf;
  int suff; //the sufficient stat we're calculating: 0->W1, 1->W2,2->W1^2,3->W1W2,4->W2^2,7->Log Lik, 5/6,-1 ->test case
};

typedef struct Param Param;

/***************************/
/** typedef functions     **/
/***************************/

//typedef void integr_fn(double *x, int n, void *ex); //is already defined in Applic.h
typedef double gsl_fn(double x, void *ex);

# endif
