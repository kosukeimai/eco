# ifndef MACROS_H
# define MACROS_H


/****************/
/** structrues **/
/****************/
/* parameters and observed data */
struct Param_old{
  double mu[2];
  double Sigma[2][2];
  double InvSigma[2][2];
  double Sigma3[3][3];
  double InvSigma3[3][3];
  int NCAR;
  double data[2]; //collect the data
  double X; //X,Y here for ease of use
  double Y;
  double normcT; //normalized const on tomog line (integrating with parameterization)
  double W[2]; //if W is known, also handy place to store E[W1] when we calculate it each step
  double Wstar[2]; //place to store E[W1*] when we calculate it each step
  double W1_lb; //lower and upper bounds for W1 and W2 (not starred)
  double W1_ub;
  double W2_lb;
  double W2_ub;
  int W1_inf; //inf: 0->(lb,ub), -1->(-inf,ub), 1->(lb,inf), 2->(-inf,inf)
  int W2_inf;
  int suff; //the sufficient stat we're calculating: 0->W1, 1->W2,2->W1^2,3->W1W2,4->W2^2,7->Log Lik, 5/6,-1 ->test case
};

typedef struct Param_old Param_old;

struct caseParam {
  double mu[2];
  double data[2]; //collect the data
  double X; //X,Y here for ease of use
  double Y;
  double normcT; //normalized const on tomog line (integrating with parameterization)
  double W[2]; //if W is known, also handy place to store E[W1] when we calculate it each step
  double Wstar[2]; //place to store E[W1*] when we calculate it each step
  double Wbounds[2][2];  //[i][j] is {j:lower,upper}-bound of W{i+1}
  int suff; //the sufficient stat we're calculating: 0->W1, 1->W2,2->W1^2,3->W1W2,4->W2^2,7->Log Lik, 5/6,-1 ->test case
  int dataType; //0=unknown, 1=(X==1),2=(X==0),3=survey
  double** Z_i; //CCAR: k x 2
};

typedef struct caseParam caseParam;

struct setParam {
  int n_samp, t_samp, s_samp,x1_samp,x0_samp,param_len,suffstat_len; //types of data sizes
  int iter, ncar, ccar, ccar_nvar, fixedRho, sem, hypTest, verbose, calcLoglik; //options
  int semDone[7]; //whether that row of the R matrix is done
  int varParam[9]; //whether the parameter is included in the R matrix
  double convergence;
  double Sigma[2][2];
  double InvSigma[2][2];
  double Sigma3[3][3];
  double InvSigma3[3][3];
  double** SigmaK; //for CCAR
  double** InvSigmaK;
  double** hypTestCoeff;
  double hypTestResult;
  double* pdTheta;
};

typedef struct setParam setParam;

struct Param {
  setParam* setP; //pointer to the singleton structure
  caseParam caseP;
};

typedef struct Param Param;

/***************************/
/** typedef functions     **/
/***************************/

//typedef void integr_fn(double *x, int n, void *ex); //is already defined in Applic.h
typedef double gsl_fn(double x, void *ex);

# endif
