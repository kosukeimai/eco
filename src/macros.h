# ifndef MACROS_H
# define MACROS_H


/****************/
/** structrues **/
/****************/

/* ENUMS
 * sufficient statistic to calculate: 0->W1*, 1->W2*, 2->(W1*)^2, 3->(W1*)(W2*), 4->(W2*)^2, 5->W1,6->W2 7->Log Lik,8->test
 * data point type: 0=general, 1= homogenous with (X==1), 2= homogenous with (X==0), 3=survey (W1 and W2 are known)
 */
 enum e_sufficient_stats {SS_W1star, SS_W2star, SS_W1star2, SS_W1W2star, SS_W2star2, SS_W1, SS_W2, SS_Loglik, SS_Test};
 typedef enum e_sufficient_stats sufficient_stat;
 enum e_datapoint_types {DPT_General,DPT_Homog_X1, DPT_Homog_X0, DPT_Survey};
 typedef enum e_datapoint_types datapoint_type;

/* parameters and observed data  -- no longer used*/
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
  sufficient_stat suff; //the sufficient stat we're calculating: 0->W1, 1->W2,2->W1^2,3->W1W2,4->W2^2,7->Log Lik, 5/6,-1 ->test case
};

typedef struct Param_old Param_old;

/**
 * The structure that holds per-record infromation
 */
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
  datapoint_type dataType;
  double** Z_i; //CCAR: k x 2
};

typedef struct caseParam caseParam;

/**
 * The structure that holds dataset infromation
 */
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
