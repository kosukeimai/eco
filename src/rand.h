double dMVN(double *Y, double *MEAN, double **SIG_INV, int dim, int give_log);
double dMVT(double *Y, double *MEAN, double **SIG_INV, int nu, int dim, int give_log);
double TruncNorm(double lb, double ub, double mu, double var);
void rMVN(double *Sample, double *mean, double **inv_Var, int size);
void rwish(double **Sample, double **S, int df, int size);
void rWish(double **Sample, double **S, int df, int size);


