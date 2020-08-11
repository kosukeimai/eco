###
### main function
###


#' Fitting Parametric Models and Quantifying Missing Information for Ecological
#' Inference in 2x2 Tables
#' 
#' \code{ecoML} is used to fit parametric models for ecological inference in
#' \eqn{2 \times 2} tables via Expectation Maximization (EM) algorithms. The
#' data is specified in proportions. At it's most basic setting, the algorithm
#' assumes that the individual-level proportions (i.e., \eqn{W_1} and
#' \eqn{W_2}) and distributed bivariate normally (after logit transformations).
#' The function calculates point estimates of the parameters for models based
#' on different assumptions. The standard errors of the point estimates are
#' also computed via Supplemented EM algorithms. Moreover, \code{ecoML}
#' quantifies the amount of missing information associated with each parameter
#' and allows researcher to examine the impact of missing information on
#' parameter estimation in ecological inference. The models and algorithms are
#' described in Imai, Lu and Strauss (2008, 2011).
#' 
#' When \code{SEM} is \code{TRUE}, \code{ecoML} computes the observed-data
#' information matrix for the parameters of interest based on Supplemented-EM
#' algorithm. The inverse of the observed-data information matrix can be used
#' to estimate the variance-covariance matrix for the parameters estimated from
#' EM algorithms. In addition, it also computes the expected complete-data
#' information matrix. Based on these two measures, one can further calculate
#' the fraction of missing information associated with each parameter. See
#' Imai, Lu and Strauss (2006) for more details about fraction of missing
#' information.
#' 
#' Moreover, when \code{hytest=TRUE}, \code{ecoML} allows to estimate the
#' parametric model under the null hypothesis that \code{mu_1=mu_2}. One can
#' then construct the likelihood ratio test to assess the hypothesis of equal
#' means. The associated fraction of missing information for the test statistic
#' can be also calculated. For details, see Imai, Lu and Strauss (2006) for
#' details.
#' 
#' @param formula A symbolic description of the model to be fit, specifying the
#' column and row margins of \eqn{2 \times 2} ecological tables. \code{Y ~ X}
#' specifies \code{Y} as the column margin (e.g., turnout) and \code{X} (e.g.,
#' percent African-American) as the row margin. Details and specific examples
#' are given below.
#' @param data An optional data frame in which to interpret the variables in
#' \code{formula}. The default is the environment in which \code{ecoML} is
#' called.
#' @param N An optional variable representing the size of the unit; e.g., the
#' total number of voters.  \code{N} needs to be a vector of same length as
#' \code{Y} and \code{X} or a scalar.
#' @param supplement An optional matrix of supplemental data. The matrix has
#' two columns, which contain additional individual-level data such as survey
#' data for \eqn{W_1} and \eqn{W_2}, respectively.  If \code{NULL}, no
#' additional individual-level data are included in the model. The default is
#' \code{NULL}.
#' @param fix.rho Logical. If \code{TRUE}, the correlation (when
#' \code{context=TRUE}) or the partial correlation (when \code{context=FALSE})
#' between \eqn{W_1} and \eqn{W_2} is fixed through the estimation. For
#' details, see Imai, Lu and Strauss(2006). The default is \code{FALSE}.
#' @param context Logical. If \code{TRUE}, the contextual effect is also
#' modeled. In this case, the row margin (i.e., X) and the individual-level
#' rates (i.e., \eqn{W_1} and \eqn{W_2}) are assumed to be distributed
#' tri-variate normally (after logit transformations). See Imai, Lu and Strauss
#' (2006) for details. The default is \code{FALSE}.
#' @param sem Logical. If \code{TRUE}, the standard errors of parameter
#' estimates are estimated via SEM algorithm, as well as the fraction of
#' missing data. The default is \code{TRUE}.
#' @param theta.start A numeric vector that specifies the starting values for
#' the mean, variance, and covariance. When \code{context = FALSE}, the
#' elements of \code{theta.start} correspond to (\eqn{E(W_1)}, \eqn{E(W_2)},
#' \eqn{var(W_1)}, \eqn{var(W_2)}, \eqn{cor(W_1,W_2)}). When \code{context =
#' TRUE}, the elements of \code{theta.start} correspond to (\eqn{E(W_1)},
#' \eqn{E(W_2)}, \eqn{var(W_1)}, \eqn{var(W_2)}, \eqn{corr(W_1, X)},
#' \eqn{corr(W_2, X)}, \eqn{corr(W_1,W_2)}). Moreover, when
#' \code{fix.rho=TRUE}, \eqn{corr(W_1,W_2)} is set to be the correlation
#' between \eqn{W_1} and \eqn{W_2} when \code{context = FALSE}, and the partial
#' correlation between \eqn{W_1} and \eqn{W_2} given \eqn{X} when \code{context
#' = FALSE}. The default is \code{c(0,0,1,1,0)}.
#' @param epsilon A positive number that specifies the convergence criterion
#' for EM algorithm. The square root of \code{epsilon} is the convergence
#' criterion for SEM algorithm. The default is \code{10^(-6)}.
#' @param maxit A positive integer specifies the maximum number of iterations
#' before the convergence criterion is met. The default is \code{1000}.
#' @param loglik Logical. If \code{TRUE}, the value of the log-likelihood
#' function at each iteration of EM is saved. The default is \code{TRUE}.
#' @param hyptest Logical. If \code{TRUE}, model is estimated under the null
#' hypothesis that means of \eqn{W1} and \eqn{W2} are the same.  The default is
#' \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, the progress of the EM and SEM
#' algorithms is printed to the screen. The default is \code{FALSE}.
#' @return An object of class \code{ecoML} containing the following elements:
#' \item{call}{The matched call.} 
#' \item{X}{The row margin, \eqn{X}.}
#' \item{Y}{The column margin, \eqn{Y}.} 
#' \item{N}{The size of each table, \eqn{N}.} 
#' \item{context}{The assumption under which model is estimated. If
#' \code{context = FALSE}, CAR assumption is adopted and no contextual effect
#' is modeled. If \code{context = TRUE}, NCAR assumption is adopted, and
#' contextual effect is modeled.} \item{sem}{Whether SEM algorithm is used to
#' estimate the standard errors and observed information matrix for the
#' parameter estimates.} 
#' \item{fix.rho}{Whether the correlation or the partial
#' correlation between \eqn{W_1} an \eqn{W_2} is fixed in the estimation.}
#' \item{r12}{If \code{fix.rho = TRUE}, the value that \eqn{corr(W_1, W_2)} is
#' fixed to.} 
#' \item{epsilon}{The precision criterion for EM convergence.
#' \eqn{\sqrt{\epsilon}} is the precision criterion for SEM convergence.}
#' \item{theta.sem}{The ML estimates of \eqn{E(W_1)},\eqn{E(W_2)},
#' \eqn{var(W_1)},\eqn{var(W_2)}, and \eqn{cov(W_1,W_2)}. If \code{context =
#' TRUE}, \eqn{E(X)},\eqn{cov(W_1,X)}, \eqn{cov(W_2,X)} are also reported.}
#' \item{W}{In-sample estimation of \eqn{W_1} and \eqn{W_2}.}
#' \item{suff.stat}{The sufficient statistics for \code{theta.em}.}
#' \item{iters.em}{Number of EM iterations before convergence is achieved.}
#' \item{iters.sem}{Number of SEM iterations before convergence is achieved.}
#' \item{loglik}{The log-likelihood of the model when convergence is achieved.}
#' \item{loglik.log.em}{A vector saving the value of the log-likelihood
#' function at each iteration of the EM algorithm.} 
#' \item{mu.log.em}{A matrix saving the unweighted mean estimation of the 
#' logit-transformed individual-level proportions (i.e., \eqn{W_1} and \eqn{W_2}) 
#' at each iteration of the EM process.} \item{Sigma.log.em}{A matrix saving the 
#' log of the variance estimation of the logit-transformed individual-level
#' proportions (i.e., \eqn{W_1} and \eqn{W_2}) at each iteration of EM process.
#' Note, non-transformed variances are displayed on the screen (when
#' \code{verbose = TRUE}).} 
#' \item{rho.fisher.em}{A matrix saving the fisher
#' transformation of the estimation of the correlations between the
#' logit-transformed individual-level proportions (i.e., \eqn{W_1} and
#' \eqn{W_2}) at each iteration of EM process.  Note, non-transformed
#' correlations are displayed on the screen (when \code{verbose = TRUE}).}
#' Moreover, when \code{sem=TRUE}, \code{ecoML} also output the following
#' values: 
#' \item{DM}{The matrix characterizing the rates of convergence of the
#' EM algorithms. Such information is also used to calculate the observed-data
#' information matrix} 
#' \item{Icom}{The (expected) complete data information
#' matrix estimated via SEM algorithm. When \code{context=FALSE, fix.rho=TRUE},
#' \code{Icom} is 4 by 4. When \code{context=FALSE, fix.rho=FALSE}, \code{Icom}
#' is 5 by 5. When \code{context=TRUE}, \code{Icom} is 9 by 9.} 
#' \item{Iobs}{The observed information matrix. The dimension of \code{Iobs} 
#' is same as \code{Icom}.} 
#' \item{Imiss}{The difference between \code{Icom} and \code{Iobs}.  
#' The dimension of \code{Imiss} is same as \code{miss}.}
#' \item{Vobs}{The (symmetrized) variance-covariance matrix of the ML parameter
#' estimates. The dimension of \code{Vobs} is same as \code{Icom}.}
#' \item{Iobs}{The (expected) complete-data variance-covariance matrix.  The
#' dimension of \code{Iobs} is same as \code{Icom}.} 
#' \item{Vobs.original}{The estimated variance-covariance matrix of the ML parameter 
#' estimates. The dimension of \code{Vobs} is same as \code{Icom}.} 
#' \item{Fmis}{The fraction of missing information associated with each parameter estimation. }
#' \item{VFmis}{The proportion of increased variance associated with each
#' parameter estimation due to observed data. } 
#' \item{Ieigen}{The largest eigen value of \code{Imiss}.} 
#' \item{Icom.trans}{The complete data information
#' matrix for the fisher transformed parameters.} 
#' \item{Iobs.trans}{The observed data information matrix for the fisher transformed parameters.}
#' \item{Fmis.trans}{The fractions of missing information associated with the
#' fisher transformed parameters.}
#' @seealso \code{eco}, \code{ecoNP}, \code{summary.ecoML}
#' @references Imai, Kosuke, Ying Lu and Aaron Strauss. (2011).  \dQuote{eco: R
#' Package for Ecological Inference in 2x2 Tables} Journal of Statistical
#' Software, Vol. 42, No. 5, pp. 1-23.
#' 
#' Imai, Kosuke, Ying Lu and Aaron Strauss. (2008).  \dQuote{Bayesian and
#' Likelihood Inference for 2 x 2 Ecological Tables: An Incomplete Data
#' Approach} Political Analysis, Vol. 16, No. 1 (Winter), pp. 41-69. 
#' @keywords models
#' @examples
#' 
#' 
#' ## load the census data
#' data(census)
#' 
#' ## NOTE: convergence has not been properly assessed for the following
#' ## examples. See Imai, Lu and Strauss (2006) for more complete analyses.
#' ## In the first example below, in the interest of time, only part of the
#' ## data set is analyzed and the convergence requirement is less stringent
#' ## than the default setting.
#' 
#' ## In the second example, the program is arbitrarily halted 100 iterations
#' ## into the simulation, before convergence.
#' 
#' ## load the Robinson's census data
#' data(census)
#' 
#' ## fit the parametric model with the default model specifications
#' \dontrun{res <- ecoML(Y ~ X, data = census[1:100,], N=census[1:100,3], 
#' 	     	  epsilon=10^(-6), verbose = TRUE)}
#' ## summarize the results
#' \dontrun{summary(res)}
#' 
#' ## fit the parametric model with some individual 
#' ## level data using the default prior specification
#' surv <- 1:600
#' \dontrun{res1 <- ecoML(Y ~ X, context = TRUE, data = census[-surv,], 
#'                    supplement = census[surv,c(4:5,1)], maxit=100, verbose = TRUE)}
#' ## summarize the results
#' \dontrun{summary(res1)}
#' 
#' @export ecoML
ecoML <- function(formula, data = parent.frame(), N=NULL, supplement = NULL, 
                  theta.start = c(0,0,1,1,0), fix.rho = FALSE,
                  context = FALSE, sem = TRUE, epsilon=10^(-6),
                  maxit = 1000, loglik = TRUE, hyptest=FALSE, verbose= FALSE) { 

  
  ## getting X and Y
  mf <- match.call()
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.matrix(eval.parent(mf$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  Y <- model.response(model.frame(tt, data=data))
  
  #n.var: total number of parameters involved in the estimation
  #n.par: number of nonstatic paramters need to estimate through EM 
  #       also need SEM
  #ndim: dimension of the multivariate normal distribution

  ndim<-2
  if (context) ndim<-3

  n.var<-2*ndim+ ndim*(ndim-1)/2
 
  n.par<-n.S<-n.var 
  if (context) {
      n.par<-n.var-2
   }

  r12<-NULL
  if (fix.rho) 
     r12<-theta.start[n.par]

  if (!context & fix.rho) n.par<-n.par-1

  flag<-as.integer(context)+2*as.integer(fix.rho)+2^2*as.integer(sem)

  ##checking data
  tmp <- checkdata(X, Y, supplement, ndim)
  bdd <- ecoBD(formula=formula, data=data)
  W1min <- bdd$Wmin[order(tmp$order.old)[1:nrow(tmp$d)],1,1]
  W1max <- bdd$Wmax[order(tmp$order.old)[1:nrow(tmp$d)],1,1]


  n <- tmp$n.samp+tmp$survey.samp+tmp$samp.X1+tmp$samp.X0
  wcol<-ndim
  if (context) {
    wcol<-wcol-1
  }
  inSample.length <- wcol*tmp$n.samp

  #if NCAR and the user did not provide a theta.start
  if (context && (length(theta.start)==5) ) 
    theta.start<-c(0,0,1,1,0,0,0)

  ## Fitting the model via EM  
  res <- .C("cEMeco", as.double(tmp$d), as.double(theta.start),
            as.integer(tmp$n.samp),  as.integer(maxit), as.double(epsilon),
            as.integer(tmp$survey.yes), as.integer(tmp$survey.samp), 
            as.double(tmp$survey.data),
            as.integer(tmp$X1type), as.integer(tmp$samp.X1), as.double(tmp$X1.W1),
            as.integer(tmp$X0type), as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
            as.double(W1min), as.double(W1max),
            as.integer(flag),as.integer(verbose),as.integer(loglik),as.integer(hyptest),
            optTheta=rep(-1.1,n.var), pdTheta=double(n.var),
            S=double(n.S+1),inSample=double(inSample.length),DMmatrix=double(n.par*n.par),
            itersUsed=as.integer(0),history=double((maxit+1)*(n.var+1)),
            PACKAGE="eco")

  ##record results from EM
  theta.em<-res$pdTheta
  theta.fisher<-param.trans(theta.em, transformation="Fisher")
  iters.em<-res$itersUsed
  mu.log.em <- matrix(rep(NA,iters.em*ndim),ncol=ndim)
  sigma.log.em <- matrix(rep(NA,iters.em*ndim),ncol=ndim)
  loglike.log.em <- as.double(rep(NA,iters.em))
  nrho<-length(theta.em)-2*ndim
  rho.fisher.em <- matrix(rep(NA,iters.em*nrho),ncol=nrho)
  for(i in 1:iters.em) {
    mu.log.em[i,1:ndim]=res$history[(i-1)*(n.var+1)+(1:ndim)]
    sigma.log.em[i,1:ndim]=res$history[(i-1)*(n.var+1)+ndim+(1:ndim)]
     if (nrho!=0)
    rho.fisher.em[i, 1:nrho]=res$history[(i-1)*(n.var+1)+2*ndim+(1:nrho)]
    loglike.log.em[i]=res$history[(i-1)*(n.var+1)+2*ndim+nrho+1]
  }

  ## In sample prediction of W
  W <- matrix(rep(NA,inSample.length),ncol=wcol)
  for (i in 1:tmp$n.samp)
    for (j in 1:wcol)
      W[i,j]=res$inSample[(i-1)*2+j]

  ## SEM step
  iters.sem<-0

   suff.stat<-res$S
  if (context)
      {
     suff.stat<-rep(0,(n.var+1))
         suff.stat[1]<-mean(logit(c(X,supplement[,3])))
         suff.stat[2:3]<-res$S[1:2]
         suff.stat[4]<-mean((logit(c(X, supplement[,3])))^2)
         suff.stat[5:6]<-res$S[3:4]
         suff.stat[7:8]<-res$S[6:7]
         suff.stat[9]<-res$S[5]
         suff.stat[10]<-res$S[8]
      }


  if (sem) 
  {

    DM <- matrix(rep(NA,n.par*n.par),ncol=n.par)

    res <- .C("cEMeco", as.double(tmp$d), as.double(theta.start),
              as.integer(tmp$n.samp),  as.integer(maxit), as.double(epsilon),
              as.integer(tmp$survey.yes), as.integer(tmp$survey.samp), 
              as.double(tmp$survey.data),
              as.integer(tmp$X1type), as.integer(tmp$samp.X1), as.double(tmp$X1.W1),
              as.integer(tmp$X0type), as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
              as.double(bdd$Wmin[,1,1]), as.double(bdd$Wmax[,1,1]),
              as.integer(flag),as.integer(verbose),as.integer(loglik),as.integer(hyptest),
              res$pdTheta, pdTheta=double(n.var), S=double(n.S+1),
              inSample=double(inSample.length),DMmatrix=double(n.par*n.par),
              itersUsed=as.integer(0),history=double((maxit+1)*(n.var+1)),
              PACKAGE="eco")     
  
    iters.sem<-res$itersUsed
    for(i in 1:n.par)
      for(j in 1:n.par)
        DM[i,j]=res$DMmatrix[(i-1)*n.par+j]


} 

 
   
  if (!context) names(theta.em)<-c("u1","u2","s1","s2","r12")
  if (context) names(theta.em)<-c("ux","u1","u2","sx","s1","s2","r1x","r2x","r12")



  ## output
  res.out<-list(call = mf, Y = Y, X = X, N = N, 
                fix.rho = fix.rho, context = context, sem=sem, epsilon=epsilon,
        theta.em=theta.em, r12=r12, 
               sigma.log = theta.fisher[(ndim+1):(2*ndim)], suff.stat = suff.stat[1:n.S],
                loglik = res$S[n.S+1], iters.em = iters.em, 
                iters.sem = iters.sem, mu.log.em = mu.log.em, 
                sigma.log.em = sigma.log.em,
                rho.fisher.em = rho.fisher.em, loglike.log.em = loglike.log.em,
                W = W)
  
  if (sem) {
    res.out$DM<-DM
#print(dim(data))
# n<-dim(data)[1]
if (!is.null(supplement)) n<-n+dim(supplement)[1]
#cat("n2=", n,"\n")

 res.info<- ecoINFO(theta.em=res.out$theta.em, suff.stat=res.out$suff.stat, DM=res.out$DM, context=context, fix.rho=fix.rho, sem=sem, r12=res.out$r12, n=n)

    res.out$DM<-res.info$DM
    res.out$Icom<-res.info$Icom
    res.out$Iobs<-res.info$Iobs
    res.out$Fmis<-res.info$Fmis
    res.out$Vobs.original<-res.info$Vobs.original
    res.out$Vobs<-res.info$Vobs
    res.out$Iobs<-res.info$Iobs
    res.out$VFmis<-res.info$VFmis
    res.out$Icom.trans<-res.info$Icom.trans
    res.out$Iobs.trans<-res.info$Iobs.trans
    res.out$Fmis.trans<-res.info$Fmis.trans
    res.out$Imiss<-res.info$Imiss
    res.out$Ieigen<-res.info$Ieigen

 res.out$Iobs<-res.info$Iobs
}

  class(res.out) <- "ecoML"
  return(res.out)
}
