#' Fitting the Parametric Bayesian Model of Ecological Inference in 2x2 Tables
#' 
#' \code{Qfun} returns the complete log-likelihood that is used to calculate
#' the fraction of missing information.
#' 
#' 
#' @param theta A vector that contains the MLE \eqn{E(W_1)},\eqn{E(W_2)},
#' \eqn{var(W_1)},\eqn{var(W_2)}, and \eqn{cov(W_1,W_2)}. Typically it is the
#' element \code{theta.em} of an object of class \code{ecoML}.
#' @param suff.stat A vector of sufficient statistics of \eqn{E(W_1)},
#' \eqn{E(W_2)}, \eqn{var(W_1)},\eqn{var(W_2)}, and \eqn{cov(W_1,W_2)}.
#' @param n A integer representing the sample size.
#' @return A single numeric value: the complete-data log-likelihood.
#' @seealso \code{ecoML}
#' @references Imai, Kosuke, Ying Lu and Aaron Strauss. (2011).  \dQuote{eco: R
#' Package for Ecological Inference in 2x2 Tables} Journal of Statistical
#' Software, Vol. 42, No. 5, pp. 1-23. 
#' 
#' Imai, Kosuke, Ying Lu and Aaron Strauss. (2008). \dQuote{Bayesian and
#' Likelihood Inference for 2 x 2 Ecological Tables: An Incomplete Data
#' Approach} Political Analysis, Vol. 16, No. 1 (Winter), pp. 41-69.
#' @keywords models
#' @export Qfun
Qfun <- function(theta, suff.stat, n) {
  mu<-rep(0,2)
  Sigma<-matrix(0, 2,2)
  Suff1<-rep(0,2)
  Suff2<-matrix(0,2,2)
  
  mu <- theta[1:2]
  Sigma[1,1]<-theta[3]
  Sigma[2,2]<-theta[4]
  Sigma[1,2]<-Sigma[2,1]<-theta[5]*sqrt(Sigma[1,1]*Sigma[2,2])

  Suff1 <- n*suff.stat[1:2]
  Suff2[1,1]<-n*suff.stat[3]
  Suff2[2,2]<-n*suff.stat[4]
  Suff2[1,2]<-n*suff.stat[5]
 
  invSigma<-solve(Sigma)

  return(-0.5*n*log(det(Sigma))-0.5*sum(diag(invSigma%*%(Suff2-mu%*%t(Suff1)-Suff1%*%t(mu)+n*mu%*%t(mu)))))

}
