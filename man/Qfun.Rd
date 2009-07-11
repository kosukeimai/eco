\name{Qfun}

\alias{Qfun}

\title{Fitting the Parametric Bayesian Model of Ecological Inference in
  2x2 Tables}

\description{ 
  \code{Qfun} returns the complete log-likelihood 
that is used to calculate the fraction of missing information.}

\usage{
Qfun(theta, suff.stat, n)
}

\arguments{
  \item{theta}{A vector that contains the MLE \eqn{E(W_1)},\eqn{E(W_2)},
    \eqn{var(W_1)},\eqn{var(W_2)}, and \eqn{cov(W_1,W_2)}. Typically it is
the element \code{theta.em} of an object of class \code{ecoML}.}
  \item{suff.stat}{A vector of sufficient statistics of \eqn{E(W_1)},
\eqn{E(W_2)}, \eqn{var(W_1)},\eqn{var(W_2)}, and \eqn{cov(W_1,W_2)}.}
\item{n}{A integer representing the sample size.}
}
\author{
  Kosuke Imai, Department of Politics, Princeton University,
  \email{kimai@Princeton.Edu}, \url{http://imai.princeton.edu};
  Ying Lu, Department of Humanities and Social Sciences in the 
Professions, Steinhardt School of Culture, Education and Human 
Development, New York University, 
  \email{yl46@Nyu.Edu}
  Aaron Strauss, Department of Politics, Princeton University,
 \email{abstraus@Princeton.Edu}.
  }

  \references{
  Imai, Kosuke, Ying Lu and Aaron Strauss. (Forthcoming).
  \dQuote{eco: R Package for Ecological Inference in 2x2 Tables}
  Journal of Statistical Software, available at
  \url{http://imai.princeton.edu/research/eco.html}

  Imai, Kosuke, Ying Lu and Aaron Strauss. (2008). \dQuote{Bayesian and
    Likelihood Inference for 2 x 2 Ecological Tables: An Incomplete Data
    Approach} Political Analysis, Vol. 16, No. 1 (Winter), pp. 41-69. available at 
  \url{http://imai.princeton.edu/research/eiall.html}
}

\seealso{\code{ecoML}}
\keyword{models}
