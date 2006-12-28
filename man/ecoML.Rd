\name{ecoML}

\alias{ecoML}

\title{Fitting Parametric Models and Quantifying Missing Information
 for Ecological Inference in 2x2 Tables}

\description{
  \code{ecoML} is used to fit parametric models for ecological 
  inference in \eqn{2 \times 2} tables via Expectation Maximization (EM)
  algorithms. It gives the point estimates of the parameters for models
  based on different assumptions. The standard errors of the point
  estimates are also computed via Supplemented EM algorithms. Moreover,
  \code{ecoML} quantifies the amount of missing information associated
  with each parameter and allows researcher to examine the impact of
  missing information on parameter estimation in ecological
  inference. The models and algorithms are described in Imai,
  Lu and Strauss(2006).
}

\usage{
   ecoML(formula, data = parent.frame(), N = NULL, supplement = NULL, 
         theta.start = c(0,0,1,1,0), fix.rho = FALSE,
         context = FALSE, sem = TRUE, epsilon = 10^(-10), 
     maxit = 1000, loglik = TRUE, hyptest = FALSE, verbose = TRUE)  
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit,
    specifying the column and row margins of \eqn{2 \times
      2} ecological tables. \code{Y ~ X} specifies \code{Y} as the
    column margin and \code{X} as the row margin. Details and specific
    examples are given below.
  } 
  \item{data}{An optional data frame in which to interpret the variables
    in \code{formula}. The default is the environment in which
    \code{ecoML} is called. 
  }
  \item{N}{An optional variable representing the size of the unit; e.g.,
    the total number of voters.
  }
  \item{supplement}{An optional matrix of supplemental data. The matrix
    has two columns, which contain additional individual-level data such
    as survey data for \eqn{W_1} and \eqn{W_2}, respectively.  If
    \code{NULL}, no additional individual-level data are included in the
    model. The default is \code{NULL}.
  }
  \item{fix.rho}{Logical. If \code{TRUE}, the correlation
    (when \code{context=TRUE}) or the partial correlation (when
    \code{context=FALSE}) between \eqn{W_1} and \eqn{W_2} 
    is fixed through the estimation. For details, see 
    Imai, Lu and Strauss(2006). The default is \code{FALSE}.
  } 
  \item{context}{Logical. If \code{TRUE}, the contextual effect is also
    modeled. See Imai, Lu and Strauss (2006) for details. The default is
    \code{FALSE}. 
  }
  \item{sem}{Logical. If \code{TRUE}, the standard errors of parameter
   estimates are estimated via SEM algorithm. The default is
    \code{TRUE}. 
  }
  \item{theta.start}{A numeric vector that specifies the starting values
    for the mean, variance, and covariance. When \code{context = FALSE},
    the elements of \code{theta.start} correspond to (\eqn{E(W_1)},
    \eqn{E(W_2)}, \eqn{var(W_1)}, \eqn{var(W_2)},
    \eqn{cor(W_1,W_2)}). When \code{context = TRUE}, the
    elements of \code{theta.start} correspond to (\eqn{E(W_1)},
    \eqn{E(W_2)}, \eqn{var(W_1)}, \eqn{var(W_2)}, \eqn{corr(W_1, X)},
    \eqn{corr(W_2, X)}, \eqn{corr(W_1,W_2)}). Moreover, when
    \code{fix.rho=TRUE}, \eqn{corr(W_1,W_2)} is set to be the
    correlation between \eqn{W_1} and \eqn{W_2} when \code{context =
      FALSE}, and the partial correlation between \eqn{W_1} and
    \eqn{W_2} given \eqn{X} when \code{context = FALSE}. The default is
    \code{c(0,0,1,1,0)}. 
  }
  \item{epsilon}{A positive number that specifies the convergence criterion
    for EM algorithm. The square root of \code{epsilon} is the convergence 
    criterion for SEM algorithm. The default is \code{10^(-10)}. 
  } 
  \item{maxit}{A positive integer specifies the maximum number of iterations
    before the convergence criterion is met. The default is \code{1000}.
  } 
  \item{loglik}{Logical. If \code{TRUE}, the value of the log-likelihood
    function at each iteration of EM is saved. The default is
    \code{TRUE}.
  } 
  \item{hyptest}{Logical. If \code{TRUE}, model is estimated under the null
    hypothesis that mean of \eqn{W1} and \eqn{W2} is the same. 
    The default is \code{FALSE}. 
  } 
  \item{verbose}{Logical. If \code{TRUE}, the progress of the EM and SEM
    algorithms is printed to the screen. The default is \code{FALSE}.
  }
}

\details{
  When \code{SEM} is \code{TRUE}, \code{ecoML} computes the observed-data 
  information matrix for the parameters of interest based on Supplemented-EM 
  algorithm. The inverse of the observed-data information matrix can be used 
  to estimate the variance-covariance matrix for the parameters estimated
  from EM algorithms. In addition, it also computes the expected complete-data 
  information matrix. Based on these two measures, one can further calculate 
  the fraction of missing information associated with each parameter. See
  Imai, Lu and Strauss (2006) for more details about fraction of missing
  information.
 
  Moreover, when \code{hytest=TRUE}, \code{ecoML} allows to estimate the 
  parametric model under the null hypothesis that \code{mu_1=mu_2}. One 
  can then construct the likelihood ratio test to assess the hypothesis of 
  equal means. The associated fraction of missing information for the test 
  statistic can be also calculated. For details, see Imai, Lu
  and Strauss (2006) for details.
}
\examples{

## load the census data
data(census)

## NOTE: convergence has not been properly assessed for the following
## examples. See Imai, Lu and Strauss (2006) for more complete analyses.
## In the first example below, in the interest of time, only part of the
## data set is analyzed and the convergence requirement is less stringent
## than the default setting.

## In the second example, the program is arbitrarily halted 100 iterations
## into the simulation, before convergence.

## load the Robinson's census data
data(census)

## fit the parametric model with the default model specifications
\dontrun{res <- ecoML(Y ~ X, data = census[1:100,],epsilon=10^(-6), verbose = TRUE)}
## summarize the results
\dontrun{summary(res)}

## obtain out-of-sample prediction
\dontrun{out <- predict(res, verbose = TRUE)}
## summarize the results
\dontrun{summary(out)}

## fit the parametric model with some individual 
## level data using the default prior specification
surv <- 1:600
\dontrun{res1 <- ecoML(Y ~ X, context = TRUE, data = census[-surv,], 
                   supplement = census[surv,c(4:5,1)], maxit=100, verbose = TRUE)}
## summarize the results
\dontrun{summary(res1)}

}

\value{
  An object of class \code{ecoML} containing the following elements:
  \item{call}{The matched call.}
  \item{X}{The row margin, \eqn{X}.}
  \item{Y}{The column margin, \eqn{Y}.}
  \item{N}{The size of each table, \eqn{N}.}
  \item{context}{The assumption under which model is estimated. If 
    \code{context = FALSE}, CAR assumption is adopted and no
    contextual effect is modeled. If \code{context = TRUE}, NCAR
    assumption is adopted, and contextual effect is modeled.}
  \item{sem}{Whether SEM algorithm is used to estimate the standard
    errors and observed information matrix for the parameter estimates.}
  \item{fix.rho}{Whether the correlation or the partial correlation between
    \eqn{W_1} an \eqn{W_2} is fixed in the estimation.}
  \item{r12}{If \code{fix.rho = TRUE}, the value that \eqn{corr(W_1,
      W_2)} is fixed to.}
  \item{epsilon}{The precision criterion for EM convergence. 
    \eqn{\sqrt{\epsilon}} is the precision criterion for SEM convergence.}
  \item{theta.sem}{The ML estimates of \eqn{E(W_1)},\eqn{E(W_2)},
    \eqn{var(W_1)},\eqn{var(W_2)}, and \eqn{cov(W_1,W_2)}. If
    \code{context = TRUE}, \eqn{E(X)},\eqn{cov(W_1,X)}, 
    \eqn{cov(W_2,X)} are also reported.}
  \item{W}{In-sample estimation of \eqn{W_1} and \eqn{W_2}.}
  \item{suff.stat}{The sufficient statistics for \code{theta.em}.}
  \item{iters.em}{Number of EM iterations before convergence is achieved.}
  \item{iters.sem}{Number of SEM iterations before convergence is achieved.}
  \item{loglik}{The log-likelihood of the model when convergence is
    achieved.}
  \item{loglik.log.em}{A vector saving the value of the log-likelihood
    function at each iteration of the EM algorithm.}
  \item{mu.log.em}{A matrix saving the mean estimation at each 
    iteration of EM.}
  \item{Sigma.log.em}{A matrix saving the variance-covariance estimation 
    at each iteration of EM.} 
  Moreover, when \code{sem=TRUE}, \code{ecoML} also output the following
  values:
  \item{DM}{The matrix characterizing the rates of convergence of the EM 
    algorithms. Such information is also used to calculate the observed-data
    information matrix}
  \item{Icom}{The (expected) complete data information matrix estimated 
    via SEM algorithm. When \code{context=FALSE, fix.rho=TRUE}, 
    \code{Icom} is 4 by 4. When \code{context=FALSE, fix.rho=FALSE}, 
    \code{Icom} is 5 by 5. When \code{context=TRUE}, \code{Icom} 
    is 9 by 9.}
  \item{Iobs}{The observed information matrix. The dimension of 
    \code{Iobs} is same as \code{Icom}.}
  \item{Imiss}{The difference between \code{Icom} and \code{Iobs}. 
    The dimension of \code{Imiss} is same as \code{miss}.}
  \item{Vobs}{The (symmetrized) variance-covariance matrix of the ML parameter
    estimates. The dimension of \code{Vobs} is same as 
    \code{Icom}.}
  \item{Iobs}{The (expected) complete-data variance-covariance matrix. 
    The dimension of \code{Iobs} is same as \code{Icom}.}
  \item{Vobs.original}{The estimated variance-covariance matrix of the 
    ML parameter  estimates. The dimension of \code{Vobs} is same as 
    \code{Icom}.}
  \item{Fmis}{The fraction of missing information associated with each 
    parameter estimation. }
  \item{VFmis}{The proportion of increased variance associated with each 
    parameter estimation due to observed data. }
  \item{Ieigen}{The largest eigen value of \code{Imiss}.}
  \item{DM}{The rate-of-convergence matrix of the SEM algorithm.}
  \item{Icom.trans}{The complete data information matrix for the fisher 
    transformed parameters.}
  \item{Iobs.trans}{The observed data information matrix for the fisher 
    transformed parameters.}
  \item{Fmis.trans}{The fractions of missing information associated with 
    the fisher transformed parameters.}
}

\author{
  Kosuke Imai, Department of Politics, Princeton University,
  \email{kimai@Princeton.Edu}, \url{http://imai.princeton.edu};
  Ying Lu, Department of Sociology, University of Colorado at Boulder, 
  \email{ying.lu@Colorado.Edu};
  Aaron Strauss, Department of Politics, Princeton University,
 \email{abstraus@Princeton.Edu}.
}

\references{
  Imai, Kosuke, Ying Lu and Aaron Strauss. (2006) \dQuote{Bayesian and
    Likelihood Inference for 2 x 2 Ecological Tables: An Incomplete Data
    Approach} Technical Report, Princeton University, available at 
  \url{http://imai.princeton.edu/research/eiall.html}
}

\seealso{\code{eco}, \code{ecoNP}, \code{summary.ecoML}}
\keyword{models}
