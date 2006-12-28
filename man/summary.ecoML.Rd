\name{summary.ecoML}

\alias{summary.ecoML}
\alias{print.summary.ecoML}

\title{Summarizing the Results for the Maximum Likelihood Parametric Model for
  Ecological Inference in 2x2 Tables}

\description{
  \code{summary} method for class \code{eco}.
}

\usage{
  \method{summary}{ecoML}(object, CI = c(2.5, 97.5),  param = TRUE, units = FALSE,
     subset = NULL, ...) 

  \method{print}{summary.ecoML}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{object}{An output object from \code{eco}.}
  \item{CI}{A vector of lower and upper bounds for the Bayesian credible
    intervals used to summarize the results. The default is the
    equal tail 95 percent credible interval.
  }
  \item{param}{Ignored.}
  \item{subset}{A numeric vector indicating the subset of the units whose 
    in-sample predications to be provided when \code{units} is 
    \code{TRUE}. The default value is \code{NULL} where the in-sample 
    predictions for each unit will be provided.
  }
  \item{units}{Logical. If \code{TRUE}, the in-sample predictions for
    each unit or for a subset of units will be provided. The default 
    value is \code{FALSE}.
  } 
  
  \item{x}{An object of class \code{summary.ecoML}.}
  \item{digits}{the number of significant digits to use when printing.}
 
  \item{...}{further arguments passed to or from other methods.}
}

\value{
  \code{summary.eco} yields an object of class \code{summary.eco}
  containing the following elements:
  \item{call}{The call from \code{eco}.}
  \item{sem}{Whether the SEM algorithm was executed, as specified by the
    user upon calling \code{ecoML}.}
  \item{fix.rho}{Whether the correlation parameter was fixed or allowed to
    vary, as specified by the user upon calling \code{ecoML}.}
  \item{epsilon}{The convergence threshold specified by the user upon
    calling \code{ecoML}.}
  \item{n.obs}{The number of units.}
  \item{iters.em}{The number iterations the EM algorithm cycled through
    before convergence or reaching the maximum number of iterations
    allowed.}
  \item{iters.sem}{The number iterations the SEM algorithm cycled through
    before convergence or reaching the maximum number of iterations
    allowed.}
  \item{loglik}{The final observed log-likelihood.}
  \item{rho}{A matrix of \code{iters.em} rows specifying the correlation
    parameters at each iteration of the EM algorithm. The number of columns
    depends on how many correlation parameters exist in the model. Column
    order is the same as the order of the parameters in \code{param.table}.}
  \item{param.table}{Final estimates of the parameter values for the model. 
    Excludes parameters fixed by the user upon calling \code{ecoML}.}
  \item{agg.table}{Aggregate estimates of the marginal means
    of \eqn{W_1} and \eqn{W_2}}
  \item{agg.wtable}{Aggregate estimates of the marginal means
    of \eqn{W_1} and \eqn{W_2} using \eqn{X} and \eqn{N} as weights.}
  If \code{units = TRUE}, the following elements are also included:
  \item{W.table}{Unit-level estimates for \eqn{W_1} and \eqn{W_2}.}

  This object can be printed by \code{print.summary.eco}
}

\seealso{\code{eco}, \code{predict.eco}}

\author{
  Kosuke Imai, Department of Politics, Princeton University,
  \email{kimai@Princeton.Edu}, \url{http://imai.princeton.edu};
  Ying Lu, Department of Sociology, University of Colorado at Boulder, 
  \email{ying.lu@Colorado.Edu}; Aaron Strauss, Department of Politics,
  Princeton University, \email{abstraus@Princeton.Edu}
  }

\keyword{methods}
