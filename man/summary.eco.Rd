\name{summary.eco}

\alias{summary.eco}
\alias{print.eco}
\alias{print.summary.eco}

\title{Summarizing the Results for the Bayesian Parametric Model for
  Ecological Inference in 2x2 Tables}

\description{
  \code{summary} method for class \code{eco}.
}

\usage{
  \method{summary}{eco}(object, CI = c(2.5, 97.5), param = TRUE,
          units = FALSE, subset = NULL, ...) 

  \method{print}{summary.eco}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{object}{An output object from \code{eco}.}
  \item{CI}{A vector of lower and upper bounds for the Bayesian credible
    intervals used to summarize the results. The default is the
    equal tail 95 percent credible interval.
  }
  \item{x}{An object of class \code{summary.eco}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{param}{Logical. If \code{TRUE}, the posterior estimates of the 
    population parameters will be provided. The default value is
    \code{TRUE}.
  }
  \item{units}{Logical. If \code{TRUE}, the in-sample predictions for
    each unit or for a subset of units will be provided. The default 
    value is \code{FALSE}.
  } 
  \item{subset}{A numeric vector indicating the subset of the units whose 
    in-sample predications to be provided when \code{units} is 
    \code{TRUE}. The default value is \code{NULL} where the in-sample 
    predictions for each unit will be provided.
  } 
  \item{...}{further arguments passed to or from other methods.}
}

\value{
  \code{summary.eco} yields an object of class \code{summary.eco}
  containing the following elements:
  \item{call}{The call from \code{eco}.}
  \item{n.obs}{The number of units.}
  \item{n.draws}{The number of Monte Carlo samples.}
  \item{agg.table}{Aggregate posterior estimates of the marginal means
    of \eqn{W_1} and \eqn{W_2} using \eqn{X} and \eqn{N} as weights.}
  If \code{param = TRUE}, the following elements are also included:
  \item{param.table}{Posterior estimates of model parameters: population
    mean estimates of \eqn{W_1} and \eqn{W_2} and their logit
    transformations.} 
  If \code{units = TRUE}, the following elements are also included:
  \item{W1.table}{Unit-level posterior estimates for \eqn{W_1}.}
  \item{W2.table}{Unit-level posterior estimates for \eqn{W_2}.}

  This object can be printed by \code{print.summary.eco}
}

\seealso{\code{eco}, \code{predict.eco}}

\author{
  Kosuke Imai, Department of Politics, Princeton University,
  \email{kimai@Princeton.Edu}, \url{http://imai.princeton.edu};
  Ying Lu, Department of Sociology, University of Colorado at Boulder, 
  \email{ying.lu@Colorado.Edu}
  }

\keyword{methods}
