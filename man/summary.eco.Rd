\name{summary.eco}

\alias{summary.eco}
\alias{print.summary.eco}

\title{Summarizing the results for the nonparametric and parametric models for ecological inference in 2x2 tables }

\description{
  \code{summary} method for class \code{eco}.
}

\usage{
  \method{summary}{eco}(object, long=FALSE, CI=c(2.5, 97.5), ...)

  \method{print}{summary.eco}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{object}{An output object from \code{eco}.}
  \item{long}{Logical. If \code{TRUE}, the in-sample predictions for
    each unit will be provided. The default value is \code{FALSE}.} 
  \item{CI}{A 2 dimensional vector of lower and upper bounds for the
    credible intervals used to summarize the results. The default is the
    equal tail 95 percent credible interval.
  }
  \item{x}{An object of class \code{summary.eco}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{...}{further arguments passed to or from other methods.}
}

\value{
  \code{summary.eco} yields an object of class \code{summary.eco}
  containing the following elements:
  \item{call}{The call from \code{eco}.}
  \item{n.obs}{The number of observations.}
  \item{coef.table}{The summary of the posterior distribution of the
    aggregate means of W1 and W2. } 
}

\seealso{\code{eco}}

\author{
  Ying Lu, Woodrow Wilson School of International and Public Affairs,
  Princeton University \email{yinglu@princeton.edu}; Kosuke Imai,
  Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}, \url{http://www.princeton.edu/~kimai}
}

\keyword{methods}
