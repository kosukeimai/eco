\name{summary.eco}

\alias{summary.eco}
\alias{print.summary.eco}

\title{Summarizing the results for the nonparametric and parametric models for ecological inference in 2 by 2 tables }

\description{
  \code{summary} method for class \code{eco}.
}

\usage{
  \method{summary}{eco}(object, CI=c(2.5, 97.5), long = FALSE, ...)

  \method{print}{summary.eco}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{object}{An output object from \code{eco}.}
  \item{CI}{A 2 dimensional vector of lower and upper bounds for the
    credible intervals used to summarize the results. The default is the
    equal tail 95 percent credible interval.
  }
  \item{x}{An object of class \code{summary.eco}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{long}{Logical. If \code{TRUE}, the in-sample predictions for
    each unit will be provided. The default value is \code{FALSE}.} 
  \item{...}{further arguments passed to or from other methods.}
}

\value{
  \code{summary.eco} yields an object of class \code{summary.eco}
  containing the following elements:
  \item{call}{The call from \code{eco}.}
  \item{nonpar}{The logical variable indicating whether the results are
    based on the parametric or nonparametric model.}
  \item{n.obs}{The number of units.}
  \item{agg.table}{Aggregate posterior estimates of the marginal means
    of \eqn{W_1} and \eqn{W_2}.}
  If \code{long = TRUE}, the following elements are also included:
  \item{W1.table}{Unit-level posterior estimates for \eqn{W_1}.}
  \item{W2.table}{Unit-level posterior estimates for \eqn{W_2}.}
  This object can be printed by \code{print.summary.eco}
}

\seealso{\code{eco}}

\author{
Kosuke Imai,
  Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}, \url{http://www.princeton.edu/~kimai};
  Ying Lu, Institute for Quantitative Social Sciences, 
  Harvard University \email{ylu@Latte.Harvard.Edu}
}

\keyword{methods}
