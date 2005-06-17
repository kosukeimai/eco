\name{predict.eco}

\alias{predict.eco}

\title{Out-of-Sample Posterior Prediction under the Parametric Bayesian
Ecological Inference Models}

\description{
  Obtains out-of-sample posterior predictions under a fitted parametric
  Bayesian model for ecological inference. \code{predict} method for
  class \code{eco}.
}

\usage{
  \method{predict}{eco}(object, newdraw = NULL, subset = NULL, verbose = FALSE, ...)
}

\arguments{
  \item{object}{An output object from \code{eco}.}
  \item{newdraw}{An optional list containing two matrices of MCMC draws
    of \eqn{\mu} and \eqn{\Sigma}. The two matrices should be named as
    \code{mu} and \code{Sigma}, respectively. The default is the
    original MCMC draws stored in \code{object}.
  }
  \item{subset}{A scalar or numerical vector specifying the row
    number(s) of \code{mu} and \code{Sigma} in the output object from
    \code{eco}. If specified, the posterior draws of parameters for
    those rows are used for posterior prediction. The default is
    \code{NULL} where all the posterior draws are used. 
  }
  \item{verbose}{logical. If \code{TRUE}, helpful messages along with a
    progress report on the Monte Carlo sampling from the posterior 
    predictive distributions are printed on the screen. The default is
    \code{FALSE}.
  }
  \item{...}{further arguments passed to or from other methods.}
}

\details{The posterior predictive values are computed using the
  Monte Carlo sample stored in the \code{eco} output (or other sample if
  \code{newdraw} is specified). Given each Monte Carlo sample of the
  parameters, we sample the vector-valued latent variable from the
  appropriate multivariate Normal distribution. Then, we apply the
  inverse logit transformation to obtain the predictive values of
  proportions, \eqn{W}. The computation may be slow if a large Monte
  Carlo sample of the model parameters is used. In either case, setting
  \code{verbose = TRUE} may be helpful in monitoring the progress of the
  code.
}

\value{
  \code{predict.eco} yields a matrix of class \code{predict.eco}
  containing the Monte Carlo sample from the posterior predictive
  distribution of inner cells of ecological tables. 
}

\seealso{\code{eco}, \code{predict.ecoNP}; eco home
page at \url{http://www.princeton.edu/~kimai/research/eco.html}}

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}; Ying Lu, Institute for Quantitative
  Social Sciences, Harvard University \email{ylu@Latte.Harvard.Edu}
}

\keyword{methods}
