\name{eco}

\alias{eco}

\title{Fitting the Parametric and Nonparametric Bayesian Models of
Ecological Inference in 2 by 2 Tables} 

\description{
  \code{eco} is used to fit the parametric and nonparametric Bayesian
  models for ecological inference in \eqn{2 \times 2} tables via Markov
  chain Monte Carlo. It gives in-sample predictions as well as
  out-of-sample predictions for population inference. The parametric
  model uses a normal/inverse-Wishart prior, while the nonparametric
  model uses a Dirichlet process prior. The models and algorithms are
  described in Imai and Lu (2004).
}

\usage{
eco(formula, data = parent.frame(), nonpar = FALSE, supplement = NULL,
    mu0 = c(0,0), tau0 = 2, nu0 = 4, S0 = diag(10,2), alpha = NULL,
    a0 = 1, b0 = 0.1, predict = FALSE, parameter = FALSE, grid = FALSE,
    n.draws = 5000, burnin = 0, thin = 0, verbose = FALSE)
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
    \code{eco} is called. 
  }  
  \item{nonpar}{Logical. If \code{TRUE}, the nonparametric
    model will be fit. Otherwise, the parametric model will be
    estimated. The default is \code{TRUE}.
  }
  \item{supplement}{A numeric matrix. The matrix has two columns, which
    contain additional individual-level data such as survey data for
    \eqn{W_1} and \eqn{W_2}, respectively.  If \code{NULL}, no
    additional individual-level data are included in the model. The
    default is \code{NULL}.
  }      
  \item{mu0}{A \eqn{2 \times 1} numeric vector. The prior mean. The
    default is \eqn{(0,0)}.
  } 
  \item{tau0}{A positive integer. The prior scale parameter. The default
    is \code{2}. 
  } 
  \item{nu0}{A positive integer. The prior degrees of freedom
    parameter. the default is \code{4}.
  } 
  \item{S0}{A \eqn{2 \times 2} numeric matrix, representing a positive
    definite prior scale matrix. The default is \code{diag(10,2)}. 
  } 
  \item{alpha}{A positive scalar. If \code{NULL}, the concentration
    parameter \eqn{\alpha} will be updated at each Gibbs draw. The prior
    parameters \code{a0} and \code{b0} need to be specified. Otherwise,
    \eqn{\alpha} is fixed at a user specified value. 
    The default is \code{NULL}.
  } 
  \item{a0}{A positive integer. The shape parameter of the gamma prior
    for \eqn{\alpha}. The default is \code{1}.
  } 
  \item{b0}{A positive integer. The scale parameter of the gamma prior
    for \eqn{\alpha}. The default is \code{0.1}.
  } 
  \item{predict}{Logical. If \code{TRUE}, out-of sample predictions will
    be returned. The default is \code{FALSE}.
  } 
  \item{parameter}{Logical. If \code{TRUE}, the Gibbs draws of the population
    parameters such as mu and sigma are returned. The default is \code{FALSE}.
  }
  \item{grid}{Logical. If \code{TRUE}, the grid method is used to sample
    \eqn{W} in the Gibbs sampler. If \code{FALSE}, the Metropolis
    algorithm is used where candidate draws are sampled from the uniform
    distribution on the tomography line for each unit. Note that the
    grid method is significantly slower than the Metropolis algorithm.
  }
  \item{n.draws}{A positive integer. The number of MCMC draws.
    The default is \code{5000}.
  }
  \item{burnin}{A positive integer. The burnin interval for the Markov
    chain; i.e. the number of initial draws that should not be stored. The
    default is \code{0}.
  }   
  \item{thin}{A positive integer. The thinning interval for the
    Markov chain; i.e. the number of Gibbs draws between the recorded
    values that are skipped. The default is \code{0}.
  } 
  \item{verbose}{Logical. If \code{TRUE}, the progress of the gibbs 
   sampler is printed to the screen. The default is \code{FALSE}.
  }
}

\details{
  An example of \eqn{2 \times 2} ecological table for racial voting is
  given below: 
  \tabular{lccc}{
    \tab black voters  \tab white voters \tab \cr
    Voted \tab \eqn{W_{1i}}  \tab \eqn{W_{2i}} \tab \eqn{Y_i} \cr
    Not voted \tab \eqn{1-W_{1i}}  \tab \eqn{1-W_{2i}} \tab \eqn{1-Y_i}  \cr
    \tab \eqn{X_i} \tab \eqn{1-X_i} \tab 
  }
  where \eqn{Y_i} and \eqn{X_i} represent the observed margins, and
  \eqn{W_1} and \eqn{W_2} are unknown variables. All variables are
  proportions and hence bounded between 0 and 1. For each \eqn{i}, the
  following deterministic relationship holds,
  \eqn{Y_i=X W_{1i}+(1-X_i)W_{2i}}.
}

\examples{

## load the registration data
data(reg)

## NOTE: convergence has not been properly assessed for the following
## examples.

## fit the parametric model to give in-sample predictions and store
## parameter estimates
res <- eco(Y ~ X, data = reg, parameter = TRUE, verbose = TRUE) 
##summarize the results
summary(res)

## fit the nonparametric model to give in-sample predictions
res1 <- eco(Y ~ X, data = reg, nonpar = TRUE, n.draws = 500, verbose = TRUE) 
##summarize the results
summary(res1)
}

\value{
  An object of class \code{eco} containing the following elements:
  \item{call}{The matched call.}
  \item{nonpar}{The logical variable indicating whether the
    nonparametric model is fit.}
  \item{X}{The row margin, \eqn{X}.}
  \item{Y}{The column margin, \eqn{Y}.}
  \item{nu0}{The prior degrees of freedom.}
  \item{tau0}{The prior scale parameter.}
  \item{mu0}{The prior mean.}
  \item{S0}{The prior scale matrix.}
  \item{burnin}{The number of initial burnin draws.}
  \item{thin}{Thinning interval.}
  \item{W1}{The posterior in-sample predictions of \eqn{W_1}.}
  \item{W2}{The posterior in-sample predictions of \eqn{W_2}.}
  The following additional elements are included in the output when
  \code{predict = TRUE}:
  \item{W1.pred}{The posterior predictive draws or out-of-sample
    predictions of \eqn{W_1}. Export only if \code{predict=TRUE}.
  }
  \item{W2.pred}{The posterior predictive draws or out-of-sample
    predictions of \eqn{W_2}. Export only if \code{predict=TRUE}.
  }
  The following additional elements are included in the output when the
  nonparametric model is fit (i.e., \code{nonpar = TRUE})
  \item{a0}{The prior shape parameter.}
  \item{b0}{The prior scale parameter.}
  The following additional elements are included in the output when
  \code{parameter = TRUE}.\n
  For the parametric model:
  \item{mu}{The posterior draws of the population mean parameter,
    \eqn{mu}.}
  \item{Sigma}{The poterior draws of the population variance matrix,
    \eqn{Sigma}.}
  For the nonparametric model:
  \item{mu1}{The posterior draws of the population mean parameter
    for \eqn{W_1}. It is an \eqn{m \times n} matrix, where \eqn{m} is
    the number of Gibbs draws saved, \eqn{n} is the number of units.}
  \item{mu2}{The posterior draws of the population mean parameter
    of \eqn{W_2}. The dimension of \code{mu2} is the same as \code{mu1}.}
  \item{Sigma11}{The posterior draws of the population variance
    parameter for \eqn{W_1}. It is an \eqn{m \times n} matrix, where \eqn{m} is the
    number of Gibbs draws saved, \eqn{n} is the number of units.}
  \item{Sigma12}{The posterior draws of the population covariance
    parameter for \eqn{W_1} and \eqn{W_2}. The dimension of
    \code{Sigma12} is the same as \code{Sigma11}.}
  \item{Sigma22}{The posterior draws of the population variance
    parameter for \eqn{W_2}. The dimension of \code{Sigma22} is same
    as \code{Sigma11}.}
  \item{alpha}{The posterior draws of \eqn{\alpha}.}
  \item{nstar}{The number of clusters at each Gibbs draw.}
}

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}, \url{http://www.princeton.edu/~kimai};
  Ying Lu, Institute for Quantitative Social Sciences, 
  Harvard University \email{ylu@Latte.Harvard.Edu}}

\references{
  Imai, Kosuke and Ying Lu. (2004) \dQuote{ Parametric and Nonparametric
  Bayesian Models for Ecological Inference in \eqn{2 \times 2}
  Tables.} Proceedings of the American Statistical Association.
  \url{http://www.princeton.edu/~kimai/research/einonpar.html} }

\seealso{\code{summary.eco}}
\keyword{models}
