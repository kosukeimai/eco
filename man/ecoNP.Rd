\name{ecoNP}

\alias{ecoNP}

\title{Fitting the Nonparametric Bayesian Models of Ecological Inference
  in 2x2 Tables}  

\description{
  \code{ecoNP} is used to fit the nonparametric Bayesian model (based on
  a Dirichlet process prior) for ecological inference in \eqn{2 \times
  2} tables via Markov chain Monte Carlo. It gives the in-sample
  predictions as well as out-of-sample predictions for population
  inference.  The models and algorithms are described in Imai, Lu and
  Strauss (2008, Forthcoming). 
}

\usage{
ecoNP(formula, data = parent.frame(), N = NULL, supplement = NULL,
      context = FALSE, mu0 = 0, tau0 = 2, nu0 = 4, S0 = 10, 
      alpha = NULL, a0 = 1, b0 = 0.1, parameter = FALSE, 
      grid = FALSE, n.draws = 5000, burnin = 0, thin = 0, 
      verbose = FALSE)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit,
    specifying the column and row margins of \eqn{2 \times
    2} ecological tables. \code{Y ~ X} specifies \code{Y} as the
    column margin (e.g., turnout) and \code{X} as the row margin
   (e.g., percent African-American). Details and specific examples 
   are given below.
  } 
  \item{data}{An optional data frame in which to interpret the variables
    in \code{formula}. The default is the environment in which
    \code{ecoNP} is called. 
  }  
  \item{N}{An optional variable representing the size of the unit; e.g.,
    the total number of voters.}
  \item{supplement}{An optional matrix of supplemental data. The matrix
    has two columns, which contain additional individual-level data such
    as survey data for \eqn{W_1} and \eqn{W_2}, respectively.  If
    \code{NULL}, no additional individual-level data are included in the
    model. The default is \code{NULL}.
  }      
  \item{context}{Logical. If \code{TRUE}, the contextual effect is also
    modeled, that is to assume the row margin \eqn{X} and the unknown 
   \eqn{W_1} and \eqn{W_2} are correlated. See Imai, Lu and Strauss
   (2008, Forthcoming) for details. The default is \code{FALSE}. 
  }
  \item{mu0}{A scalar or a numeric vector that specifies the prior mean 
    for the mean parameter \eqn{\mu} of the base prior distribution \eqn{G_0}
    (see Imai, Lu and Strauss (2008, Forthcoming) for detailed 
    descriptions of Dirichlete prior and the normal base prior distribution) . 
    If it is a scalar, then its value will be repeated to yield a vector 
    of the length of \eqn{\mu}, otherwise,
    it needs to be a vector of same length as \eqn{\mu}.
    When \code{context=TRUE }, the length of \eqn{\mu} is 3, 
    otherwise it is 2. The default is \code{0}. 
  } 
  \item{tau0}{A positive integer representing the scale parameter of the 
  Normal-Inverse Wishart prior for the mean and variance parameter
  \eqn{(\mu_i, \Sigma_i)} of each observation. The default is \code{2}.} 

  \item{nu0}{A positive integer representing the prior degrees of
  freedom of the variance matrix \eqn{\Sigma_i}. the default is \code{4}.
  } 
  \item{S0}{A positive scalar or a positive definite matrix that specifies 
  the prior scale matrix for the variance matrix \eqn{\Sigma_i}. If it is 
  a scalar, then the prior scale matrix will be a diagonal matrix with 
  the same dimensions as \eqn{\Sigma_i} and the diagonal elements all
  take value of \code{S0}, otherwise \code{S0} needs to have same 
  dimensions as \eqn{\Sigma_i}. When \code{context=TRUE}, \eqn{\Sigma} is a 
  \eqn{3 \times 3} matrix, otherwise, it is \eqn{2 \times 2}. 
   The default is \code{10}.
  } 
  \item{alpha}{A positive scalar representing a user-specified fixed
    value of the concentration parameter, \eqn{\alpha}. If \code{NULL},
    \eqn{\alpha} will be updated at each Gibbs draw, and its prior
    parameters \code{a0} and \code{b0} need to be specified. The default
    is \code{NULL}. 
  } 
  \item{a0}{A positive integer representing the value of shape parameter
    of the gamma prior distribution for \eqn{\alpha}. The default is \code{1}.
  } 
  \item{b0}{A positive integer representing the value of the scale
    parameter of the gamma prior distribution for \eqn{\alpha}. The
    default is \code{0.1}.
  }
  \item{parameter}{Logical. If \code{TRUE}, the Gibbs draws of the population
    parameters, \eqn{\mu} and \eqn{\Sigma}, are returned in addition to
    the in-sample predictions of the missing internal cells,
    \eqn{W}. The default is \code{FALSE}. This needs to be set to
    \code{TRUE} if one wishes to make population inferences through
    \code{predict.eco}. See an example below.
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
  \item{verbose}{Logical. If \code{TRUE}, the progress of the Gibbs 
   sampler is printed to the screen. The default is \code{FALSE}.
  }
}

\examples{

## load the registration data
data(reg)

## NOTE: We set the number of MCMC draws to be a very small number in
## the following examples; i.e., convergence has not been properly
## assessed. See Imai, Lu and Strauss (2006) for more complete examples.

## fit the nonparametric model to give in-sample predictions
## store the parameters to make population inference later
res <- ecoNP(Y ~ X, data = reg, n.draws = 50, param = TRUE, verbose = TRUE) 

##summarize the results
summary(res)

## obtain out-of-sample prediction
out <- predict(res, verbose = TRUE)

## summarize the results
summary(out)

## density plots of the out-of-sample predictions
par(mfrow=c(2,1))
plot(density(out[,1]), main = "W1")
plot(density(out[,2]), main = "W2")


## load the Robinson's census data
data(census)

## fit the parametric model with contextual effects and N 
## using the default prior specification

res1 <- ecoNP(Y ~ X, N = N, context = TRUE, param = TRUE, data = census, 
n.draws = 25, verbose = TRUE)

## summarize the results
summary(res1)

## out-of sample prediction 
pres1 <- predict(res1)
summary(pres1)
}

\value{
  An object of class \code{ecoNP} containing the following elements:
  \item{call}{The matched call.}
  \item{X}{The row margin, \eqn{X}.}
  \item{Y}{The column margin, \eqn{Y}.}
  \item{burnin}{The number of initial burnin draws.}
  \item{thin}{The thinning interval.}
  \item{nu0}{The prior degrees of freedom.}
  \item{tau0}{The prior scale parameter.}
  \item{mu0}{The prior mean.}
  \item{S0}{The prior scale matrix.}
  \item{a0}{The prior shape parameter.}
  \item{b0}{The prior scale parameter.}
  \item{W}{A three dimensional array storing the posterior in-sample
  predictions of \eqn{W}. The first dimension indexes the Monte Carlo
  draws, the second dimension indexes the columns of the table, and the
  third dimension represents the observations.}
  \item{Wmin}{A numeric matrix storing the lower bounds of \eqn{W}.}
  \item{Wmax}{A numeric matrix storing the upper bounds of \eqn{W}.}
  The following additional elements are included in the output when
  \code{parameter = TRUE}.
  \item{mu}{A three dimensional array storing the posterior draws of the
  population mean parameter, \eqn{\mu}. The first dimension indexes the
  Monte Carlo draws, the second dimension indexes the columns of the
  table, and the third dimension represents the observations.} 
  \item{Sigma}{A three dimensional array storing the posterior draws of
  the population variance matrix, \eqn{\Sigma}. The first dimension
  indexes the Monte Carlo draws, the second dimension indexes the
  parameters, and the third dimension represents the observations. }
  \item{alpha}{The posterior draws of \eqn{\alpha}.}
  \item{nstar}{The number of clusters at each Gibbs draw.}
}

\author{
  Kosuke Imai, Department of Politics, Princeton University,
  \email{kimai@Princeton.Edu}, \url{http://imai.princeton.edu};
  Ying Lu, Department of Sociology, University of Colorado at Boulder, 
  \email{ying.lu@Colorado.Edu}
}

\references{
  Imai, Kosuke, Ying Lu and Aaron Strauss. (Forthcoming).
  \dQuote{eco: R Package for Ecological Inference in 2x2 Tables}
  Journal of Statistical Software, available at
  \url{http://imai.princeton.edu/research/eco.html}

  Imai, Kosuke, Ying Lu and Aaron Strauss. (2008).
  \dQuote{Bayesian and Likelihood Inference for 2 x 2 Ecological Tables:
    An Incomplete Data Approach} Political Analysis, Vol. 16, No. 1 (Winter), pp. 41-69. available at
 \url{http://imai.princeton.edu/research/eiall.html}
}

\seealso{\code{eco}, \code{ecoML}, \code{predict.eco}, \code{summary.ecoNP}}
\keyword{models}
