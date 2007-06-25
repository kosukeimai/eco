\name{eco}

\alias{eco}

\title{Fitting the Parametric Bayesian Model of Ecological Inference in
  2x2 Tables}

\description{
  \code{eco} is used to fit the parametric Bayesian model (based on a
  Normal/Inverse-Wishart prior) for ecological inference in \eqn{2
  \times 2} tables via Markov chain Monte Carlo. It gives the in-sample
  predictions as well as the estimates of the model parameters. The
  model and algorithm are described in Imai, Lu and Strauss (Forthcoming 
  a and b). 
  }

\usage{
eco(formula, data = parent.frame(), N = NULL, supplement = NULL, 
    context = FALSE, mu0 = 0, tau0 = 2, nu0 = 4, S0 = 10,
    mu.start = 0, Sigma.start = 10, parameter = TRUE,
    grid = FALSE, n.draws = 5000, burnin = 0, thin = 0, 
    verbose = FALSE)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit,
    specifying the column and row margins of \eqn{2 \times
    2} ecological tables. \code{Y ~ X} specifies \code{Y} as the
    column margin  (e.g., turnout) and \code{X} as the row margin 
    (e.g., percent African-American). Details and specific
    examples are given below.
  } 
  \item{data}{An optional data frame in which to interpret the variables
    in \code{formula}. The default is the environment in which
    \code{eco} is called. 
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
   (Forthcoming b) for details. The default is \code{FALSE}. 
  }
  \item{mu0}{A scalar or a numeric vector that specifies the prior mean 
    for the mean parameter \eqn{\mu} for \eqn{(W_1,W_2)} (or for 
    \eqn{(W_1, W_2, X)} if \code{context=TRUE}). When the input of 
    \code{mu0} is a scalar, its value will be repeated to yield a vector of 
     the length of \eqn{\mu}, otherwise, it needs to be a vector of same 
    length as \eqn{\mu}. When \code{context=TRUE}, the length of \eqn{\mu} 
    is 3, otherwise it is 2. The default is \code{0}. 
  } 
  \item{tau0}{A positive integer representing the scale parameter of the 
  Normal-Inverse Wishart prior for the mean and variance parameter
  \eqn{(\mu, \Sigma)}. The default is \code{2}.} 
  \item{nu0}{A positive integer representing the prior degrees of
  freedom of the Normal-Inverse Wishart prior for the mean and 
  variance parameter \eqn{(\mu, \Sigma)}. The default is \code{4}.
  } 
  \item{S0}{A positive scalar or a positive definite matrix that specifies 
  the prior scale matrix of the Normal-Inverse Wishart prior for the 
  mean and variance parameter \eqn{(\mu, \Sigma)} . If it is 
  a scalar, then the prior scale matrix will be a diagonal matrix with 
  the same dimensions as \eqn{\Sigma} and the diagonal elements all take 
  value of \code{S0}, otherwise \code{S0} needs to have same dimensions as 
  \eqn{\Sigma}. When \code{context=TRUE}, \eqn{\Sigma} is a 
  \eqn{3 \times 3} matrix, otherwise, it is \eqn{2 \times 2}. 
   The default is \code{10}.
  } 
  \item{mu.start}{A scalar or a numeric vector that specifies the 
    starting values of the mean parameter \eqn{\mu}. 
    If it is a scalar, then its value will be repeated to 
    yield a vector of the length of \eqn{\mu}, otherwise,
    it needs to be a vector of same length as \eqn{\mu}.
    When \code{context=FALSE}, the length of \eqn{\mu} is 2, 
    otherwise it is 3. The default is \code{0}. 
  }
  \item{Sigma.start}{A scalar or a positive definite matrix
   that specified the starting value of the variance matrix 
   \eqn{\Sigma}. If it is  a scalar, then the prior scale 
  matrix will be a diagonal matrix with the same dimensions 
  as \eqn{\Sigma} and the diagonal elements all take value 
  of \code{S0}, otherwise \code{S0} needs to have same dimensions as 
  \eqn{\Sigma}. When \code{context=TRUE}, \eqn{\Sigma} is a 
  \eqn{3 \times 3} matrix, otherwise, it is \eqn{2 \times 2}. 
   The default is \code{10}.
  } 
  \item{parameter}{Logical. If \code{TRUE}, the Gibbs draws of the population
    parameters, \eqn{\mu} and \eqn{\Sigma}, are returned in addition to
    the in-sample predictions of the missing internal cells,
    \eqn{W}. The default is \code{TRUE}.
  }
  \item{grid}{Logical. If \code{TRUE}, the grid method is used to sample
    \eqn{W} in the Gibbs sampler. If \code{FALSE}, the Metropolis
    algorithm is used where candidate draws are sampled from the uniform
    distribution on the tomography line for each unit. Note that the
    grid method is significantly slower than the Metropolis algorithm.
    The default is \code{FALSE}.
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

\details{
  An example of \eqn{2 \times 2} ecological table for racial voting is
  given below: 
  \tabular{llccc}{
    \tab \tab black voters  \tab  white voters \tab \cr
    \tab vote \tab \eqn{W_{1i}}  \tab \eqn{W_{2i}} \tab \eqn{Y_i} \cr
    \tab not vote \tab \eqn{1-W_{1i}}  \tab \eqn{1-W_{2i}} \tab \eqn{1-Y_i}  \cr
    \tab \tab \eqn{X_i} \tab \eqn{1-X_i} \tab 
  }
  where \eqn{Y_i} and \eqn{X_i} represent the observed margins, and
  \eqn{W_1} and \eqn{W_2} are unknown variables. In this exmaple, 
  \eqn{Y_i} is the turnout rate in the ith precint, \eqn{X_i} is the
  proproption of African American in the ith precinct. The unknowns 
  \eqn{W_{1i}} an d\eqn{W_{2i}} are the black and white turnout, 
  respectively. All variables are proportions and hence bounded 
  between 0 and 1. For each \eqn{i}, the following deterministic 
  relationship holds, \eqn{Y_i=X_i W_{1i}+(1-X_i)W_{2i}}.
}

\examples{

## load the registration data
data(reg)

## NOTE: convergence has not been properly assessed for the following
## examples. See Imai, Lu and Strauss (Forthcoming a and b) for more
## complete analyses.

## fit the parametric model with the default prior specification
res <- eco(Y ~ X, data = reg, verbose = TRUE)
## summarize the results
summary(res)

## obtain out-of-sample prediction
out <- predict(res, verbose = TRUE)
## summarize the results
summary(out)

## load the Robinson's census data
data(census)

## fit the parametric model with contextual effects and N 
## using the default prior specification
res1 <- eco(Y ~ X, N = N, context = TRUE, data = census, verbose = TRUE)
## summarize the results
summary(res1)

## obtain out-of-sample prediction
out1 <- predict(res1, verbose = TRUE)
## summarize the results
summary(out1)

}

\value{
  An object of class \code{eco} containing the following elements:
  \item{call}{The matched call.}
  \item{X}{The row margin, \eqn{X}.}
  \item{Y}{The column margin, \eqn{Y}.}
  \item{N}{The size of each table, \eqn{N}.}
  \item{burnin}{The number of initial burnin draws.}
  \item{thin}{The thinning interval.}
  \item{nu0}{The prior degrees of freedom.}
  \item{tau0}{The prior scale parameter.}
  \item{mu0}{The prior mean.}
  \item{S0}{The prior scale matrix.}
  \item{W}{A three dimensional array storing the posterior in-sample
  predictions of \eqn{W}. The first dimension indexes the Monte Carlo
  draws, the second dimension indexes the columns of the table, and the
  third dimension represents the observations.}
  \item{Wmin}{A numeric matrix storing the lower bounds of \eqn{W}.}
  \item{Wmax}{A numeric matrix storing the upper bounds of \eqn{W}.}
  The following additional elements are included in the output when
  \code{parameter = TRUE}.
  \item{mu}{The posterior draws of the population mean parameter,
    \eqn{\mu}.}
  \item{Sigma}{The posterior draws of the population variance matrix,
    \eqn{\Sigma}.}
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

  Imai, Kosuke, Ying Lu and Aaron Strauss. (Forthcoming). \dQuote{Bayesian and
    Likelihood Inference for 2 x 2 Ecological Tables: An Incomplete Data
    Approach} Political Analysis, available at 
  \url{http://imai.princeton.edu/research/eiall.html}
}

\seealso{\code{ecoML}, \code{ecoNP}, \code{predict.eco}, \code{summary.eco}}
\keyword{models}
