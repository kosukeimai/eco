\name{eco}

\alias{ECO}
\alias{eco}

\title{Fitting the Parametric and Nonparametric Bayesian Models of
Ecological Inference in 2 x 2 Tables} 

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
eco <- function(Y, X, data = parent.frame(), n.draws = 5000, burnin = 0,
                thin = 5, verbose = FALSE, nonpar =  TRUE, nu0 = 4, tau0 = 1,
                mu0 = c(0,0), S0 = diag(8,2), supplement=NULL, alpha = NULL,
                a0 = 1, b0 = 0.1, predict = FALSE, parameter = FALSE) 
}

\arguments{
  \item{Y}{A numeric vector of proportions, representing the
    weighted average of the missing internal cells of a \eqn{2 \times  
      2} ecological table.
  } 
  \item{X}{A numeric vector of proportions, representing the weights.
  }
  \item{data}{An optional data frame in which to interpret the variables
    in \code{Y} and \code{X}. The default is the environment in which
    \code{eco} is called. 
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
    values that are skipped. The default is \code{5}.
  } 
  \item{verbose}{Logical. If \code{TRUE}, the progress of the gibbs 
   sampler is printed to the screen. The default is \code{FALSE}.
  }
  \item{nonpar}{Logical. If \code{TRUE}, the nonparametric
    model will be fit. Otherwise, the parametric model will be
    estimated. The default is \code{TRUE}.
  } 
  \item{nu0}{A positive integer. The prior degrees of freedom
    parameter. the default is \code{4}.
  } 
  \item{tau0}{A positive integer. The prior scale parameter. The default
    is \code{2}. 
  } 
  \item{mu0}{A \eqn{2 \times 1} numeric vector. The prior mean. The
    default is \eqn{(0,0)}.
  } 
  \item{S0}{A \eqn{2 \times 2} numeric matrix, representing a positive
    definite prior scale matrix. The default is \code{diag(10,2)}. 
  } 
  \item{supplement}{A numeric matrix. The matrix has two columns, which
    contain additional individual-level data such as survey data for
    \eqn{W_1} and \eqn{W_2}, respectively.  If \code{NULL}, no
    additional individual-level data are included in the model. The
    default is \code{NULL}.
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
}

\details{
  An example of \eqn{2 \times 2} ecological table for racial voting is
  as following: 
  \tabular{lccc}{
    \tab black voters  \tab white voters \tab \cr
    Voted \tab \eqn{W_{1i}}  \tab \eqn{W_{2i}} \tab \eqn{Y_i} \cr
    Not voted \tab \eqn{1-W_{1i}}  \tab \eqn{1-W_{2i}} \tab \eqn{1-Y_i}  \cr
    \tab \eqn{X_i} \tab \eqn{1-X_i} \tab 
  }
  where \eqn{Y_i} and \eqn{X_i} represent the observed margins, and
  \eqn{W_1} and \eqn{W_2} are unknown variables. The following
  deterministic relationship holds for each \eqn{i}:
  \eqn{Y_i=X W_{1i}+(1-X_i)W_{2i}}
}

\examples{

## load the registration data
data(reg)

## run the nonparametric model to give in-sample & out-of sample predictions
res <- eco(Y = Y, X = X, data = reg, n.draws = 50, verbose = TRUE) 

##summarize the results
summary(res)
}

\value{
  An object of class \code{eco} containing the following elements:
  \item{model}{The name of the model is used to produce the
    predictions. If the nonparametric model is used, \code{
      model=``Dirichlet Process prior''}; if the parametric model is used,
    \code{ model=``Normal prior''}.
  } 
  \item{X}{The vector of data \eqn{X}.}
  \item{Y}{The vector of data \eqn{Y}.}
  \item{nu0}{The prior degrees of freedom.}
  \item{tau0}{The prior scale parameter.}
  \item{mu0}{The prior means.}
  \item{S0}{The prior scale matrix.}
  \item{burnin}{The number of initial burnin draws.}
  \item{thin}{Thinning interval.}
  \item{mu1.post}{The posterior draws of the population mean parameter
    of \eqn{W_1}. In nonparametric model, \code{mu1.post} is a m by n
    matrix, where m is the number of Gibbs draws saved, n is the number of
    observations. In parametric model, \code{mu1.post} is a vector of length
    m. Export only if \code{parameter=TRUE}. }
  \item{mu2.post}{The posterior draws of the population mean parameter
    of \eqn{W_2}. The dimension of \code{mu2.post} is same as
    \code{mu1.post}. Export only if \code{parameter=TRUE}.}
  \item{Sigma11.post}{The posterior draws of the population variance
    parameter of \eqn{W_1}. In nonparametric model, \code{Sigma11.post}
    is a m by n matrix, where m is the number of Gibbs draws saved, n is the
    number of observations. In parametric model, \code{Sigma11.post} is a
    vector of length m. Export only if \code{parameter=TRUE}.}
  \item{Sigma12.post}{The posterior draws of the population covariance
    parameter between \eqn{W_1} and \eqn{W_2}. The dimension of
    \code{Sigma12.post} is same as \code{Sigma11.post}. Export only if
    \code{parameter=TRUE}.}
  \item{Sigma22.post}{The posterior draws of the population variance
    parameter of \eqn{W_2}. The dimension of \code{Sigma22.post} is same
    as \code{Sigma11.post}. Export only if \code{parameter=TRUE}.}
  \item{W1.post}{The posterior draws or in-sample predictions of
    \eqn{W_1}.}
  \item{W2.post}{The posterior draws or in-sample predictions of
    \eqn{W_2}.
  }
  \item{W1.pred}{The posterior predictive draws or out-of-sample
    predictions of \eqn{W_1}. Export only if \code{predict=TRUE}.
  }
  \item{W2.pred}{The posterior predictive draws or out-of-sample
    predictions of \eqn{W_2}. Export only if \code{predict=TRUE}.
  }
  The following additional elements are included in the output when the
  nonparametric model is fit:
  \item{alpha}{Whether \eqn{\alpha} is being updated at each Gibbs
    draw.}
  \item{a0}{The prior shape parameter.}
  \item{b0}{The prior scale parameter.}
  \item{a.post}{The Gibbs draws of \eqn{\alpha}.}
  \item{nstar}{The number of clusters at each Gibbs draw.}
}

\author{
  Ying Lu, Woodrow Wilson School of International and Public Affairs,
  Princeton University \email{yinglu@Princeton.Edu};
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}, \url{http://www.princeton.edu/~kimai}}

\references{
  Imai, Kosuke and Ying Lu. (2004) \dQuote{ Parametric and Nonparametric
  Bayesian Models for Ecological Inference in \eqn{2 \times 2}
  Tables.} Proceedings of the American Statistical Association.
  \url{http://www.princeton.edu/~kimai/research/einonpar.html} }

\seealso{\code{summary.eco}}
\keyword{models}
