\name{eco}

\alias{eco}

\title{Fitting the Parametric Bayesian Model of Ecological Inference in
  2 by 2 Tables}  

\description{
  \code{eco} is used to fit the parametric Bayesian model (based on a
  Normal/Inverse-Wishart prior) for ecological inference in \eqn{2
  \times 2} tables via Markov chain Monte Carlo. It gives in-sample
  predictions as well as the estimates of the model parameters. The
  model and algorithm are described in Imai and Lu (2004).
}

\usage{
eco(formula, data = parent.frame(), supplement = NULL, mu0 = c(0,0),
    tau0 = 2, nu0 = 4, S0 = diag(10,2), parameter = TRUE,
    grid = FALSE, n.draws = 5000, burnin = 0, thin = 0, verbose = FALSE)
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
  \item{parameter}{Logical. If \code{TRUE}, the Gibbs draws of the population
    parameters such as mu and sigma are returned. The default is \code{TRUE}.
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

## fit the parametric model with the default prior specification
res <- eco(Y ~ X, data = reg, verbose = TRUE)
## summarize the results
summary(res)

## obtain out-of-sample prediction
out <- predict(res, verbose=TRUE)
## summarize the results
summary(out)

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
  \code{parameter = TRUE}.
  \item{mu}{The posterior draws of the population mean parameter,
    \eqn{mu}.}
  \item{Sigma}{The poterior draws of the population variance matrix,
    \eqn{Sigma}.}
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

\seealso{\code{ecoNP}, \code{predict.eco}, \code{summary.eco}}
\keyword{models}
