\name{ecoML}

\alias{ecoML}

\title{Fitting Parametric Models and Quantifying Missing Information
 for Ecological Inference in 2x2 Tables}

\description{
  \code{ecoML} is used to fit parametric models for ecological 
  inference in \eqn{2 \times 2} tables via EM algorithms. It gives the point 
  estimates of the parameters for models based on different assumptions. 
  The standard errors of the point estimates are also computed via SEM 
  algorithms. Moreover, \code{ecoML} quantifies the amount of missing 
  information associated with each parameter and allows researcher to 
  examine the impact of missing information on parameter estimation in 
  ecological inference. The models and detailed algorithms 
  are described in Imai, Lu and Strauss(2006). 
}

\usage{
   ecoML(formula, data = parent.frame(), N=NULL, supplement = NULL, 
         theta.start = c(0,0,1,1,0), fix.rho = TRUE,
         context = FALSE, sem = TRUE, epsilon=10^(-10), 
	maxit = 1000, loglik = TRUE, verbose= TRUE)  
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
  \item{N}{An optional variable representing the size of the unit; e.g.,
    the total number of voters.}
  \item{supplement}{An optional matrix of supplemental data. The matrix
    has two columns, which contain additional individual-level data such
    as survey data for \eqn{W_1} and \eqn{W_2}, respectively.  If
    \code{NULL}, no additional individual-level data are included in the
    model. The default is \code{NULL}.
  }
  \item{fix.rho}{Logical. If \code{TRUE}, the correlation
  (when \code{context=TRUE}) or the partial correlation(when
   \code{context=FALSE}) between \eqn{W_1} and \eqn{W_2} 
    is fixed through the estimation. For details, see 
     Imai,Lu and Strauss(2006). The default is \code{TRUE}.
  } 
  \item{context}{Logical. If \code{TRUE}, the contextual effect is also
    modeled. See Imai and Lu (2005) for details. The default is
    \code{FALSE}. 
  }
  \item{sem}{Logical. If \code{TRUE}, the standard errors of parameter
   estimates are estimated via SEM alogirthm. The default is
    \code{TRUE}. 
  }
  \item{theta.start}{A numeric vector that spcifies the starting values
    for the mean, variance and covariance. When \code{context=FALSE},
    the elements of \code{theta.start} correspond to (\eqn{E(W_1), E(W_2), 
      V(W_1), V(W_2),\rho(W_1,W_2)}). When \code{context=TRUE}, the 
    elements of \code{theta.start} correspond to (\eqn{E(W_1), E(W_2), 
      V(W_1),V(W_2),\rho(W_1X),\rho(W_2X),\rho(W_1,W_2)}). Moreover, when
    \code{fix.rho=TRUE}, \eqn{\rho(W_1,W_2)} is set to be the correlation
    between \eqn{W_1} and \eqn{W_2} when \code{context=FALSE}, and 
    the partial correlation between \eqn{W_1} and \eqn{W_2} given \eqn{X}
    when \code{context=FALSE}. The default is \code{c(0,0,1,1,0)}.
  }
  
  \item{epsilon}{A positive number that specifies the convergence criterion
    for EM algorithm. The square root of \code{epsilon} is the convergence 
    criterion for SEM algorithm. The default is \code{10^(-10)}. 
  } 
  \item{maxit}{A positive integer specifies the maximum number of iterations
    before the convergence criterion is met. The default is \code{1000}.
  } 
  \item{loglik}{Logical. If \code{TRUE}, the loglikelihood at each iteration
    of EM is saved. The default is \code{TRUE}. 
  } 
  \item{verbose}{Logical. If \code{TRUE}, the progress of the Gibbs 
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
  \eqn{Y_i=X_i W_{1i}+(1-X_i)W_{2i}}.

 In particular, we consider two assumptions, CAR when 
\eqn{W_1} and \eqn{W_2} are independent from the 
contextual variable \eqn{X}, and NCAR when they are 
correlated with \eqn{X}.
}
\examples{

## load the census data
data(census)

## NOTE: convergence has not been properly assessed for the following
## examples. See Imai, Lu and Strauss (2006) for more complete analyses.

## fit the parametric model with the default prior specification
#res <- ecoML(Y ~ X, data = census, verbose = TRUE)
## summarize the results
##summary(res)  not run yet

## obtain out-of-sample prediction
##out <- predict(res, verbose = TRUE)
## summarize the results
##summary(out)

## load the Robinson's census data
data(census)

## fit the parametric model with contextual effects and some individual 
## level data
## using the default prior specification
surv<-1:600
res1 <- ecoML(Y ~ X, context = TRUE, data = census[-surv,], 
   	supplement=census[surv,c(4:5,1)], maxit=100,verbose = TRUE)
## summarize the results--not run yet
##summary(res1)

}

\value{
  An object of class \code{ecoML} containing the following elements:
  \item{call}{The matched call.}
  \item{X}{The row margin, \eqn{X}.}
  \item{Y}{The column margin, \eqn{Y}.}
  \item{N}{The size of each table, \eqn{N}.}
  \item{context}{The assumption under which model is estimated. If 
    \code{context=FALSE}, CAR assumption is adopted and no
    contextual effect is modelled. If \code{context=TRUE}, NCAR
    assumtion is adopted, and contextual effect is modelled.}
  \item{sem}{Whether SEM alogorithm is used to estimate the standard
    errors and observed information matrix for the parameter estimates.}
  \item{rho.fix}{Whether the correlation or the partial correlation between
    \eqn{W_1} an \eqn{W_2} is fixed in the estimation.}
  \item{r12}{If \code{rho.fix=TRUE}, the value that \eqn{\r12} is fixed to.}
  \item{epsilon}{The precision criterion for EM convergence. 
    \eqn{\sqrt{\epsilon}} is the precision criterion for SEM convergence.}
  \item{theta.sem}{The EM estimates of \eqn{E(W_1)},\eqn{E(W_2)},
    \eqn{V(W_1)},\eqn{V(W_2)}, and \eqn{COV(W_1,W_2)}. If 
    \code{context=TRUE}, \eqn{E(X)},\eqn{COV(W_1,X)}, 
    \eqn{COV(W_2,X)} are also reported.}
  \item{W}{Insample estimation of \eqn{W_1} and \eqn{W_2}.}
  \item{suff.stat}{The sufficient statistics for \code{theta.em}.}
  \item{iters.em}{Number of EM iterations before convergence is achieved.}
  \item{iters.sem}{Number of SEM iterations before convergence is achieved.}
  \item{loglik}{The loglihood of the model when convergence is achieved.}
  \item{loglik.log.em}{A vector saving the loglikelihood at each 
    iteration of EM.}
  \item{mu.log.em}{A matrix saving the mean estimation at each 
    iteration of EM.}
  \item{Sigma.log.em}{A matrix saving the variance-covariance estimation 
    at each iteration of EM.} 
  Moreover, when \code{sem=TRUE}, \code{ecoML} also output the following
  values:
  \item{Icom}{The (expected) complete data information matrix estimated 
    via SEM algorithm. When \code{context=FALSE, fix.rho=TRUE}, 
    \code{Icom} is 4 by 4. When \code{context=FALSE, fix.rho=FALSE}, 
    \code{Icom} is 5 by 5. When \code{context=TRUE}, \code{Icom} 
    is 9 by 9.}
  \item{Iobs}{The observed information matrix. The dimension of 
    \code{Iobs} is same as \code{Icom}.}
  \item{Imiss}{The difference between \code{Icom} and \code{Iobs}. 
    The dimension of \code{Imiss} is same as \code{miss}.}
  \item{Vobs}{The (symmetrized) variance-covariance matrix of the EM paramter
    estimates. The dimension of \code{Vobs} is same as 
    \code{Icom}.}
  \item{Vobs.original}{The estimated variance-covariance matrix of the 
    EM paramter  estimates. The dimension of \code{Vobs} is same as 
    \code{Icom}.}
  \item{Fmis}{The fraction of missing information associated with each 
    parameter estimation. }
  \item{Ieigen}{The largest eigen value of \code{Imiss}.}
  \item{DM}{The rate-of-converence matrix of the SEM algorithm. For details,
    see Meng and Rubin(1999).}
  \item{Icom.trans}{The complete data information matrix for the fisher 
    transformed parameters.}
  \item{Iobs.trans}{The observed data information matrix for the fisher 
    transformed parameters.}
  \item{Fmis.trans}{The fractions of missing information associated with 
    the fisher transformed parameters.}
}

\author{
  Kosuke Imai, Department of Politics, Princeton University,
  \email{kimai@Princeton.Edu}, \url{http://imai.princeton.edu};
  Ying Lu, Institute for Quantitative Social Sciences, 
  Harvard University \email{ylu@Latte.Harvard.Edu};
  Aaron Strauss, Department of Politics, Princeton University,
 \email{abstraus@Princeton.Edu}.
}

\references{
  Imai, Kosuke, Ying Lu and Aaron Strauss. (2006) \dQuote{Quantifying
  Missing Information in Ecological Inference.} Working Paper, Princeton
  University, available at 
 \url{http://imai.princeton.edu/research/eiinfo.html}
}

\seealso{\code{ecoNP}, \code{eco}}
\keyword{models}
