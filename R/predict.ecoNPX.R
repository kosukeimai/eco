#' Out-of-Sample Posterior Prediction under the Nonparametric Bayesian Model
#' for Ecological Inference in 2x2 Tables
#' 
#' Obtains out-of-sample posterior predictions under the fitted nonparametric
#' Bayesian model for ecological inference. \code{predict} method for class
#' \code{ecoNP} and \code{ecoNPX}.
#' 
#' The posterior predictive values are computed using the Monte Carlo sample
#' stored in the \code{eco} or \code{ecoNP} output (or other sample if
#' \code{newdraw} is specified). Given each Monte Carlo sample of the
#' parameters, we sample the vector-valued latent variable from the appropriate
#' multivariate Normal distribution. Then, we apply the inverse logit
#' transformation to obtain the predictive values of proportions, \eqn{W}. The
#' computation may be slow (especially for the nonparametric model) if a large
#' Monte Carlo sample of the model parameters is used. In either case, setting
#' \code{verbose = TRUE} may be helpful in monitoring the progress of the code.
#' 
#' @aliases predict.ecoNPX
#' @param object An output object from \code{ecoNP}.
#' @param newdraw An optional list containing two matrices (or three
#' dimensional arrays for the nonparametric model) of MCMC draws of \eqn{\mu}
#' and \eqn{\Sigma}. Those elements should be named as \code{mu} and
#' \code{Sigma}, respectively. The default is the original MCMC draws stored in
#' \code{object}.
#' @param subset A scalar or numerical vector specifying the row number(s) of
#' \code{mu} and \code{Sigma} in the output object from \code{eco}. If
#' specified, the posterior draws of parameters for those rows are used for
#' posterior prediction. The default is \code{NULL} where all the posterior
#' draws are used.
#' @param obs An integer or vector of integers specifying the observation
#' number(s) whose posterior draws will be used for predictions. The default is
#' \code{NULL} where all the observations in the data set are selected.
#' @param cond logical. If \code{TRUE}, then the conditional prediction will
#' made for the parametric model with contextual effects. The default is
#' \code{FALSE}.
#' @param verbose logical. If \code{TRUE}, helpful messages along with a
#' progress report on the Monte Carlo sampling from the posterior predictive
#' distributions are printed on the screen. The default is \code{FALSE}.
#' @param ... further arguments passed to or from other methods.
#' @return \code{predict.eco} yields a matrix of class \code{predict.eco}
#' containing the Monte Carlo sample from the posterior predictive distribution
#' of inner cells of ecological tables. \code{summary.predict.eco} will
#' summarize the output, and \code{print.summary.predict.eco} will print the
#' summary.
#' @author Kosuke Imai, Department of Politics, Princeton University,
#' \email{kimai@@Princeton.Edu}, \url{http://imai.princeton.edu}; Ying Lu,
#' Center for Promoting Research Involving Innovative Statistical Methodology
#' (PRIISM), New York University \email{ying.lu@@nyu.Edu}
#' @seealso \code{eco}, \code{ecoNP}, \code{summary.eco}, \code{summary.ecoNP}
#' @keywords methods
predict.ecoNPX <- function(object, newdraw = NULL, subset = NULL,
                           obs = NULL, cond = FALSE, verbose = FALSE, ...){

  if (is.null(newdraw) && is.null(object$mu))
    stop("Posterior draws of mu and Sigma must be supplied")
  else if (!is.null(newdraw)){
    if (is.null(newdraw$mu) && is.null(newdraw$Sigma))
      stop("Posterior draws of both mu and Sigma must be supplied.")
    object <- newdraw
  }

  n.draws <- dim(object$mu)[1]
  n <- dim(object$mu)[3]
  mu <- aperm(coef(object, subset = subset, obs = obs), c(2,3,1))
  
  if (is.null(subset))
    subset <- 1:n.draws
  if (is.null(obs))
    obs <- 1:n
  Sigma <- aperm(object$Sigma[subset,,obs], c(2,3,1))

  if (cond) { # conditional prediction
    X <- object$X
    res <- .C("preDPX", as.double(mu), as.double(Sigma), as.double(X),
              as.integer(n), as.integer(n.draws), as.integer(2),
              as.integer(verbose), pdStore = double(n.draws*2*n),
              PACKAGE="eco")$pdStore
    res <- matrix(res, ncol=2, nrow=n.draws*n, byrow=TRUE)
    colnames(res) <- c("W1", "W2")
  }
  else { # unconditional prediction
    res <- .C("preDP", as.double(mu), as.double(Sigma), as.integer(n),
              as.integer(n.draws), as.integer(3), as.integer(verbose),
              pdStore = double(n.draws*3*n), PACKAGE="eco")$pdStore
    
    res <- matrix(res, ncol=3, nrow=n.draws*n, byrow=TRUE)
    colnames(res) <- c("W1", "W2", "X")
  }
  
  class(res) <- c("predict.eco", "matrix")
  return(res)
}
