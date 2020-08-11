#' Out-of-Sample Posterior Prediction under the Parametric Bayesian Model for
#' Ecological Inference in 2x2 Tables
#' 
#' Obtains out-of-sample posterior predictions under the fitted parametric
#' Bayesian model for ecological inference. \code{predict} method for class
#' \code{eco} and \code{ecoX}.
#' 
#' The posterior predictive values are computed using the Monte Carlo sample
#' stored in the \code{eco} output (or other sample if \code{newdraw} is
#' specified). Given each Monte Carlo sample of the parameters, we sample the
#' vector-valued latent variable from the appropriate multivariate Normal
#' distribution. Then, we apply the inverse logit transformation to obtain the
#' predictive values of proportions, \eqn{W}. The computation may be slow
#' (especially for the nonparametric model) if a large Monte Carlo sample of
#' the model parameters is used. In either case, setting \code{verbose = TRUE}
#' may be helpful in monitoring the progress of the code.
#' 
#' @aliases predict.ecoX
#' @param object An output object from \code{eco} or \code{ecoNP}.
#' @param newdraw An optional list containing two matrices (or three
#' dimensional arrays for the nonparametric model) of MCMC draws of \eqn{\mu}
#' and \eqn{\Sigma}. Those elements should be named as \code{mu} and
#' \code{Sigma}, respectively. The default is the original MCMC draws stored in
#' \code{object}.
#' @param newdata An optional data frame containing a new data set for which
#' posterior predictions will be made. The new data set must have the same
#' variable names as those in the original data.
#' @param subset A scalar or numerical vector specifying the row number(s) of
#' \code{mu} and \code{Sigma} in the output object from \code{eco}. If
#' specified, the posterior draws of parameters for those rows are used for
#' posterior prediction. The default is \code{NULL} where all the posterior
#' draws are used.
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
#' @seealso \code{eco}, \code{predict.ecoNP}
#' @keywords methods
#' @export
predict.ecoX <- function(object, newdraw = NULL, subset = NULL,
                         newdata = NULL, cond = FALSE, verbose = FALSE, ...){

  if (is.null(newdraw) && is.null(object$mu))
    stop("Posterior draws of mu and Sigma must be supplied")
  else if (!is.null(newdraw)){
    if (is.null(newdraw$mu) && is.null(newdraw$Sigma))
      stop("Posterior draws of both mu and Sigma must be supplied.")
    object <- newdraw
  }

  if (cond) { ## conditional prediction
    mu <- coef(object, subset = subset)
    n.draws <- nrow(mu)
    if (is.null(subset))
      subset <- 1:n.draws
    Sigma <- object$Sigma[subset,]
    if (is.null(newdata))
      X <- object$X
    else {
      mf <- match.call()
      if (is.matrix(eval.parent(mf$newdata)))
        data <- as.data.frame(data)
      tt <- terms(object)
      attr(tt, "intercept") <- 0
      X <- model.matrix(tt, newdata)
    }
    n <- nrow(X)
    res <- .C("preBaseX", as.double(X), as.double(mu), as.double(t(Sigma)),
              as.integer(length(c(X))), as.integer(nrow(mu)),
              as.integer(verbose),
              pdStore = double(n.draws*n*2), PACKAGE="eco")$pdStore
    res <- array(res, c(2, n, n.draws), dimnames=list(c("W1", "W2"),
                                          rownames(X), 1:n.draws))
    class(res) <- c("predict.ecoX", "array")
  }
  else {
    res <- predict.eco(object, newdraw = newdraw, subset = subset,
                       newdata = newdata, verbose = verbose, ...)
    colnames(res) <- c("W1", "W2", "X")
  }
  return(res)
}
