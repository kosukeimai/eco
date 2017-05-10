#' Summarizing the Results for the Bayesian Parametric Model for Ecological
#' Inference in 2x2 Tables
#' 
#' \code{summary} method for class \code{eco}.
#' 
#' 
#' @aliases summary.eco print.eco
#' @param object An output object from \code{eco}.
#' @param CI A vector of lower and upper bounds for the Bayesian credible
#' intervals used to summarize the results. The default is the equal tail 95
#' percent credible interval.
#' @param param Logical. If \code{TRUE}, the posterior estimates of the
#' population parameters will be provided. The default value is \code{TRUE}.
#' @param units Logical. If \code{TRUE}, the in-sample predictions for each
#' unit or for a subset of units will be provided. The default value is
#' \code{FALSE}.
#' @param subset A numeric vector indicating the subset of the units whose
#' in-sample predications to be provided when \code{units} is \code{TRUE}. The
#' default value is \code{NULL} where the in-sample predictions for each unit
#' will be provided.
#' @param ... further arguments passed to or from other methods.
#' @return \code{summary.eco} yields an object of class \code{summary.eco}
#' containing the following elements: 
#' \item{call}{The call from \code{eco}.}
#' \item{n.obs}{The number of units.} 
#' \item{n.draws}{The number of Monte Carlo samples.} 
#' \item{agg.table}{Aggregate posterior estimates of the marginal
#' means of \eqn{W_1} and \eqn{W_2} using \eqn{X} and \eqn{N} as weights.} If
#' \code{param = TRUE}, the following elements are also included:
#' \item{param.table}{Posterior estimates of model parameters: population mean
#' estimates of \eqn{W_1} and \eqn{W_2} and their logit transformations.} If
#' \code{units = TRUE}, the following elements are also included:
#' \item{W1.table}{Unit-level posterior estimates for \eqn{W_1}.}
#' \item{W2.table}{Unit-level posterior estimates for \eqn{W_2}.}
#' 
#' This object can be printed by \code{print.summary.eco}
#' @author Kosuke Imai, Department of Politics, Princeton University,
#' \email{kimai@@Princeton.Edu}, \url{http://imai.princeton.edu}; Ying Lu,
#' Center for Promoting Research Involving Innovative Statistical Methodology
#' (PRIISM), New York University \email{ying.lu@@nyu.Edu}
#' @seealso \code{eco}, \code{predict.eco}
#' @keywords methods
summary.eco <- function(object, CI = c(2.5, 97.5), param = TRUE,
                        units = FALSE, subset = NULL,...) { 

  n.obs <- ncol(object$W[,1,])
  n.draws <- nrow(object$W[,1,])
      
  if (is.null(subset)) subset <- 1:n.obs 
  else if (!is.numeric(subset))
    stop("Subset should be a numeric vector.")
  else if (!all(subset %in% c(1:n.obs)))
    stop("Subset should be any numbers in 1:obs.")
  
  table.names<-c("mean", "std.dev", paste(min(CI), "%", sep=" "),
                 paste(max(CI), "%", sep=" ")) 

  
  agg.table <-agg.wtable <-NULL
  N<-rep(1, length(object$X))
  W1.agg.mean <- as.vector(object$W[,1,]%*% (object$X*N/sum(object$X*N)))
  W2.agg.mean <- as.vector(object$W[,2,]%*% ((1-object$X)*N/sum((1-object$X)*N)))

  agg.table <- rbind(cbind(mean(W1.agg.mean), sd(W1.agg.mean), 
                           quantile(W1.agg.mean, min(CI)/100), 
                           quantile(W1.agg.mean, max(CI)/100)),
                     cbind(mean(W2.agg.mean), sd(W2.agg.mean), 
                           quantile(W2.agg.mean, min(CI)/100), 
                           quantile(W2.agg.mean, max(CI)/100)))
  colnames(agg.table) <- table.names
  rownames(agg.table) <- c("W1", "W2")

    
  if (!is.null(object$N)) {
    N <- object$N

    W1.agg.wmean <- as.vector(object$W[,1,] %*% (object$X*N/sum(object$X*N)))
    W2.agg.wmean <- as.vector(object$W[,2,] %*% ((1-object$X)*N/sum((1-object$X)*N)))
    agg.wtable <- rbind(cbind(mean(W1.agg.wmean), sd(W1.agg.wmean),	
                           quantile(W1.agg.wmean, min(CI)/100),			
                           quantile(W1.agg.wmean, max(CI)/100)),	
                     cbind(mean(W2.agg.wmean), sd(W2.agg.wmean), 
                           quantile(W2.agg.wmean, min(CI)/100), 
                           quantile(W2.agg.wmean, max(CI)/100)))
    colnames(agg.wtable) <- table.names
    rownames(agg.wtable) <- c("W1", "W2")
  }

  
  if (units) {
     W1.table <- cbind(apply(object$W[,1,subset], 2, mean), 
                       apply(object$W[,1,subset], 2, sd),
                       apply(object$W[,1,subset], 2, quantile, min(CI)/100),
                       apply(object$W[,1,subset], 2, quantile, max(CI)/100))
     W2.table <- cbind(apply(object$W[,2,subset], 2, mean), 
                       apply(object$W[,2,subset], 2, sd),
                       apply(object$W[,2,subset], 2, quantile, min(CI)/100),
                       apply(object$W[,2,subset], 2, quantile, max(CI)/100))
     colnames(W2.table) <- colnames(W1.table) <- table.names
     rownames(W1.table) <- rownames(W2.table) <- row.names(object$X[subset])
   }
  else
     W1.table <- W2.table <- NULL
  
  if (param) {
    if (is.null(object$mu) || is.null(object$Sigma))
      stop("Parameters are missing values.")
    else {
      param <- cbind(object$mu, object$Sigma)
      param.table <- cbind(apply(param, 2, mean), apply(param, 2, sd),
                           apply(param, 2, quantile, min(CI)/100),
                           apply(param, 2, quantile, max(CI)/100))
      colnames(param.table) <- table.names
    }
  }
  else
    param.table <- NULL
  
  ans <- list(call = object$call, W1.table = W1.table, W2.table = W2.table,
              agg.table = agg.table, agg.wtable=agg.wtable, param.table = param.table,
              n.draws = n.draws, n.obs = n.obs) 
  
  class(ans) <-"summary.eco"
  return(ans)
}
