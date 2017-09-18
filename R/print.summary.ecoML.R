## for simlicity, this summary function only reports parameters related to W_1 and W_2

#' Print the Summary of the Results for the Maximum Likelihood Parametric Model for
#' Ecological Inference in 2x2 Tables
#' 
#' \code{summary} method for class \code{eco}.
#' 
#' 
#' @aliases print.summary.ecoML
#' @param x An object of class \code{summary.ecoML}.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#' @return \code{summary.eco} yields an object of class \code{summary.eco}
#' containing the following elements: 
#' \item{call}{The call from \code{eco}.}
#' \item{sem}{Whether the SEM algorithm was executed, as specified by the user
#' upon calling \code{ecoML}.} 
#' \item{fix.rho}{Whether the correlation parameter was fixed or allowed to vary, 
#' as specified by the user upon calling \code{ecoML}.} 
#' \item{epsilon}{The convergence threshold specified by the
#' user upon calling \code{ecoML}.} 
#' \item{n.obs}{The number of units.}
#' \item{iters.em}{The number iterations the EM algorithm cycled through before
#' convergence or reaching the maximum number of iterations allowed.}
#' \item{iters.sem}{The number iterations the SEM algorithm cycled through
#' before convergence or reaching the maximum number of iterations allowed.}
#' \item{loglik}{The final observed log-likelihood.} 
#' \item{rho}{A matrix of \code{iters.em} rows specifying the correlation parameters 
#' at each iteration of the EM algorithm. The number of columns depends on how many 
#' correlation parameters exist in the model. Column order is the same as the order of the
#' parameters in \code{param.table}.} 
#' \item{param.table}{Final estimates of the parameter values for the model.  
#' Excludes parameters fixed by the user upon calling \code{ecoML}.  
#' See \code{ecoML} documentation for order of parameters.} 
#' \item{agg.table}{Aggregate estimates of the marginal means of \eqn{W_1} and \eqn{W_2}} 
#' \item{agg.wtable}{Aggregate estimates of the marginal means of \eqn{W_1} and \eqn{W_2} 
#' using \eqn{X} and \eqn{N} as weights.} If \code{units = TRUE}, the following elements 
#' are also included:
#' \item{W.table}{Unit-level estimates for \eqn{W_1} and \eqn{W_2}.}
#' 
#' This object can be printed by \code{print.summary.eco}
#' @author Kosuke Imai, Department of Politics, Princeton University,
#' \email{kimai@@Princeton.Edu}, \url{http://imai.princeton.edu}; Ying Lu,
#' Center for Promoting Research Involving Innovative Statistical Methodology
#' (PRIISM), New York University \email{ying.lu@@nyu.Edu}; Aaron Strauss,
#' Department of Politics, Princeton University,
#' \email{abstraus@@Princeton.Edu}
#' @seealso \code{ecoML}
#' @keywords methods
print.summary.ecoML <- function(x, digits=max(3,
                                     getOption("digits")-3), ...) {

  cat("\nCall: ", paste(deparse(x$call), sep="\n", collapse="\n"))
  cat("\n")
  if (!is.null(x$param.table)) {
    cat("\n*** Parameter Estimates ***\n")
    if (x$fix.rho)
      cat("\nOriginal Model Parameters (rho is fixed at ", x$rho, "):\n", sep="")   
    else
      cat("\nOriginal Model Parameters:\n")
    print(x$param.table, digits=digits, na.print="NA",...)
  }

  cat("\n*** Insample Predictions ***\n")
  cat("\nUnweighted:\n")
  print(x$agg.table, digits=digits, na.print="NA",...)
  
  if (!is.null(x$agg.wtable)) {
  cat("\nWeighted:\n")
  print(x$agg.wtable, digits=digits, na.print="NA",...)
  }
  if (!is.null(x$W.table)) {
    cat("\n\nUnit-level Estimates of W:\n")
    print(x$W.table, digits=digits, na.print="NA",...)
  }

  cat("\n\nLog-likelihood:", x$loglik)
  cat("\nNumber of Observations:", x$n.obs)
  cat("\nNumber of EM iterations:", x$iters.em)
  if (x$sem)
    cat("\nNumber of SEM iterations:", x$iters.sem)
  cat("\nConvergence threshold for EM:", x$epsilon)
  
  cat("\n\n")
  invisible(x)
}
