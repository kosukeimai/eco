#' Print the Summary of the Results for the Bayesian Nonparametric Model for Ecological
#' Inference in 2x2 Tables
#' 
#' \code{summary} method for class \code{ecoNP}.
#' 
#' 
#' @aliases print.summary.ecoNP
#' @param x An object of class \code{summary.ecoNP}.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#' @return \code{summary.ecoNP} yields an object of class \code{summary.ecoNP}
#' containing the following elements: 
#' \item{call}{The call from \code{ecoNP}.}
#' \item{n.obs}{The number of units.} 
#' \item{n.draws}{The number of Monte Carlo samples.} 
#' \item{agg.table}{Aggregate posterior estimates of the marginal
#' means of \eqn{W_1} and \eqn{W_2} using \eqn{X} and \eqn{N} as weights.} If
#' \code{param = TRUE}, the following elements are also included:
#' \item{param.table}{Posterior estimates of model parameters: population mean
#' estimates of \eqn{W_1} and \eqn{W_2}. If \code{subset} is specified, only a
#' subset of the population parameters are included.} If \code{unit = TRUE},
#' the following elements are also included: 
#' \item{W1.table}{Unit-level posterior estimates for \eqn{W_1}.} 
#' \item{W2.table}{Unit-level posterior estimates for \eqn{W_2}.}
#' 
#' This object can be printed by \code{print.summary.ecoNP}
#' @seealso \code{ecoNP}, \code{predict.eco}
#' @keywords methods
#' @export
print.summary.ecoNP <- function(x, digits=max(3, getOption("digits")-3), ...) 
     {
    cat("\nCall: ") 
    cat(paste(deparse(x$call), sep="\n", collapse="\n"))

    cat("\n\nIn-sample Predictions:\n")
    cat("\nUnweighted:\n")
        print(x$agg.table, digits=digits, na.print="NA",...)
        
    if (!is.null(x$agg.wtable)) {
    cat("\nWeighted:\n")
        print(x$agg.wtable, digits=digits, na.print="NA",...)
    
    }
        cat("\nNumber of Units:", x$n.obs)
        cat("\nNumber of Monte Carlo Draws:", x$n.draws)
        if (!is.null(x$param.table)) {
          tt <- x$param.table
          cat("\nParameter Estimates of mu1:\n")
          print(tt$mu1.table, digits=digits, na.print="NA",...)
          cat("\nParameter Estimates of mu2:\n")
          print(tt$mu2.table, digits=digits, na.print="NA",...)
          cat("\nParameter Estimates of Sigma11:\n")
          print(tt$Sigma11.table, digits=digits, na.print="NA",...)
          cat("\nParameter Estimates of Sigma12:\n")
          print(tt$Sigma12.table, digits=digits, na.print="NA",...)
          cat("\nParameter Estimates of Sigma22:\n")
          print(tt$Sigma22.table, digits=digits, na.print="NA",...)
    }

    if (!is.null(x$W1.table)) {
          cat("\n\nUnit-level Estimates of W1:\n")
          print(x$W1.table, digits=digits, na.print="NA",...)
          cat("\n\nUnit-level Estimates of W2:\n")
          print(x$W2.table, digits=digits, na.print="NA",...)
        }

     cat("\n")
     invisible(x)
}
