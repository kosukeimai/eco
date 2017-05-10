#' Calculating the Bounds for Ecological Inference in RxC Tables
#' 
#' \code{ecoBD} is used to calculate the bounds for missing internal cells of
#' \eqn{R \times C} ecological table. The data can be entered either in the
#' form of counts or proportions.
#' 
#' The data may be entered either in the form of counts or proportions.  If
#' proportions are used, \code{formula} may omit the last row and/or column of
#' tables, which can be calculated from the remaining margins.  For example,
#' \code{Y ~ X} specifies \code{Y} as the first column margin and \code{X} as
#' the first row margin in \eqn{2 \times 2} tables.  If counts are used,
#' \code{formula} may omit the last row and/or column margin of the table only
#' if \code{N} is supplied. In this example, the columns will be labeled as
#' \code{X} and \code{not X}, and the rows will be labeled as \code{Y} and
#' \code{not Y}.
#' 
#' For larger tables, one can use \code{cbind()} and \code{+}. For example,
#' \code{cbind(Y1, Y2, Y3) ~ X1 + X2 + X3 + X4)} specifies \eqn{3 \times 4}
#' tables.
#' 
#' An \eqn{R \times C} ecological table in the form of counts: \tabular{lcccc}{
#' \eqn{n_{i11}} \tab \eqn{n_{i12}} \tab \dots{} \tab \eqn{n_{i1C}} \tab
#' \eqn{n_{i1.}} \cr \eqn{n_{i21}} \tab \eqn{n_{i22}} \tab \dots{} \tab
#' \eqn{n_{i2C}} \tab \eqn{n_{i2.}} \cr \dots{} \tab \dots{} \tab \dots{} \tab
#' \dots{} \tab \dots{}\cr \eqn{n_{iR1}} \tab \eqn{n_{iR2}} \tab \dots{} \tab
#' \eqn{n_{iRC}} \tab \eqn{n_{iR.}} \cr \eqn{n_{i.1}} \tab \eqn{n_{i.2}} \tab
#' \dots{} \tab \eqn{n_{i.C}} \tab \eqn{N_i} } where \eqn{n_{nr.}} and
#' \eqn{n_{i.c}} represent the observed margins, \eqn{N_i} represents the size
#' of the table, and \eqn{n_{irc}} are unknown variables. Note that for each
#' \eqn{i}, the following deterministic relationships hold; \eqn{n_{ir.} =
#' \sum_{c=1}^C n_{irc}} for \eqn{r=1,\dots,R}, and \eqn{n_{i.c}=\sum_{r=1}^R
#' n_{irc}} for \eqn{c=1,\dots,C}. Then, each of the unknown inner cells can be
#' bounded in the following manner, \deqn{\max(0, n_{ir.}+n_{i.c}-N_i) \le
#' n_{irc} \le \min(n_{ir.}, n_{i.c}).} If the size of tables, \code{N}, is
#' provided,
#' 
#' An \eqn{R \times C} ecological table in the form of proportions:
#' \tabular{lcccc}{ \eqn{W_{i11}} \tab \eqn{W_{i12}} \tab \dots{} \tab
#' \eqn{W_{i1C}} \tab \eqn{Y_{i1}} \cr \eqn{W_{i21}} \tab \eqn{W_{i22}} \tab
#' \dots{} \tab \eqn{W_{i2C}} \tab \eqn{Y_{i2}} \cr \dots{} \tab \dots{} \tab
#' \dots{} \tab \dots{} \tab \dots{} \cr \eqn{W_{iR1}} \tab \eqn{W_{iR2}} \tab
#' \dots{} \tab \eqn{W_{iRC}} \tab \eqn{Y_{iR}} \cr \eqn{X_{i1}} \tab
#' \eqn{X_{i2}} \tab \dots{} \tab \eqn{X_{iC}} \tab } where \eqn{Y_{ir}} and
#' \eqn{X_{ic}} represent the observed margins, and \eqn{W_{irc}} are unknown
#' variables. Note that for each \eqn{i}, the following deterministic
#' relationships hold; \eqn{Y_{ir} = \sum_{c=1}^C X_{ic} W_{irc}} for
#' \eqn{r=1,\dots,R}, and \eqn{\sum_{r=1}^R W_{irc}=1} for \eqn{c=1,\dots,C}.
#' Then, each of the inner cells of the table can be bounded in the following
#' manner, \deqn{\max(0, (X_{ic} + Y_{ir}-1)/X_{ic}) \le W_{irc} \le \min(1,
#' Y_{ir}/X_{ir}).}
#' 
#' @param formula A symbolic description of ecological table to be used,
#' specifying the column and row margins of \eqn{R \times C} ecological tables.
#' Details and specific examples are given below.
#' @param data An optional data frame in which to interpret the variables in
#' \code{formula}. The default is the environment in which \code{ecoBD} is
#' called.
#' @param N An optional variable representing the size of the unit; e.g., the
#' total number of voters. If \code{formula} is entered as counts and the last
#' row and/or column is omitted, this input is necessary.
#' @return An object of class \code{ecoBD} containing the following elements
#' (When three dimensional arrays are used, the first dimension indexes the
#' observations, the second dimension indexes the row numbers, and the third
#' dimension indexes the column numbers): 
#' \item{call}{The matched call.}
#' \item{X}{A matrix of the observed row margin, \eqn{X}.} 
#' \item{Y}{A matrix of the observed column margin, \eqn{Y}.} 
#' \item{N}{A vector of the size of ecological tables, \eqn{N}.} 
#' \item{aggWmin}{A three dimensional array of
#' aggregate lower bounds for proportions.} 
#' \item{aggWmax}{A three dimensional array of aggregate upper bounds for proportions.} 
#' \item{Wmin}{A three dimensional array of lower bounds for proportions.} 
#' \item{Wmax}{A three dimensional array of upper bounds for proportions.} 
#' \item{Nmin}{A three dimensional array of lower bounds for counts.} 
#' \item{Nmax}{A three dimensional array of upper bounds for counts.} The object 
#' can be printed through \code{print.ecoBD}.
#' @author Kosuke Imai, Department of Politics, Princeton University
#' \email{kimai@@Princeton.Edu}, \url{http://imai.princeton.edu/}; Ying Lu,
#' Center for Promoting Research Involving Innovative Statistical Methodology
#' (PRIISM), New York University \email{ying.lu@@nyu.Edu}
#' @seealso \code{eco}, \code{ecoNP}
#' @references Imai, Kosuke, Ying Lu and Aaron Strauss. (2011) \dQuote{eco: R
#' Package for Ecological Inference in 2x2 Tables} Journal of Statistical
#' Software, Vol. 42, No. 5, pp. 1-23. available at
#' \url{http://imai.princeton.edu/software/eco.html}
#' 
#' Imai, Kosuke, Ying Lu and Aaron Strauss. (2008) \dQuote{Bayesian and
#' Likelihood Inference for 2 x 2 Ecological Tables: An Incomplete Data
#' Approach} Political Analysis, Vol. 16, No. 1, (Winter), pp. 41-69.
#' available at \url{http://imai.princeton.edu/research/eiall.html}
#' @keywords models
#' @examples
#' 
#' 
#' ## load the registration data
#' data(reg)
#' 
#' ## calculate the bounds
#' res <- ecoBD(Y ~ X, N = N, data = reg)
#' ## print the results
#' print(res)
#' 
ecoBD <- function(formula, data = parent.frame(), N=NULL){
  mf <- match.call()
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  vnames <- attr(tt, "variables")
  vnamesR <- vnames[[2]]
  
  if (is.matrix(eval.parent(mf$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  Y <- as.matrix(model.response(model.frame(tt, data = data)))
  N <- eval(mf$N, data)
  n.obs <- nrow(X)

  ## counts
  if (all(X>1) & all(Y>1)) {
    if (!is.null(N)) {
      if (!all(apply(X, 1, sum) == N))
        X <- cbind(X, N-apply(X, 1, sum))
      if (!all(apply(Y, 1, sum) == N))
        Y <- cbind(Y, N-apply(Y, 1, sum))
      if(any(X<0) || any(Y<0))
        stop("Invalid inputs for X, Y, or/and N")
    }
    else {
      if (!all(apply(X, 1, sum) == apply(Y, 1, sum)))
        stop("X and Y do not sum to the same number. Input N.")
      N <- apply(X, 1, sum)
    }
    C <- ncol(X)
    R <- ncol(Y)
    Wmin <- Wmax <- Nmin <- Nmax <- array(NA, c(n.obs, R, C))
    clab <- rlab <- NULL
    if (length(vnames) == 3)
      clab <- c(vnames[[3]], paste("not",vnames[[3]]))
    else {
      for (j in 1:C) {
        if ((j == C) & (length(vnames) < j+2))
          clab <- c(clab, "other")
        else
          clab <- c(clab, vnames[[j+2]])
      }
    }
    if (length(vnamesR) == 1)
      rlab <- c(vnamesR, paste("not",vnamesR))
    else {
      for (i in 1:R) {
        if ((i == R) & (length(vnamesR) < i+1))
          rlab <- c(rlab, "other")
        else
          rlab <- c(rlab, vnamesR[[i]])
      }
    }
    for (i in 1:R) {
      for (j in 1:C) {
        Nmin[,i,j] <- apply(cbind(0, X[,j]+Y[,i]-N), 1, max)
        Nmax[,i,j] <- apply(cbind(Y[,i], X[,j]), 1, min)
        Wmin[,i,j] <- Nmin[,i,j]/X[,j]
        Wmax[,i,j] <- Nmax[,i,j]/X[,j]
      }
    }
    dimnames(Wmin) <- dimnames(Wmax) <- dimnames(Nmin) <-
      dimnames(Nmax) <-
        list(if (is.null(rownames(X))) 1:n.obs else rownames(X),
             rlab, clab)
  }
  else { ## proportions
    if (any(apply(X, 1, sum) > 1.000000001))
      stop("invalid input for X")
    if (any(apply(X, 1, sum) < 0.9999999999))
      X <- cbind(X, 1-X)
    if (any(apply(Y, 1, sum) > 1.0000000001))
      stop("invalid input for Y")
    if (any(apply(Y, 1, sum) < 0.9999999999))
      Y <- cbind(Y, 1-Y)
    C <- ncol(X)
    R <- ncol(Y)
    Wmin <- Wmax <- array(NA, c(n.obs, R, C))
    clab <- rlab <- NULL
    if (length(vnames) == 3)
      clab <- c(vnames[[3]], paste("not",vnames[[3]]))
    else {
      for (j in 1:C) {
        if ((j == C) & (length(vnames) < j+2))
          clab <- c(clab, "other")
        else
          clab <- c(clab, vnames[[j+2]])
      }
    }
    if (length(vnamesR) == 1)
      rlab <- c(vnamesR, paste("not",vnamesR))
    else {
      for (i in 1:R) {
        if ((i == R) & (length(vnamesR) < i+1))
          rlab <- c(rlab, "other")
        else
          rlab <- c(rlab, vnamesR[[i]])
      }
    }
    for (i in 1:R) {
      for (j in 1:C) {
        Wmin[,i,j] <- apply(cbind(0, (X[,j]+Y[,i]-1)/X[,j]), 1, max)
        Wmax[,i,j] <- apply(cbind(1, Y[,i]/X[,j]), 1, min)
      }
    }
    dimnames(Wmin) <- dimnames(Wmax) <-
      list(if (is.null(rownames(X))) 1:n.obs else rownames(X),
           rlab, clab)
    colnames(X) <- clab
    colnames(Y) <- rlab
    if (!is.null(N)) {
      Nmin <- Nmax <- array(NA, c(n.obs, R, C), dimnames =
                            dimnames(Wmin))
      for (i in 1:R) 
        for (j in 1:C) {
          Nmin[,i,j] <- Wmin[,i,j]*X[,j]*N
          Nmax[,i,j] <- Wmax[,i,j]*X[,j]*N
        }
    }
    else
      Nmin <- Nmax <- NULL
  }

  ## aggregate bounds
  aggWmin <- aggWmax <- matrix(NA, R, C, dimnames =
                               list(dimnames(Wmin)[[2]], dimnames(Wmin)[[3]]))
  if (is.null(N))
    for (j in 1:C) {
      aggWmin[,j] <- apply(Wmin[,,j], 2, weighted.mean, X[,j])
      aggWmax[,j] <- apply(Wmax[,,j], 2, weighted.mean, X[,j])
    }
  else
    for (j in 1:C) {
      aggWmin[,j] <- apply(Wmin[,,j], 2, weighted.mean, X[,j]*N)
      aggWmax[,j] <- apply(Wmax[,,j], 2, weighted.mean, X[,j]*N)
    }

  if (!is.null(Nmin) & !is.null(Nmax)) {
    aggNmin <- aggNmax <- matrix(NA, R, C, dimnames =
                                 list(dimnames(Nmin)[[2]], dimnames(Nmin)[[3]]))
    for (j in 1:C) {
      aggNmin[,j] <- apply(Nmin[,,j], 2, sum)
      aggNmax[,j] <- apply(Nmax[,,j], 2, sum)
    }
  }
  else
    aggNmin <- aggNmax <- NULL
    
  ## output
  res <- list(call = mf, X = X, Y = Y, N = N, aggWmin = aggWmin,
              aggWmax = aggWmax, aggNmin = aggNmin, aggNmax = aggNmax,
              Wmin = Wmin, Wmax = Wmax, Nmin = Nmin, Nmax = Nmax)
  class(res) <- c("ecoBD", "eco")
  return(res)
}
