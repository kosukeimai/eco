\name{ecoBD}

\alias{ecoBD}

\title{Calculating the Bounds for Ecological Inference in RxC Tables}

\description{
  \code{ecoBD} is used to calculate the bounds for missing internal cells
  of \eqn{R \times C} ecological table. The data can be entered either
  in the form of counts or proportions.
}

\usage{
ecoBD(formula, data = parent.frame(), N = NULL)
}

\arguments{
  \item{formula}{A symbolic description of ecological table to be used,
    specifying the column and row margins of \eqn{R \times
      C} ecological tables.
    Details and specific
    examples are given below. 
  } 
  \item{data}{An optional data frame in which to interpret the variables
    in \code{formula}. The default is the environment in which
    \code{ecoBD} is called. 
  }
  \item{N}{An optional variable representing the size of the unit; e.g.,
    the total number of voters. If \code{formula} is entered as counts
    and the last row and/or column is omitted, this input is necessary.}
}

\details{
  The data may be entered either in the form of counts or proportions.
  If proportions are used, \code{formula} may omit the last row and/or
  column of tables, which can be calculated from the remaining margins.
  For example, \code{Y ~ X} specifies \code{Y} as the first column
  margin and \code{X} as the first row margin in \eqn{2 \times 2} tables.
  If counts are used, \code{formula} may omit the last row and/or column
  margin of the table only if \code{N} is supplied. In this example, the
  columns will be labeled as \code{X} and \code{not X}, and the rows
  will be labeled as \code{Y} and \code{not Y}.

  For larger tables, one can use \code{cbind()} and \code{+}. For
  example, \code{cbind(Y1,
    Y2, Y3) ~ X1 + X2 + X3 + X4)} specifies \eqn{3 \times 4} tables.
 
  An \eqn{R \times C} ecological table in the form of counts:
  \tabular{lcccc}{
    \eqn{n_{i11}} \tab \eqn{n_{i12}}  \tab \dots \tab \eqn{n_{i1C}} \tab \eqn{n_{i1.}} \cr
    \eqn{n_{i21}} \tab \eqn{n_{i22}}  \tab \dots \tab \eqn{n_{i2C}} \tab \eqn{n_{i2.}}  \cr
    \dots         \tab \dots          \tab \dots \tab \dots         \tab \dots\cr
    \eqn{n_{iR1}} \tab \eqn{n_{iR2}}  \tab \dots \tab \eqn{n_{iRC}} \tab \eqn{n_{iR.}}  \cr
    \eqn{n_{i.1}}  \tab \eqn{n_{i.2}} \tab \dots \tab \eqn{n_{i.C}} \tab \eqn{N_i}
  }
  where \eqn{n_{nr.}} and \eqn{n_{i.c}} represent the observed margins,
  \eqn{N_i} represents the size of the table, and \eqn{n_{irc}} are
  unknown variables. Note that for each \eqn{i}, the following
  deterministic relationships hold; \eqn{n_{ir.} = \sum_{c=1}^C
    n_{irc}} for \eqn{r=1,\dots,R}, and \eqn{n_{i.c}=\sum_{r=1}^R n_{irc}}
  for \eqn{c=1,\dots,C}. Then, each of the unknown inner cells can be
  bounded in the following manner,
  \deqn{\max(0, n_{ir.}+n_{i.c}-N_i) \le n_{irc} \le \min(n_{ir.}, n_{i.c}).}
  If the size of tables, \code{N}, is provided,
  
  An \eqn{R \times C} ecological table in the form of proportions:
  \tabular{lcccc}{
    \eqn{W_{i11}} \tab \eqn{W_{i12}} \tab \dots \tab \eqn{W_{i1C}} \tab \eqn{Y_{i1}} \cr
    \eqn{W_{i21}} \tab \eqn{W_{i22}} \tab \dots \tab \eqn{W_{i2C}} \tab \eqn{Y_{i2}}  \cr
    \dots         \tab \dots  \tab \dots \tab \dots \tab \dots \cr
    \eqn{W_{iR1}} \tab \eqn{W_{iR2}} \tab \dots \tab \eqn{W_{iRC}} \tab \eqn{Y_{iR}}  \cr
    \eqn{X_{i1}}  \tab \eqn{X_{i2}}  \tab \dots \tab \eqn{X_{iC}} \tab 
  }
  where \eqn{Y_{ir}} and \eqn{X_{ic}} represent the observed margins,
  and \eqn{W_{irc}} are unknown variables. Note that for each \eqn{i},
  the following deterministic relationships hold; \eqn{Y_{ir} =
    \sum_{c=1}^C X_{ic} W_{irc}} for \eqn{r=1,\dots,R}, and \eqn{\sum_{r=1}^R
    W_{irc}=1} for \eqn{c=1,\dots,C}. Then, each of
  the inner cells of the table can be bounded in the following manner,
  \deqn{\max(0, (X_{ic} + Y_{ir}-1)/X_{ic}) \le W_{irc}
    \le \min(1, Y_{ir}/X_{ir}).}
}

\examples{

## load the registration data
data(reg)

## calculate the bounds
res <- ecoBD(Y ~ X, N = N, data = reg)
## print the results
print(res)
}

\value{
  An object of class \code{ecoBD} containing the following elements
  (When three dimensional arrays are used, the first dimension indexes
  the observations, the second dimension indexes the row numbers, and
  the third dimension indexes the column numbers): 
  \item{call}{The matched call.}
  \item{X}{A matrix of the observed row margin, \eqn{X}.}
  \item{Y}{A matrix of the observed column margin, \eqn{Y}.}
  \item{N}{A vector of the size of ecological tables, \eqn{N}.}  
  \item{aggWmin}{A three dimensional array of aggregate lower bounds for
    proportions.}
  \item{aggWmax}{A three dimensional array of aggregate upper bounds for
    proportions.}
  \item{Wmin}{A three dimensional array of lower bounds for
    proportions.}
  \item{Wmax}{A three dimensional array of upper bounds for
    proportions.}
  \item{Nmin}{A three dimensional array of lower bounds for
    counts.}
  \item{Nmax}{A three dimensional array of upper bounds for
    counts.}
  The object can be printed through \code{print.ecoBD}.
}

\references{ Imai, Kosuke, Ying Lu and Aaron Strauss. (2011)
  \dQuote{eco: R Package for Ecological Inference in 2x2 Tables}
  Journal of Statistical Software, Vol. 42, No. 5, pp. 1-23. available
  at \url{http://imai.princeton.edu/software/eco.html}

  Imai, Kosuke, Ying Lu and Aaron Strauss. (2008)
  \dQuote{Bayesian and Likelihood Inference for 2 x 2 Ecological
  Tables: An Incomplete Data Approach} Political Analysis, Vol. 16,
  No. 1, (Winter), pp. 41-69.  available at
  \url{http://imai.princeton.edu/research/eiall.html} }

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}, \url{http://imai.princeton.edu/};
  Ying Lu, Center for Promoting Research Involving Innovative Statistical Methodology (PRIISM), New York University
   \email{ying.lu@nyu.Edu}
}

\seealso{\code{eco}, \code{ecoNP}}
\keyword{models}
