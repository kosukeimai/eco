\name{reg}

\docType{data}

\alias{reg}

\title{Voter Registration in US Southern States} 

\description{
  This data set contains the racial composition (\eqn{X}), the
  registration rate (\eqn{Y}), the number of eligible voters
  (\eqn{N}) as well as the actual observed racial registration rates
  (\eqn{BB} and \eqn{BW}) for every county in four US southern states: 
  Florida, Louisiana, North Carolina, and South Carolina.
}

\usage{data(reg)}

\format{A data frame containing 5 variables and 275 observations 
 \tabular{lll}{
    X \tab numeric \tab the fraction of Black people \cr
    Y \tab numeric \tab the fraction of perople who register \cr
    N \tab numeric \tab the total number of people in each county \cr
    BB \tab numeric \tab the actual fraction of Black who register \cr
    BW \tab numeric \tab the actual fraction of White who register 
  }
}
\references{Gary King (1997) \dQuote{A Solution to the Ecological
    Inference Problem: Reconstructing Individual Behavior from Aggregate
    Data}. Princeton University Press, Princeton, NJ.
}

\keyword{datasets}
