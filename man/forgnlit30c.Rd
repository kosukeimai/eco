\name{forgnlit30c}

\docType{data}

\alias{forgnlit30c}

\title{Foreign-born literacy in 1930, County Level}

\description{
  This data set contains, on a county level, the proportion of
  white residents ten years and older who are foreign born, and
  the proportion of those residents who are literate.  Data come from
  the 1930 census and were first analyzed by Robinson (1950). Counties
  with fewer than 100 foreign born residents are dropped.
}

\usage{data(forgnlit30c)}

\format{A data frame containing 6 variables and 1976 observations
 \tabular{lll}{
    X \tab numeric \tab proportion of the white population at least 10
     years of age that is foreign born \cr
    Y \tab numeric \tab proportion of the white population at least 10
     years of age that is illiterate \cr
    W1 \tab numeric \tab proportion of the foreign-born white population
     at least 10 years of age that is illiterate \cr
    W2 \tab numeric \tab proportion of the native-born white population
     at least 10 years of age that is illiterate \cr
    state \tab numeric \tab the ICPSR state code \cr
    county \tab numeric \tab the ICPSR (within state) county code 
  }
}
\references{
 Robinson, W.S. (1950). ``Ecological Correlations and the Behavior
 of Individuals.'' \emph{American Sociological Review}, vol. 15,
 pp.351-357.
}

\keyword{datasets}
