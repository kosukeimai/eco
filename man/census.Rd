\name{census}

\docType{data}

\alias{census}

\title{Black Illiteracy Rates in 1910 US Census} 

\description{
  This data set contains the proportion of the residents over 10 years
  of age who are black (\eqn{X}), the proportion of those who can read
  (\eqn{Y}), the total population over 10 years old (\eqn{N}) 
  as well as the actual black literacy rate(\eqn{TrueBB}) and white literacy
  rate(\eqn{TrueBW}) for 1040 counties in the US. The dataset was originally
  analyzed by Robison(1950) at the state level. King(1997) recode 
  the 1910 census at county level.
 
}

\usage{data(census)}

\format{A data frame containing 5 variables and 1040 observations 
 \tabular{lll}{
    X \tab numeric \tab the proportion of Black(\eqn{10+})  \cr
    Y \tab numeric \tab the fraction of perople(\eqn{10+}) who can read \cr
    N \tab numeric \tab the total number of people(\eqn{10+}) in 
		each county \cr
    TrueBB \tab numeric \tab the actual Black literacy rate \cr
    TrueBW \tab numeric \tab the actual White literacy rate
  }
}
\references{
 Robinson, W.S. (1950) Ecological Correlations and the Behavior 
 of Individuals. American Sociologic Review, 15:351-357. \cr
 Gary King (1997) \dQuote{A Solution to the Ecological
 Inference Problem: Reconstructing Individual Behavior from Aggregate
 Data}. Princeton University Press, Princeton, NJ.
}

\keyword{datasets}
