

#' Black Illiteracy Rates in 1910 US Census
#' 
#' This data set contains the proportion of the residents who are black, the
#' proportion of those who can read, the total population as well as the actual
#' black literacy rate and white literacy rate for 1040 counties in the US. The
#' dataset was originally analyzed by Robinson (1950) at the state level. King
#' (1997) recoded the 1910 census at county level. The data set only includes
#' those who are older than 10 years of age.
#' 
#' 
#' @name census
#' @docType data
#' @format A data frame containing 5 variables and 1040 observations
#' \tabular{lll}{ X \tab numeric \tab the proportion of Black residents in each
#' county\cr Y \tab numeric \tab the overall literacy rates in each county\cr N
#' \tab numeric \tab the total number of residents in each county \cr W1 \tab
#' numeric \tab the actual Black literacy rate \cr W2 \tab numeric \tab the
#' actual White literacy rate }
#' @references Robinson, W.S. (1950). ``Ecological Correlations and the
#' Behavior of Individuals.'' \emph{American Sociological Review}, vol. 15,
#' pp.351-357. \cr \cr King, G. (1997). \dQuote{A Solution to the Ecological
#' Inference Problem: Reconstructing Individual Behavior from Aggregate Data}.
#' Princeton University Press, Princeton, NJ.
#' @keywords datasets
NULL





#' Foreign-born literacy in 1930
#' 
#' This data set contains, on a state level, the proportion of white residents
#' ten years and older who are foreign born, and the proportion of those
#' residents who are literate.  Data come from the 1930 census and were first
#' analyzed by Robinson (1950).
#' 
#' 
#' @name forgnlit30
#' @docType data
#' @format A data frame containing 5 variables and 48 observations
#' \tabular{lll}{ X \tab numeric \tab proportion of the white population at
#' least 10 years of age that is foreign born \cr Y \tab numeric \tab
#' proportion of the white population at least 10 years of age that is
#' illiterate \cr W1 \tab numeric \tab proportion of the foreign-born white
#' population at least 10 years of age that is illiterate \cr W2 \tab numeric
#' \tab proportion of the native-born white population at least 10 years of age
#' that is illiterate \cr ICPSR \tab numeric \tab the ICPSR state code }
#' @references Robinson, W.S. (1950). ``Ecological Correlations and the
#' Behavior of Individuals.'' \emph{American Sociological Review}, vol. 15,
#' pp.351-357.
#' @keywords datasets
NULL





#' Foreign-born literacy in 1930, County Level
#' 
#' This data set contains, on a county level, the proportion of white residents
#' ten years and older who are foreign born, and the proportion of those
#' residents who are literate.  Data come from the 1930 census and were first
#' analyzed by Robinson (1950). Counties with fewer than 100 foreign born
#' residents are dropped.
#' 
#' 
#' @name forgnlit30c
#' @docType data
#' @format A data frame containing 6 variables and 1976 observations
#' \tabular{lll}{ X \tab numeric \tab proportion of the white population at
#' least 10 years of age that is foreign born \cr Y \tab numeric \tab
#' proportion of the white population at least 10 years of age that is
#' illiterate \cr W1 \tab numeric \tab proportion of the foreign-born white
#' population at least 10 years of age that is illiterate \cr W2 \tab numeric
#' \tab proportion of the native-born white population at least 10 years of age
#' that is illiterate \cr state \tab numeric \tab the ICPSR state code \cr
#' county \tab numeric \tab the ICPSR (within state) county code }
#' @references Robinson, W.S. (1950). ``Ecological Correlations and the
#' Behavior of Individuals.'' \emph{American Sociological Review}, vol. 15,
#' pp.351-357.
#' @keywords datasets
NULL





#' Electoral Results for the House and Presidential Races in 1988
#' 
#' This data set contains, on a House district level, the percentage of the
#' vote for the Democratic House candidate, the percentage of the vote for the
#' Democratic presidential candidate (Dukakis), the number of voters who voted
#' for a major party candidate in the presidential race, and the ratio of
#' voters in the House race versus the number who cast a ballot for President.
#' Eleven (11) uncontested races are not included.  Dataset compiled and
#' analyzed by Burden and Kimball (1988). Complete dataset and documentation
#' available at ICSPR study number 1140.
#' 
#' 
#' @name housep88
#' @docType data
#' @format A data frame containing 5 variables and 424 observations
#' \tabular{lll}{ X \tab numeric \tab proportion voting for the Democrat in the
#' presidential race \cr Y \tab numeric \tab proportion voting for the Democrat
#' in the House race \cr N \tab numeric \tab number of major party voters in
#' the presidential contest \cr HPCT \tab numeric \tab House election turnout
#' divided by presidential election turnout (set to 1 if House turnout exceeds
#' presidential turnout) \cr DIST \tab numeric \tab 4-digit ICPSR state and
#' district code: first 2 digits for the state code, last two digits for the
#' district number (e.g., 2106=IL 6th) }
#' @references Burden, Barry C. and David C. Kimball (1988). ``A New Approach
#' To Ticket- Splitting.'' The American Political Science Review. vol 92., no.
#' 3, pp. 553-544.
#' @keywords datasets
NULL





#' Voter Registration in US Southern States
#' 
#' This data set contains the racial composition, the registration rate, the
#' number of eligible voters as well as the actual observed racial registration
#' rates for every county in four US southern states: Florida, Louisiana, North
#' Carolina, and South Carolina.
#' 
#' 
#' @name reg
#' @docType data
#' @format A data frame containing 5 variables and 275 observations
#' \tabular{lll}{ X \tab numeric \tab the fraction of Black voters \cr Y \tab
#' numeric \tab the fraction of voters who registered themselves\cr N \tab
#' numeric \tab the total number of voters in each county \cr W1 \tab numeric
#' \tab the actual fraction of Black voters who registered themselves \cr W2
#' \tab numeric \tab the actual fraction of White voters who registered
#' themselves }
#' @references King, G. (1997). \dQuote{A Solution to the Ecological Inference
#' Problem: Reconstructing Individual Behavior from Aggregate Data}. Princeton
#' University Press, Princeton, NJ.
#' @keywords datasets
NULL





#' Black voting rates for Wallace for President, 1968
#' 
#' This data set contains, on a county level, the proportion of county
#' residents who are Black and the proportion of presidential votes cast for
#' Wallace.  Demographic data is based on the 1960 census. Presidential returns
#' are from ICPSR study 13.  County data from 10 southern states (Alabama,
#' Arkansas, Georgia, Florida, Louisiana, Mississippi, North Carolina, South
#' Carolina, Tennessee, Texas) are included. (Virginia is excluded due to the
#' difficulty of matching counties between the datasets.)  This data is
#' analyzed in Wallace and Segal (1973).
#' 
#' 
#' @name wallace
#' @docType data
#' @format A data frame containing 3 variables and 1009 observations
#' \tabular{lll}{ X \tab numeric \tab proportion of the population that is
#' Black \cr Y \tab numeric \tab proportion presidential votes cast for Wallace
#' \cr FIPS \tab numeric \tab the FIPS county code }
#' @references Wasserman, Ira M. and David R. Segal (1973). ``Aggregation
#' Effects in the Ecological Study of Presidential Voting.'' American Journal
#' of Political Science. vol. 17, pp. 177-81.
#' @keywords datasets
NULL



