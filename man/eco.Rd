\name{eco}
\alias{ECO}
\alias{eco}
\title{Parametric and Nonparametric Models for Ecological Inference in 2 by 2 Tables}

\description{
\code{eco} is used to produce the in-sample predictions of the
unknown internal cells for \eqn{2 \times 2} ecological tables and/or
out-of-sample predictions. Two models based on data augmentation.
The parametric model uses a Normal/Inverse Wishart prior while the
nonparametric model uses a Dirichlet Process Prior with a base
prior distribution of Normal/Inverse Wishart. Both models are
fitted via Markov chain Monte Carlo. For details, refer to Imai
and Lu (2004).}

\usage{
eco(n.gen, prior=TRUE, link=1, Y, X, nu0=4, tau0=1, mu0=c(0,0),
S0=diag(8,2), a.update=FALSE, burn.in=0, nth.draw=5,
predict=TRUE, a0=0.1, b0=0.1, survey=FALSE, survey.W1=NA,
survey.W2=NA)
}

\arguments{
\item{n.gen}{ A positive integer. The number of MCMC draws.
The default is 5000. }
\item{prior}{ Logical. If \code{TRUE}, the module for nonparametric model is called. Otherwise the module for parametric model is called. The default is \code{TRUE}. }
\item{link}{ The choice of link function. \code{link} equals to 1 corresponds to ``logit'' transformation; equals to 2 corresponds to ``probit'' transformation; equals to 3 corresponds to ``c-log-log'' transformation. The default is 1. }
\item{Y}{ A numeric vector. \code{Y} contains the weighted sums of the internal cells. The weights are \code{X}.}
\item{X} A numeric vector. \code{X} contains the weights.}
\item{nu0}{ A positive integer. The prior degrees of freedom parameter. the default is 4. For details see the descriptions of the models.}
\item{tau0}{ A positive integer. A scale parameter. The defualt is 1. For details see the descriptions of the models.}
\item{mu0}{ A \eqn{2 \times 1} numeric vector. The prior mean. The defualt is \eqn{(0,0)}. For details see the descriptions of the models.}
\item{S0}{ A \eqn{2 \times 2} numeric matrix. A positive definite prior scale matrix. The default is \code{diag(8,2)}. For details see the descriptions of the models.}  
\item{a.update}{Logical. If \code{TRUE}, the concentration parameter \eqn{\alpha} will be updated at each Gibbs draw. Otherwise, \eqn{\alpha} is fixed at given value. }
\item{burn.in}{ A positive integer. The burnin interval for the Markov chain; i.e. the number of initial draws that should not be stored. The default is 0.}  
\item{nth.draw}{ A positive integer. The thinning interval for the Markov chain; i.e. the number of Gibbs draws between the recorded values that are skipped. Thge default is 5.}
\item{predict}{ Logical. If \code{TRUE}, out-of sample predictions will be produced. The default is \code{TRUE}.}
\item{a0}{ A positive integer. The shape parameter of the Gamma prior for \eqn{\alpha}. The default is 0.1.}
\item{b0}{ A positive integer. The scale parameter of the Gamma prior for \eqn{\alpha}. The default is 0.1.}
\item{survey}{ Logical. If \code{TRUE}, the additional individual-level data such as survey data is included in the model.} 
\item{survey.W1}{ A numeric vector. The additional individual-level data of \eqn{W_1}.}
\item{survey.W2} {A numeric vector. The additional individual-level data of \eqn{W_2}.}
%\item[\code{x1type}] 
%\item[\code{x1.W1}] 
%\item[\code{x0type}] 
%\item[\code{x0.W2}] 
}

\details{
An example of ecological table for racial voting is as following:
    \tabular{lccc}{
      \tab black voters  \tab white voters \tab \cr
       Voted \tab \eqn{W_{1i}}  \tab \eqn{W_{2i}} \tab \eqn{Y_i} \cr
      Not Voted \tab \eqn{1-W_{1i}}  \tab \eqn{1-W_{2i}} \tab \eqn{1-Y_i}  \cr
      \tab \eqn{X_i} \tab \eqn{1-X_i} \tab 
  }

In particular, the data vectors of \eqn{Y} and \eqn{X} and the unknown \eqn{W_1} and \eqn{W_2} have the following deterministic relationship:
\eqn{Y=X W_1+(1-X)W_2}

Use \eqn{1} for the choice of \code{link} in the beta version.

Always specify \code{a.update=TRUE} in the beta version.
}

\examples{
## load the matproii data data(matproii) 
## run the nonparametric model to produce in-sample and out-of sample predictions.

data(matproii)
data<-cbind(matproii[,1:2], matproii[, 4:5])
names(data)<-c(``X'', ``Y'', ``BB'', ``BW'')

Y <- data[,2]
X <- data[,1]

test <-eco(n.gen=20000, prior=TRUE, link=1, predict=TRUE, burn.in=5000,
            nth.draw=5,S0=S,Y, X, nu0=nu, tau0=tau,
            a.update=alpha, a0=1, b0=0.1)

W1.post<-colMeans(test$w1.post) #create in-sample predicitions of \eqn{W_1}
W2.post<-colMeans(test$w2.post) #create in-sample predicitions of \eqn{W_2}
}


\value{
An object of class \code{eco} containing the following elements:
\item{model}{ The name of the model is used to produce the predictions. If the nonparametric model is used, \code{ model=''Dirichlet Process prior''}; if the parametric model is used, \code{ model=''Normal prior''}. }
\item{ALPHA}{ Whether \eqn{\alpha} is being updated at each Gibbs draw. Only applies when \code{ model=''Dirichlet Process prior''}.}
\item{burn.in}{ The number of initial burnin draws.}
\item{X}{ The vector of data \eqn{X}.}
\item{Y}{ The vector of data \eqn{Y}.}
\item{nu0}{ The prior degrees of freedom.}
\item{tau0}{ The prior scale parameter.}
\item {mu0}{The prior means.}
\item{S0}{ The prior scale matrix.}
\item {a0}{ The prior shape parameter. Only applies when \code{ model=''Dirichlet Process prior''}.}
\item{b0}{ The prior scale parameter. Only applies when \code{ model=''Dirichlet Process prior''}.}
\item{mu1.post}{ The posterior draws of the population mean parameter. }
\item{mu2.post}{ The posterior draws of the population mean parameter. }
\item{Sigma11.post}{ The posterior draws of the population variance parameter. }
\item{Sigma12.post}{ The posterior draws of the population covariance parameter. }
\item{Sigma22.post}{ The posterior draws of the population variance parameter.} 
\item{w1.post}{ The posterior draws of \eqn{W_1}, i.e. the in-sample predicitions of \eqn{W_1}.}
\item{w2.post}{ The posterior draws of \eqn{W_2}, i.e. the in-sample predicitions of \eqn{W_2}.}
\item{w1.pred}{ The posterior predicitive draws of \eqn{W_1}, i.e. the out-of-sample predicitions of \eqn{W_1}.}
\item{w2.pred}{ The posterior predicitive draws of \eqn{W_2}, i.e. the out-of-sample predicitions of \eqn{W_2}.}
\item{Y.pred}{ The posterior predicitive draws of \eqn{Y}.}
\item{a.post}{ The Gibbs draws of \eqn{\alpha}. Only applies when \code{ model=''Dirichlet Process prior''}.}
\item{nstar}{ The number of clusters at each Gibbs draw. Only applies when \code{ model=''Dirichlet Process prior''}.}
}
\author{
Kosuke Imai, Department of Politics, Princeton University
\email{kimai@Princeton.Edu},
\url{http://www.princeton.edu/~kimai}; Ying Lu, Woodrow Wilson School of International and Public Affairs, Princeton University \email{yinglu@Princeton.Edu}.
}
\references{
Imai, Kosuke and Ying Lu. (2004) \dQuote{ Parametric and Nonparametric Bayesian Models for Ecological Inference in \eqn{2 \times 2} Tables. }
\url{http://www.princeton.edu/~kimai/research/einonpar.html} }

\seealso{}
\keyword{ecological inferences}
