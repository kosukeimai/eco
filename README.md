# R package eco [![Build Status](https://travis-ci.org/kosukeimai/eco.svg?branch=master)](https://travis-ci.org/kosukeimai/eco)  [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/eco)](https://cran.r-project.org/package=eco)
## eco: Ecological inference in 2 x 2 Tables

We implement the Bayesian and likelihood methods proposed 
  in Imai, Lu, and Strauss (2008, 2011) for ecological inference in 2 
  by 2 tables as well as the method of bounds introduced by Duncan and 
  Davis (1953).  The package fits both parametric and nonparametric 
  models using either the Expectation-Maximization algorithms (for 
  likelihood models) or the Markov chain Monte Carlo algorithms (for 
  Bayesian models).  For all models, the individual-level data can be 
  directly incorporated into the estimation whenever such data are available.
  Along with in-sample and out-of-sample predictions, the package also
  provides a functionality which allows one to quantify the effect of data
  aggregation on parameter estimation and hypothesis testing under the
  parametric likelihood models.
