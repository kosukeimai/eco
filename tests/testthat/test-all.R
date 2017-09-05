rm(list=ls())
library(eco)
library(testthat)
context("tests eco")


# set random seed
set.seed(12345)

test_that("tests eco on registration data", {
  ## load the data
  data(reg)

  # fit the parametric model with the default prior specification
  res <- eco(Y ~ X, data = reg, verbose = TRUE)

  # summarize the results
  x <- summary(res)
  expect_that(length(x), is_equivalent_to(8))
  expect_true("W2.table" %in% names(x))
  expect_that(round(x$param.table[2,1], 3), is_equivalent_to(2.976))
  expect_that(round(x$param.table[3,4], 3), is_equivalent_to(8.238))
  
  # obtain out-of-sample prediction
  out <- predict(res, verbose = TRUE)
  # summarize the results
  x <- summary(out)
  expect_that(length(x), is_equivalent_to(2))
  expect_true("n.draws" %in% names(x))
  expect_that(round(x$W.table[1,1], 3), is_equivalent_to(0.490))
  expect_that(round(x$W.table[2,3], 3), is_equivalent_to(0.307))
})  
  
test_that("tests eco on Robinson census", {
  # load the Robinson's census data
  data(census)

  # fit the parametric model with contextual effects and N using the default prior specification
  res1 <- eco(Y ~ X, N = N, context = TRUE, data = census, verbose = TRUE)

  # summarize the results
  x <- summary(res1)
  expect_that(length(x), is_equivalent_to(8))
  expect_true("W2.table" %in% names(x))
  expect_that(round(x$param.table[2,3], 3), is_equivalent_to(2.068))
  expect_that(round(x$agg.wtable[1,3], 3), is_equivalent_to(0.663))
  
  # obtain out-of-sample prediction
  out1 <- predict(res1, verbose = TRUE)
  # summarize the results
  x <- summary(out1)
  expect_that(length(x), is_equivalent_to(2))
  expect_true("n.draws" %in% names(x))
  expect_that(round(x$W.table[1,3], 2), is_equivalent_to(0.45))
  expect_that(round(x$W.table[3,1], 2), is_equivalent_to(0.34))
})




# ecoBD
test_that("tests ecoBD on registration data", {
  # load the registration data
  data(reg)
  # calculate the bounds
  x <- ecoBD(Y ~ X, N = N, data = reg)
  expect_that(length(x), is_equivalent_to(12))
  expect_true("aggWmin" %in% names(x))
  expect_that(round(x$aggWmin[1,1], 3), is_equivalent_to(0.217))
  expect_that(x$aggNmax[2,2], is_equivalent_to(2046800))
})





# ecoML
test_that("tests ecoML on census data", {
  # load the census data
  data(census)

  # fit the parametric model with the default model specifications
  res <- ecoML(Y ~ X, data = census[1:100,], N=census[1:100,3], epsilon=10^(-6), verbose = TRUE)
  # summarize the results
  x <- summary(res)
  expect_that(length(x), is_equivalent_to(13))
  expect_true("iters.sem" %in% names(x))
  expect_that(round(x$loglik, 1), is_equivalent_to(-70.8))
  expect_that(round(x$param.table[2,3], 3), is_equivalent_to(0.072))
  
  #####################################################################################
  # NOTE: this example does not work! There is no predict.ecoML defined in the package. 
  #
  # obtain out-of-sample prediction
  # out <- predict(res, verbose = TRUE)
  # summarize the results
  # summary(out)
  #####################################################################################
  
  # fit the parametric model with some individual 
  # level data using the default prior specification
  surv <- 1:600
  res1 <- ecoML(Y ~ X, context = TRUE, data = census[-surv,], supplement = census[surv,c(4:5,1)], maxit=100, verbose = TRUE)
  # summarize the results
  x <- summary(res1)
  expect_that(length(x), is_equivalent_to(13))
  expect_true("iters.sem" %in% names(x))
  expect_that(round(x$loglik, 1), is_equivalent_to(-3481.9))
  expect_that(round(x$param.table[2,3], 3), is_equivalent_to(0.006))
  expect_true(is.na(x$param.table[2,2]))
  expect_true(is.na(x$param.table[2,5]))
  expect_false(is.na(x$param.table[2,6]))
})





# set random seed
set.seed(12345)

# ecoNP
test_that("tests ecoNP on census data", {
  # load the registration data
  data(reg)

  # NOTE: We set the number of MCMC draws to be a very small number in
  # the following examples; i.e., convergence has not been properly
  # assessed. See Imai, Lu and Strauss (2006) for more complete examples.

  # fit the nonparametric model to give in-sample predictions
  # store the parameters to make population inference later
  res <- ecoNP(Y ~ X, data = reg, n.draws = 50, param = TRUE, verbose = TRUE)
  #summarize the results
  x <- summary(res)
  expect_that(length(x), is_equivalent_to(8))
  expect_true(is.null(x$agg.wtable))
  expect_that(round(x$agg.table[1,2], 4), is_equivalent_to(0.0406))
  expect_that(round(x$agg.table[2,3], 3), is_equivalent_to(0.813))
  
  # obtain out-of-sample prediction
  out <- predict(res, verbose = TRUE)
  # summarize the results
  x <- summary(out)
  expect_that(length(x), is_equivalent_to(2))
  expect_true("n.draws" %in% names(x))
  expect_that(round(x$W.table[1,3], 3), is_equivalent_to(0.026))
  expect_that(round(x$W.table[2,1], 3), is_equivalent_to(0.814))
  
  # density plots of the out-of-sample predictions
  # par(mfrow=c(2,1))
  # plot(density(out[,1]), main = "W1")
  # plot(density(out[,2]), main = "W2")
})

test_that("tests ecoNP on Robinson census data", {
  # load the Robinson's census data
  data(census)

  # fit the parametric model with contextual effects and N 
  # using the default prior specification

  res1 <- ecoNP(Y ~ X, N = N, context = TRUE, param = TRUE, data = census, n.draws = 25, verbose = TRUE)
  # summarize the results
  x <- summary(res1)
  expect_that(length(x), is_equivalent_to(8))
  expect_false(is.null(x$agg.wtable))
  expect_that(round(x$agg.table[1,2], 3), is_equivalent_to(0.010))
  expect_that(round(x$agg.table[2,3], 3), is_equivalent_to(0.869))
  expect_that(round(x$agg.wtable[1,2], 4), is_equivalent_to(0.0095))
  expect_that(round(x$agg.wtable[2,3], 4), is_equivalent_to(0.9005))
  
  # out-of sample prediction 
  pres1 <- predict(res1)
  x <- summary(pres1)
  expect_that(length(x), is_equivalent_to(2))
  expect_true("n.draws" %in% names(x))
  expect_that(round(x$W.table[1,3], 3), is_equivalent_to(0.133))
  expect_that(round(x$W.table[2,1], 3), is_equivalent_to(0.843))
})


