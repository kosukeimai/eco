rm(list=ls())
library(eco)
library(testthat)
context("tests eco")

accuracy1 <- ifelse(capabilities("long.double"), 0.002, 0.005)
accuracy2 <- ifelse(capabilities("long.double"), 0.2, 0.5)

# set random seed
set.seed(12345)

# for the tests that may take a long time to finish, skip them
donotrun = 1

test_that("tests eco on registration data", {
  ## load the data
  data(reg)

  # fit the parametric model with the default prior specification
  res <- eco(Y ~ X, data = reg, verbose = TRUE)

  # summarize the results
  x <- summary(res)
  expect_that(length(x), is_equivalent_to(8))
  expect_true("W2.table" %in% names(x))
  expect_equal(x$param.table[2,1], 2.976172, tolerance = accuracy1)
  expect_equal(x$param.table[3,4], 8.238363, tolerance = accuracy1)

  # obtain out-of-sample prediction
  out <- predict(res, verbose = TRUE)
  # summarize the results
  x <- summary(out)
  expect_that(length(x), is_equivalent_to(2))
  expect_true("n.draws" %in% names(x))
  expect_equal(x$W.table[1,1], 0.4896912, tolerance = accuracy1)
  expect_equal(x$W.table[2,3], 0.3071461, tolerance = accuracy1)
})  
  
if (!donotrun) test_that("tests eco on Robinson census", {
  # load the Robinson's census data
  data(census)

  # fit the parametric model with contextual effects and N using the default prior specification
  res1 <- eco(Y ~ X, N = N, context = TRUE, data = census, verbose = TRUE)

  # summarize the results
  x <- summary(res1)
  expect_that(length(x), is_equivalent_to(8))
  expect_true("W2.table" %in% names(x))
  expect_equal(x$param.table[2,3], 2.068228, tolerance = accuracy1)
  expect_equal(x$agg.wtable[1,3], 0.6631372, tolerance = accuracy1)

  # obtain out-of-sample prediction
  out1 <- predict(res1, verbose = TRUE)
  # summarize the results
  x <- summary(out1)
  expect_that(length(x), is_equivalent_to(2))
  expect_true("n.draws" %in% names(x))
  expect_equal(x$W.table[1,3], 0.4499757, tolerance = accuracy1)
  expect_equal(x$W.table[3,1], 0.3400277, tolerance = accuracy1)
})

# ecoBD
test_that("tests ecoBD on registration data", {
  # load the registration data
  data(reg)
  # calculate the bounds
  x <- ecoBD(Y ~ X, N = N, data = reg)
  expect_that(length(x), is_equivalent_to(12))
  expect_true("aggWmin" %in% names(x))
  expect_equal(x$aggWmin[1,1], 0.216785, tolerance = accuracy1)
  expect_that(x$aggNmax[2,2], is_equivalent_to(2046800))
})

# ecoML
if (!donotrun) test_that("tests ecoML on census data", {
  # load the census data
  data(census)

  # fit the parametric model with the default model specifications
  res <- ecoML(Y ~ X, data = census[1:100,], N=census[1:100,3], epsilon=10^(-6), verbose = TRUE)
  # summarize the results
  x <- summary(res)
  expect_that(length(x), is_equivalent_to(13))
  expect_true("iters.sem" %in% names(x))
  expect_equal(x$loglik, -70.80674, tolerance = accuracy2)
  expect_equal(x$param.table[2,3], 0.07192878, tolerance = accuracy1)

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
  expect_equal(x$loglik, -3481.877, tolerance = accuracy2)
  expect_equal(x$param.table[2,3], 0.006055498, tolerance = accuracy1)
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
  expect_equal(x$agg.table[1,2], 0.04059766, tolerance = accuracy1)
  expect_equal(x$agg.table[2,3], 0.8129786, tolerance = accuracy1)

  # obtain out-of-sample prediction
  out <- predict(res, verbose = TRUE)
  # summarize the results
  x <- summary(out)
  expect_that(length(x), is_equivalent_to(2))
  expect_true("n.draws" %in% names(x))
  expect_equal(x$W.table[1,3], 0.02617743, tolerance = accuracy1)
  expect_equal(x$W.table[2,1], 0.8137116, tolerance = accuracy1)

  # density plots of the out-of-sample predictions
  # par(mfrow=c(2,1))
  # plot(density(out[,1]), main = "W1")
  # plot(density(out[,2]), main = "W2")
})

if (!donotrun) test_that("tests ecoNP on Robinson census data", {
  # load the Robinson's census data
  data(census)

  # fit the parametric model with contextual effects and N 
  # using the default prior specification

  res1 <- ecoNP(Y ~ X, N = N, context = TRUE, param = TRUE, data = census, n.draws = 25, verbose = TRUE)
  # summarize the results
  x <- summary(res1)
  expect_that(length(x), is_equivalent_to(8))
  expect_false(is.null(x$agg.wtable))
  expect_equal(x$agg.table[1,2], 0.0134030, tolerance = accuracy1)
  expect_equal(x$agg.table[2,3], 0.8709344, tolerance = accuracy1)
  expect_equal(x$agg.wtable[1,2], 0.0113289, tolerance = accuracy1)
  expect_equal(x$agg.wtable[2,3], 0.9033144, tolerance = accuracy1)
 
  # out-of sample prediction 
  pres1 <- predict(res1)
  x <- summary(pres1)
  expect_that(length(x), is_equivalent_to(2))
  expect_true("n.draws" %in% names(x))
  expect_equal(x$W.table[1,3], 0.1333375, tolerance = accuracy1)
  expect_equal(x$W.table[2,1], 0.8434944, tolerance = accuracy1)
})


