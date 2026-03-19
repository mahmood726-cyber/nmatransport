library(testthat)
library(nmatransport)

test_that("assess_transitivity works", {
  set.seed(123)
  df <- data.frame(
    treat1 = c("A", "A", "B", "C"),
    treat2 = c("B", "C", "C", "D"),
    mod1 = c(10, 12, 11, 15),
    mod2 = c(5, 6, 5.5, 8)
  )
  
  res <- assess_transitivity(df, modifiers = c("mod1", "mod2"))
  
  expect_true(is.list(res))
  expect_true(!is.null(res$global_transitivity_violation))
  expect_equal(res$n_comparisons, 4)
})

test_that("compute_nma_weights produces weights summing to 1", {
  set.seed(123)
  df <- data.frame(
    mod1 = rnorm(10, 50, 5),
    mod2 = rnorm(10, 25, 2)
  )
  
  target <- list(mod1 = 55, mod2 = 28)
  weights <- compute_nma_weights(df, target)
  
  expect_equal(sum(weights), 1, tolerance = 1e-4)
})

test_that("transported_nma runs", {
  # Mocking a small network
  df <- data.frame(
    studlab = paste0("S", 1:4),
    treat1 = c("A", "A", "B", "B"),
    treat2 = c("B", "C", "C", "D"),
    TE = c(0.1, 0.2, 0.1, 0.3),
    seTE = c(0.05, 0.05, 0.05, 0.05)
  )
  
  weights <- rep(0.25, 4)
  
  # Check if netmeta is available
  if (requireNamespace("netmeta", quietly = TRUE)) {
    res <- transported_nma(df, weights = weights)
    expect_s3_class(res, "transported_nma")
    expect_true(res$weights_used)
  }
})
