test_that("plot.bsimGp: dependencies are available", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("tidyr")
  ok_pipe <- requireNamespace("magrittr", quietly = TRUE)
  if (!ok_pipe) skip("magrittr not installed (pipe needed in plot.bsimGp)")
  succeed()
})

# ---- bsimGp object ----
make_bsimGp <- function(modelName = "gpFisher", chains = 1L, n = 25L, p = 3L) {
  stopifnot(chains >= 1, p >= 1, n >= 2)
  set.seed(1)
  x <- matrix(rnorm(20 * p), 20, p)

  if (identical(modelName, "gpSpike")) {
    nm_theta <- paste0("indexstar[", seq_len(p), "]")
  } else {
    nm_theta <- paste0("index[", seq_len(p), "]")
  }
  nm_sigma <- "sigma2"

  mk_mat <- function(nr = n) {
    th  <- matrix(rnorm(nr * p), nr, p); colnames(th) <- nm_theta
    s2  <- matrix(rexp(nr, rate = 2), nr, 1); colnames(s2) <- nm_sigma
    cbind(th, s2)
  }

  sampling <- if (chains == 1L) mk_mat(n) else lapply(seq_len(chains), function(i) mk_mat(n))
  structure(
    list(modelName = modelName, input = list(data = list(x = x)), sampling = sampling),
    class = "bsimGp"
  )
}


make_bsimSpline <- function(modelName = "bsFisher", chains = 1L, n = 25L, p = 3L) {
  stopifnot(chains >= 1, p >= 1, n >= 2)
  set.seed(1)
  x <- matrix(rnorm(20 * p), 20, p)

  nm_theta <- paste0("index[", seq_len(p), "]")
  nm_sigma <- "sigma2"

  mk_mat <- function(nr = n) {
    th  <- matrix(rnorm(nr * p), nr, p); colnames(th) <- nm_theta
    s2  <- matrix(rexp(nr, rate = 2), nr, 1); colnames(s2) <- nm_sigma
    cbind(th, s2)
  }

  sampling <- if (chains == 1L) mk_mat(n) else lapply(seq_len(chains), function(i) mk_mat(n))
  structure(
    list(modelName = modelName, input = list(data = list(x = x)), sampling = sampling),
    class = "bsimSpline"
  )
}


test_that("plot: single chain (non-spike) returns ggplot with correct mappings", {
  skip_if_not_installed("ggplot2"); skip_if_not_installed("tidyr")
  objgp <- make_bsimGp(modelName = "gpFisher", chains = 1L, n = 30L, p = 2L)
  objbs <- make_bsimSpline(modelName = "bsFisher", chains = 1L, n = 30L, p = 2L)

  gp <- expect_s3_class(BayesSIM:::plot.bsimGp(objgp), "ggplot")
  bs <- expect_s3_class(BayesSIM:::plot.bsimSpline(objbs), "ggplot")

  expect_true(all(c("Iteration","Parameter","Value") %in% colnames(gp$data)))
  expect_equal(levels(gp$data$Parameter),
               c(paste0("index[", 1:2, "]"), "sigma2"))
  expect_true("x" %in% names(gp$mapping))
  expect_true("y" %in% names(gp$mapping))
  expect_false("colour" %in% names(gp$mapping))
  expect_true(inherits(gp$facet, "FacetWrap"))

  expect_true(all(c("Iteration","Parameter","Value") %in% colnames(bs$data)))
  expect_equal(levels(bs$data$Parameter),
               c(paste0("index[", 1:2, "]"), "sigma2"))
  expect_true("x" %in% names(bs$mapping))
  expect_true("y" %in% names(bs$mapping))
  expect_false("colour" %in% names(bs$mapping))
  expect_true(inherits(bs$facet, "FacetWrap"))

})

test_that("plot: multi-chain adds chain color and iteration repetition", {
  skip_if_not_installed("ggplot2"); skip_if_not_installed("tidyr")
  objgp <- make_bsimGp(modelName = "gpFisher", chains = 3L, n = 15L, p = 3L)
  objbs <- make_bsimSpline(modelName = "bsFisher", chains = 3L, n = 15L, p = 3L)

  gp <- expect_s3_class(BayesSIM:::plot.bsimGp(objgp), "ggplot")
  bs <- expect_s3_class(BayesSIM:::plot.bsimSpline(objbs), "ggplot")

  expect_true("colour" %in% names(gp$mapping))
  expect_true("num.chain" %in% colnames(gp$data))
  expect_true(is.factor(gp$data$num.chain))
  expect_identical(levels(gp$data$num.chain), as.character(1:3))
  iters_by_chain <- tapply(gp$data$Iteration, gp$data$num.chain, function(v) length(unique(v)))
  expect_true(all(iters_by_chain == length(unique(gp$data$Iteration[gp$data$num.chain == levels(gp$data$num.chain)[1]]))))

  expect_true("colour" %in% names(bs$mapping))
  expect_true("num.chain" %in% colnames(bs$data))
  expect_true(is.factor(bs$data$num.chain))
  expect_identical(levels(bs$data$num.chain), as.character(1:3))
  iters_by_chain <- tapply(bs$data$Iteration, bs$data$num.chain, function(v) length(unique(v)))
  expect_true(all(iters_by_chain == length(unique(bs$data$Iteration[bs$data$num.chain == levels(bs$data$num.chain)[1]]))))
})

test_that("plot: gpSpike normalizes indexstar to index and facets correctly", {
  skip_if_not_installed("ggplot2"); skip_if_not_installed("tidyr")

  objgp <- make_bsimGp(modelName = "gpSpike", chains = 1L, n = 20L, p = 4L)
  gp <- expect_s3_class(BayesSIM:::plot.bsimGp(objgp), "ggplot")

  expect_equal(levels(gp$data$Parameter),
               c(paste0("index[", 1:4, "]"), "sigma2"))
  expect_true(is.numeric(gp$data$Value))
  expect_false(anyNA(gp$data$Value))
})

