# tests/testthat/test-summary-bsimSpline.R
set.seed(1)

# ---- bsimSpline object ----
make_bsim_obj <- function(modelName = "bsSphere", chains = 1L, n = 30L, p = 3L) {
  stopifnot(chains >= 1, p >= 1)
  x <- matrix(rnorm(n * p), n, p)

  namesIndex     <- paste0("index[", 1:p, "]")
  namesIndicator <- paste0("nu[", 1:p, "]")
  namesSigma     <- "sigma2"

  mk_mat <- function() {
    idx <- matrix(rnorm(n * p), n, p); colnames(idx) <- namesIndex
    nu  <- matrix(sample(0:1, n * p, replace = TRUE), n, p); colnames(nu) <- namesIndicator
    s2  <- matrix(rexp(n, rate = 2), n, 1); colnames(s2) <- namesSigma
    if (modelName %in% c("bsSphere", "bsSpike")) {
      cbind(s2, nu, idx)
    } else {
      cbind(idx, s2)
    }
  }

  if (chains == 1L) {
    mat <- mk_mat()
    sampling <- coda::mcmc(mat)
  } else {
    lst <- replicate(chains, coda::mcmc(mk_mat()), simplify = FALSE)
    sampling <- coda::mcmc.list(lst)
  }

  structure(
    list(
      modelName = modelName,
      input = list(data = list(x = x)),
      sampling = sampling
    ),
    class = "bsimSpline"
  )
}

make_gp_obj <- function(modelName = "bsSphere", chains = 1L, n = 30L, p = 3L) {
  stopifnot(chains >= 1, p >= 1)
  x <- matrix(rnorm(n * p), n, p)
  if (modelName == "gpSpike"){
    namesIndex     <- paste0("indexstar[", 1:p, "]")
  } else{
    namesIndex     <- paste0("index[", 1:p, "]")
  }
  namesIndicator <- paste0("nu[", 1:p, "]")
  namesSigma     <- "sigma2"

  mk_mat <- function() {
    idx <- matrix(rnorm(n * p), n, p); colnames(idx) <- namesIndex
    nu  <- matrix(sample(0:1, n * p, replace = TRUE), n, p); colnames(nu) <- namesIndicator
    s2  <- matrix(rexp(n, rate = 2), n, 1); colnames(s2) <- namesSigma
    if (modelName == "gpSpike") {
      cbind(s2, nu, idx)
    } else {
      cbind(idx, s2)
    }
  }

  if (chains == 1L) {
    mat <- mk_mat()
    sampling <- coda::mcmc(mat)
  } else {
    lst <- replicate(chains, coda::mcmc(mk_mat()), simplify = FALSE)
    sampling <- coda::mcmc.list(lst)
  }

  structure(
    list(
      modelName = modelName,
      input = list(data = list(x = x)),
      sampling = sampling
    ),
    class = "bsimGp"
  )
}


# ------------------------------------------------------------
test_that("summary | variable selection: single chain", {
  objbs <- make_bsim_obj(modelName = "bsSphere", chains = 1L, n = 40L, p = 3L)
  objgp <- make_gp_obj(modelName = "gpSpike", chains = 1L, n = 40L, p = 3L)

  resbs <- expect_invisible(summary.bsimSpline(objbs, verbose = FALSE))
  expect_type(resbs, "list")
  expect_true(is.data.frame(resbs$all.chain))

  resgp <- expect_invisible(summary.bsimGp(objgp, verbose = FALSE))
  expect_type(resgp, "list")
  expect_true(is.data.frame(resgp$all.chain))

  p <- ncol(objbs$input$data$x)
  namesIndex     <- paste0("index[", 1:p, "]")
  namesIndicator <- paste0("nu[", 1:p, "]")
  namesSigma     <- "sigma2"
  expect_identical(rownames(resbs$all.chain),
                   c(namesIndicator, namesIndex, namesSigma))
  expect_identical(rownames(resgp$all.chain),
                   c(namesIndicator, namesIndex, namesSigma))

  # expect_true(all(c("mean","median","sd","LB","UB","ESS") %in% colnames(res$all.chain)))
#
#   expect_output(summary.bsimSpline(objbs, verbose = TRUE), "Summary statistics of all chains")
#   expect_output(summary.bsimGp(objgp, verbose = TRUE), "Summary statistics of all chains")
})

test_that("summary | variable selection: multiple chains", {
  objbs <- make_bsim_obj(modelName = "bsSphere", chains = 2L, n = 30L, p = 4L)
  resbs <- expect_invisible(summary.bsimSpline(objbs, verbose = FALSE))
  objgp <- make_gp_obj(modelName = "gpSpike", chains = 2L, n = 30L, p = 4L)
  resgp <- expect_invisible(summary.bsimGp(objgp, verbose = FALSE))

  expect_true(is.list(resbs$chain))
  expect_length(resbs$chain, 2L)
  expect_identical(names(resbs$chain), c("chain1","chain2"))
  expect_true(all(vapply(resbs$chain, is.data.frame, logical(1))))
  expect_true(all(c("gelman","ESS") %in% colnames(resbs$all.chain)))

  expect_true(is.list(resgp$chain))
  expect_length(resgp$chain, 2L)
  expect_identical(names(resgp$chain), c("chain1","chain2"))
  expect_true(all(vapply(resgp$chain, is.data.frame, logical(1))))
  expect_true(all(c("gelman","ESS") %in% colnames(resgp$all.chain)))

  p <- ncol(objbs$input$data$x)
  expect_identical(
    rownames(resbs$all.chain),
    c(paste0("nu[", 1:p, "]"), paste0("index[", 1:p, "]"), "sigma2")
  )
  expect_identical(
    rownames(resgp$all.chain),
    c(paste0("nu[", 1:p, "]"), paste0("index[", 1:p, "]"), "sigma2")
  )
})

test_that("summary | others: multiple chain", {
  objbs <- make_bsim_obj(modelName = "bsFisher", chains = 1L, n = 35L, p = 2L)
  resbs <- expect_invisible(summary.bsimSpline(objbs, verbose = FALSE))
  objgp <- make_gp_obj(modelName = "gpFisher", chains = 1L, n = 35L, p = 2L)
  resgp <- expect_invisible(summary.bsimGp(objgp, verbose = FALSE))

  expect_true(is.data.frame(resbs$all.chain))
  expect_true("ESS" %in% colnames(resbs$all.chain))
  expect_identical(
    rownames(resbs$all.chain),
    c(paste0("index[", 1:2, "]"), "sigma2")
  )

  expect_true(is.data.frame(resgp$all.chain))
  expect_true("ESS" %in% colnames(resgp$all.chain))
  expect_identical(
    rownames(resgp$all.chain),
    c(paste0("index[", 1:2, "]"), "sigma2")
  )
})

test_that("summary | others: single chain", {
  objbs <- make_bsim_obj(modelName = "bsFisher", chains = 3L, n = 35L, p = 2L)
  resbs <- expect_invisible(summary.bsimSpline(objbs, verbose = FALSE))
  objgp <- make_gp_obj(modelName = "gpFisher", chains = 3L, n = 35L, p = 2L)
  resgp <- expect_invisible(summary.bsimGp(objgp, verbose = FALSE))

  expect_true(all(c("gelman","ESS") %in% colnames(resbs$all.chain)))
  expect_length(resbs$chain, 3L)
  expect_true(all(c("gelman","ESS") %in% colnames(resgp$all.chain)))
  expect_length(resgp$chain, 3L)
})

