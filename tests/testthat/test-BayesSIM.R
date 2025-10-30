# BayesSIM
# Test for input
# tests/testthat/test-input-args.R
test_that("check input", {
  base_args <- list(
    x = matrix(rnorm(300), 100, 3),
    y = rnorm(100),
    sampling = FALSE,
    fitted   = FALSE,
    niter    = 10,
    nburnin  = 5,
    thin     = 1,
    thin2    = NULL,
    nchain   = 2
  )

  run_call <- function(overrides = list()) {
    args <- base_args
    for (nm in names(overrides)) args[[nm]] <- overrides[[nm]]
    do.call(BayesSIM, args)
  }

  expect_no_error(run_call())

  # --- sampling / fitted: logical ---
  expect_error(run_call(list(sampling = 1L)),           "'sampling' argument should be logical.", fixed = TRUE)
  expect_error(run_call(list(sampling = c(TRUE, FALSE))),"'sampling' argument should be scalar.", fixed = TRUE)
  expect_error(run_call(list(fitted = "yes")),           "'fitted' argument should be logical.",   fixed = TRUE)
  expect_error(run_call(list(fitted = c(TRUE, TRUE))),   "'fitted' argument should be scalar.",   fixed = TRUE)

  # --- niter: numeric scalar & not NA ---
  expect_error(run_call(list(niter = "10")),    "'niter' argument should be numeric scalar.", fixed = TRUE)
  expect_error(run_call(list(niter = c(10, 20))),"'niter' argument should be numeric scalar.", fixed = TRUE)
  expect_error(run_call(list(niter = NA_real_)),"'niter' argument should be numeric scalar.", fixed = TRUE)

  # --- nburnin ---
  expect_error(run_call(list(nburnin = "5")),      "'nburnin' argument should be numeric scalar.", fixed = TRUE)
  expect_error(run_call(list(nburnin = c(1, 2))),  "'nburnin' argument should be numeric scalar.", fixed = TRUE)
  expect_error(run_call(list(nburnin = NA_real_)), "'nburnin' argument should be numeric scalar.", fixed = TRUE)

  # --- thin ---
  expect_error(run_call(list(thin = "1")),       "'thin' argument should be numeric scalar.", fixed = TRUE)
  expect_error(run_call(list(thin = c(1, 2))),   "'thin' argument should be numeric scalar.", fixed = TRUE)
  expect_error(run_call(list(thin = NA_real_)),  "'thin' argument should be numeric scalar.", fixed = TRUE)

  # --- thin2 ---
  expect_no_error(run_call(list(thin2 = NULL)))
  expect_error(run_call(list(thin2 = "1")),        "'thin2' argument should be numeric scalar.", fixed = TRUE)
  expect_error(run_call(list(thin2 = c(1, 2))),    "'thin2' argument should be numeric scalar.", fixed = TRUE)
  expect_error(run_call(list(thin2 = NA_real_)),   "'thin2' argument should be numeric scalar.", fixed = TRUE)

  # --- nchain ---
  expect_error(run_call(list(nchain = "2")),      "'nchain' argument should be numeric scalar.", fixed = TRUE)
  expect_error(run_call(list(nchain = c(1, 2))),  "'nchain' argument should be numeric scalar.", fixed = TRUE)
  expect_error(run_call(list(nchain = NA_real_)), "'nchain' argument should be numeric scalar.", fixed = TRUE)
})


test_that("parameters/connect to correct models", {
  X2 <- matrix(rnorm(20), 10, 2); y <- rnorm(10)
  X6 <- matrix(rnorm(60), 10, 6)

  with_mocked_bindings(
    {
      f1 <- BayesSIM(X2, y, indexprior="fisher", link="bspline",
                     sampling=TRUE, fitted=FALSE, niter=10, nburnin=2,
                     thin=1, thin2=NULL, nchain=2, setSeed=FALSE)
      expect_equal(f1$name, "bsFisher"); expect_identical(f1$p, 2L)

      f2 <- BayesSIM(X2, y, indexprior="sphere", link="gp",
                     sampling=TRUE, fitted=TRUE, method="Nystrom",
                     lowerB = NULL, upperB = NULL,
                     niter=5, nburnin=1, thin=1, nchain=1, setSeed=TRUE)
      expect_equal(f2$name, "gpSphere"); expect_equal(f2$method, "Nystrom")

      f3l <- BayesSIM(X2, y, indexprior="polar", link="gp",
                      sampling=FALSE, fitted=TRUE, niter=3, nburnin=1, thin=1, nchain=1, setSeed=FALSE)
      f3h <- BayesSIM(X6, y, indexprior="polar", link="gp",
                      sampling=FALSE, fitted=TRUE, niter=3, nburnin=1, thin=1, nchain=1, setSeed=FALSE)
      expect_equal(f3l$name, "gpPolar")
      expect_equal(f3h$name, "gpPolarHigh")

      expect_error(BayesSIM(X2, y, indexprior="fisher", link="wrong",
                            sampling=TRUE, fitted=TRUE, niter=2, nburnin=1, thin=1, nchain=1, setSeed=FALSE),
                   "Wrong link function name!")
      expect_error(BayesSIM(X2, y, indexprior="wrong", link="bspline",
                            sampling=TRUE, fitted=TRUE, niter=2, nburnin=1, thin=1, nchain=1, setSeed=FALSE),
                   "Wrong index prior name!")
    },

    ## ---- functions ----
    prior.param.default = function(indexprior, link)
      list(alpha = 1, tag = paste(indexprior, link)),

    init.param.default  = function(indexprior, link)
      list(theta = 0, tag = paste(indexprior, link)),

    param.check = function(user, template) TRUE,

    bsFisher = function(x, y, prior, init, sampling, fitted, monitors2, niter, nburnin, thin, thin2, nchain, setSeed) {
      list(name = "bsFisher", p = ncol(x), prior = prior, init = init, thin2 = thin2)
    },
    bsSphere = function(...)   list(name = "bsSphere"),
    bsPolar  = function(x, y, ...)      list(name = "bsPolar", p = ncol(x)),
    bsSpike  = function(...)   list(name = "bsSpike"),

    gpFisher = function(...)   list(name = "gpFisher"),
    gpSphere = function(x, y, prior, init, sampling, fitted, method,
                        lowerB, upperB,
                        monitors2, niter, nburnin, thin, thin2, nchain, setSeed) {
      list(name = "gpSphere", method = method, p = ncol(x))
    },
    gpPolar  = function(x, y, ...)      list(name = "gpPolar", p = ncol(x)),
    gpPolarHigh = function(x, y, ...)   list(name = "gpPolarHigh", p = ncol(x)),
    gpSpike  = function(...)   list(name = "gpSpike"),

    .package = "BayesSIM"
  )
})
