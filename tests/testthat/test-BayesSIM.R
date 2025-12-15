# BayesSIM
# Test for input
# tests/testthat/test-input-args.R
test_that("check input", {
  dataset <- data.frame(matrix(rnorm(300), 100, 3), y = rnorm(100))
  base_args <- list(
    formula = y ~ .,
    data = dataset,
    niter    = 10,
    nburnin  = 5,
    thin     = 1,
    nchain   = 2
  )

  run_call <- function(overrides = list()) {
    args <- base_args
    for (nm in names(overrides)) args[[nm]] <- overrides[[nm]]
    do.call(BayesSIM, args)
    do.call(BayesSIM.setup, args)
  }

  # expect_no_error(run_call())

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

  # --- nchain ---
  expect_error(run_call(list(nchain = "2")),      "'nchain' argument should be numeric scalar.", fixed = TRUE)
  expect_error(run_call(list(nchain = c(1, 2))),  "'nchain' argument should be numeric scalar.", fixed = TRUE)
  expect_error(run_call(list(nchain = NA_real_)), "'nchain' argument should be numeric scalar.", fixed = TRUE)
})


test_that("BayesSIM errors on wrong indexprior", {
  data <- data.frame(y = rnorm(10), x1 = rnorm(10))
  formula <- y ~ x1

  expect_error(
    BayesSIM(formula, data, indexprior = "wrong", link = "bspline"),
    "Wrong index prior name!"
  )
})

test_that("BayesSIM errors on wrong link", {
  data <- data.frame(y = rnorm(10), x1 = rnorm(10))
  formula <- y ~ x1

  expect_error(
    BayesSIM(formula, data, indexprior = "fisher", link = "wrong"),
    "Wrong link function name!"
  )

  expect_error(
    BayesSIM_setup(formula, data, indexprior = "fisher", link = "wrong"),
    "Wrong link function name!"
  )
})

test_that("BayesSIM have wrong output structure", {
  data <- data.frame(y = rnorm(10), x1 = rnorm(10))
  formula <- y ~ x1

  expect_error(
    BayesSIM(formula, data, indexprior = "fisher", link = "wrong"),
    "Wrong link function name!"
  )

  expect_error(
    BayesSIM_setup(formula, data, indexprior = "fisher", link = "wrong"),
    "Wrong link function name!"
  )
})

# Specific model
test_that("bsFisher_setup runs and returns bsimSetup", {
  skip_if_not_installed("nimble")
  skip_on_cran()

  set.seed(123)
  n <- 50
  dat <- data.frame(
    y  = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n)
  )

  fit <- bsFisher_setup(
    formula  = y ~ .,
    data     = dat,
    prior    = NULL,
    init     = NULL,
    niter    = 50,
    nburnin  = 10,
    thin     = 2,
    nchain   = 1,
    setSeed  = FALSE
  )

  # check class
  expect_equal(class(fit), "bsimSetup")

  # check structure
  expect_true(is.list(fit$input))
  expect_equal(fit$modelName, "bsFisher")

  # save original data
  expect_true(is.matrix(fit$input$origdata$x))
  expect_true(is.matrix(fit$input$origdata$y))
  expect_equal(nrow(fit$input$origdata$x), n)
  expect_equal(nrow(fit$input$origdata$y), n)

  # check sampling options
  expect_false(fit$input$samplingOptions$sampling)
  expect_equal(fit$input$samplingOptions$niter, 50)
  expect_equal(fit$input$samplingOptions$nburnin, 10)
  expect_equal(fit$input$samplingOptions$thin, 2)
  expect_equal(fit$input$samplingOptions$nchain, 1)

  expect_error(
    bsFisher_setup(
      formula  = y ~ x1 + x2,
      data     = dat,
      prior    = NULL,
      init     = NULL,
      nchain   = 2,
      setSeed  = "wrong_type"
    ),
    "'setSeed' argument should be logical or numeric vector."
  )

  y  <- rnorm(10)
  x1 <- rnorm(10)
  x2 <- rnorm(10)
  mat_data <- cbind(y, x1, x2)
  expect_error(
    bsFisher_setup(
      formula  = y ~ x1 + x2,
      data     = mat_data,
      prior    = NULL,
      init     = NULL
    ),
    "data should be an data.frame."
  )
})

# prior/initial parameter
test_that("prior/initial parameter", {
  skip_if_not_installed("nimble")
  skip_on_cran()

  set.seed(123)
  n <- 50
  dat <- data.frame(
    y  = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n)
  )

  prior <- prior_param(indexprior = "sphere", link = "bspline")
  init <- init_param(indexprior = "fisher", link = "bspline")



  expect_error(BayesSIM_setup(
    formula  = y ~ ., data = dat,
    prior    = prior, init     = init,
    niter    = 50, nburnin  = 10,
    thin     = 2, nchain   = 1, setSeed  = FALSE))

})

