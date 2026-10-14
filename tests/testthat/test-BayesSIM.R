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
    "'setSeed' must be logical or numeric."
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

# Control global environment
# ---- Shared helpers -----------------------------------------------------
.fns <- c(
  # a_common
  "quickSortOrderIndexOnly", "nimOrder", "nimSort",
  "sampleQuantile_nim", "quantile_nimble",

  # aa_bspline_ver3
  "SplineState", "any_duplicated", "mat_wo_col1", "update_spline_df", "update_x_index",
  "update_knot_sequence", "get_basis_simple", "simplify_knots", "get_inside_x",
  "gen_default_internal_knots", "SplineBase1", "SplineBase2", "basis", "bsNimble", "bsBasis",

  # bsFisher
  "postll_bspline_fisher","nimNorm","rW","besselI_nimble","Stheta",
  "estBeta_fisher","gvcCV","transX_fisher",
  "pred_bsplineFisher","indexSampler_bspline_fisher","betaSampler_bspline_fisher",

  # bsSphere
  "transX_sp","estBetaInit_sp","computeSig","ComputeS","logdet_nim",
  "postll_bspline_sphere","postll_knots","betaFunction","computeA1_nu1",
  "computeA1_nu0","computeA2_add","computeA2_delete","rnormTrun","newKnots",
  "prop_add","prop_delete","rtruncnorm","pred_bsplineSphere",
  "nuSampler_bspline_sphere","indexSampler_bspline_sphere","knotsSampler_bspline_sphere",
  "betaSampler_bspline_sphere","sigma2Sampler_bspline_sphere",

  # gpPolar
  "alphaTheta","Xlinear","invcov","expcov_gpPolar","expcovTest_gpPolar",
  "obj_btt_theta","thetaPrior","pred_gpPolar","gibbsSampler_sigma2",
  "gibbsSampler_kappa","MH_thetaeta",

  # gpspike
  "expcov_gpSpike","expcovnn_gpSpike","computeA1","computeA2",
  "llFunLambda","llFunThetaV1","llFunThetaV2","transitionTheta",
  "multiplyMatrixByConstant","pred_gpSpike",
  "SamplingLambda_gp_spike","SamplingThetaV_gp_spike",
  "gibbsSigma2_gp_spike","gibbsGam_gp_spike",

  # gpSphere
  "expcov_gpSphere","expcovTest_gpSphere","conBeta","obj_btt",
  "obj_btt_EB","pred_gpSphere","indexSampler_gpSphere","sigma2Sampler_gpSphere",
  "optSampler",

  # utils
  "pred_fitted",
  # distribution
  "dvMFnim", "rvMFnim", "dunitSphere", "runitSphere","dKnotsSimple", "rKnotsSimple"

)

make_sentinel <- function(nm) {
  structure(list(token = paste0("SENTINEL::", nm)),
            class = "bayessim_sentinel")
}

check_globalenv <- function(expr) {
  sentinels <- setNames(lapply(.fns, make_sentinel), .fns)
  list2env(sentinels, envir = .GlobalEnv)
  before <- ls(envir = .GlobalEnv, all.names = TRUE)

  on.exit({
    after     <- ls(envir = .GlobalEnv, all.names = TRUE)
    new_names <- setdiff(after, before)
    non_nimble_new <- new_names[!grepl("^\\.__C__", new_names)]
    to_clean  <- intersect(c(.fns, non_nimble_new), after)
    if (length(to_clean)) rm(list = to_clean, envir = .GlobalEnv)
  }, add = TRUE)

  result  <- tryCatch(expr, error = function(e) e)
  errored <- inherits(result, "error")

  intact <- vapply(.fns, function(nm) {
    exists(nm, envir = .GlobalEnv, inherits = FALSE) &&
      identical(get(nm, envir = .GlobalEnv, inherits = FALSE),
                sentinels[[nm]])
  }, logical(1))

  after      <- ls(envir = .GlobalEnv, all.names = TRUE)
  all_leaked <- setdiff(after, before)
  is_nimble  <- grepl("^\\.__|^\\.requireCachedGenerics$", all_leaked)

  list(
    broken       = .fns[!intact],
    leaked       = setdiff(all_leaked[!is_nimble], .fns),
    nimble_leaks = all_leaked[is_nimble],
    errored      = errored,
    result       = result
  )
}

# Simulated dataset used by the success-path tests
make_simdata <- function(n = 200, d = 4, seed = 123) {
  set.seed(seed)
  theta <- c(2, 1, 1, 1); theta <- theta / sqrt(sum(theta^2))
  f     <- function(u) u^2 * exp(u)
  X     <- matrix(runif(n * d, -1, 1), nrow = n)
  y     <- as.vector(f(X %*% theta)) + rnorm(n, 0, 0.5)
  dat   <- data.frame(X, y = y)
  colnames(dat) <- c(paste0("X", 1:d), "y")
  dat
}

# =========================================================================
#  TEST 1 — Fitting: 8 model combinations, ERROR path
#           (BayesSIM must reject non-data-frame input without touching
#            .GlobalEnv, for every advertised link × indexprior combo)
# =========================================================================

test_that("BayesSIM() leaves .GlobalEnv clean on invalid input (all 8 combos)", {
  bad_data     <- matrix(1:8, nrow = 4)
  index_priors <- c("fisher", "sphere", "polar", "spike")
  links        <- c("bspline", "gp")

  for (link in links) {
    for (ip in index_priors) {
      combo <- sprintf("link = %s, indexprior = %s", link, ip)

      res <- check_globalenv(
        BayesSIM(y ~ x, data = bad_data,
                 indexprior = ip, link = link)
      )

      expect_true(res$errored,
                  info = sprintf("Expected an error for %s", combo))
      expect_length(res$broken, 0L)
      expect_length(res$leaked, 0L)
    }
  }
})

# =========================================================================
#  TEST 2 — Compilation: BayesSIM_setup() + compileModelAndMCMC()
# =========================================================================

test_that("compileModelAndMCMC() leaves .GlobalEnv clean", {
  # skip_on_cran()
  simdata <- make_simdata()
  models <- BayesSIM_setup(y ~ ., data = simdata)

  res <- check_globalenv({
    compileModelAndMCMC(models)
  })

  expect_false(res$errored,
               info = if (res$errored) conditionMessage(res$result) else "")
  expect_length(res$broken, 0L)
  expect_length(res$leaked, 0L)
})

# =========================================================================
#  TEST 3 — Success path: fit + predict
# =========================================================================

test_that("BayesSIM() + predict() leave .GlobalEnv clean (success path)", {
  # skip_on_cran()
  simdata <- make_simdata()
  newx    <- data.frame(X1 = rnorm(10), X2 = rnorm(10),
                        X3 = rnorm(10), X4 = rnorm(10))

  # --- 1) fit ---
  fit_res <- check_globalenv(
    BayesSIM(y ~ ., data = simdata,
             niter = 5000, nburnin = 1000, nchain = 1)
  )
  expect_false(fit_res$errored,
               info = if (fit_res$errored) conditionMessage(fit_res$result) else "")
  expect_length(fit_res$broken, 0L)
  expect_length(fit_res$leaked, 0L)

  # --- 2) predict on the fit from step 1 ---
  fit <- fit_res$result
  pred_res <- check_globalenv(
    predict(fit, newdata = newx, interval = "credible", level = 0.95)
  )
  expect_false(pred_res$errored,
               info = if (pred_res$errored) conditionMessage(pred_res$result) else "")
  expect_length(pred_res$broken, 0L)
  expect_length(pred_res$leaked, 0L)
})
