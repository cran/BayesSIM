# tests/testthat/test-bayesim-generics.R
test_that("bsim and bsimSetup constructors create correct classes", {
  dummy <- bsim(
    coefficients = c(a = 1, b = 2),
    ses          = c(0.1, 0.2),
    residuals    = c(0.5, -0.5),
    fitted.values = c(1.5, 1.5),
    linear.predictors = c(0.2, 0.3),
    gof          = 0.25,
    input        = list(),
    samples      = NULL,
    model        = NULL,
    sampler      = NULL,
    modelName    = "bsFisher"
  )

  expect_s3_class(dummy, "bsim")
  expect_equal(dummy$modelName, "bsFisher")
  expect_true("coefficients" %in% names(dummy))

  setup <- bsimSetup(
    coefficients = NULL,
    ses          = NULL,
    residuals    = NULL,
    fitted.values = NULL,
    linear.predictors = NULL,
    gof          = NULL,
    input        = list(foo = 1),
    samples      = NULL,
    model        = "model",
    sampler      = "sampler",
    modelName    = "bsFisher"
  )

  expect_s3_class(setup, "bsimSetup")
  expect_equal(setup$input$foo, 1)
  expect_equal(setup$modelName, "bsFisher")
})

test_that("bsim.predict constructor builds bsimPred object", {
  obj <- bsim.predict(
    fitted   = data.frame(mean = 1:3),
    truey    = 1:3,
    idxValue = c(0.1, 0.2, 0.3),
    level    = 0.95
  )

  expect_s3_class(obj, "bsimPred")
  expect_equal(obj$level, 0.95)
  expect_equal(length(obj$idxValue), 3)
})

test_that("as_bsim combines setup and samples using sampleBind and bsimFit_pointest", {
  skip_if_not_installed("coda")

  # dummy setup object (bsimSetup-like)
  X <- matrix(rnorm(4), ncol = 2)
  colnames(X) <- c("x1", "x2")
  Y <- matrix(rnorm(2), ncol = 1)

  setup <- list(
    input = list(
      origdata = list(x = X, y = Y),
      samplingOptions = list(nchain = 1),
      formula = c("y", "~", "x1 + x2")
    ),
    defModel   = "model_stub",
    defSampler = "sampler_stub",
    modelName  = "bsFisher"
  )
  class(setup) <- "bsimSetup"

  # mcmc.out with "samples" element
  samp_mat <- cbind(
    "index[1]" = rnorm(10),
    "index[2]" = rnorm(10),
    "sigma2"   = rexp(10)
  )
  mcmc_samples <- coda::mcmc(samp_mat)
  mcmc_out <- list(samples = mcmc_samples)

  env_as <- environment(as_bsim)

  testthat::local_mocked_bindings(
    sampleBind = function(x, nchain) x,
    bsimFit_pointest = function(samples, X, Y) {
      list(
        coefficients       = c(index1 = 0.5, index2 = -0.5),
        ses_coef           = c(0.1, 0.2),
        se                 = c(0.1, 0.2),
        residuals          = as.numeric(Y) - mean(as.numeric(Y)),
        fitted.values      = rep(mean(as.numeric(Y)), length(Y)),
        linear.predictors  = X %*% c(0.5, -0.5),
        gof                = 0.123
      )
    },
    .env = env_as
  )

  fit <- as_bsim(setup, mcmc_out)

  expect_s3_class(fit, "bsim")
  expect_equal(fit$modelName, "bsFisher")
  expect_true(all(c("coefficients", "residuals", "fitted.values") %in% names(fit)))
  expect_true(inherits(fit$samples, "mcmc"))
})

test_that("print.bsim prints header and coefficients when samples exist", {
  skip_if_not_installed("coda")

  X <- matrix(rnorm(4), ncol = 2)
  colnames(X) <- c("x1", "x2")
  Y <- matrix(rnorm(2), ncol = 1)

  obj <- list(
    coefficients = c(index1 = 0.5, index2 = -0.5),
    ses_coef     = c(0.1, 0.2),
    residuals    = c(0.1, -0.1),
    fitted.values = c(1, 2),
    linear.predictors = c(0.3, 0.4),
    gof          = 0.25,
    input = list(
      origdata = list(x = X, y = Y),
      formula  = c("y", "~", "x1 + x2")
    ),
    samples   = coda::mcmc(matrix(rnorm(20), ncol = 2)),
    modelName = "bsFisher"
  )
  class(obj) <- "bsim"

  expect_output(
    print(obj),
    "BayesSIM model"
  )
  expect_output(
    print(obj),
    "Single-index estimates:"
  )
})

test_that("print.bsim prints message when samples are NULL", {
  X <- matrix(rnorm(4), ncol = 2)
  colnames(X) <- c("x1", "x2")
  Y <- matrix(rnorm(2), ncol = 1)

  obj <- list(
    coefficients = c(index1 = 0.5, index2 = -0.5),
    ses_coef     = c(0.1, 0.2),
    residuals    = c(0.1, -0.1),
    fitted.values = c(1, 2),
    linear.predictors = c(0.3, 0.4),
    gof          = 0.25,
    input = list(
      origdata = list(x = X, y = Y),
      formula  = c("y", "~", "x1 + x2")
    ),
    samples   = NULL,
    modelName = "bsFisher"
  )
  class(obj) <- "bsim"

  expect_output(
    print(obj),
    "Try compiling the object and draw posterior samples"
  )
})

test_that("coef.bsim returns mean coefficients and SE when requested", {
  obj <- list(
    coefficients = c(index1 = 0.5, index2 = -0.5),
    ses_coef     = c(0.1, 0.2),
    input = list(
      origdata = list(x = matrix(rnorm(4), ncol = 2)),
      samplingOptions = list(nchain = 1)
    ),
    samples = NULL
  )
  class(obj) <- "bsim"

  # method = "mean"
  co <- coef.bsim(obj, method = "mean", se = FALSE)
  expect_equal(co, obj$coefficients)

  co_se <- coef.bsim(obj, method = "mean", se = TRUE)
  expect_true(is.matrix(co_se))
  expect_equal(rownames(co_se), c("est", "Std.error"))
})

test_that("coef.bsim with median uses sampleBind and MAD", {
  skip_if_not_installed("coda")

  X <- matrix(rnorm(4), ncol = 2)
  colnames(X) <- c("x1", "x2")

  samp_mat <- cbind(
    "index[1]" = rnorm(20, 1),
    "index[2]" = rnorm(20, -1),
    "sigma2"   = rexp(20)
  )
  samples <- coda::mcmc(samp_mat)

  obj <- list(
    coefficients = c(index1 = 0.5, index2 = -0.5),
    ses_coef     = c(0.1, 0.2),
    input = list(
      origdata = list(x = X),
      samplingOptions = list(nchain = 1)
    ),
    samples = samples
  )
  class(obj) <- "bsim"

  env_coef <- environment(coef.bsim)
  testthat::local_mocked_bindings(
    sampleBind = function(x, nchain) x,
    .env = env_coef
  )

  med <- coef.bsim(obj, method = "median", se = FALSE)
  expect_equal(length(med), 2)
  expect_equal(names(med), colnames(X))
})


test_that("residuals.bsim returns mean and median residuals", {
  skip_if_not_installed("coda")

  X <- matrix(rnorm(4), ncol = 2)
  colnames(X) <- c("x1", "x2")
  Y <- matrix(rnorm(2), ncol = 1)

  samp_mat <- cbind(
    "linkFunction[1, 1]" = rnorm(20),
    "linkFunction[2, 1]" = rnorm(20)
  )
  samples <- coda::mcmc(samp_mat)

  obj <- list(
    residuals = c(0.1, -0.1),
    input = list(
      origdata = list(x = X, y = Y),
      samplingOptions = list(nchain = 1)
    ),
    samples = samples
  )
  class(obj) <- "bsim"

  env_res <- environment(residuals.bsim)
  testthat::local_mocked_bindings(
    sampleBind = function(x, nchain) x,
    .env = env_res
  )

  # mean
  r_mean <- residuals(obj, method = "mean")
  expect_equal(r_mean, obj$residuals)

  # median
  r_med <- residuals(obj, method = "median")
  expect_length(r_med, 2)
})

test_that("fitted.bsim works for response/linpred and mean/median", {
  skip_if_not_installed("coda")

  X <- matrix(rnorm(4), ncol = 2)
  colnames(X) <- c("x1", "x2")
  Y <- matrix(rnorm(2), ncol = 1)

  samp_mat <- cbind(
    "linkFunction[1, 1]" = rnorm(20),
    "linkFunction[2, 1]" = rnorm(20),
    "Xlin[1]"            = rnorm(20),
    "Xlin[2]"            = rnorm(20)
  )
  samples <- coda::mcmc(samp_mat)

  obj <- list(
    fitted.values     = c(1, 2),
    linear.predictors = c(0.3, 0.4),
    input = list(
      origdata = list(x = X, y = Y),
      samplingOptions = list(nchain = 1)
    ),
    samples   = samples
  )
  class(obj) <- "bsim"

  env_fit <- environment(fitted.bsim)
  testthat::local_mocked_bindings(
    sampleBind = function(x, nchain) x,
    .env = env_fit
  )

  # response, mean
  fr <- fitted(obj, type = "response", method = "mean")
  expect_equal(fr, obj$fitted.values)

  # response, median
  fr_med <- fitted(obj, type = "response", method = "median")
  expect_length(fr_med, 2)

  # linpred, mean
  fl <- fitted(obj, type = "linpred", method = "mean")
  expect_length(fl, 2)

  # linpred, median
  fl_med <- fitted(obj, type = "linpred", method = "median")
  expect_length(fl_med, 2)
})

test_that("summary.bsim produces summary.bsim object and print works", {
  skip_if_not_installed("coda")

  X <- matrix(rnorm(4), ncol = 2)
  colnames(X) <- c("x1", "x2")
  Y <- matrix(rnorm(2), ncol = 1)

  samp_mat <- cbind(
    "index[1]" = rnorm(40),
    "index[2]" = rnorm(40),
    "sigma2"   = rexp(40)
  )
  samples <- coda::mcmc(samp_mat)

  obj <- list(
    input = list(
      origdata = list(x = X, y = Y),
      samplingOptions = list(nchain = 1, monitors = NULL)
    ),
    samples   = samples,
    modelName = "bsFisher"
  )
  class(obj) <- "bsim"

  s <- summary(obj)
  expect_s3_class(s, "summary.bsim")
  expect_true("all.chain" %in% names(s))
  expect_true(is.data.frame(s$all.chain))

  expect_output(
    print(s),
    "Summary statistics of all chains"
  )
})

test_that("nimTraceplot returns a ggplot object", {
  skip_if_not_installed("coda")
  skip_if_not_installed("ggplot2")

  X <- matrix(rnorm(4), ncol = 2)
  colnames(X) <- c("x1", "x2")
  Y <- matrix(rnorm(2), ncol = 1)

  samp_mat <- cbind(
    "index[1]" = rnorm(40),
    "index[2]" = rnorm(40),
    "sigma2"   = rexp(40)
  )
  samples <- coda::mcmc(samp_mat)

  obj <- list(
    input = list(
      origdata = list(x = X, y = Y),
      samplingOptions = list(nchain = 1, monitors = NULL)
    ),
    samples   = samples,
    modelName = "bsFisher"
  )
  class(obj) <- "bsim"

  p <- nimTraceplot(obj)
  expect_s3_class(p, "ggplot")
})

test_that("plot.bsim calls predict.bsim and plot_bsim_fitted with correct interval flag", {
  skip_if_not_installed("ggplot2")

  X <- matrix(rnorm(4), ncol = 2)
  colnames(X) <- c("x1", "x2")
  Y <- matrix(rnorm(2), ncol = 1)

  obj <- list(
    input = list(
      origdata = list(x = X, y = Y),
      samplingOptions = list(nchain = 1, monitors = NULL),
      formula = c("y", "~", "x1 + x2")
    ),
    samples   = NULL,
    modelName = "bsFisher"
  )
  class(obj) <- "bsim"

  env_plot <- environment(plot.bsim)

  testthat::local_mocked_bindings(
    predict.bsim = function(object, ...) {
      bsim.predict(
        fitted   = data.frame(mean = c(1, 2)),
        truey    = c(1, 2),
        idxValue = c(0.1, 0.2),
        level    = 0.95
      )
    },
    plot_bsim_fitted = function(pred, interval) {
      list(pred = pred, interval = interval)
    },
    .env = env_plot
  )

  res1 <- plot(obj, method = "mean", interval = TRUE)
  expect_true(is.list(res1))
  expect_true(res1$interval)

  res2 <- plot(obj, method = "mean", interval = FALSE)
  expect_false(res2$interval)
})

test_that("predict.bsim with type = 'index' returns bsimPred object", {
  skip_if_not_installed("coda")
  skip_if_not_installed("nimble")

  X <- matrix(rnorm(4), ncol = 2)
  colnames(X) <- c("x1", "x2")
  Y <- matrix(rnorm(2), ncol = 1)

  samp_mat <- cbind(
    "index[1]" = rnorm(20),
    "index[2]" = rnorm(20),
    "Xlin[1]"  = rnorm(20),
    "Xlin[2]"  = rnorm(20),
    "sigma2"   = rexp(20)
  )
  samples <- coda::mcmc(samp_mat)

  obj <- list(
    input = list(
      origdata = list(x = X, y = Y),
      samplingOptions = list(nchain = 1, monitors = NULL),
      formula = c("y", "~", "x1 + x2"),
      prior   = list()
    ),
    samples   = samples,
    modelName = "bsFisher"
  )
  class(obj) <- "bsim"

  env_pred <- environment(predict.bsim)
  testthat::local_mocked_bindings(
    sampleBind = function(x, nchain) x,
    .env = env_pred
  )

  res <- predict.bsim(obj, type = "index", method = "mean",
                      interval = "none", level = 0.9)
  expect_s3_class(res, "bsimPred")
  expect_equal(length(res$idxValue), nrow(X))
})

# test_that("predict.bsim validates arguments", {
#   skip_if_not_installed("nimble")
#   skip_if_not_installed("coda")
#
#   X <- matrix(rnorm(4), ncol = 2)
#   colnames(X) <- c("x1", "x2")
#   Y <- matrix(rnorm(2), ncol = 1)
#
#   samp_mat <- cbind(
#     "index[1]" = rnorm(20),
#     "index[2]" = rnorm(20),
#     "Xlin[1]"  = rnorm(20),
#     "Xlin[2]"  = rnorm(20),
#     "sigma2"   = rexp(20)
#   )
#   samples <- coda::mcmc(samp_mat)
#
#   obj <- list(
#     input = list(
#       origdata = list(x = X, y = Y),
#       samplingOptions = list(nchain = 1, monitors = NULL),
#       formula = c("y", "~", "x1 + x2")
#     ),
#     samples   = samples,
#     modelName = "bsFisher"
#   )
#   class(obj) <- "bsim"
#
#   env_pred <- environment(predict.bsim)
#   testthat::local_mocked_bindings(
#     sampleBind = function(x, nchain) x,
#     .env = env_pred
#   )
#
#   expect_error(
#     predict.bsim(obj, newdata = matrix(1, 2, 2)),
#     "newdata should be data.frame."
#   )
# })

test_that("getVarMonitor returns correct name/list from expand_monitor", {
  obj <- list(
    defModel = list(
      getNodeNames = function(stochOnly = TRUE, includeData = FALSE) {
        c("index[1]", "index[2]", "sigma2")
      }
    )
  )

  env_getVar <- environment(getVarMonitor)
  testthat::local_mocked_bindings(
    expand_monitor = function(possibleList) {
      list(
        baseName = paste0(possibleList, "_base"),
        out      = paste0(possibleList, "_out")
      )
    },
    .env = env_getVar
  )

  # type = "name"
  nm <- getVarMonitor(obj, type = "name")
  expect_equal(nm, c("index[1]_base", "index[2]_base", "sigma2_base"))

  # type = "list"
  ls <- getVarMonitor(obj, type = "list")
  expect_equal(ls, c("index[1]_out", "index[2]_out", "sigma2_out"))
})

test_that("getModelDef returns BUGScode", {
  obj <- list(
    defModel = list(
      modelDef = list(
        BUGScode = "model { y ~ dnorm(0,1) }"
      )
    )
  )
  expect_equal(getModelDef(obj), "model { y ~ dnorm(0,1) }")
})

test_that("getInit returns initial list", {
  obj <- list(input = list(init = list(a = 1, b = 2)))
  expect_equal(getInit(obj), list(a = 1, b = 2))
})

test_that("compileModelAndMCMC compiles nimble model and sampler", {
  skip_if_not_installed("nimble")
  skip_on_cran()

  code <- nimble::nimbleCode({
    for (i in 1:N) {
      y[i] ~ dnorm(mu, 1)
    }
    mu ~ dnorm(0, 1)
  })

  const <- list(N = 3)
  dat   <- list(y = rnorm(3))
  inits <- list(mu = 0)

  model <- nimble::nimbleModel(code, constants = const, data = dat, inits = inits)
  conf  <- nimble::configureMCMC(model)
  mcmc  <- nimble::buildMCMC(conf)

  fullmodel <- list(defModel = model, defSampler = mcmc)
  res <- compileModelAndMCMC(fullmodel)

  expect_true("model" %in% names(res))
  expect_true("sampler" %in% names(res))
})

test_that("get_model and get_sampler extract compiled objects", {
  obj <- list(
    model   = "compiled_model_stub",
    sampler = "compiled_sampler_stub"
  )

  expect_equal(get_model(obj), "compiled_model_stub")
  expect_equal(get_sampler(obj), "compiled_sampler_stub")
})
