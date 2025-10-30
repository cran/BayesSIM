#' Bayesian single-index regression with B-spline link and half-unit sphere prior
#'
#'
#' @description
#' Fits a single-index model \eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n}
#' where the link \eqn{f(\cdot)} is represented by B-spline link and the
#' index vector \eqn{\theta} is on half-unit sphere.
#'
#' @param x Numeric data.frame/matrix of predictors. Each row is an observation.
#' @param y Numeric response vector/matrix.
#' @param prior Optional named list of prior settings with sublists:
#' \describe{
#' \item{\code{index}}{ \code{nu} is binary inclusion indicators prior for variable selection in the index:
#'     \code{list(r1, r2)} gives the Beta hyperprior on the Bernoulli–inclusion probability \eqn{w},
#'     inducing \eqn{p(\nu) \propto \mathrm{Beta}(r_1 + n_\nu, r_2 + p - n_\nu)} (default \code{r1 = 1, r2 = 2}).}
#' \item{\code{link}}{ B-spline knots, basis and coefficient setup.
#' \enumerate{
#'   \item{\code{knots} Free-knot prior for the spline. \code{lambda_k} is the Poisson mean for the number of
#'     interior knot \eqn{k} (default \code{5}). \code{maxknots} is the maximum number of interior knots.
#'     If \code{maxknots} is \code{NULL}, the number of interior knots is randomly drawn from a Poisson distribution.}
#'   \item{\code{basis} For the basis of B-spline,  \code{degree} is the spline
#'     degree (default \code{2}).}
#'   \item{\code{beta} For the coefficient of B-spline, conjugate normal prior on \eqn{\beta} with covariance \eqn{\tau \Sigma_0} is assigned.
#'    By default, \code{mu} is a zero vector, \code{tau} is set to the sample size,
#'     and \code{Sigma0} is the identity matrix of dimension \eqn{1 + k + m},
#'     where \eqn{k} is the number of interior knots and \eqn{m} is the spline order (degree + 1).}}
#' }
#'
#' \item{\code{sigma2}}{Error-variance prior hyperparameters. An Inverse-Gamma prior is assigned to \eqn{\sigma^2}
#'  where \code{shape} is shape parameter and \code{rate} is rate parameter of inverse gamma distribution. (default \code{shape = 0.001, rate = 0.001}).
#'     Small values for shape and rate parameters yield a weakly-informative prior on \eqn{\sigma^{2}}.}
#' }
#'
#' @param init Optional named list of initial values. If the values are not assigned, they are randomly sampled from prior.
#' \describe{
#' \item{\code{index}}{\code{nu} is binary vector indicating active predictors for the index.
#' \code{index} is initial unit-norm index vector \eqn{\theta} (automatically normalized if necessary, with the first nonzero element made positive for identifiability).
#'     The vector length must match the number of columns in data \code{x}.
#'     Ideally, positions where \code{nu} has a value of 1 should correspond to nonzero elements in \eqn{\theta}; elements corresponding to \code{nu} = 0 will be set to zero.}
#'  \item{\code{link}}{\code{k} is initial number of interior knots. \code{knots} is initial vector of interior knot positions in \eqn{[0, 1]}, automatically scaled to the true boundary.
#'     Length of this vector should be equal to the initial value of \code{k}.
#'     \code{beta} is initial vector of spline coefficients. Length should be equal to the initial number of basis functions with intercept (\eqn{1 + k + m}).}
#'   \item{\code{sigma2}}{Initial scalar error variance. (default \code{0.01})}
#' }
#'
#' @param sampling Logical. If \code{TRUE} (default), run MCMC; otherwise return prepared nimble model objects without sampling.
#' @param fitted Logical. If \code{TRUE} (default), fitted values drawn from posterior distribution are included in the output and
#' \code{c("linkFunction", "beta", "k", "knots", "numBasis", "a_alpha", "b_alpha", "Xlin")} is monitored for prediction.
#' @param monitors2 Optional character vector of additional monitor nodes. To check the names of the nodes, set \code{fit <- bsSphere(x, y, sampling = FALSE)} and then inspect the variable names stored in the model object using \code{fit$model$getVarNames()}.
#' @param niter Integer. Total MCMC iterations (default \code{10000}).
#' @param nburnin Integer. Burn-in iterations (default \code{1000}).
#' @param thin Integer. Thinning for monitors1 (default \code{1}).
#' @param thin2 Integer. Optional thinning for \code{monitors2} (default \code{1}).
#' @param nchain Integer. Number of MCMC chains (default \code{1}).
#' @param setSeed Logical or numeric argument.  Further details are provided in \link[nimble]{runMCMC}.
#'
#' @details
#' \strong{Model} The single–index model uses a \eqn{m}-order polynomial spline with \eqn{k} interior knots and intercept as follows:
#' \eqn{f(t) = \sum_{j=1}^{1+m+k} B_j(t)\,\beta_j} on \eqn{[a, b]} with \eqn{t_i = x_i' \theta, i = 1,\cdots, n}
#' and \eqn{\|\theta\|_2 = 1}. \eqn{\{\beta_j\}_{j=1}^{m+k+1}} are spline coefficient and \eqn{a_\theta} and \eqn{ b_\theta} are boundary knots where \eqn{a = min(t_i, i = 1, \cdots, n)},
#' and \eqn{b = max(t_i, i = 1,\cdots, n)}. Variable selection is encoded by a binary vector \eqn{\nu}, equivalently
#' setting components of \eqn{\theta} to zero when \eqn{\nu_i = 0}.
#'
#' \strong{Priors}
#' \itemize{
#'   \item Free‑knot prior: \eqn{k \sim \mathrm{Poisson}(\lambda_k)}, knot locations \eqn{\xi_i, i = 1,...,k} via a Dirichlet prior on the scaled interval.
#'   \item Beta–Bernoulli hierarchy for \eqn{\nu}, yielding a Beta–Binomial prior.
#'   \item Spherical prior on the index: uniform on the half‑sphere of dimension \eqn{n_\nu}with first nonzero positive.
#'   \item Conjugate normal–inverse-gamma on \eqn{(\beta, \sigma^2)} enables analytic integration and a lower‑dimensional marginal target for RJMCMC.
#' }
#'
#' \strong{Sampling} Posterior exploration follows a hybrid RJMCMC with six move types:
#' add/remove predictor \eqn{\nu}, update \eqn{\theta}, add/delete/relocate a knot. The \eqn{\theta} update is a random‑walk
#' Metropolis via local rotations in a two‑coordinate subspace; knot moves use simple proposals with tractable acceptance ratios.
#' Further sampling method is described in Wang(2009).
#'
#' @return A \code{list} typically containing:
#' \describe{
#'   \item{\code{model}}{Nimble model}
#'   \item{\code{sampler}}{Nimble sampler}
#'   \item{\code{sampling}}{Posterior draws of \eqn{\nu}, \eqn{\theta}, \eqn{\sigma^2}, and nodes for fitted values by default. Variables specified in \code{monitors2} will be added if provided.}
#'   \item{\code{fitted}}{If \code{fitted = TRUE}, in-sample fitted values is given.}
#'   \item{\code{input}}{List of data,input values for prior and initial values, and computation time without compiling.}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 200; d <- 4
#' theta <- c(2, 1, 1, 1); theta <- theta / sqrt(sum(theta^2))
#' f <- function(u) u^2 * exp(u)
#' sigma <- 0.5
#' X <- matrix(runif(n * d, -1, 1), nrow = n)
#' index_vals <- as.vector(X %*% theta)
#' y <- f(index_vals) + rnorm(n, 0, sigma)
#'
#' # One-call version
#' fit <- bsSphere(X, y)
#'
#' # Split version
#' models <- bsSphere(X, y, sampling = FALSE)
#' Ccompile <- compileModelAndMCMC(models)
#' mcmc.out <- runMCMC(Ccompile$mcmc, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, init = models$input$init,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#' }
#'
#' @references
#' Wang, H.-B. (2009). Bayesian estimation and variable selection for single index models.
#' \emph{Computational Statistics & Data Analysis}, 53, 2617–2627.
#'
#' Hornik, K., & Grün, B. (2014). movMF: an R package for fitting mixtures of von Mises-Fisher distributions.
#' \emph{Journal of Statistical Software}, 58, 1-31.
#'
#'
#' @export

bsSphere <- function(x, y,
                     prior = list(
                       index = list(nu = list(r1 = 1, r2 = 1)),
                       link = list(knots = list(lambda_k = 5, maxknots = NULL),
                                   basis = list(degree = 2),
                                   beta = list(mu = NULL, tau = NULL, Sigma0 = NULL)),
                       sigma2 = list(shape = 0.001, rate = 0.001)),
                     init = list(index = list(nu = NULL, index = NULL),
                                 link = list(k = NULL, knots = NULL, beta = NULL),
                                 sigma2 = 0.01),
                     sampling = TRUE, fitted = TRUE,
                     monitors2 = NULL, niter = 10000, nburnin=1000,
                     thin = 1, thin2 = NULL, nchain = 1, setSeed = FALSE
){

  start1 <- Sys.time()
  index <- 0; k <- 0; knots <- 0; sigma2 <- 0

  # check sampling, prior, init parameters for independent execution
  checkOutput <- validate_and_finalize_args(
    sampling, fitted, niter, nburnin, thin, thin2, nchain,
    prior, init, "sphere", "bspline"
  )
  prior <- checkOutput$priorlist_final
  init <- checkOutput$initlist_final

  # environment
  envobj <- ls(envir=.GlobalEnv)
  on.exit(rm(list=ls(envir=.GlobalEnv)[which(!ls(envir=.GlobalEnv)%in%envobj)],envir=.GlobalEnv))
  .fns <- c(
    # a_common
    "quickSortOrderIndexOnly", "nimOrder", "nimSort",
    "sampleQuantile_nim", "quantile_nimble",

    # aa_bspline_ver3
    "SplineState", "any_duplicated", "mat_wo_col1", "update_spline_df", "update_x_index",
    "update_knot_sequence", "get_basis_simple", "simplify_knots", "get_inside_x",
    "gen_default_internal_knots", "SplineBase1", "SplineBase2", "basis", "bsNimble", "bsBasis",

    # bsSphere
    "transX_sp","estBetaInit_sp","computeSig","ComputeS","logdet_nim",
    "postll_bspline_sphere","postll_knots","betaFunction","computeA1_nu1",
    "computeA1_nu0","computeA2_add","computeA2_delete","rnormTrun","newKnots",
    "prop_add","prop_delete","rtruncnorm","pred_bsplineSphere",
    "nuSampler_bspline_sphere","indexSampler_bspline_sphere","knotsSampler_bspline_sphere",
    "betaSampler_bspline_sphere","sigma2Sampler_bspline_sphere",

    # utils
    "pred_fitted"
  )

  pkg <- "BayesSIM"
  ns <- asNamespace(pkg)
  list2env(mget(.fns, envir = ns, inherits = FALSE), envir = globalenv())



  # check data dimension
  if (!is.matrix(x) & !is.data.frame(x)){stop("x is not matrix/data.frame.")}
  if (!is.vector(y) & !is.matrix(y)){stop("y is not vector or matrix.")}
  if (is.matrix(y)){
    if ((ncol(y) != 1)){
      stop("y should be scalar vector or matrix.")
    }
  }

  X <- as.matrix(x)
  Y <- matrix(y, ncol = 1)

  if (nrow(X) != nrow(Y)){
    stop("x and y have different dimension.")
  }

  # data dimension
  N <- length(Y)
  p <- ncol(X)


  # Model code
  Rmodelcode <- nimbleCode({
    for (i in 1:p) {
      nu[i] ~ dbern(r1/(r1+r2))
    }

    # index
    index[1:p] ~ dunitSphere(p) ## half unit sphere

    # Linear predictor
    for (i in 1:N){
      Xlin[i] <- sum(X[i, 1:p] * index[1:p])
    }

    # boundary knots
    a_alpha <- min(Xlin[1:N])
    b_alpha <- max(Xlin[1:N])

    # knots
    k ~ dpois(lambda_k)
    knots[1:maxknots] ~ dKnotsSimple(a = a_alpha, b = b_alpha, k = k, alpha = alpha[1:maxknots])

    # Design matrix - b spline basis
    numBasis <- degree + k + 1 # real dimension
    Xmat[1:N, 1:maxBasis] <- transX_sp(Xlin[1:N], degree, knots[1:maxknots], k, maxBasis, a_alpha, b_alpha)

    # likelihood
    sigma2 ~ dinvgamma(a_sig, b_sig)
    covBeta[1:maxBasis, 1:maxBasis] <- (sigma2 * tau) * Sigma0[1:maxBasis, 1:maxBasis]
    beta[1:maxBasis, 1] ~ dmnorm(mubeta[1:maxBasis],
                                 cov = covBeta[1:maxBasis, 1:maxBasis])

    linkFunction[1:N, 1] <- Xmat[1:N, 1:maxBasis] %*% beta[1:maxBasis, 1]
    covY[1:N, 1:N] <- sigma2 * Sigma[1:N, 1:N]
    Y[1:N, 1] ~ dmnorm(linkFunction[1:N, 1], cov = covY[1:N, 1:N])

  })


  # Prior parameters
  ## check data dimension and save
  ## index - nu
  if (is.null(prior$index$nu$r1) || is.null(prior$index$nu$r2)){
    stop("nu prior is not available.")
  }
  if (prior$index$nu$r1 <= 0 || prior$index$nu$r2 <= 0){
    stop("nu prior should be positive.")
  }

  r1 <- prior$index$nu$r1
  r2 <- prior$index$nu$r2

  ## link - knots
  if (is.null(prior$link$knots$lambda_k) || is.null(prior$link$knots$lambda_k)){
    stop("knots prior is not available.")
  }
  if (prior$link$knots$lambda_k <= 0){
    stop("knots prior should be positive.")
  }

  lambda_k <- prior$link$knots$lambda_k

  maxknots <- ifelse(is.null(prior$link$knots$maxknots),
                     qpois(1-0.01, lambda = lambda_k), prior$link$knots$maxknots)

  if (length(maxknots) >= 2 || maxknots < 0){
    stop("Prior link (maxknots) has incorrect value.")
  }

  ## link - basis
 if (is.null(prior$link$basis$degree)||length(prior$link$basis$degree) >= 2 || prior$link$basis$degree< 0){
    stop("Prior basis (degree) has incorrect value.")
  } else{
    degree <- prior$link$basis$degree
  }

  alpha <- rep(1/maxknots, maxknots)
  maxBasis <- degree + 1 + maxknots


  ## beta
  ## tau
  tau <- ifelse(is.null(prior$link$beta$tau), N, prior$link$beta$tau)
  if (length(tau) >= 2 || tau < 0){
    stop("Prior beta (tau) has incorrect value.")
  }

  ## mu
  if (!is.vector(prior$link$beta$mu) & !is.null(prior$link$beta$mu)){
    stop("Prior beta (mu) should be vector.")
  }
  if (!is.null(prior$link$beta$mu) & length(prior$link$beta$mu) < maxBasis ){
    stop("Incorrect dimention on prior beta (mu).")
  }

  if (is.null(prior$link$beta$mu)){
    mubeta <- rep(0, maxBasis)
  } else{
    mubeta <- c(prior$link$beta$mu, rep(0, maxBasis - length(prior$link$beta$mu)))
  }



  ## Sigma0

  if (is.null(prior$link$beta$Sigma0)){
    Sigma0 <- diag(1, maxBasis)
  } else{
    Sigma0 <- prior$link$beta$Sigma0
  }
  if (!is.null(prior$link$beta$Sigma0) & !is.matrix(prior$link$beta$Sigma0)){
    stop("Prior beta (Sigma0) should be matrix.")
  }

  Sigma <- diag(1, N)


  ## sigma
  if (is.null(prior$sigma2$shape)||length(prior$sigma2$shape) >= 2 || prior$sigma2$shape < 0){
    stop("Prior sigma2 (shape) has incorrect value.")
  } else{
    a_sig <- prior$sigma2$shape
  }

  if (is.null(prior$sigma2$rate)||length(prior$sigma2$rate) >= 2 || prior$sigma2$rate < 0){
    stop("Prior sigma2 (rate) has incorrect value.")
  } else{
    b_sig <- prior$sigma2$rate
  }


  # Initialize
  ## nu
  if (!is.null(init$index$nu) & (length(init$index$nu) != p)){
    stop("Incorrect dimension on initial value of nu.")
  }

  if (!is.null(init$index$nu) & (!all(init$index$nu %in% c(0, 1)))){
    stop("nu should consist of 0 or 1.")
  }

  ## index
  if (!is.null(init$index$index) & (length(init$index$index) != p)){
    stop("Incorrect dimension on initial value of index")
  }

  # seed
  seedNum <- rep(FALSE, nchain)
  if (!is.logical(setSeed) & !is.numeric(setSeed)){
    stop("'setSeed' argument should be logical or numeric vector.")
  }
  if (is.logical(setSeed) & (setSeed == TRUE)){
    seedNum <- seq(1, nchain, 1)
  }
  if (is.numeric(setSeed)){
    if (length(setSeed) == nchain){
      seedNum <- setSeed
    } else if(length(setSeed) !=  nchain){
      stop("The length of 'setSeed' should be equal to the number of chain.")
    }
  }





  inits_list <- lapply(seq_len(nchain),
                       function(j) initFunction_bS(nu = init$index$nu, index = init$index$index,
                                                   k = init$link$k, beta = init$link$beta,
                                                   sigma2 = init$sigma2, knots = init$link$knots,
                                                   p = p, X = X, Y = Y, lambda_k = lambda_k,
                                                   maxknots = maxknots, degree = degree,
                                                   maxBasis = maxBasis, tau = tau,
                                                   Sigma0 = Sigma0,setSeed = seedNum[j]))
  firstInit <- inits_list[[1]]


  # Build model
  message("Build Model")
  suppressMessages(simpleModel <- nimbleModel(Rmodelcode,
                             data = list(X = X, Y = Y),
                             constants = list(p = p, N = N, r1 = r1, r2 = r2,
                                              lambda_k = lambda_k,maxknots = maxknots,
                                              a_sig = a_sig, b_sig = b_sig, alpha= alpha,
                                              tau = tau, degree = degree, maxBasis = maxBasis,
                                              mubeta = mubeta, Sigma0 = Sigma0,
                                              Sigma = Sigma),
                             inits = firstInit))

  # Assign samplers
  message("Assign samplers")
  monitorsList <- c("nu" ,"index", "sigma2")
  if (fitted){
    monitorsList <- c(monitorsList, "linkFunction", "beta",
                      "k", "knots", "numBasis", "a_alpha", "b_alpha", "Xlin")
  }
  if (is.null(monitors2)){
    suppressMessages(mcmcConf <- configureMCMC(simpleModel,
                              monitors = monitorsList,
                              print = FALSE))
  } else{
    suppressMessages(mcmcConf <- configureMCMC(simpleModel,
                              monitors = monitorsList, monitors2 = monitors2,
                              print = FALSE))
  }


  # mcmcConf1$printSamplers(executionOrder= TRUE)
  mcmcConf$removeSamplers(c("nu"))
  mcmcConf$addSampler(target = c("nu"),
                       type   = nuSampler_bspline_sphere)

  mcmcConf$removeSamplers(c("index"))
  mcmcConf$addSampler(target = c("index"),
                       type   = indexSampler_bspline_sphere)

  mcmcConf$removeSamplers(c("k", "knots"))
  mcmcConf$addSampler(target = c("k","knots"),
                       type   = knotsSampler_bspline_sphere)

  mcmcConf$removeSamplers(c("beta"))
  mcmcConf$addSampler(target = c("beta"),
                       type   = betaSampler_bspline_sphere)

  mcmcConf$removeSamplers(c("sigma2"))
  mcmcConf$addSampler(target = c("sigma2"),
                       type   = sigma2Sampler_bspline_sphere)


  mcmc1 <- buildMCMC(mcmcConf)
  end1 <- Sys.time()


  if (!sampling){ # not mcmc sampling
    mcmc.out <- NULL
    fittedResult <- NULL
    sampMCMC <- NULL

  } else{
    # Compile
    message("Compile Model")
    suppressMessages(CsimpleModel <- compileNimble(simpleModel))
    message("Compile MCMC")
    suppressMessages(Cmcmc <- compileNimble(mcmc1,
                                            project = simpleModel,
                                            resetFunctions = TRUE))

    # Sampling
    start2 <- Sys.time()
    message("Run MCMC")
    mcmc.out <- NULL
    if (setSeed == FALSE){
      seedNum <- setSeed
    }
    if (is.null(monitors2)){
      mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                          thin = thin,
                          nchains = nchain, setSeed = seedNum,
                          inits = inits_list,
                          summary = FALSE, samplesAsCodaMCMC = TRUE)
    } else{
      mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                          thin = thin, thin2 = thin2,
                          nchains = nchain, setSeed = seedNum,
                          inits = inits_list,
                          summary = FALSE, samplesAsCodaMCMC = TRUE)
    }
    end2 <- Sys.time()

    # output
    start3 <- Sys.time()
    ## combine all chaines
    samples <- NULL
    sampMCMC <- mcmc.out
    # if (enableWAIC){
    #   sampMCMC <- mcmc.out$samples
    # } else{
    #   sampMCMC <- mcmc.out
    # }
    if (nchain > 1){
      for (i in 1:nchain){
        samples <- rbind(samples, mcmc.out[[i]])
      }
    } else if (nchain == 1){
      samples <- mcmc.out
    }

    if (fitted){ # posterior fitted value output (mean, median, sd)
      # fittedResult <- NULL
      message("Compute posterior fitted value")
      # namesBeta <- paste0("index", 1:p)
      namesLink <- paste0("linkFunction[", 1:N, ", 1]")
      namesSigma <- "sigma2"
      LinkFunction_samples <- samples[, namesLink]
      sigma2_samples <- samples[, namesSigma]

      message("Compile function..")
      suppressMessages(cpred_fitted <- compileNimble(pred_fitted))
      message("Computing predicted value..")
      fittedValue <- cpred_fitted(LinkFunction_samples,
                                   sigma2_samples)

      fittedResult <- fittedValue
      } else{
      fittedResult <- NULL
    }

  }

  end3 <- Sys.time()

  # inputOptions <- NULL
  if (!sampling){
    time <- NULL
  } else if (!fitted){
    samp_time <- difftime(end2, start2, units = "secs") + difftime(end1, start1, units = "secs")
    time <- list(samp = samp_time)
  } else{
    samp_time <- difftime(end2, start2, units = "secs") + difftime(end1, start1, units = "secs")
    fitted_time <- difftime(end3, start3, units = "secs")
    time <- list(samp = samp_time, fitted = fitted_time)
  }

  inputOptions <- list(data = list(x = X, y = Y),
                       prior = list(index = list(nu = list(r1 = r1, r2 = r2)),
                                    link = list(knots = list(lambda_k = lambda_k, maxknots = maxknots),
                                                basis = list(degree = degree),
                                                beta = list(mu = mubeta, tau = tau, sigma0 = Sigma0)),
                                    sigma2 = list(shape = a_sig, rate = b_sig)),
                       init = inits_list,
                       time = time)


  out <- list(model = simpleModel,
              sampler = mcmc1,
              sampling = sampMCMC,
              fitted = fittedResult,
              input = inputOptions,
              modelName = "bsSphere")
  class(out) = "bsimSpline"


  return(out)
}
