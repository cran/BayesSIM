#' Bayesian single-index regression with Gaussian process link and spike-and-slab prior
#' @description Fits a single-index model \eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n}
#' where index vector \eqn{\theta} has a spike and slab prior and
#' the link \eqn{f(\cdot)} is represented by Gaussian process and the
#'
#'
#' @param x Numeric data.frame/matrix of predictors. Each row is an observation.
#' @param y Numeric response numeric vector/matrix. Other types  are not available
#' @param prior Optional named list of prior settings with sublists:
#' \describe{
#'     \item{\code{index}}{Spike and slab prior hyperparameters: Beta-binomial for variable selection (default \code{r1 = 1, r2 = 1}),
#'     and normal distribution for selected variables (default: \eqn{N(0, \sigma_{\theta}^{2}}))}
#'     \item{\code{link}}{Gaussian process prior hyperparameters \code{lambda}: Inverse-Gamma prior is assigned for \eqn{\lambda^{-1}}
#'     (default \code{inv_lambda_shape = 1, inv_lambda_rate = 0.1})}
#'     \item{\code{sigma2}}{Error-variance prior hyperparameters. An Inverse-Gamma prior is assigned to \eqn{\sigma^2}
#'         where \code{shape} is shape parameter and \code{rate} is rate parameter of inverse gamma distribution.
#'         (default \code{shape = 0.001, rate = 100})}
#'   }
#' @param init Optional named list of initial values. If the values are not assigned, they are randomly sampled from prior.
#' \describe{
#' \item{\code{index}}{
#' \enumerate{
#'      \item{\code{pi}: Initial selecting variable probability. (default: \code{0.5})}
#'      \item{\code{nu}: Initial vector of inclusion indicators . By default, each \code{nu} is randomly drawn by  \eqn{Bernoulli(1/2)}}
#'      \item{\code{index}: Initial vector of index. By default, each element of index vector, which is chosen by nu, is proposed by normal distribution.}
#' }}
#' \item{\code{link}}{Initial scalar of lambda (\code{inv_lambda}) for covariance of Gaussian process.}
#' \item{\code{sigma2}}{Initial scalar error variance. (default: \code{0.01})}
#' }
#' @param sampling Logical. If \code{TRUE} (default), run MCMC; otherwise return prepared nimble model objects without sampling.
#' @param fitted Logical. If \code{fitted = FALSE}, fitted values are not drawn and only \code{c("nu", "indexstar", "sigma2")} are monitored.
#' If \code{fitted = TRUE} (default), fitted values drawn from posterior distribution are included in the output and \code{c("Xlin", "invlambda")} is additionally monitored for prediction.
#' @param monitors2 Optional character vector of additional monitor nodes. To check the names of the nodes, set \code{fit <- gpSpike(x, y, sampling = FALSE)} and then inspect the variable names stored in the model object using \code{fit$model$getVarNames()}.
#' @param niter Integer. Total MCMC iterations (default \code{10000}).
#' @param nburnin Integer. Burn-in iterations (default \code{1000}).
#' @param thin Integer. Thinning for monitors1 (default \code{1}).
#' @param thin2 Integer. Optional thinning for \code{monitors2} (default \code{1}).
#' @param nchain Integer. Number of MCMC chains (default \code{1}).
#' @param setSeed Logical or numeric argument.  Further details are provided in \link[nimble]{runMCMC}.
#'
#' @details
#' \strong{Model} The single–index model is specified as \eqn{Y_i = f(X_i' \theta) + \epsilon_i},
#' where \eqn{\theta} is a p-dimensional index vector subject to a spike-and-slab
#' prior for variable selection. The link function \eqn{f(\cdot)} is modeled
#' using a Gaussian process prior with zero mean and squared exponential covariance
#' kernel \eqn{K(x_1, x_2) = \exp\{-\rho {(x_1 - x_2)^{T}\theta}^2\}},
#' where \eqn{\rho} determines the smoothness of \eqn{f}.
#' The covariance kernel is re-parameterized to \eqn{\exp\{-{(x_1 - x_2)^{T}\theta^{*}}^2\}} where
#' \eqn{\rho = ||\theta^{*}||} and
#' \eqn{\theta = ||\theta||^{-1}\theta^{*}}.
#' Therefore, \eqn{\theta^{*}} is sampled in MCMC.
#'
#' \strong{Priors}
#' \itemize{
#' \item Inclusion indicators \eqn{\nu_l}: Bernoulli(\eqn{\pi}).
#'   \item Inclusion probability \eqn{\pi}: Beta(\eqn{r_1, r_2}).
#'   \item Slab coefficients \eqn{\theta_l^*}: Gaussian \eqn{N(0, \sigma_\theta^2)}.
#'   \item GP precision \eqn{\lambda^{-1}}: Gamma(\eqn{a_\lambda, b_\lambda}).
#'   \item Error precision \eqn{(\sigma^2)^{-1}}: Gamma(\eqn{a_\sigma, b_\sigma}).
#' }
#'
#' \strong{Sampling} A random walk Metropolis algorithm is used to sample \eqn{\lambda^{-1}}
#' and a Metropolis-Hastings algorithm is used for the main parameters \eqn{(\theta^{*}, \nu)}.
#' The variance \eqn{\sigma^2} is directly sampled from posterior distribution.
#' \eqn{f} is not directly sampled by MCMC.
#'
#'
#' @return A \code{list} typically containing:
#' \describe{
#'   \item{\code{model}}{Nimble model}
#'   \item{\code{sampler}}{Nimble sampler}
#'   \item{\code{sampling}}{Posterior draws of \eqn{\nu}, \eqn{\theta^*}, \eqn{\sigma^2}, and nodes for fitted values by default. Variables specified in \code{monitors2} will be added if provided.}
#'   \item{\code{fitted}}{If \code{fitted = TRUE}, in-sample fitted values.}
#'   \item{\code{input}}{List of input values for prior, initial values and execution time without compiling.}
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
#' # One tool version
#' fit <- gpSpike(X, as.vector(y))
#'
#' # Split version
#' models <- gpSpike(X, as.vector(y), sampling = FALSE)
#' Ccompile <- compileModelAndMCMC(models)
#' mcmc.out <- runMCMC(Ccompile$mcmc, niter = 5000, nburnin = 1000, thin = 1,
#'                    nchains = 1, setSeed = TRUE, init = models$input$init,
#'                    summary = TRUE, samplesAsCodaMCMC = TRUE)
#' }
#'
#' @references
#' McGee, G., Wilson, A., Webster, T. F., & Coull, B. A. (2023).
#' Bayesian multiple index models for environmental mixtures.
#' \emph{Biometrics}, 79(1), 462-474.
#'
#'
#' @export
#'


gpSpike <- function(x, y,
                    prior = list(
                            index = list(r1 = 1, r2 = 1, sigma_theta = 0.25),
                            link = list(inv_lambda_shape = 1, inv_lambda_rate = 0.1),
                            sigma2 = list(shape = 0.001, rate = 0.001)),
                    init = list(index = list(pi = 0.5, nu = NULL, index = NULL),
                                link = list(inv_lambda = NULL), sigma2 = NULL),
                    sampling = TRUE, fitted = TRUE,
                    monitors2 = NULL, niter = 10000, nburnin=1000,
                    thin = 1, thin2 = NULL, nchain = 1, setSeed = FALSE
){
  start1 <- Sys.time()
  index_raw <- 0; nu <- 0; invlambda <- 0; sigma2 <- 0

  # check sampling, prior, init parameters for independent execution
  checkOutput <- validate_and_finalize_args(
    sampling, fitted, niter, nburnin, thin, thin2, nchain,
    prior, init, "spike", "gp"
  )
  prior <- checkOutput$priorlist_final
  init <- checkOutput$initlist_final

  # environment
  envobj <- ls(envir=.GlobalEnv)
  on.exit(rm(list=ls(envir=.GlobalEnv)[which(!ls(envir=.GlobalEnv)%in%envobj)],envir=.GlobalEnv))


  # aa_bspline_ver3
  .fns <- c(
    # a_common
    "quickSortOrderIndexOnly", "nimOrder", "nimSort",
    "sampleQuantile_nim", "quantile_nimble",

    # gpspike
    "expcov_gpSpike","expcovnn_gpSpike","computeA1","computeA2",
    "llFunLambda","llFunThetaV1","llFunThetaV2","transitionTheta",
    "multiplyMatrixByConstant","pred_gpSpike",
    "SamplingLambda_gp_spike","SamplingThetaV_gp_spike",
    "gibbsSigma2_gp_spike","gibbsGam_gp_spike",


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


  # standardized X
  meanX <- apply(X, 2, mean); sdX <- apply(X, 2, sd)
  scaleX <- matrix(0, ncol = p, nrow = N)
  for (i in 1:p){
    scaleX[,i] <- (X[,i] - meanX[i])/sdX[i]
  }

  ## zero mean
  # z <- 1
  # Z <- matrix(rep(0, N), nrow = N)
  # init_gam <- rep(0, z)
  # init_Zlin <- as.vector(Z %*% matrix(init_gam, nrow = z))


  # Prior parameters
  ## check data dimension and save
  ## theta
  if (is.null(prior$index$r1)||
      length(prior$index$r1) >= 2 ||
      prior$index$r1 < 0){
    stop("Prior index (r1) has incorrect value.")
  } else{
    a0 <- prior$index$r1
  }

  if (is.null(prior$index$r2)||
      length(prior$index$r2) >= 2 ||
      prior$index$r2 < 0){
    stop("Prior index (r2) has incorrect value.")
  } else{
    b0 <- prior$index$r2
  }

  if (is.null(prior$index$sigma_theta)||
      length(prior$index$sigma_theta) >= 2 ||
      prior$index$sigma_theta < 0){
    stop("Prior index (sigma_theta) has incorrect value.")
  } else{
    sigma_theta <- prior$index$sigma_theta
  }

  ## link
  if ((!is.null(prior$link$inv_lambda_shape))&
      (length(prior$link$inv_lambda_shape) >= 2 ||
      prior$link$inv_lambda_shape < 0)){
    stop("Prior inverse lambda (a) has incorrect value.")
  }
  a_lam <- prior$link$inv_lambda_shape


  if ((!is.null(prior$link$inv_lambda_rate))&
      (length(prior$link$inv_lambda_rate) >= 2 ||
      prior$link$inv_lambda_rate < 0)){
    stop("Prior inverse lambda (b) has incorrect value.")
  }
  b_lam <- prior$link$inv_lambda_rate


  ## sigma
  if (is.null(prior$sigma2$shape)||length(prior$sigma2$shape) >= 2 || prior$sigma2$shape < 0){
    stop("Prior sigma2 (a) has incorrect value.")
  } else{
    a_sig <- prior$sigma2$shape
  }

  if (is.null(prior$sigma2$rate)||length(prior$sigma2$rate) >= 2 || prior$sigma2$rate < 0){
    stop("Prior sigma2 (b) has incorrect value.")
  } else{
    b_sig <- prior$sigma2$rate
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


  # initial value
  inits_list <- lapply(seq_len(nchain),
                       function(j) initfunction_gpSpike(x = scaleX, y = Y, p = p,
                                                        pi = init$index$pi, nu = init$index$nu,
                                                        index = init$index$index, inv_lambda = init$link$inv_lambda,
                                                        sigma2 = init$sigma2,
                                                        setSeed = seedNum[j]))
  firstInit <- inits_list[[1]]



  # Model code
  Rmodelcode <- nimbleCode({
    # thetastar: spike and slab
    pi ~ dbeta(a0, b0)
    for (j in 1:p) {
      nu[j] ~ dbern(pi)
      index_raw[j] ~ dnorm(0, sd = sigma_theta)
      indexstar[j] <- index_raw[j] * nu[j]
    }

    # Linear predictor
    for (i in 1:N){
      Xlin[i] <- sum(X[i,1:p] * indexstar[1:p])
    }

    # likelihood
    invlambda ~ dgamma(a_lam, rate = b_lam)
    sigma2 ~ dinvgamma(a_sig, rate = b_sig)
    Ki[1:N, 1:N] <- expcov_gpSpike(Xlin[1:N], invlambda) # invlambda * (I + lambda^(-1)*K)
    cov[1:N, 1:N] <- sigma2 * Ki[1:N, 1:N]

    Y[1:N, 1] ~ dmnorm(mu0[1:N], cov = cov[1:N, 1:N])

  })

  # Build model
  message("Build Model")
  suppressMessages(simpleModel <- nimbleModel(Rmodelcode,
                        data = list(X = scaleX, Y = Y),
                        constants = list(p = p, N = N,
                                         a0 = a0, b0 = b0,
                                         a_lam = a_lam, b_lam = b_lam,
                                         a_sig = a_sig, b_sig = b_sig,
                                         sigma_theta = sigma_theta,
                                         mu0 = rep(0, N)),
                        inits = firstInit))


  # Assign samplers
  message("Assign samplers")
  monitorsList <- c("nu", "indexstar", "sigma2")
  if (fitted){
    monitorsList <- c(monitorsList, "Xlin", "invlambda")
  }

  if (is.null(monitors2)){
    suppressMessages(mcmcConf <- configureMCMC(simpleModel,
                              monitors = monitorsList,
                              print = FALSE))
  } else{
    suppressMessages(mcmcConf <- configureMCMC(simpleModel,
                              monitors = monitorsList,
                              monitors2 = monitors2,
                              print = FALSE))
  }

  # 1. invlambda
  mcmcConf$removeSamplers(c("invlambda"))
  mcmcConf$addSampler(target = c("invlambda"),
                       type   = SamplingLambda_gp_spike)

  # 2. index_raw, v, (indexstar)
  mcmcConf$removeSamplers(c("index_raw", "nu"))
  mcmcConf$addSampler(target = c("index_raw", "nu"),
                       type   = SamplingThetaV_gp_spike)

  # 3. sigma^2
  mcmcConf$removeSamplers(c("sigma2"))
  mcmcConf$addSampler(target = c("sigma2"),
                       type   = gibbsSigma2_gp_spike)
  # mcmcConf1$printSamplers(executionOrder= TRUE)

  # 4. gam
  # mcmcConf$removeSamplers(c("gam"))
  # mcmcConf$addSampler(target = c("gam"),
  #                      type   = "gibbsGam_gp_spike")


  mcmc1 <- buildMCMC(mcmcConf)
  end1 <- Sys.time()

  if (!sampling){
    mcmc.out <- NULL
    fittedResult <- NULL
    sampMCMC <- NULL

  } else{
    # Compile
    message("Compile Model")
    suppressMessages(CsimpleModel <- compileNimble(simpleModel))
    message("Compile MCMC")
    suppressMessages(Cmcmc <- compileNimble(mcmc1, project = simpleModel,
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
    end2 <- Sys.time()
    start3 <- Sys.time()

    # output
    # fittedResult <- NULL
    if (fitted){ # posterior fitted value output (mean, median, sd)
      message("Compute posterior fitted value")

      nsamp <- nrow(samples)
      namesXlin <- paste0("Xlin[", 1:N, "]")
      XlinSample <- samples[, namesXlin]
      sigma2_samples <- samples[, "sigma2"]
      namesIndex <- paste0("indexstar[", 1:p, "]")
      indexstarSample <- samples[, namesIndex]
      invlambdaSample <- samples[,"invlambda"]
      newdataMat <- scaleX

      message("Compile function..")
      suppressMessages(cpred_gpSpike <- compileNimble(pred_gpSpike))
      message("Computing predicted value..")
      fittedValue <- cpred_gpSpike(newdataMat, nsamp, as.vector(Y), indexstarSample,
                                XlinSample, sigma2_samples, invlambdaSample)

      # fittedResult <- data.frame(mean = apply(fittedValue, 2, mean),
      #                            median = apply(fittedValue, 2, median),
      #                            sd = apply(fittedValue, 2, sd),
      #                            LB = apply(fittedValue, 2, quantile, probs = 0.025, na.rm = TRUE),
      #                            UB = apply(fittedValue, 2, quantile, probs = 0.975, na.rm = TRUE))

      fittedResult <- fittedValue
    } else{
      fittedResult <- NULL
    }
  }
  end3 <- Sys.time()

  ## Input options
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
                       prior = list(index = list(r1 = a0, r2 = b0, sigma_theta = sigma_theta),
                                    link = list(inv_lambda_shape = a_lam, inv_lambda_rate = b_lam),
                                    sigma2 = list(shape = a_sig, rate = b_sig)),
                       init = inits_list,
                       time = time)



  out <- list(model = simpleModel, sampler = mcmc1, sampling = sampMCMC,
              fitted = fittedResult, input = inputOptions,
              modelName = "gpSpike")
  class(out) = "bsimGp"
  return(out)
}

