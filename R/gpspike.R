#' Bayesian Single-Index Regression with Gaussian Process Link and Spike-and-Slab Prior
#'
#' @description Fits a single-index model \eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n},
#' where index vector \eqn{\theta} has a spike and slab prior and
#' the link \eqn{f(\cdot)} is represented by Gaussian process and the
#'
#' @inheritParams bsFisher
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
#' \item The variable selection indicator \eqn{\nu} has Beta–Bernoulli hierarchy prior. Beta hyper-prior on the Bernoulli–inclusion probability \eqn{w},
#'     inducing \eqn{p(\nu) \propto \mathrm{Beta}(r_1 + n_\nu, r_2 + p - n_\nu)} where  \eqn{n_\nu = \Sigma_{i=1}^{p}I(\nu_i = 1)}.
#'     \eqn{r_1, r_2} are shape and rate parameter of beta distribution.
#' \item Slab coefficients \eqn{\theta} have normal distribution with zero mean and \eqn{\sigma_\theta^2} variance.
#'   \item GP precision \eqn{\lambda^{-1}} follows gamma distribution with shape parameter \eqn{a_\lambda}, and rate parameter \eqn{b_\lambda}.
#'   \item Error precision \eqn{(\sigma^2)^{-1}} follows gamma distribution with shape parameter \eqn{a_\sigma}, and rate parameter \eqn{b_\sigma}.
#' }
#'
#' \strong{Sampling} A random walk Metropolis algorithm is used to sample \eqn{\lambda^{-1}}
#' and a Metropolis-Hastings algorithm is used for the main parameters \eqn{(\theta^{*}, \nu)}.
#' The variance \eqn{\sigma^2} is directly sampled from posterior distribution.
#' \eqn{f} is not directly sampled by MCMC.
#'
#' \strong{Prior hyper-parameters}
#' These are the prior hyper-parameters set in the function. You can define new values for each parameter in \link{prior_param}.
#' \enumerate{
#' \item Index vector: \code{index_nu_r1, index_nu_r2} gives the shape and rate parameter of beta-binomial prior, respectively.
#'     For slab prior, normal distribution with zero mean is assigned for selected variables \eqn{\theta}. \code{index_sigma_theta} is for variance of \eqn{\theta}, and it is assigned 0.25 by default.
#'     \item Link function: Inverse gamma prior is assigned for hyper-parameters \eqn{\lambda^{-1}}
#'     `link_inv_lambda_shape` is shape parameter and `link_inv_lambda_rate` is rate parameter of inverse gamma distribution.
#'     (default \code{link_inv_lambda_shape = 1, link_inv_lambda_rate = 0.1})
#'     \item Error variance (\code{sigma2}): An Inverse gamma prior is assigned to \eqn{\sigma^2}
#'         where \code{sigma2_shape} is shape parameter and \code{sigma2_rate} is rate parameter of inverse gamma distribution.
#'         (default \code{sigma2_shape = 0.001, sigma2_rate = 100})
#' }
#'
#' \strong{Initial values}
#' These are the initial values set in the function. You can define new values for each initial value in \link{init_param}
#' \enumerate{
#' \item Index vector:
#' \itemize{
#'      \item{\code{index_pi}: Initial selecting variable probability. (default: \code{0.5})}
#'      \item{\code{index_nu}: Initial vector of inclusion indicators . By default, each \code{index_nu} is randomly drawn by  \eqn{Bernoulli(1/2)}}
#'      \item{\code{index}: Initial vector of index. By default, each element of index vector, which is chosen by indicator, is proposed by normal distribution.}
#' }
#' \item Link function: Initial scalar of lambda (\code{link_inv_lambda}) for covariance kernel of Gaussian process.
#' \item Error variance (\code{sigma2}): Initial scalar error variance. (default: \code{0.01})
#'
#' }
#'
#'
#'@inherit bsFisher return
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
#' simdata <- data.frame(x = X, y = y)
#' colnames(simdata) <- c(paste0("X", 1:4), "y")
#'
#' # One tool version
#' fit1 <- gpSpike(y ~ ., data = simdata,
#'                niter = 5000, nburnin = 1000)
#'
#' # Split version
#' models <- gpSpike_setup(y ~ ., data = simdata)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, inits = initList,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#' fit2 <- as_bsim(models, mcmc.out)
#' summary(fit2)
#' }
#'
#' @references
#' McGee, G., Wilson, A., Webster, T. F., & Coull, B. A. (2023).
#' Bayesian multiple index models for environmental mixtures.
#' \emph{Biometrics}, 79(1), 462-474.
#'
#' @name gpSpike
#' @export
gpSpike <- function(formula, data,
                    prior = NULL,
                    init = NULL,
                    monitors = NULL, niter = 10000, nburnin=1000,
                    thin = 1, nchain = 1, setSeed = FALSE
){
  return(gpSpike.default(formula = formula, data = data,
                         prior = prior,
                         init = init,
                         sampling = TRUE, monitors = monitors,
                         niter = niter, nburnin=nburnin,
                         thin = thin, nchain = nchain, setSeed = setSeed))
}

#' @rdname gpSpike
#' @export
gpSpike_setup <- function(formula, data,
                          prior = NULL,
                          init = NULL,
                          monitors = NULL, niter = 10000, nburnin=1000,
                          thin = 1, nchain = 1, setSeed = FALSE
){
  return(gpSpike.default(formula = formula, data = data,
                         prior = prior,
                         init = init,
                         sampling = FALSE, monitors = monitors,
                         niter = niter, nburnin=nburnin,
                         thin = thin, nchain = nchain, setSeed = setSeed))
}

gpSpike.default <- function(formula, data,
                    prior = NULL,
                    init = NULL,
                    sampling = TRUE,
                    monitors = NULL, niter = 10000, nburnin=1000,
                    thin = 1, nchain = 1, setSeed = FALSE
){
  start1 <- Sys.time()
  index_raw <- 0; nu <- 0; invlambda <- 0; sigma2 <- 0

  # check sampling, prior, init parameters for independent execution
  checkOutput <- validate_and_finalize_args(
    sampling, niter, nburnin, thin, nchain,
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



  # Data
  if (!is.data.frame(data)){
    stop("data should be an data.frame.")
  }

  Call <- match.call()
  indx <- match(c("formula","data"), names(Call), nomatch = 0L)
  if (indx[1] == 0L)
    stop("a 'formula' argument is required")
  temp <- Call[c(1L,indx)]
  temp[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(temp)
  Terms <- attr(m,"terms")
  formula <- as.character(formula)
  response.name <- formula[2]
  data.name <- strsplit(formula[3]," \\+ ")[[1]]
  int.flag <- any(strsplit(formula[3]," \\* ")[[1]] == formula[3])
  if(data.name[1]=="."){
    tot.name <- response.name
  } else{
    tot.name <- c(response.name ,data.name)
  }
  if(!int.flag){
    stop("BayesSIM cannot treat interaction terms")
  }else if(!sum(duplicated(c(colnames(data),tot.name))[-c(1:ncol(data))])==
           length(tot.name)){
    stop(paste(paste(tot.name[duplicated(c(colnames(data),
                                           tot.name))[-c(1:ncol(data))]],collapse=","),
               " is/are not in your data"))
  }else{
    origY <- data[ ,response.name]
    if(data.name[1]=="."){
      origX <- data[,colnames(data) != response.name]
    }else {
      origX <- data[ ,data.name,drop=FALSE]
    }
  }

  # X = origX, Y = origY
  X <- as.matrix(origX)
  Y <- as.matrix(origY)

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
    index[1:p] <- indexstar[1:p]/sqrt(sum(indexstar[1:p]^2))

    # Linear predictor
    for (i in 1:N){
      Xlin[i] <- sum(X[i,1:p] * index[1:p])
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
  monitorsList <- c("nu", "index", "sigma2", "Xlin", "invlambda", "pi",
                    "index_raw")

  suppressMessages(mcmcConf <- configureMCMC(simpleModel,
                                             monitors = monitorsList,
                                             print = FALSE))


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


  if (!sampling){
    mcmc.out <- NULL
    samples <- NULL

  } else{
    # Compile
    start2 <- Sys.time()
    message("Compile Model")
    suppressMessages(CsimpleModel <- compileNimble(simpleModel))
    message("Compile MCMC")
    suppressMessages(Cmcmc <- compileNimble(mcmc1, project = simpleModel,
                                            resetFunctions = TRUE))
    end2 <- Sys.time()

    # Sampling
    message("Run MCMC")
    if (setSeed == FALSE){
      seedNum <- setSeed
    }
    mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                        thin = thin,
                        nchains = nchain, setSeed = seedNum,
                        inits = inits_list,
                        summary = FALSE, samplesAsCodaMCMC = TRUE)


    samples <- NULL

    if (nchain > 1){
      for (i in 1:nchain){
        samples <- rbind(samples, mcmc.out[[i]])
      }
    } else if (nchain == 1){
      samples <- mcmc.out
    }

    # # standardize index (make length = 1)
    # namesIndex <- paste0("indexstar", 1:p)
    # samples[, namesIndex] <- t(apply(samples[, namesIndex, drop = FALSE], 1, function(x) {
    #   x / sqrt(sum(x^2))
    # }))


  }
  end1 <- Sys.time()

  ## Input options
  ## Input options
  if (!sampling){
    time <- NULL
  } else{
    time <- difftime(end1, start1, units = "secs") - difftime(end2, start2, units = "secs")
  }


  # results ---------------------------------------------------------------------
  #Model point estimation
  modelESTLIST <- bsimFit_pointest(samples, X, Y)

  inputOptions <- list(origdata = list(x = X, y = Y), formula = formula,
                       prior = list(index = list(r1 = a0, r2 = b0, sigma_theta = sigma_theta),
                                    link = list(inv_lambda_shape = a_lam, inv_lambda_rate = b_lam),
                                    sigma2 = list(shape = a_sig, rate = b_sig)),
                       samplingOptions = list(sampling = sampling,
                                              monitors = monitors, niter = niter,
                                              nburnin = nburnin, thin = thin,
                                              nchain = nchain, setSeed = setSeed),
                       init = inits_list, # initial value for MCMC
                       time = time)


  if (sampling){
    out <- list(coefficients = modelESTLIST$coefficients,
                ses_coef = modelESTLIST$ses_coef, se = modelESTLIST$se,
                residuals = modelESTLIST$residuals,
                fitted.values = modelESTLIST$fitted.values,
                linear.predictors = modelESTLIST$linear.predictors,
                gof = modelESTLIST$gof,
                samples = mcmc.out, # return of runMCMC
                input = inputOptions,
                defModel = simpleModel, defSampler = mcmc1,
                modelName = "gpSpike")

    class(out) = "bsim"


  } else{
    out <- list(input = inputOptions,
                defModel = simpleModel,
                defSampler = mcmc1,
                modelName = "gpSpike")

    class(out) = "bsimSetup"

  }
  return(out)

}

