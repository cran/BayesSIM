#' Bayesian single-index regression with Gaussian process link and one-to-one polar transformation
#' @name gpPolar
#'
#' @description
#' Fits a single–index model \eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n} where
#' the index \eqn{\theta} is specified and computed via a one-to-one polar
#' transformation, and the link \eqn{f(\cdot)} is represented with a Gaussian
#' process.
#'
#' @param x Numeric data.frame/matrix of predictors. Each row is an observation.
#' @param y Numeric response numeric vector/matrix. Other types  are not available.
#'
#' @param prior Optional named list of prior settings with sublists:
#' \describe{
#' \item{\code{index}}{ \code{psi} is polar angle and rescaled Beta distribution on \eqn{[0, \pi]} is assigned.
#'     Only shape parameter \code{alpha} of \eqn{p-1} dimension vector is needed since rate parameters is computed to satisfy \eqn{\theta_{j, \text{MAP}}}.
#'     By default, the shape parameter for each element of the polar vector is set to \code{5000}.}
#'
#' \item{\code{link}}{{Prior for the smoothness parameter \code{kappa} in the Gaussian process kernel: Prior for \eqn{\kappa} is discrete uniform of equally spaced grid points
#'     in \eqn{[\kappa_{\text{min}}, \kappa_{\text{max}}}].
#'      \code{min_kappa} is minimum value of kappa (default \code{0.5}), \code{max_kappa} is maximum value of kappa (default \code{4}),
#'      and \code{grid.width} is space (default \code{0.1}).}}
#' \item{\code{sigma2}}{Error-variance prior hyperparameters. An Inverse-Gamma prior is assigned to \eqn{\sigma^2}
#'         where \code{shape} is shape parameter and \code{rate} is rate parameter of inverse gamma distribution.
#'         (default \code{shape = 2, rate = 0.01})}
#'
#'   }
#' @param init Optional named list of initial values. If the values are not assigned, they are randomly sampled from prior.
#' \describe{
#'      \item{\code{index}}{Initial vector of polar angle \code{psi} with \eqn{p-1} dimension. Each element of angle is between 0 and \eqn{\pi}.}
#'      \item{\code{link}}{Initial scalar scale parameter of covariance kernel \code{kappa}. (default: \code{2})}
#'      \item{\code{sigma2}}{Initial scalar error variance. (default: \code{0.01})}
#' }
#'
#' @param sampling Logical. If \code{TRUE} (default), run MCMC; otherwise return prepared nimble model objects without sampling.
#' @param fitted Logical. If \code{TRUE} (default), fitted values drawn from posterior distribution are included in the output and \code{c("linkFunction", "kappa", "Xlin")} is monitored for prediction.
#' @param monitors2 Optional character vector of additional monitor nodes. To check the names of the nodes, set \code{fit <- gpPolar(x, y, sampling = FALSE)} and then inspect the variable names stored in the model object using \code{fit$model$getVarNames()}.
#' @param niter Integer. Total MCMC iterations (default \code{10000}).
#' @param nburnin Integer. Burn-in iterations (default \code{1000}).
#' @param thin Integer. Thinning for monitors1 (default \code{1}).
#' @param thin2 Integer. Optional thinning for \code{monitors2} (default \code{1}).
#' @param nchain Integer. Number of MCMC chains (default \code{1}).
#' @param setSeed Logical or numeric argument.  Further details are provided in \link[nimble]{runMCMC}.
#'
#' @details
#' \strong{Model} The single–index model is specified as \eqn{Y_i = f(X_i'{\theta}) + \epsilon_i},
#' where the index vector \eqn{\theta} lies on the unit sphere with (\eqn{\|\theta\|_2=1}) with non-zero first component
#' to ensure identifiability and is parameterized via a one-to-one polar transformation with angle \eqn{\psi_1,...,\psi_{p-1}}.
#' Sampling is  performed on the angular parameters \eqn{\theta} defining
#' the index vector. The link function \eqn{f(\cdot)} is modeled by a Gaussian process
#' prior with zero mean and an Ornstein–Uhlenbeck (OU) covariance kernel
#' \eqn{\exp(-\kappa |t_i - t_j|), i, j = 1,\ldots, N}, where \eqn{\kappa} is a bandwidth (smoothness)
#' parameter and \eqn{t_i, t_j} is index value(\eqn{t_i = X_i'\theta}).

#'
#' \strong{Priors}
#' \itemize{
#'   \item Index vector: Uniform prior on the unit sphere (\eqn{\|\theta\|_2=1}).
#'   \item Bandwidth parameter \eqn{\kappa}: discrete uniform prior over a fixed grid.
#'   \item Error variance \eqn{\sigma^2}: Inverse–Gamma prior.
#'
#' }
#'
#' \strong{Sampling} For \code{gpPolar()}, \eqn{\theta} is sampled by Metropolis-Hastings and updated with \eqn{f},
#' \eqn{\kappa} is chosen by grid search method that maximizes likelihood,
#' \eqn{\sigma^2} are sampled from their full conditional
#' distributions using Gibbs sampling.
#' Since \eqn{\kappa} is sampled by grid search, more than 5 dimension of variables \code{gpPolarHigh()} is recommended.
#' For \code{gpPolarHigh()}, all sampling parameters' samplers are assigned by nimble.
#'
#'
#' @return A \code{list} typically containing:
#' \describe{
#'   \item{\code{model}}{Nimble model}
#'   \item{\code{sampler}}{Nimble sampler}
#'   \item{\code{sampling}}{Posterior draws of \eqn{\theta}, \eqn{\sigma^2}, and nodes for fitted values by default. Variables specified in \code{monitors2} will be added if provided.}
#'   \item{\code{fitted}}{If \code{fitted = TRUE}, in-sample fitted values is given.}
#'   \item{\code{input}}{List of data,input values for prior and initial values, and computation time without compiling.}
#' }
#'
#' @examples
#' \donttest{
#' library(MASS)
#' N <- 100    # Sample Size
#' p <- 3
#' mu <- c(0,0,0)
#' rho <- 0
#' Cx <- rbind(c(1,rho,rho), c(rho,1,rho), c(rho, rho,1))
#' X <- mvrnorm(n = N, mu=mu, Sigma=Cx, tol=1e-8)
#' alpha <- c(1,1,1)
#' alpha <- alpha/sqrt(sum(alpha^2))
#' z <- matrix(0,N)
#' z <- X %*% alpha
#' z <- z[,1]
#' sigma <- 0.3
#' f <- exp(z)
#' y <- f + rnorm(N, 0, sd=sigma) # adding noise
#' y <- y-mean(y)
#' f_all <- f
#' x_all <- X
#' y_all <- y
#' data_frame <- cbind(x_all, y, f)
#' colnames(data_frame) = c('x1', 'x2', 'x3', 'y','f')
#'
#' # One-call version
#' fit1 <- gpPolar(X, y)
#' fit2 <- gpPolarHigh(X, y)
#'
#' # Split version
#' models1 <- gpPolar(X, y, sampling = FALSE)
#' models2 <- gpPolarHigh(X, y, sampling = FALSE)
#' Ccompile1 <- compileModelAndMCMC(models1)
#' Ccompile2 <- compileModelAndMCMC(models2)
#' mcmc.out1 <- runMCMC(Ccompile1$mcmc, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, init = models1$input$init,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#' mcmc.out2 <- runMCMC(Ccompile2$mcmc, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, init = models2$input$init,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#' }
#'
#' @references
#' Dhara, K., Lipsitz, S., Pati, D., & Sinha, D. (2019). A new Bayesian single index model with or without covariates missing at random.
#' \emph{Bayesian analysis}, 15(3), 759.
#'
#'
#' @export

gpPolar <- function(x, y,
                    prior = list(index = list(psi = list(alpha = NULL)),
                                 link = list(kappa = list(min_kappa = 0.5, max_kappa = 4, grid.width = 0.1)),
                                 sigma2 = list(shape = 2, rate = 0.01)),
                    init = list(index = list(psi = NULL),
                                 link = list(kappa = 2),
                                 sigma2 = 0.01),
                    sampling = TRUE, fitted = TRUE,
                    monitors2 = NULL, niter = 10000, nburnin=1000,
                    thin = 1, thin2 = NULL, nchain = 1, setSeed = FALSE
                    ){

  start1 <- Sys.time()
  sigma2 <- 0; psi <- 0

  # check sampling, prior, init parameters for independent execution
  checkOutput <- validate_and_finalize_args(
    sampling, fitted, niter, nburnin, thin, thin2, nchain,
    prior, init, "polar", "gp"
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

    # gpPolar
    "alphaTheta","Xlinear","invcov","expcov_gpPolar","expcovTest_gpPolar",
    "obj_btt_theta","thetaPrior","pred_gpPolar","gibbsSampler_sigma2",
    "gibbsSampler_kappa","MH_thetaeta",


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
  Rmodel <- nimbleCode({
    # 0. Likelihood(y)
    Sigma[1:N,1:N] <-  diag(rep(sigma2, N))
    y[1:N] ~ dmnorm(linkFunction[1:N], cov = Sigma[1:N,1:N])

    # 1. linkFunction(f) - gibbs
    cov[1:N, 1:N] <- expcov_gpPolar(Xlin[1:N], kappa)
    linkFunction[1:N] ~ dmnorm(mu0[1:N], cov = cov[1:N, 1:N])

    # 2. sigma2 - gibbs
    sigma2 ~ dinvgamma(a, b)

    # 3. hyperprior-kappa(MAP)
    kappa ~ dunif(kappa_a, kappa_b)

    # 4. index, psi
    for(j in 1:(p-1)) {
      d[j] ~ dunif(1e-10, 1e10)
      psi[j] ~ dbeta(c[j], d[j])
    }
    index[1:p] <- alphaTheta(psi[1:(p-1)]*pi)
    Xlin[1:N] <- Xlinear(index[1:p], x[1:N, 1:p])
  })


  # Prior parameters
  ## psi
  if (!is.null(prior$index$psi$alpha) & length(prior$index$psi$alpha) != (p-1))
  {stop("Prior psi has incorrect dimension")}
  if (is.null(prior$index$psi$alpha)){
    psi_c <- c(rep(5000,(p-1)))
  } else{
    psi_c <- prior$index$psi$alpha
  }

  ## sigma2 - scalar, > 0
  if (is.null(prior$sigma2$shape)||length(prior$sigma2$shape) >= 2 || prior$sigma2$shape < 0){
    stop("Prior sigma2 (a) has incorrect value.")
  } else{
    sigma2_shape <- prior$sigma2$shape
  }

  if (is.null(prior$sigma2$rate)||length(prior$sigma2$rate) >= 2 || prior$sigma2$rate < 0){
    stop("Prior sigma2 (b) has incorrect value.")
  } else{
    sigma2_rate <- prior$sigma2$rate
  }

  # kappa
  if (is.null(prior$link$kappa$min_kappa)||length(prior$link$kappa$min_kappa) >= 2 || prior$link$kappa$min_kappa < 0){
    stop("Prior kappa (min_kappa) has incorrect value.")
  } else{
    kappa_min <- prior$link$kappa$min_kappa
  }

  if (is.null(prior$link$kappa$max_kappa)||length(prior$link$kappa$max_kappa) >= 2 || prior$link$kappa$max_kappa < 0){
    stop("Prior kappa (max_kappa) has incorrect value.")
  } else{
    kappa_max <- prior$link$kappa$max_kappa
  }

  if (is.null(prior$link$kappa$grid.width)||length(prior$link$kappa$grid.width) >= 2 || prior$link$kappa$grid.width < 0){
    stop("Prior kappa (grid.width) has incorrect value.")
  } else{
    kappa_grid_width <- prior$link$kappa$grid.width
  }

  # Initialize
  init_psi <- init$link$psi
  if (!is.null(init_psi) & length(init_psi) != (p-1)){
    stop("Initial psi has incorrect dimension")
  }

  if (is.null(init$link$kappa)||length(init$link$kappa) >= 2 || init$link$kappa < 0){
    stop("Initial kappa has incorrect value.")
  } else{
    init_kappa <- init$link$kappa
  }

  if (is.null(init$sigma2)||length(init$sigma2) >= 2 || init$sigma2 < 0){
    stop("Initial sigma2 has incorrect value.")
  } else{
    init_sigma2 <- init$sigma2
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


  inits_list <- lapply(seq_len(nchain), function(j) initfunction_gpPolar(x = X, y = as.vector(Y), kappa_init = init_kappa,
                                                                         sigma2_init = init_sigma2, psi_init = init_psi,
                                                                         grid.with = kappa_grid_width,
                                                                         sig_a = sigma2_shape, sig_b = sigma2_rate,
                                                                         setSeed = seedNum[j]))
  firstInit <- inits_list[[1]]




  message("Build Model")
  suppressMessages(simpleModel <- nimbleModel(Rmodel,
                             data = list(x = X,
                                         y = as.vector(Y)),
                             constants = list(c = psi_c,
                                              kappa_a = kappa_min, kappa_b = kappa_max,
                                              a = sigma2_shape, b = sigma2_rate,
                                              p = p, N = N, mu0 = rep(0, N), pi = pi),
                             inits = firstInit))

  # Assign samplers
  message("Assign samplers")
  # monitorsList <-  c("linkFunction","index", "psi", "kappa", "sigma", "d")
  monitorsList <- c("index", "sigma2")
  if (fitted){
    monitorsList <- c(monitorsList, "linkFunction", "kappa", "Xlin")
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


  mcmcConf$removeSamplers(c("sigma2"))
  mcmcConf$addSampler(target = c("sigma2"),
                      type   = gibbsSampler_sigma2)

  mcmcConf$removeSamplers(c("kappa"))
  mcmcConf$addSampler(target = c("kappa"),
                      type   = gibbsSampler_kappa,
                      control = list(grid.width = 0.1))

  mcmcConf$removeSamplers(c("psi", "d", "linkFunction"))
  mcmcConf$addSampler(target = c("psi", "d", "linkFunction"),
                      type   = MH_thetaeta)


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
                          nchains = nchain, setSeed = seedNum, inits = inits_list,
                          summary = FALSE, samplesAsCodaMCMC = TRUE)
    } else{
      mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                          thin = thin, thin2 = thin2,
                          nchains = nchain, setSeed = seedNum,
                          inits = inits_list,
                          summary = FALSE, samplesAsCodaMCMC = TRUE)
    }
    # output
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
    if (fitted){ # posterior fitted value output (mean, median, sd)
      message("Compute posterior fitted value")


      # namesBeta <- paste0("theta", 1:p)
      namesLink <- paste0("linkFunction[", 1:N, "]")
      namesSigma <- "sigma2"
      LinkFunction_samples <- samples[, namesLink]
      sigma2_samples <- samples[, namesSigma]
      n <- nrow(LinkFunction_samples)
      p <- ncol(LinkFunction_samples)

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
                       prior = list(index = list(psi = list(alpha = psi_c)),
                                    link = list(kappa = list(min_kappa = kappa_min, max_kappa = kappa_max,
                                                             grid.width = kappa_grid_width)),
                                    sigma2 = list(shape = sigma2_shape, rate = sigma2_rate)),
                       # initial value for MCMC
                       init = inits_list,
                       time = time)



  out <- list(model = simpleModel, sampler = mcmc1, sampling = sampMCMC,
              fitted = fittedResult, input = inputOptions,
              modelName = "gpPolar")
  class(out) = "bsimGp"
  return(out)
}


# Operation
# library(MASS)
# N = 100    # Sample Size
# p = 3
# mu=c(0,0,0)
# rho=0
# Cx<-rbind(c(1,rho,rho), c(rho,1,rho), c(rho, rho,1))
# x<-mvrnorm(n = N, mu=mu, Sigma=Cx, tol=1e-8)
# alpha <- c(1,1,1)
# alpha <- alpha/sqrt(sum(alpha^2))
# z <- matrix(0,N)
# z <- x %*% alpha
# z <- z[,1]
# sigma <- .3
# f <- exp(z)
# y <- f + rnorm(N,0,sd=sigma) # adding noise
# y <- y-mean(y)
#
# f_all <- f
# x_all <- x
# z_all <- z
# y_all <- y
# data_frame <- cbind(x_all, y, f)
# colnames(data_frame) = c('x1', 'x2', 'x3', 'y','f')
#
# # One tool version
# modelPolar <- gpPolar(x, y, fitted = TRUE)
# # #
# # Split version
# models <- gpPolar(x, y, sampling = FALSE)
# Ccompile <- compileModelAndMCMC(models)
# mcmc.out <- runMCMC(Ccompile$mcmc, niter = 10000, nburnin = 2000, thin = 1,
#                     nchains = 1, setSeed = TRUE,
#                     summary = TRUE, samplesAsCodaMCMC = TRUE)
