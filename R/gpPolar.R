#' Bayesian Single-Index Regression with Gaussian Process Link and One-to-One Polar Transformation
#'
#' @description
#' Fits a single–index model \eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n}, where
#' the index \eqn{\theta} is specified and computed via a one-to-one polar
#' transformation, and the link \eqn{f(\cdot)} is represented with a Gaussian process.
#' @inheritParams bsFisher
#'
#' @details
#' \strong{Model} The single–index model is specified as \eqn{Y_i = f(X_i'{\theta}) + \epsilon_i},
#' where the index vector \eqn{\theta} lies on the unit sphere with (\eqn{\|\theta\|_2=1}) with non-zero first component
#' to ensure identifiability and is parameterized via a one-to-one polar transformation with angle \eqn{\psi_1,...,\psi_{p-1}}.
#'
#' The mapping is
#' \deqn{
#' \begin{aligned}
#' \theta_1 &= \sin(\psi_1),\\
#' \theta_i &= \Big(\prod_{j=1}^{i-1}\cos(\psi_j)\Big)\sin(\psi_i), \quad i=2,\dots,p-1,\\
#' \theta_p &= \prod_{j=1}^{p-1}\cos(\psi_j).
#' \end{aligned}
#' }
#' The vector is then scaled to unit length.
#'
#' Sampling is  performed on the angular parameters \eqn{\theta} defining
#' the index vector. The link function \eqn{f(\cdot)} is modeled by a Gaussian process
#' prior with zero mean and an Ornstein–Uhlenbeck (OU) covariance kernel
#' \eqn{\exp(-\kappa \cdot |t_i - t_j|), i, j = 1,\ldots, n}, where \eqn{\kappa} is a bandwidth (smoothness)
#' parameter and \eqn{t_i, t_j} is index value (\eqn{t_i = X_i'\theta}).

#'
#' \strong{Priors}
#' \itemize{
#'   \item \eqn{\psi} is \eqn{p-1} dimension of polar angle of index vector and re-scaled Beta distribution on \eqn{[0, \pi]} is assigned for prior.
#'   \item Prior for \eqn{\kappa} (bandwidth parameter) is discrete uniform of equally spaced grid points in \eqn{[\kappa_{\text{min}}, \kappa_{\text{max}}}].
#'   \item Inverse gamma prior on \eqn{\sigma^2} with shape parameter \eqn{a_\sigma} and rate parameter \eqn{b_\sigma}.
#'
#' }
#'
#' \strong{Sampling} For \code{gpPolar}, \eqn{\theta} is sampled by Metropolis-Hastings and updated with \eqn{f},
#' \eqn{\kappa} is chosen by grid search method that maximizes likelihood,
#' \eqn{\sigma^2} are sampled from their full conditional
#' distributions using Gibbs sampling.
#' Since \eqn{\kappa} is sampled by grid search, more than 5 dimension of variables \code{gpPolarHigh} is recommended.
#' For \code{gpPolarHigh}, all sampling parameters' samplers are assigned by nimble.
#'
#' \strong{Prior hyper-parameters}
#' These are the prior hyper-parameters set in the function. You can define new values for each parameter in \link{prior_param}.
#' \enumerate{
#' \item Index vector: Only shape parameter \code{index_psi_alpha} of \eqn{p-1} dimension vector is needed since rate parameters is computed to satisfy \eqn{\theta_{j, \text{MAP}}}.
#'     By default, the shape parameter for each element of the polar vector is set to \code{5000}.
#'
#' \item Link function:
#'      \code{link_kappa_min} is minimum value of kappa (default \code{0.5}), \code{link_kappa_max} is maximum value of kappa (default \code{4}),
#'      and \code{link_kappa_grid_width} is space (default \code{0.1}).
#' \item Error variance (\code{sigma2}): An Inverse gamma prior is assigned to \eqn{\sigma^2}
#'         where \code{sigma2_shape} is shape parameter and \code{sigma2_rate} is rate parameter of inverse gamma distribution.
#'         (default \code{sigma2_shape = 2, sigma2_rate = 0.01})
#' }
#'
#' \strong{Initial values}
#' These are the initial values set in the function. You can define new values for each initial value in \link{init_param}.
#' \enumerate{
#'
#' \item Index vector: Initial vector of polar angle \code{index_psi} with \eqn{p-1} dimension. Each element of angle is between 0 and \eqn{\pi}.
#'      \item Link function: Initial scalar scale parameter of covariance kernel \code{link_kappa}. (default: \code{2})
#'      \item Error variance (\code{sigma2}): Initial scalar error variance. (default: \code{0.01})
#' }
#'
#'
#'
#' @inherit bsFisher return
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
#' simdata <- cbind(x_all, y, f)
#' simdata <- as.data.frame(simdata)
#' colnames(simdata) = c('x1', 'x2', 'x3', 'y','f')
#'
#' # One tool version
#' fit1 <- gpPolar(y ~ x1 + x2 + x3, data = simdata,
#'                 niter = 5000, nburnin = 1000, nchain = 1)
#' fit2 <- gpPolarHigh(y ~ x1 + x2 + x3, data = simdata,
#'                     niter = 5000, nburnin = 1000, nchain = 1)
#'
#' # Split version
#' models1 <- gpPolar_setup(y ~ x1 + x2 + x3, data = simdata)
#' models2 <- gpPolarHigh_setup(y ~ x1 + x2 + x3, data = simdata)
#' Ccompile1 <- compileModelAndMCMC(models1)
#' Ccompile2 <- compileModelAndMCMC(models2)
#' sampler1 <- get_sampler(Ccompile1)
#' sampler2 <- get_sampler(Ccompile2)
#' initList1 <- getInit(models1)
#' initList2 <- getInit(models2)
#' mcmc.out1 <- runMCMC(sampler1, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, init = initList1,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#' mcmc.out2 <- runMCMC(sampler2, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, init = initList2,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#' fit1_split <- as_bsim(models1, mcmc.out1)
#' fit2_split <- as_bsim(models2, mcmc.out2)
#' }
#'
#' @references
#' Dhara, K., Lipsitz, S., Pati, D., & Sinha, D. (2019). A new Bayesian single index model with or without covariates missing at random.
#' \emph{Bayesian analysis}, 15(3), 759.
#'
#' @name gpPolar
#' @export
gpPolar <- function(formula, data,
                   prior = NULL,
                   init = NULL,
                   monitors = NULL, niter = 10000, nburnin=1000,
                   thin = 1, nchain = 1, setSeed = FALSE){

  return(gpPolar.default(formula = formula, data = data,
                          prior = prior,
                          init = init,
                          sampling = TRUE, monitors = monitors,
                         niter = niter, nburnin=nburnin,
                         thin = thin, nchain = nchain, setSeed = setSeed))
}

#' @rdname gpPolar
#' @export
gpPolar_setup <- function(formula, data,
                          prior = NULL,
                          init = NULL,
                          monitors = NULL, niter = 10000, nburnin=1000,
                          thin = 1, nchain = 1, setSeed = FALSE){
  return(gpPolar.default(formula = formula, data = data,
                         prior = prior,
                         init = init,
                         sampling = FALSE, monitors = monitors,
                         niter = niter, nburnin=nburnin,
                         thin = thin, nchain = nchain, setSeed = setSeed))
}

gpPolar.default <- function(formula, data,
                    prior = NULL,
                    init = NULL,
                    sampling = TRUE,
                    monitors = NULL, niter = 10000, nburnin=1000,
                    thin = 1, nchain = 1, setSeed = FALSE){

  start1 <- Sys.time()
  sigma2 <- 0; psi <- 0

  # check sampling, prior, init parameters for independent execution
  checkOutput <- validate_and_finalize_args(
    sampling, niter, nburnin, thin, nchain,
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



  # data
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
  if (is.null(prior$link$kappa$min)||length(prior$link$kappa$min) >= 2 || prior$link$kappa$min < 0){
    stop("Prior kappa (min) has incorrect value.")
  } else{
    kappa_min <- prior$link$kappa$min
  }

  if (is.null(prior$link$kappa$max)||length(prior$link$kappa$max) >= 2 || prior$link$kappa$max < 0){
    stop("Prior kappa (max) has incorrect value.")
  } else{
    kappa_max <- prior$link$kappa$max
  }

  if (is.null(prior$link$kappa$grid_width)||length(prior$link$kappa$grid_width) >= 2 || prior$link$kappa$grid_width < 0){
    stop("Prior kappa (grid_width) has incorrect value.")
  } else{
    kappa_grid_width <- prior$link$kappa$grid_width
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
  monitorsList <- c("index", "sigma2", "linkFunction", "kappa", "Xlin", "d", "psi")
  suppressMessages(mcmcConf <- configureMCMC(simpleModel,
                                             monitors = monitorsList,
                                             print = FALSE))


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
    samples <- NULL
    mcmc.out <- NULL

  } else{
    # Compile
    start2 <- Sys.time()
    message("Compile Model")
    suppressMessages(CsimpleModel <- compileNimble(simpleModel))
    message("Compile MCMC")
    suppressMessages(Cmcmc <- compileNimble(mcmc1,
                                            project = simpleModel,
                                            resetFunctions = TRUE))
    end2 <- Sys.time()

    # Sampling
    message("Run MCMC")
    if (setSeed == FALSE){
      seedNum <- setSeed
    }
    mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                        thin = thin,
                        nchains = nchain, setSeed = seedNum, inits = inits_list,
                        summary = FALSE, samplesAsCodaMCMC = TRUE)
    # output
    if (nchain > 1){
      for (i in 1:nchain){
        samples <- rbind(samples, mcmc.out[[i]])
      }
    } else{
      samples <- mcmc.out
    }
  }

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
                       prior = list(index = list(psi = list(alpha = psi_c)),
                                    link = list(kappa = list(min = kappa_min, max = kappa_max,
                                                             grid_width = kappa_grid_width)),
                                    sigma2 = list(shape = sigma2_shape, rate = sigma2_rate)),
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
                modelName = "gpPolar")

    class(out) = "bsim"


  } else{
    out <- list(input = inputOptions,
                defModel = simpleModel,
                defSampler = mcmc1,
                modelName = "gpPolar")

    class(out) = "bsimSetup"

  }
  return(out)

}

