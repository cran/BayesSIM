#' @rdname gpPolar
#' @export
gpPolarHigh <- function(x, y,
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
    # 1. index, psi
    for(j in 1:(p-1)) {
      d[j] ~ dunif(1e-10, 1e10)
      psi[j] ~ dbeta(c[j], d[j])
    }
    index[1:p] <- alphaTheta(psi[1:(p-1)]*pi)
    Xlin[1:N] <- Xlinear(index[1:p], x[1:N, 1:p])

    # 2. sigma2 - gibbs
    sigma2 ~ dinvgamma(a, b)

    # 2. hyperprior-kappa
    kappa ~ dunif(kappa_a, kappa_b)

    # 1. linkFunction(f)
    cov[1:N, 1:N] <- expcov_gpPolar(Xlin[1:N], kappa)
    linkFunction[1:N] ~ dmnorm(mu0[1:N], cov = cov[1:N, 1:N])

    # 0. Likelihood(y)
    Sigma[1:N,1:N] <-  diag(rep(sigma2, N))
    y[1:N] ~ dmnorm(linkFunction[1:N], cov = Sigma[1:N,1:N])


  })


  # Prior parameters
  ## psi
  if (!is.null(prior$index$psi$alpha) & length(prior$index$psi$alpha) != (p-1))
  {stop("Prior psi has incorrect dimension")}
  if (is.null(prior$index$psi$alpha)){
    psi_c <- rep(5000,(p-1))
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

  inits_list <- lapply(seq_len(nchain),
                       function(j) initfunction_gpPolar2(X = X, kappa = init_kappa,
                                                         sigma2 = init_sigma2, psi = init_psi,
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
    if (setSeed == FALSE){
      seedNum <- setSeed
    }
    mcmc.out <- NULL
    if (is.null(monitors2)){
      mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                          thin = thin,
                          nchains = nchain, setSeed = seedNum, inits = inits_list,
                          summary = FALSE, samplesAsCodaMCMC = TRUE)
    } else{
      mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                          thin = thin, thin2 = thin2,
                          nchains = nchain, setSeed = seedNum, inits = inits_list,
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

