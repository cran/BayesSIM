#' @rdname gpPolar
#' @export
gpPolarHigh <- function(formula, data,
                        prior = NULL,
                        init = NULL,
                        monitors = NULL, niter = 10000, nburnin=1000,
                        thin = 1, nchain = 1, setSeed = FALSE
){
  return(gpPolarHigh.default(formula = formula, data = data,
                             prior = prior,
                             init = init,
                             sampling = TRUE, monitors = monitors,
                             niter = niter, nburnin=nburnin,
                             thin = thin, nchain = nchain, setSeed = setSeed))
}

#' @rdname gpPolar
#' @export
gpPolarHigh_setup <-function(formula, data,
                             prior = NULL,
                             init = NULL,
                             monitors = NULL, niter = 10000, nburnin=1000,
                             thin = 1, nchain = 1, setSeed = FALSE
){
  return(gpPolarHigh.default(formula = formula, data = data,
                             prior = prior,
                             init = init,
                             sampling = FALSE, monitors = monitors,
                             niter = niter, nburnin=nburnin,
                             thin = thin, nchain = nchain, setSeed = setSeed))
}


gpPolarHigh.default <- function(formula, data,
                    prior = NULL,
                    init = NULL,
                    sampling = TRUE,
                    monitors = NULL, niter = 10000, nburnin=1000,
                    thin = 1, nchain = 1, setSeed = FALSE
){
  start1 <- Sys.time()
  sigma2 <- 0; psi <- 0;

  # check sampling, prior, init parameters for independent execution
  checkOutput <- validate_and_finalize_args(
    sampling, niter, nburnin, thin, nchain,
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
  monitorsList <- c("index", "sigma2", "linkFunction", "kappa", "Xlin", "d", "psi")
  suppressMessages(mcmcConf <- configureMCMC(simpleModel,
                                             monitors = monitorsList,
                                             print = FALSE))



  mcmcConf$removeSamplers(c("sigma2"))
  mcmcConf$addSampler(target = c("sigma2"),
                      type   = gibbsSampler_sigma2)


  mcmc1 <- buildMCMC(mcmcConf)
  end1 <- Sys.time()

  if (!sampling){
    mcmc.out <- NULL
    samples <- NULL

  } else{
    start2 <- Sys.time()
    # Compile
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
    samples <- NULL

    if (nchain > 1){
      for (i in 1:nchain){
        samples <- rbind(samples, mcmc.out[[i]])
      }
    } else if (nchain == 1){
      samples <- mcmc.out
    }
  }
  end1 <- 0
  ## Input options
  if (!sampling){
    time <- NULL
  } else{
    time <- difftime(end1, start1, units = "secs") - difftime(end2, start2, units = "secs")
  }


  # Results ------------------------------------------------------------------
  # Model point estimation
  modelESTLIST <- bsimFit_pointest(samples, X, Y)

  # input options
  inputOptions <- list(origdata = list(x = X, y = Y), formula = formula,
                       prior = list(index = list(psi = list(alpha = psi_c)),
                                    link = list(kappa = list(min = kappa_min, max = kappa_max,
                                                             grid_width = kappa_grid_width)),
                                    sigma2 = list(shape = sigma2_shape, rate = sigma2_rate)),
                       init = inits_list,
                       samplingOptions = list(sampling = sampling,
                                              monitors = monitors, niter = niter,
                                              nburnin = nburnin, thin = thin,
                                              nchain = nchain, setSeed = setSeed),
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

