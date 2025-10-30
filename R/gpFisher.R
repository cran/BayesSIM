#' Bayesian single-index regression with Gaussian process link and von Mises-Fisher prior
#'
#' @description
#' Fits a single–index model\eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n} where
#' the index \eqn{\theta} lies on the unit sphere with von Mises-Fisher prior, and the link \eqn{f(\cdot)} is represented
#' with Gaussian process.
#'
#' @param x Numeric data.frame/matrix of predictors. Each row is an observation.
#' @param y Numeric response numeric vector/matrix. Other types  are not available.
#' @param prior Optional named list of prior settings with sublists:
#' \describe{
#'   \item{\code{index}}{von Mises--Fisher prior for the projection vector \eqn{\theta}.
#'     \code{direction} is a unit direction vector of the von Mises--Fisher distribution.
#'     If \code{direction} is \code{NULL}, the estimated vector from projection pursuit regression is assigned.
#'     \code{dispersion} is the concentration parameter \eqn{c_{\mathrm{prior}} > 0}. (default \code{150})
#'   }
#' \item{\code{link}}{
#' \enumerate{
#'    \item{\code{lenghscale}: Prior of length-scale parameter for covariance kernel. Gamma distribution is assigned for \eqn{l} (\eqn{\text{G}(\alpha_l, \beta_l))}.
#'     \code{shape} is shape parameter (default \code{1/8}) and \code{rate} is rate parameter of \code{lenghscale} (default \code{1/8})}
#'     \item{\code{amp}: Prior of amplitude parameter for covariance kernel. Log-normal distribution is assigned for \eqn{\eta}: \eqn{\log(\eta) \sim \mathrm{N}(a_\eta, b_\eta)}
#'      \code{a_amp} is mean(default \code{-1}), and \code{b_amp} is standard deviation(default \code{1})}
#' }
#' }
#' \item{\code{sigma2}}{Error-variance prior hyperparameters. An inverse-gamma prior is assigned to \eqn{\sigma^2}
#'         where \code{shape} is shape parameter and \code{rate} is rate parameter of inverse gamma distribution.
#'         (default \code{shape = 1, rate = 1})}
#'   }
#'
#' @param init Optional named list of initial values. If the values are not assigned, they are randomly sampled from prior.
#' \describe{
#'      \item{\code{index}}{Initial unit index vector \eqn{\theta}. By default, the vector is sampled from the von Mises--Fisher prior.}
#'      \item{\code{link}}{\code{lenghscale} is initial scalar range parameter. (default: \code{0.1})
#'      \code{amp} is initial scalar scale parameter. (default: \code{1})}
#'      \item{\code{sigma2}}{Initial scalar error variance. (default: \code{1})}
#' }
#'
#' @param sampling Logical. If \code{TRUE} (default), run MCMC; otherwise return
#'   prepared nimble model objects without sampling.
#' @param fitted Logical. If \code{TRUE} (default), posterior fitted values are included in the output.
#' Also, if \code{"sampling = FALSE"}, parameters for prediction(\code{c("linkFunction", "Xlin", "lengthscale", "amp")}) is additionally monitored.
#' @param monitors2 Optional character vector of additional monitor nodes. To check the names of the nodes, set \code{fit <- gpFisher(x, y, sampling = FALSE)} and then inspect the variable names stored in the model object using \code{fit$model$getVarNames()}.
#' @param niter Integer. Total MCMC iterations (default \code{10000}).
#' @param nburnin Integer. Burn-in iterations (default \code{1000}).
#' @param thin Integer. Thinning for primary monitors (default \code{1}).
#' @param thin2 Integer. Optional thinning for \code{monitors2} (default \code{1}).
#' @param nchain Integer. Number of MCMC chains (default \code{1}).
#' @param setSeed Logical or numeric argument.  Further details are provided in \link[nimble]{runMCMC}.
#'
#' @details
#' \strong{Model} The single-index model uses Gaussian process with zero mean and and covariance kernel \eqn{\eta \text{exp}(-\frac{(t_i-t_j)^2}{l})} as a link function, where \eqn{t_i, t_j, j = 1, \ldots, n} is index value.
#' Index vector should be length 1.
#'
#' \strong{Priors}
#' \itemize{
#'   \item von Mises–Fisher prior on the index \eqn{\theta}: direction \code{prior$index$direction}, concentration \code{prior$index$dispersion}.
#'   \item Covariance kernel: \eqn{\eta \sim \text{lognormal}(a_\eta, b_\eta)} , \eqn{l \sim \text{G}(\alpha_l, \beta_l)}
#'   \item Error variance \eqn{\sigma^2}: \eqn{IG(a_\sigma, b_\sigma)}.
#'
#' }
#'
#' \strong{Sampling} All sampling parameters' samplers are assigned by nimble.
#'
#' @return A \code{list} typically containing:
#' \describe{
#'   \item{\code{model}}{Nimble model}
#'   \item{\code{sampler}}{Nimble sampler}
#'   \item{\code{sampling}}{Posterior draws of \eqn{\theta}, \eqn{\sigma^2}, and nodes for fitted values by default. Variables specified in \code{monitors2} will be added if provided.}
#'   \item{\code{fitted}}{If \code{fitted = TRUE}, summary values of in-sample fitted values are included.}
#'   \item{\code{input}}{List of data,input values for prior and initial values, and computation time without compiling.}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(20250818)
#' N <- 60; p <- 2
#' x1 <- runif(N, -3, 5)
#' x2 <- runif(N, -3, 5)
#' beta1 <- 0.45; beta2 <- sqrt(1-beta1^2)
#' sigma <- sqrt(0.0036)
#' xlin <- x1*beta1 + x2*beta2
#' eta <- 0.1*xlin + sin(0.5*xlin)^2
#' y <- rnorm(N, eta, sigma)
#' x <- matrix(c(x1, x2), ncol = 2)
#'
#' # One-call version
#' fit <- gpFisher(x = x, y = y, nchain = 3, fitted = TRUE)
#'
#' # Split version
#' models <- gpFisher(x = x, y = y, nchain = 1, sampling = FALSE)
#' Ccompile <- compileModelAndMCMC(models)
#' mcmc.out <- runMCMC(Ccompile$mcmc, niter = 5000, nburnin = 1000, thin = 1,
#'                    nchains = 1, setSeed = TRUE, inits = models$input$init,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#' }
#'
#' @references
#'
#' Antoniadis, A., Grégoire, G., & McKeague, I. W. (2004).
#' Bayesian estimation in single-index models. \emph{Statistica Sinica}, 1147-1164.
#'
#' Choi, T., Shi, J. Q., & Wang, B. (2011).
#' A Gaussian process regression approach to a single-index model.
#' \emph{Journal of Nonparametric Statistics}, 23(1), 21-36.
#'
#' Hornik, K., & Grün, B. (2014). movMF: an R package for fitting mixtures of von Mises-Fisher distributions.
#' \emph{Journal of Statistical Software}, 58, 1-31.
#'
#' @export
#'
gpFisher <- function(x, y,
                     prior = list(
                       index = list(direction = NULL, dispersion = 150),
                       link = list(lengthscale = list(shape = 1/8, rate = 1/8),
                                   amp = list(a_amp = -1, b_amp = 1)),
                       sigma2 = list(shape = 1, rate = 1)),
                     init = list(index = NULL,
                                 link = list(lengthscale = 0.1, amp = 1),
                                 sigma2 = 1),
                     sampling = TRUE, fitted = FALSE,
                     monitors2 = NULL, niter = 10000, nburnin=1000,
                     thin = 1, thin2 = NULL, nchain = 1, setSeed = FALSE){
  start1 <- Sys.time()
  index_temp <- 0; log_amp <- 0; lengthscale <- 0; sigma2 <- 0; Id <- 0
  amp <- 0

  # check sampling, prior, init parameters for independent execution
  checkOutput <- validate_and_finalize_args(
    sampling, fitted, niter, nburnin, thin, thin2, nchain,
    prior, init, "fisher", "gp"
  )
  prior <- checkOutput$priorlist_final
  init <- checkOutput$initlist_final

  # environment
  envobj <- ls(envir=.GlobalEnv)
  on.exit(rm(list=ls(envir=.GlobalEnv)[which(!ls(envir=.GlobalEnv)%in%envobj)],envir=.GlobalEnv))

  .fns <- c(
    "nimNorm","postll_bspline_fisher","nimNorm","rW","besselI_nimble","Stheta",
   "estBeta_fisher","gvcCV",
    # gpSphere
    "expcov_gpSphere","expcovTest_gpSphere","sigma2Sampler_gpSphere",

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

  # model code
  modelCode <- nimbleCode({
    # index
    index_temp[1:p] ~ dvMFnim(theta = theta_prior[1:p])
    index[1:p] <- index_temp[1:p]/sqrt(sum(index_temp[1:p]^2))

    # Linear predictor
    for (i in 1:N){
      Xlin[i] <- sum(X[i,1:p] * index[1:p])
    }

    # hyperprior
    lengthscale ~ dgamma(shape, rate)
    amp ~ dlnorm(a_amp, b_amp)

    # linkFunction
    cov[1:N, 1:N] <- expcov_gpSphere(Xlin[1:N], lengthscale, amp)
    linkFunction[1:N] ~ dmnorm(mu0[1:N], cov = cov[1:N, 1:N])

    # likelihood - gibbs
    sigma2 ~ dinvgamma(a_sig, b_sig)
    Sigma[1:N,1:N] <- sigma2 * Id[1:N,1:N]
    Y[1:N] ~ dmnorm(linkFunction[1:N], cov = Sigma[1:N,1:N])
  })

  # Check hyper-parameters
  ## index
  if (length(prior$index$dispersion) >= 2|| !is.numeric(prior$index$dispersion)){
    stop("Dispersion of index vector should be scalar.")
  }
  if (prior$index$dispersion < 0){
    stop("Dispersion of index vector should be postive")
  }

  dispersion_prior <- prior$index$dispersion

  if (!is.null(prior$index$direction) & length(prior$index$direction) != (p)){
    stop("Prior index vector has incorrect dimension.")
  }
  if (!is.null(prior$index$direction) & sum((prior$index$direction)^2) != 1){
    stop("Prior index vector should be unit vector.")
  }

  if (is.null(prior$index$direction)){
    # index_direct_input <- c(rep(0, p-1), 1)
    index_direct_input <- as.vector(ppr(X, y, nterms = 1)$alpha)
    index_direct_input <- index_direct_input/sqrt(sum(index_direct_input^2))
  } else{
    index_direct_input <- prior$index$direction
  }
  index_direction <- index_direct_input * dispersion_prior


  # lengthscale
  ## shape
  if (is.null(prior$link$lengthscale$shape)||
      length(prior$link$lengthscale$shape) >= 2 ||
      prior$link$lengthscale$shape < 0){
    stop("Prior lengthscale (shape) has incorrect value.")
  } else{
    shape <- prior$link$lengthscale$shape
  }

  ## rate
  if (is.null(prior$link$lengthscale$rate)||
      length(prior$link$lengthscale$rate) >= 2 ||
      prior$link$lengthscale$rate < 0){
    stop("Prior lengthscale (rate) has incorrect value.")
  } else{
    rate <- prior$link$lengthscale$rate
  }

  # amp - log normal
  ## a_amp: mean
  if (is.null(prior$link$amp$a_amp)||
      length(prior$link$amp$a_amp) >= 2 ){
    stop("Prior amp (a_amp) has incorrect value.")
  } else{
    a_amp <- prior$link$amp$a_amp
  }

  ## b_amp: sd
  if (is.null(prior$link$amp$b_amp)||
      length(prior$link$amp$b_amp) >= 2 ||
      prior$link$amp$b_amp < 0){
    stop("Prior amp (b_amp) has incorrect value.")
  } else{
    b_amp <- prior$link$amp$b_amp
  }

  # sigma2: inverse gamma
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

  # initialize - function

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
                       function(j) initFunction_gpFisher(X = X, Y = as.vector(Y),
                                                                        index = init$index,
                                                                        lengthscale = init$link$lengthscale,
                                                                        amp  = init$link$amp,
                                                                        sigma2  = init$sigma2,
                                                                        setSeed = seedNum[j]))
  firstInit <- inits_list[[1]]

  message("Build Model")
  suppressMessages(simpleModel <- nimbleModel(modelCode,
                             data = list(X = X, Y = as.vector(Y)),
                             constants = list(theta_prior = index_direction,
                                              shape = shape, rate = rate,
                                              a_sig = a_sig, b_sig = b_sig,
                                              a_amp = a_amp, b_amp = b_amp,
                                              p = p, N = N,
                                              Id = diag(1, N), mu0 = rep(0, N)),
                             inits = firstInit))

  # sampling
  monitorsList <- c("index", "sigma2")
  if (fitted){
    monitorsList <- c(monitorsList, "linkFunction", "Xlin", "lengthscale", "amp")
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

  # mcmcConf$removeSamplers('lengthscale')
  # mcmcConf$addSampler(target = "lengthscale", type = "RW",
  #                     control = list(log = TRUE))
  #
  mcmcConf$removeSamplers('sigma2')
  mcmcConf$addSampler(target = "sigma2",
                      type = sigma2Sampler_gpSphere)
  #
  # mcmcConf$setSamplerExecutionOrder(c(2, 1, 4, 3, 5))

  message("Build MCMC")
  mcmc1 <- buildMCMC(mcmcConf)
  end1 <- Sys.time()

  if (!sampling){
    mcmc.out <- NULL
    fittedResult <- NULL
    sampMCMC <- NULL


  } else{
    message("Compile Model")
    suppressMessages(CModel <- compileNimble(simpleModel,
                                             resetFunctions = TRUE))

    message("Compile MCMC")
    suppressMessages(Cmcmc <- compileNimble(mcmc1))

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
                          thin = thin, thin2 = thin2, inits = inits_list,
                          nchains = nchain, setSeed = seedNum,
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

    if (fitted){ # posterior fitted value output (mean, median, sd)
      message("Compute posterior fitted value")

      # namesBeta <- paste0("lengthscale", 1:p)
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

    } else{ # 그냥 NULL
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
                       prior = list(index = list(direction = index_direct_input, dispersion = dispersion_prior),
                                    link= list(lengthscale = list(shape = shape, rate = rate),
                                               amp = list(a_amp = a_amp, b_amp = b_amp)),
                                    sigma2 = list(shape = a_sig, rate = b_sig)),
                       init = inits_list,
                       time = time)

  out <- list(model = simpleModel, sampler = mcmc1, sampling = sampMCMC,
              fitted = fittedResult, input = inputOptions,
              modelName = "gpFisher")
  class(out) = "bsimGp"
  return(out)
}


