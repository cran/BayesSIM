#' Bayesian Single-Index Regression with Gaussian Process Link and von Mises-Fisher Prior
#'
#' @description
#' Fits a single–index model \eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n}, where
#' the index \eqn{\theta} lies on the unit sphere with von Mises-Fisher prior, and the link \eqn{f(\cdot)} is represented
#' with Gaussian process.
#'
#' @inheritParams bsFisher
#'
#' @details
#' \strong{Model} The single-index model uses Gaussian process with zero mean and and covariance kernel \eqn{\eta \cdot \text{exp}(-\frac{(t_i-t_j)^2}{l})} as a link function, where \eqn{t_i, t_j, j = 1, \ldots, n} is index value.
#' Index vector should be length 1.
#'
#' \strong{Priors}
#' \itemize{
#' \item von Mises–Fisher prior on the index \eqn{\theta} with direction and concentration.
#' \item Covariance kernel: Amplitude, \eqn{\eta}, follows log normal distribution with mean \eqn{a_\eta} and variance \eqn{b_\eta}.
#' Length-scale parameter follows gamma distribution with shape parameter \eqn{\alpha_l} and rate parameter \eqn{\beta_l}.
#' \item Inverse gamma prior on \eqn{\sigma^2} with shape parameter \eqn{a_\sigma} and rate parameter \eqn{b_\sigma}.
#'
#' }
#'
#' \strong{Sampling} All sampling parameters' samplers are assigned by nimble.
#'
#' \strong{Prior hyper-parameters}
#' These are the prior hyper-parameters set in the function. You can define new values for each parameter in \link{prior_param}.
#' \enumerate{
#'   \item Index vector: von Mises--Fisher prior for the projection vector \eqn{\theta} (\code{index}).
#'         \code{index_direction} is a unit direction vector of the von Mises--Fisher distribution.
#'         By default, the estimated vector from projection pursuit regression is assigned.
#'         \code{index_dispersion} is the positive concentration parameter. By default, \code{150} is assigned.
#'
#' \item Link function:
#' \itemize{
#'    \item{Length-scale:Gamma distribution is assigned for length-scale parameter, \eqn{l}.
#'     \code{link_lengthscale_shape} is shape parameter (default \code{1/8}) and \code{link_lengthscale_rate} is rate parameter of \code{lengthscale}. (default \code{1/8})}
#'     \item{Amplitude: Log-normal distribution is assigned for amplitude parameter, \eqn{\eta}.
#'      \code{link_amp_a} is mean (default \code{-1}), and \code{link_amp_b} is variance. (default \code{1})}
#' }
#'
#' \item Error variance (\code{sigma2}): An inverse-gamma prior is assigned to \eqn{\sigma^2}
#'         where \code{sigma2_shape} is shape parameter and \code{sigma2_rate} is rate parameter of inverse gamma distribution.
#'         (default \code{sigma2_shape = 1, sigma2_rate = 1})
#'         }
#'
#'
#' \strong{Initial values}
#' These are the initial values set in the function. You can define new values for each initial value in \link{init_param}.
#' \enumerate{
#'       \item Index vector (\code{index}): Initial unit index vector \eqn{\theta}. By default, the vector is sampled from the von Mises--Fisher prior.
#'      \item Link function: \code{link_lengthscale} is initial scalar length-scale parameter. (default: \code{0.1})
#'      \code{link_amp} is initial scalar amplitude parameter. (default: \code{1})
#'      \item Error variance (\code{sigma2}): Initial scalar error variance. (default: \code{1})
#' }
#'
#' @inherit bsFisher return
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' N <- 60; p <- 2
#' x1 <- runif(N, -3, 5)
#' x2 <- runif(N, -3, 5)
#' beta1 <- 0.45; beta2 <- sqrt(1-beta1^2)
#' sigma <- sqrt(0.0036)
#' xlin <- x1*beta1 + x2*beta2
#' eta <- 0.1*xlin + sin(0.5*xlin)^2
#' y <- rnorm(N, eta, sigma)
#' x <- matrix(c(x1, x2), ncol = 2)
#' simdata <- data.frame(x = x, y = y)
#' colnames(simdata) <- c("X1", "X2", "y")
#'
#' # One tool version
#' fit1 <- gpFisher(y ~ ., data = simdata, nchain = 1, niter = 1000, nburnin = 100)
#'
#' # Split version
#' models <- gpFisher_setup(y ~ ., data = simdata, nchain = 1)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 1000, nburnin = 100, thin = 1,
#'                    nchains = 1, setSeed = TRUE, inits = initList,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#' fit2 <- as_bsim(models, mcmc.out)
#' summary(fit2)
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
#' @name gpFisher
#' @export
gpFisher <- function(formula, data,
                     prior = NULL,
                     init = NULL,
                     monitors = NULL, niter = 10000, nburnin=1000,
                     thin = 1, nchain = 1, setSeed = FALSE){
  return(gpFisher.default(formula = formula, data = data,
                          prior = prior,
                          init = init, sampling = TRUE,
                          monitors = monitors, niter = niter, nburnin = nburnin,
                          thin = thin, nchain = nchain, setSeed = setSeed))
}

#' @rdname gpFisher
#' @export
gpFisher_setup <- function(formula, data,
                           prior = NULL,
                           init = NULL,
                           monitors = NULL, niter = 10000, nburnin=1000,
                           thin = 1, nchain = 1, setSeed = FALSE){
  return(gpFisher.default(formula = formula, data = data,
                          prior = prior,
                          init = init, sampling = FALSE,
                          monitors = monitors, niter = niter, nburnin = nburnin,
                          thin = thin, nchain = nchain, setSeed = setSeed))
}

gpFisher.default <- function(formula, data,
                             prior = NULL,
                             init = NULL,
                             sampling = TRUE, monitors = NULL, niter = 10000, nburnin=1000,
                             thin = 1, nchain = 1, setSeed = FALSE){
  start1 <- Sys.time()
  index_temp <- 0; log_amp <- 0; lengthscale <- 0; sigma2 <- 0; Id <- 0
  amp <- 0

  # check sampling, prior, init parameters for independent execution
  checkOutput <- validate_and_finalize_args(
    sampling, niter, nburnin, thin, nchain,
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

   # distribution
   "dvMFnim", "rvMFnim",

    # utils
    "pred_fitted"
  )

  pkg <- "BayesSIM"
  ns <- asNamespace(pkg)
  list2env(mget(.fns, envir = ns, inherits = FALSE), envir = globalenv())

  suppressMessages(
    nimble::registerDistributions(list(
      dvMFnim = list(
        BUGSdist = "dvMFnim(theta)",
        types    = c("value = double(1)", "theta = double(1)"),
        discrete = FALSE
      )
    ), verbose = FALSE)
  )


  # check data dimension
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
    index_direct_input <- as.vector(ppr(X, Y, nterms = 1)$alpha)
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
  monitorsList <- c("index", "sigma2", "linkFunction", "Xlin", "lengthscale", "amp",
                    "index_temp")

  suppressMessages(mcmcConf <- configureMCMC(simpleModel,
                                             monitors = monitorsList,
                                             print = FALSE))

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


  if (!sampling){
    samples <- NULL

  } else{
    start2 <- Sys.time()
    message("Compile Model")
    suppressMessages(CModel <- compileNimble(simpleModel,
                                             resetFunctions = TRUE))

    message("Compile MCMC")
    suppressMessages(Cmcmc <- compileNimble(mcmc1))
    end2 <- Sys.time()

    message("Run MCMC")
    mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                        thin = thin,
                        nchains = nchain, setSeed = seedNum, inits = inits_list,
                        summary = FALSE, samplesAsCodaMCMC = TRUE)
    # sampMCMC <- mcmc.out
    samples <- sampleBind(mcmc.out, nchain)


  }
  end1 <- Sys.time()

  # inputOptions <- NULL
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
                       prior = list(index = list(direction = index_direct_input, dispersion = dispersion_prior),
                                    link= list(lengthscale = list(shape = shape, rate = rate),
                                               amp = list(a_amp = a_amp, b_amp = b_amp)),
                                    sigma2 = list(shape = a_sig, rate = b_sig)),
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
                modelName = "gpFisher")

    class(out) = "bsim"


  } else{
    out <- list(input = inputOptions,
                defModel = simpleModel,
                defSampler = mcmc1,
                modelName = "gpFisher")

    class(out) = "bsimSetup"

  }
  return(out)



}


