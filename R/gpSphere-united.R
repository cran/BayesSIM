#' Bayesian Single-Index Regression with Gaussian Process Link and Unit Sphere Prior
#'
#' @description
#' Fits a single–index model \eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n}, where
#' the index \eqn{\theta} lies on the unit sphere, and the link \eqn{f(\cdot)} is represented
#' with Gaussian process.
#'
#' @inheritParams bsFisher
#' @param method Character, `gpSphere` model has 3 different types of sampling method, fully Bayesian method (\code{"FB"}), empirical Bayes approach (\code{"EB"}), and empirical Gibbs sampler (\code{"EG"}).
#' Assign one sampler method. Empirical sampling approach is recommended in high-dimensional data. By default, fully Bayesian approach is assigned.
#' @param lowerB This parameter is only for `gpSphere` model. Numeric vector of element-wise lower bounds for the \code{"L-BFGS-B"} method.
#' When the empirical Bayes or Gibbs sampler method is used, the marginal likelihood is optimized via \code{optim(method = "L-BFGS-B")}.
#' The vector must be ordered as \code{c(index vector, lengthscale, amp, sigma2)}. Note that \code{sigma2} is included only for the empirical Bayes method (omit it for Gibbs).
#' By default, the lower bounds are \code{-1} for the index vector and \code{-1e2} for logarithm of \code{lengthscale}, \code{amp}, and (if present) \code{sigma2}.
#' @param upperB This parameter is only for `gpSphere` model. Numeric vector of element-wise upper bounds for the \code{"L-BFGS-B"} method.
#' When the empirical Bayes or Gibbs sampler method is used, the marginal likelihood is optimized via \code{optim(method = "L-BFGS-B")}.
#' The vector must be ordered as \code{c(index vector, lengthscale, amp, sigma2)}. Note that \code{sigma2} is included only for the empirical Bayes method (omit it for Gibbs).
#' By default, the upper bounds are \code{1} for the index vector and \code{1e2} for logarithm of \code{lengthscale}, \code{amp}, and (if present) \code{sigma2}.
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
#' \item Inverse-Gamma prior on \eqn{\sigma^2} with shape parameter \eqn{a_\sigma} and rate parameter \eqn{b_\sigma}.
#'
#' }
#'
#' \strong{Sampling} In the fully Bayesian approach, \eqn{\theta}, \eqn{l}, and \eqn{\eta}
#' are updated via the Metropolis–Hastings algorithm, while \eqn{f} and
#' \eqn{\sigma^2} are sampled using Gibbs sampling.
#'
#' In the empirical Bayes approach, \eqn{\theta}, \eqn{l}, \eqn{\eta},
#' and \eqn{\sigma^2} are estimated by maximum a posteriori (MAP), and
#' \eqn{f} is sampled from its full conditional posterior distribution.
#'
#' In the empirical Gibbs sampler, \eqn{\theta}, \eqn{l}, and \eqn{\eta}
#' are estimated by MAP, whereas \eqn{f} and \eqn{\sigma^2} are sampled
#' via Gibbs sampling.
#'
#' For estimation via MAP, effective sample size or potential scale reduction factor is meaningless.
#'
#' \strong{Prior hyper-parameters}
#' These are the prior hyper-parameters set in the function. You can define new values for each parameter in \link{prior_param}.
#' \enumerate{
#' \item Index vector: Nothing to assign.
#' \item Link function:
#' \itemize{
#'    \item{Length-scale:Gamma distribution is assigned for length-scale parameter, \eqn{l}.
#'     \code{link_lengthscale_shape} is shape parameter (default \code{1/8}) and \code{link_lengthscale_rate} is rate parameter of \code{lengthscale}. (default \code{1/8})}
#'     \item{Amplitude: Log-normal distribution is assigned for amplitude parameter, \eqn{\eta}.
#'      \code{link_amp_a} is mean (default \code{-1}), and \code{link_amp_b} is variance. (default \code{1})}
#' }
#'
#' \item Error variance (\code{sigma2}): inverse gamma prior is assigned to \eqn{\sigma^2}
#'         where \code{sigma2_shape} is shape parameter and \code{sigma2_rate} is rate parameter of inverse gamma distribution.
#'         (default \code{sigma2_shape = 1, sigma2_rate = 1})
#'
#' }
#'
#' \strong{Initial values}
#' These are the initial values set in the function. You can define new values for each initial value in \link{init_param}.
#' \enumerate{
#' \item Index vector (\code{index}): Initial unit index vector. By default, vector is randomly drawn from normal distribution and standardized.
#'      \item Link function: \code{link_lengthscale} is initial scalar length-scale parameter. (default: \code{0.1})
#'      \code{link_amp} is initial scalar amplitude parameter. (default: \code{1})
#'      \item Error variance (\code{sigma2}): Initial scalar error variance. (default: \code{1})
#'
#' }
#'
#'
#' @inherit bsFisher return
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
#' fit1 <- gpSphere(y ~ ., method = "EB", data = simdata,
#'                  niter = 1000, nburnin = 100)
#'
#' # Split version
#' models <- gpSphere_setup(y ~ ., method = "EB", data = simdata)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 1000, nburnin = 100, thin = 1,
#'                     nchains = 1, setSeed = TRUE, inits = initList,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#' fit2 <- as_bsim(models, mcmc.out)
#' # The estimates computed by MAP - standard error of the esitmate is meaningless.
#' summary(fit2)
#' }
#' @references
#' Choi, T., Shi, J. Q., & Wang, B. (2011).
#' A Gaussian process regression approach to a single-index model.
#' \emph{Journal of Nonparametric Statistics}, 23(1), 21-36.
#' @name gpSphere
#' @export
gpSphere <- function(formula, data,
                     prior = NULL,
                     init = NULL,
                     method = "FB",
                     lowerB = NULL, upperB = NULL, monitors = NULL, niter = 10000,
                     nburnin=1000, thin = 1, nchain = 1, setSeed = FALSE){
  return(gpSphere.default(formula = formula, data = data,
                          prior = prior, init = init,
                          sampling = TRUE, method = method,
                          lowerB = lowerB, upperB = upperB, monitors = monitors,
                          niter = niter,
                          nburnin=nburnin, thin = thin, nchain = nchain, setSeed = setSeed))

}

#' @rdname gpSphere
#' @export
gpSphere_setup <- function(formula, data,
                           prior = NULL,
                           init = NULL,
                           method = "FB",
                           lowerB = NULL, upperB = NULL, monitors = NULL, niter = 10000,
                           nburnin=1000, thin = 1, nchain = 1, setSeed = FALSE){

  return(gpSphere.default(formula = formula, data = data,
                          prior = prior, init = init,
                          sampling = FALSE, method = method,
                          lowerB = lowerB, upperB = upperB, monitors = monitors,
                          niter = niter,
                          nburnin=nburnin, thin = thin, nchain = nchain, setSeed = setSeed))
}


gpSphere.default <- function(formula, data,
                             prior = NULL,
                             init = NULL,
                             sampling = TRUE, method = "FB",
                             lowerB = NULL, upperB = NULL, monitors = NULL, niter = 10000,
                             nburnin=1000, thin = 1, nchain = 1, setSeed = FALSE){

  start1 <- Sys.time()
  index <- 0; log_amp <- 0; lengthscale <- 0; sigma2 <- 0; Id <- 0
  index0 <- 0; amp <- 0

  # check sampling, prior, init parameters for independent execution
  checkOutput <- validate_and_finalize_args(
    sampling, niter, nburnin, thin, nchain,
    prior, init, "sphere", "gp"
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

    # gpSphere
    "expcov_gpSphere","expcovTest_gpSphere","conBeta","obj_btt",
    "obj_btt_EB","pred_gpSphere","indexSampler_gpSphere","sigma2Sampler_gpSphere",
    "optSampler",

    # distribution
    "dunitSphere", "runitSphere",

    # utils
    "pred_fitted"
  )

  pkg <- "BayesSIM"
  ns <- asNamespace(pkg)
  list2env(mget(.fns, envir = ns, inherits = FALSE), envir = globalenv())

  suppressMessages(
    nimble::registerDistributions(list(
      dunitSphere = list(
        BUGSdist     = "dunitSphere(dim)",
        types        = c("value = double(1)",
                         "dim   = double(0)"),
        discrete     = FALSE
      )
    ), verbose = FALSE)

  )



  # check data dimension
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

  if (method == "EB" & nchain > 1){
    stop("Only single chain is allowed.")
  }

  # data dimension
  N <- length(Y)
  p <- ncol(X)

  # model code
  modelCode_FB <- nimbleCode({
    # beta - MH
    index0[1:p] ~ dunitSphere(p)
    index[1:p] <- index0[1:p]/sqrt(sum(index0[1:p]^2))

    # hyperprior - MH
    lengthscale ~ dgamma(shape, rate)
    amp ~ dlnorm(a_amp, b_amp)

    # Linear predictor
    for (i in 1:N){
      Xlin[i] <- sum(X[i,1:p] * index[1:p])
    }

    sigma2 ~ dinvgamma(a_sig, b_sig)

    # linkFunction
    cov[1:N, 1:N] <- expcov_gpSphere(Xlin[1:N], lengthscale, amp)
    linkFunction[1:N] ~ dmnorm(mu0[1:N], cov = cov[1:N, 1:N])

    # likelihood - gibbs
    Sigma[1:N,1:N] <- sigma2 * Id[1:N,1:N]
    Y[1:N] ~ dmnorm(linkFunction[1:N], cov = Sigma[1:N,1:N])

  })

  modelCode_EG <- nimbleCode({
    # index - MH
    index[1:p] ~ dunitSphere(p)

    # hyperprior - MH
    lengthscale ~ dgamma(shape, rate)
    amp ~ dlnorm(a_amp, b_amp)

    # Linear predictor
    for (i in 1:N){
      Xlin[i] <- sum(X[i,1:p] * index[1:p])
    }

    sigma2 ~ dinvgamma(a_sig, b_sig)

    # linkFunction
    cov[1:N, 1:N] <- expcov_gpSphere(Xlin[1:N], lengthscale, amp)
    linkFunction[1:N] ~ dmnorm(mu0[1:N], cov = cov[1:N, 1:N])

    # likelihood - gibbs
    Sigma[1:N,1:N] <- sigma2 * Id[1:N,1:N]
    Y[1:N] ~ dmnorm(linkFunction[1:N], cov = Sigma[1:N,1:N])

  })

  modelCode_EB <- nimbleCode({
    # MAP method is used
    # Linear predictor
    for (i in 1:N){
      Xlin[i] <- sum(X[i,1:p] * index[1:p])
    }
    # linkFunction
    cov[1:N, 1:N] <- expcov_gpSphere(Xlin[1:N], lengthscale, amp)
    linkFunction[1:N] ~ dmnorm(mu0[1:N], cov = cov[1:N, 1:N])

    # likelihood - gibbs
    Sigma[1:N,1:N] <- sigma2 * Id[1:N,1:N]
    Y[1:N] ~ dmnorm(linkFunction[1:N], cov = Sigma[1:N,1:N])

  })


  # Check hyper-parameters
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



  # For all initial values
  inits_list <- lapply(seq_len(nchain),
                       function(j) initFunction_gpSphere(N = N, p = p, X = X, Y = as.vector(Y),
                                                         index = init$index,
                                                         lengthscale = init$link$lengthscale,
                                                         amp  = init$link$amp,
                                                         sigma2  = init$sigma2, method = method,
                                                         setSeed = seedNum[j]))

  firstInit <- inits_list[[1]]

  # boundarys for optimizing (using MAP) ---------------------------------------
  if (method %in% c("EB", "EG")){
    # lower bound
    if (!is.null(lowerB) & method == "EB" &  length(lowerB) != (p+3)){
      stop("The length of lower bound (lowerB) is not appropriate!")
    }

    if (!is.null(lowerB) & method == "EG" &  length(lowerB) != (p+2)){
      stop("The length of lower bound (lowerB) is not appropriate!")
    }


    if (is.null(lowerB) & method == "EB"){ # index, lengthscale, amp, sigma2
      lowerB <- c(rep(-1, p),-1e3, -1e3, -1e3)
    } else if (is.null(lowerB) & method == "EG"){  # index, lengthscale, amp
      lowerB <- c(rep(-1, p),-1e3, -1e3)
    }


    # Upper bound
    if (!is.null(upperB) & method == "EB" & length(upperB) != (p+3)){
      stop("The length of upper bound (upperB) is not appropriate!")
    }

    if (!is.null(upperB) & method == "EG" & length(upperB) != (p+2)){
      stop("The length of upper bound (upperB) is not appropriate!")
    }


    if (is.null(upperB) & method == "EB"){ # index, lengthscale, amp, sigma2
      upperB <- c(rep(1, p),1e3, 1e3, 1e3)
    } else if (is.null(upperB) & method == "EG"){  # index, lengthscale, amp
      upperB <- c(rep(1, p),1e3, 1e3)
    }
  }

  # Initial values -------------------------------------------------------
  if (method == "EB"){ # only nchain = 1
    obj_fn <- obj_btt_EB(X, Y, p, shape, rate,
                   a_amp, b_amp, a_sig, b_sig)
    # lowerB <- c(rep(-1, p),-1e2, -1e2, -1e2)
    # upperB <- c(rep(1, p),1e2, 1e2, 1e2)
    tempIndex <- as.vector(ppr(X,Y, nterms = 1)$alpha) # initial value for optimizing

    current <- c(tempIndex, log(firstInit$lengthscale), log(firstInit$amp),
                 firstInit$sigma2)

    message("Find MAP..")

    tryCatch(
      {
        optRes  <- nimOptim(par    = current,
                            fn     = obj_fn$run,
                            method = "L-BFGS-B",
                            lower  = lowerB,
                            upper  = upperB)
      },
      error = function(e) {
        stop(
          paste0(
            "Optimization was not successful.\n",
            "Likely cause: the specified lower/upper bounds are inappropriate for this dataset.\n",
            "Recommendation: adjust the boundarys and verify initial values to enable convergence.\n",
          "Specific error message in optim:"
            ), conditionMessage(e), call. = FALSE)
      }
    )

    proposedIndex <- optRes$par[1:p]
    if (proposedIndex[1] < 0){
      proposedIndex <- (proposedIndex) * (-1)
    }
    map_index  <- proposedIndex/sqrt(sum((proposedIndex)^2))
    map_lengthscale <- exp(optRes$par[(p+1)])
    map_amp   <- exp(optRes$par[(p+2)])
    map_sigma2 <- exp(optRes$par[(p+3)])

    # EB initlist
    init_Xlin <- as.vector(X %*% matrix(map_index, nrow = p))
    init_cov <-  expcov_gpSphere(init_Xlin, map_lengthscale, map_amp)
    firstInit <- list(Xlin = init_Xlin,
                      cov = init_cov,
                      linkFunction = mvtnorm::rmvnorm(1, rep(0, N), sigma = init_cov)[1,],
                      Sigma = map_sigma2 * diag(1, N))
    inits_list <- firstInit
    modelCode <- modelCode_EB

  } else{ # FB, EG

    if (method == "FB"){
      modelCode <- modelCode_FB
    } else if (method == "EG"){
      modelCode <- modelCode_EG
    }
  }

  message("Build Model")

  if (method %in% c("FB", "EG")){
    suppressMessages(
      simpleModel <- nimbleModel(modelCode,
                               data = list(X = X, Y = as.vector(Y)),
                               constants = list(shape = shape, rate = rate,
                                                a_sig = a_sig, b_sig = b_sig,
                                                a_amp = a_amp, b_amp = b_amp,
                                                p = p, N = N,
                                                Id = diag(1, N), mu0 = rep(0, N)),
                               inits = firstInit)
    )
    monitorsList <- c("index", "sigma2", "linkFunction", "Xlin", "lengthscale", "amp")

  } else{ # EB
    suppressMessages(
      simpleModel <- nimbleModel(modelCode,
                                 data = list(X = X, Y = as.vector(Y)),
                                 constants = list(p = p, N = N, index = map_index,
                                                  lengthscale = map_lengthscale,
                                                  amp = map_amp, sigma2 = map_sigma2,
                                                  Id = diag(1, N), mu0 = rep(0, N)),
                                 inits = firstInit)
    )
    monitorsList <- c("linkFunction", "Xlin")
  }

  # sampling
  suppressMessages(mcmcConf <- configureMCMC(simpleModel,
                                             monitors = monitorsList,
                                             print = FALSE))

  if (method == "FB"){
    mcmcConf$removeSamplers('lengthscale')
    mcmcConf$addSampler(target = "lengthscale", type = "RW",
                        control = list(log = TRUE))

    mcmcConf$removeSamplers('sigma2')
    mcmcConf$addSampler(target = "sigma2",
                        type = sigma2Sampler_gpSphere)
    mcmcConf$setSamplerExecutionOrder(c(1, 2, 4, 3, 5))

  } else if (method == "EG"){ # Empirical gibbs
    mcmcConf$removeSamplers(c('index','lengthscale','amp'))
    mcmcConf$addSampler(target = c('index','lengthscale','amp'),
                        type   = optSampler,
                        control = list(lowerB = lowerB,
                                       upperB = upperB))
    mcmcConf$setSamplerExecutionOrder(c(3, 2, 1))

  } else if (method == "EB"){ # Empirical bayes
    # mcmcConf$removeSamplers(c('index','lengthscale','amp','sigma2'))
    # mcmcConf$addSampler(target = c('index','lengthscale','amp','sigma2'),
    #                     type   = optSampler,
    #                     control = list(sigmaTrue = TRUE))
    # mcmcConf$setSamplerExecutionOrder(c(2, 1))

  } else{
    stop("Wrong sampling method!")
  }

  message("Build MCMC")
  mcmc1 <- buildMCMC(mcmcConf)


  if (!sampling){
    samples <- NULL

  } else{
    start2 <- Sys.time()
    message("Compile Model")
    suppressMessages(
      CModel <- compileNimble(simpleModel, resetFunctions = TRUE)
    )
    message("Compile MCMC")
    suppressMessages(Cmcmc <- compileNimble(mcmc1))
    end2 <- Sys.time()

    message("Run MCMC")
    mcmc.out <- NULL

    if (setSeed == FALSE){
      seedNum <- setSeed
    }
    if (method == "EB"){
      mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                          thin = thin,nchains = nchain, setSeed = seedNum,
                          inits = inits_list,
                          summary = FALSE, samplesAsCodaMCMC = TRUE)
    } else if (method == "EG"){
      tryCatch(
        {
          mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                              thin = thin,nchains = nchain, setSeed = seedNum,
                              inits = inits_list, summary = FALSE, samplesAsCodaMCMC = TRUE)
          },
          error = function(e) {
            stop(paste0(
                "Optimization was not successful.\n",
                "Likely cause: the specified lower/upper bounds are inappropriate for this dataset.\n",
                "Recommendation: adjust the boundarys and verify initial values to enable convergence.\n",
                "Specific error message in optim:"
              ), conditionMessage(e), call. = FALSE)
          }
        )
    } else{ # FB
      mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                          thin = thin, nchains = nchain, setSeed = seedNum,
                          inits = inits_list, summary = FALSE, samplesAsCodaMCMC = TRUE)
      }

    sampMCMC <- mcmc.out
    if (method == "EB"){
      nmcmcsamp <- nrow(sampMCMC)
      sampMCMC <- as.matrix(mcmc.out)
      index_mat <- matrix(rep(map_index, nmcmcsamp), ncol = p, byrow = TRUE)
      colnames(index_mat) <- paste0("index[",1:p, "]")
      hype_mat <- matrix(rep(c(map_lengthscale, map_amp, map_sigma2), nmcmcsamp),
                         ncol = 3, byrow = TRUE)
      colnames(hype_mat) <- c("lengthscale", "amp", "sigma2")
      sampMCMC <- cbind(sampMCMC, index_mat, hype_mat)
    }

    # Make single matrix for estimation
    samples <- NULL
    if (nchain > 1){
      for (i in 1:nchain){
        samples <- rbind(samples, sampMCMC[[i]])
      }
    } else if (nchain == 1){
      samples <- sampMCMC
    }
  }
  end1 <- Sys.time()


  ## Input options
  if (!sampling){
    time <- NULL
  } else{
    time <- difftime(end1, start1, units = "secs") - difftime(end2, start2, units = "secs")
  }

  # Results ------------------------------------------------------------------
  #Model point estimation
  modelESTLIST <- bsimFit_pointest(samples, X, Y)

  inputOptions <- list(origdata = list(x = X, y = Y), formula = formula,
                       prior = list(index = NULL,
                                    link = list(lengthscale = list(shape = shape, rate = rate),
                                                amp = list(a_amp = a_amp, b_amp = b_amp)),
                                    sigma2 = list(shape = a_sig, rate = b_sig)),
                       init = inits_list,
                       samplingOptions = list(lowerB = lowerB, upperB = upperB,
                                              sampling = sampling, method = method,
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
                samples = sampMCMC, # return of runMCMC
                input = inputOptions,
                defModel = simpleModel, defSampler = mcmc1,
                modelName = "gpSphere")

    class(out) = "bsim"


  } else{
    out <- list(input = inputOptions,
                defModel = simpleModel,
                defSampler = mcmc1,
                modelName = "gpSphere")

    class(out) = "bsimSetup"

  }
  return(out)


}
