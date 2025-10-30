#' Bayesian single-index regression with Gaussian process link and unit sphere prior
#'
#'
#' @description
#' Fits a single–index model \eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n} where
#' the index \eqn{\theta} lies on the unit sphere, and the link \eqn{f(\cdot)} is represented
#' with Gaussian process.
#'
#' @param x Numeric data.frame/matrix of predictors. Each row is an observation.
#' @param y Numeric response numeric vector/matrix. Other types  are not available.
#' @param prior Optional named list of prior settings with sublists:
#' \describe{
#' \item{\code{index}}{Nothing to assign.}
#' \item{\code{link}}{
#' \enumerate{
#'    \item{\code{lenghscale}: Prior of length-scale parameter for covariance kernel. Gamma distribution is assigned for \eqn{l}: \eqn{\text{G}(\alpha_l, \beta_l)}
#'     \code{shape} is shape parameter (default \code{1/8}) and \code{rate} is rate parameter of lengthscale \eqn{l}. (default \code{1/8})}
#'     \item{\code{amp}: Prior of amplitude parameter for covariance kernel. Log-normal distribution is assigned for \eqn{\eta}: \eqn{\log(\eta) \sim \mathrm{N}(a_\eta, b_\eta)}
#'      \code{a_amp} is mean (default \code{-1}), and \code{b_amp} is standard deviation (default \code{1})}
#' }
#' }
#' \item{\code{sigma2}}{Error-variance prior hyperparameters. An inverse-gamma prior is assigned to \eqn{\sigma^2}
#'         where \code{shape} is shape parameter and \code{rate} is rate parameter of inverse gamma distribution.
#'         (default \code{shape = 1, rate = 1})}
#'
#'   }
#' @param init Optional named list of initial values. If the values are not assigned, they are randomly sampled from prior.
#' \describe{
#'      \item{\code{index}}{Initial unit index vector. By default, vector is randomly drawn from normal distribution and standardized.}
#'      \item{\code{link}}{\code{lenghscale} is initial scalar range parameter. (default: \code{0.1})
#'      \code{amp} is initial scalar scale parameter. (default: \code{1})
#'      }
#'      \item{\code{sigma2}}{Initial scalar error variance. (default: \code{1})}
#' }
#'
#' @param sampling Logical. If \code{TRUE} (default), run MCMC; otherwise return
#'   prepared nimble model objects without sampling.
#' @param fitted Logical. If \code{TRUE} (default), posterior fitted values are included in the output.
#' Also, if \code{"sampling = FALSE"}, parameters for prediction(\code{c("linkFunction", "Xlin", "lengthscale", "amp")}) is additionally monitored.
#' @param method Character, Gp-uniform model has 3 different types of sampling method, fully Bayesian method (\code{"FB"}), empirical Bayes approach (\code{"EB"}), and empirical Gibbs sampler (\code{"EG"}).
#' Assign one sampler method. Empirical sampling approach is recommended in high-dimensional data. By default, fully Bayesian approach is assigned.
#' @param lowerB Numeric vector of element-wise lower bounds for the \code{"L-BFGS-B"} method.
#' When the empirical Bayes or Gibbs sampler method is used, the marginal likelihood is optimized via \code{optim(method = "L-BFGS-B")}.
#' The vector must be ordered as \code{c(index vector, lengthscale, amp, sigma2)}; note that \code{sigma2} is included only for the empirical Bayes method (omit it for Gibbs).
#' By default, the lower bounds are \code{-1} for the index vector and \code{-1e2} for logarithm of \code{lengthscale}, \code{amp}, and (if present) \code{sigma2}.
#' @param upperB Numeric vector of element-wise upper bounds for the \code{"L-BFGS-B"} method.
#' When the empirical Bayes or Gibbs sampler method is used, the marginal likelihood is optimized via \code{optim(method = "L-BFGS-B")}.
#' The vector must be ordered as \code{c(index vector, lengthscale, amp, sigma2)}; note that \code{sigma2} is included only for the empirical Bayes method (omit it for Gibbs).
#' By default, the upper bounds are \code{1} for the index vector and \code{1e2} for logarithm of \code{lengthscale}, \code{amp}, and (if present) \code{sigma2}.
#' @param monitors2 Optional character vector of additional monitor nodes. To check the names of the nodes, set \code{fit <- gpSphere(x, y, sampling = FALSE)} and then inspect the variable names stored in the model object using \code{fit$model$getVarNames()}.
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
#'   \item Index vector: Uniform prior with \eqn{||\theta|| =1}
#'   \item Covariance kernel: \eqn{\eta \sim \text{lognormal}(a_\eta, b_\eta)} , \eqn{l \sim \text{G}(\alpha_l, \beta_l)}
#'   \item Error variance \eqn{\sigma^2}: \eqn{IG(a_\sigma, b_\sigma)}
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
#' @return A \code{list} typically containing:
#' \describe{
#'   \item{\code{model}}{Nimble model}
#'   \item{\code{sampler}}{Nimble sampler}
#'   \item{\code{sampling}}{Posterior draws of \eqn{\theta}, \eqn{\sigma^2}, and nodes for fitted values by default. Variables specified in \code{monitors2} will be added if provided.}
#'   \item{\code{fitted}}{If \code{fitted = TRUE}, summary values of in-sample fitted values are included.}
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
#' # One-call version
#' fit <- gpSphere(X, y, method = "EB")
#'
#' # Split version
#' model <- gpSphere(X, y, method = "EB", sampling = FALSE)
#' Ccompile <- compileModelAndMCMC(model)
#' mcmc.out <- runMCMC(Ccompile$mcmc, niter = 5000, nburnin = 1000, thin = 1,
#'                    nchains = 1, setSeed = TRUE, inits = model$input$init,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#' }
#'
#' @references
#' Choi, T., Shi, J. Q., & Wang, B. (2011).
#' A Gaussian process regression approach to a single-index model.
#' \emph{Journal of Nonparametric Statistics}, 23(1), 21-36.
#'
#' @export


gpSphere <- function(x, y,
                      prior = list(index = NULL,
                                   link = list(lengthscale = list(shape = 1/8, rate = 1/8),
                                               amp = list(a_amp = -1, b_amp = 1)),
                                   sigma2 = list(shape = 1, rate = 1)),
                      init = list(index = list(index = NULL),
                                  link = list(lengthscale = 0.1, amp = 1),
                                  sigma2 = 1),
                      sampling = TRUE, fitted = TRUE, method = "FB",
                     lowerB = NULL, upperB = NULL, monitors2 = NULL, niter = 10000,
                     nburnin=1000, thin = 1, thin2 = NULL, nchain = 1, setSeed = FALSE){

  start1 <- Sys.time()
  index <- 0; log_amp <- 0; lengthscale <- 0; sigma2 <- 0; Id <- 0
  index0 <- 0; amp <- 0

  # check sampling, prior, init parameters for independent execution
  checkOutput <- validate_and_finalize_args(
    sampling, fitted, niter, nburnin, thin, thin2, nchain,
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

  if (method == "EB" & nchain > 1){
    stop("Only single chain is allowed.")
  }

  # data dimension
  N <- length(Y)
  p <- ncol(X)

  # model code
  modelCode_fb <- nimbleCode({
    # beta - MH
    index0[1:p] ~ dunitSphere(p)
    index[1:p] <- index0[1:p]/sqrt(sum(index0[1:p]^2))

    # hyperprior - MH
    lengthscale ~ dgamma(shape, rate)
    amp ~ dlnorm(a_amp, b_amp)
#
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

  # model code
  modelCode_EG <- nimbleCode({
    # beta - MH
    index[1:p] ~ dunitSphere(p)

    # hyperprior - MH
    lengthscale ~ dgamma(shape, rate)
    amp ~ dlnorm(a_amp, b_amp)
    # log_amp ~ dnorm(a_amp, b_amp)
    # amp <- exp(log_amp)

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

  inits_list <- lapply(seq_len(nchain),
                       function(j) initFunction_gpSphere(N = N, p = p, X = X, Y = as.vector(Y),
                                                       index = init$index$index,
                                                       lengthscale = init$link$lengthscale,
                                                       amp  = init$link$amp,
                                                       sigma2  = init$sigma2, method = method,
                                                       setSeed = seedNum[j]))

  firstInit <- inits_list[[1]]

  # boundarys for optimizing
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

  if (method == "EB"){ # only nchain = 1
    obj_fn <- obj_btt_EB(X, Y, p, shape, rate,
                   a_amp, b_amp, a_sig, b_sig)
    # lowerB <- c(rep(-1, p),-1e2, -1e2, -1e2)
    # upperB <- c(rep(1, p),1e2, 1e2, 1e2)
    tempIndex <- as.vector(ppr(X,Y, nterms = 1)$alpha)

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
                      linkFunction =mvtnorm::rmvnorm(1, rep(0, N), sigma = init_cov)[1,],
                      Sigma = map_sigma2 * diag(1, N))

  }


  if (method == "FB"){
    modelCode <- modelCode_fb
  } else if (method == "EG"){
    modelCode <- modelCode_EG
  } else{
    modelCode <- modelCode_EB
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
  }

  # sampling
  if (method %in% c("FB", "EG")){
    monitorsList <- c("index", "sigma2")
    if (fitted){
      monitorsList <- c(monitorsList, "linkFunction", "Xlin", "lengthscale", "amp")
    }

  } else{ # EB
    monitorsList <- c("linkFunction")
    if (fitted){
      monitorsList <- c(monitorsList, "Xlin")
    }
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
  end1 <- Sys.time()

  if (!sampling){
    mcmc.out <- NULL
    fittedResult <- NULL
    sampMCMC <- NULL


  } else{
    message("Compile Model")
    suppressMessages(
      CModel <- compileNimble(simpleModel, resetFunctions = TRUE)
    )

    message("Compile MCMC")
    suppressMessages(Cmcmc <- compileNimble(mcmc1))

    start2 <- Sys.time()
    message("Run MCMC")
    mcmc.out <- NULL
    if (setSeed == FALSE){
      seedNum <- setSeed
    }
    if (is.null(monitors2)){
      if (method == "EB"){
        mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                            thin = thin,
                            nchains = nchain, setSeed = seedNum,
                            summary = FALSE, samplesAsCodaMCMC = TRUE)
      } else if (method == "EG"){
        tryCatch(
          {
            mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                                thin = thin,
                                nchains = nchain, setSeed = seedNum, inits = inits_list,
                                summary = FALSE, samplesAsCodaMCMC = TRUE)
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
      } else{
        mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                            thin = thin,
                            nchains = nchain, setSeed = seedNum, inits = inits_list,
                            summary = FALSE, samplesAsCodaMCMC = TRUE)
      }

    } else{
      if (method == "EB"){
        mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                            thin = thin, thin2 = thin2,
                            nchains = nchain, setSeed = seedNum,
                            summary = FALSE, samplesAsCodaMCMC = TRUE)
      } else if (method == "EG"){
        tryCatch(
          {
            mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                                thin = thin, thin2 = thin2, inits = inits_list,
                                nchains = nchain, setSeed = seedNum,
                                summary = FALSE, samplesAsCodaMCMC = TRUE)
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
      } else{
        mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                            thin = thin, thin2 = thin2, inits = inits_list,
                            nchains = nchain, setSeed = seedNum,
                            summary = FALSE, samplesAsCodaMCMC = TRUE)
      }

    }
    samples <- NULL
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
      if (method == "EB"){
        sigma2_samples <- rep(map_sigma2, nrow(LinkFunction_samples))
      } else{
        sigma2_samples <- samples[, namesSigma]
      }
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
                       prior = list(index = NULL,
                                    link = list(lengthscale = list(shape = shape, rate = rate),
                                                amp = list(a_amp = a_amp, b_amp = b_amp)),
                                    sigma2 = list(shape = a_sig, rate = b_sig)),
                       init = inits_list,
                       time = time)

  out <- list(model = simpleModel, sampler = mcmc1, sampling = sampMCMC,
              fitted = fittedResult, input = inputOptions,
              modelName = "gpSphere")
  class(out) = "bsimGp"
  return(out)
}
