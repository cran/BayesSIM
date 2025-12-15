#' Bayesian Single-Index Regression with B-Spline Link and Spike-and-Slab Prior
#'
#' @description Fits a single-index model \eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n}
#' where the link \eqn{f(\cdot)} is represented by B-spline link function and the
#' index vector \eqn{\theta} has spike-and-slab prior.
#' @inheritParams bsFisher
#'
#' @details
#' \strong{Model} The single–index model uses a \eqn{m}-order polynomial spline with \eqn{k} interior knots as follows:
#' \deqn{f(t) = \sum_{j=1}^{m+k} B_j(t)\,\beta_j} on \eqn{[a, b]} with \eqn{t_i = X_i' \theta, i = 1,\cdots, n}
#' and \eqn{\|\theta\|_2 = 1}. \eqn{\{\beta_j\}_{j=1}^{m+k}} are spline coefficient and \eqn{a_\theta} and \eqn{ b_\theta} are boundary knots where \eqn{a_\theta = min(t_i, i = 1, \cdots, n) - \delta},
#' and \eqn{b_\theta = max(t_i, i = 1,\cdots, n) + \delta}. \eqn{\theta} is a p-dimensional index vector subject to a spike-and-slab
#' prior for variable selection with binary indicator variable \eqn{\nu}.
#'
#' \strong{Priors}
#' \itemize{
#' \item The variable selection indicator \eqn{\nu} has Beta–Bernoulli hierarchy prior. Beta hyper-prior on the Bernoulli–inclusion probability \eqn{w},
#'     inducing \eqn{p(\nu) \propto \mathrm{Beta}(r_1 + n_\nu, r_2 + p - n_\nu)} where  \eqn{n_\nu = \Sigma_{i=1}^{p}I(\nu_i = 1)}.
#'     \eqn{r_1, r_2} are shape and rate parameter of beta distribution.
#' \item Slab coefficients \eqn{\theta} have normal distribution with zero mean and \eqn{\sigma_\theta^2} variance.
#' \item Conditioned on \eqn{\theta} and \eqn{\sigma^2}, the link coefficients \eqn{\beta = (\beta_1,\ldots,\beta_{m+k})^\top} follow
#'        normal distribution with estimated mean vector \eqn{\hat{\beta}_{\theta} = (X_{\theta}'X_{\theta})^{-1}X_{\theta}'Y} and
#'        covariance \eqn{\sigma^2 (X_{\theta}^\top X_{\theta})^{-1}}, where \eqn{X_{\theta}} is the B-spline basis design matrix.
#' \item Inverse gamma prior on \eqn{\sigma^2} with shape parameter \eqn{a_\sigma} and rate parameter \eqn{b_\sigma}.
#' }
#'
#' \strong{Sampling}
#' Samplers are automatically assigned by nimble.
#'
#' \strong{Prior hyper-parameters}
#' These are the prior hyper-parameters set in the function. You can define new values for each parameter in \link{prior_param}.
#' \enumerate{
#' \item Index vector: \code{index_nu_r1, index_nu_r2} gives the shape and rate parameter of beta-binomial prior, respectively.
#'     For slab prior, normal distribution with zero mean is assigned for selected variables \eqn{\theta}. \code{index_sigma_theta} is for variance of \eqn{\theta}, and it is assigned 0.25 by default.
#'     \item Link function: B-spline basis and coefficient of B-spline setup.
#'          \itemize{
#'          \item{basis: For the basis of B-spline, \code{link_basis_df} is the number of basis functions (default \code{21}), \code{link_basis_degree} is the spline degree (default \code{2}) and \code{link_basis_delta} is a small jitter for boundary-knot spacing control (default \code{0.01}).}
#'          \item{beta: For the coefficient of B-spline, multivariate normal prior is assigned with mean \code{link_beta_mu}, and covariance \code{link_beta_cov}. By default, \eqn{\mathcal{N}_p\!\big(0, \mathrm{I}_p\big)}}
#' }
#'   \item Error variance (\code{sigma2}):  Inverse gamma prior is assigned to \eqn{\sigma^2}
#'         where \code{sigma2_shape} is shape parameter and \code{sigma2_rate} is rate parameter of inverse gamma distribution.
#'         (default \code{sigma2_shape = 0.001, sigma2_rate = 100})
#'         }
#'
#' \strong{Initial values}
#' These are the initial values set in the function. You can define new values for each initial value in \link{init_param}.
#' \enumerate{
#' \item Index vector:
#'      \itemize{
#'      \item{\code{index_pi} Initial selecting variable probability. (default: \code{0.5})}
#'      \item{\code{index_nu} Initial vector of inclusion indicators . By default, each nu is randomly drawn by  \eqn{Bernoulli(1/2)}}
#'      \item{\code{index} Initial vector of index. By default, each element of index vector, which is chosen by \eqn{\nu}, is proposed by normal distribution.}
#'      }
#'      \item Link function: Initial spline coefficients (\code{link_beta}). By default,
#'          \eqn{\big(X_{\theta}^\top X_{\theta} + \rho\, \mathrm{I}\big)^{-1} X_{\theta}^\top Y} is computed,
#'     where \eqn{X_{\theta}} is the B-spline basis design matrix.
#'   \item Error variance (\code{sigma2}): Initial scalar error variance (default \code{0.01}).
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
#' fit1 <- bsSpike(y ~ ., data = simdata,
#'                 niter = 5000, nburnin = 1000, nchain = 1)
#'
#' # Split version
#' models <- bsSpike_setup(y ~ ., data = simdata)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                    nchains = 1, setSeed = TRUE, inits = initList,
#'                    summary = TRUE, samplesAsCodaMCMC = TRUE)
#' fit2 <- as_bsim(models, mcmc.out)
#' summary(fit2)
#' }
#'
#' @references
#' Antoniadis, A., Grégoire, G., & McKeague, I. W. (2004).
#' Bayesian estimation in single-index models. \emph{Statistica Sinica}, 1147-1164.
#'
#' Hornik, K., & Grün, B. (2014). movMF: an R package for fitting mixtures of von Mises-Fisher distributions.
#' \emph{Journal of Statistical Software}, 58, 1-31.
#'
#' McGee, G., Wilson, A., Webster, T. F., & Coull, B. A. (2023).
#' Bayesian multiple index models for environmental mixtures.
#' \emph{Biometrics}, 79(1), 462-474.
#'
#' @name bsSpike
#' @export
bsSpike <- function(formula, data,
                    prior = NULL,
                    init = NULL,
                    monitors = NULL, niter = 10000, nburnin=1000,
                    thin = 1, nchain = 1, setSeed = FALSE
){
  return(
    bsSpike.default(formula = formula, data  = data,
                    prior = prior,
                    init = init,
                    sampling = TRUE,
                    monitors = monitors, niter = niter, nburnin = nburnin,
                    thin = thin, nchain = nchain, setSeed = setSeed)
    )

}

#' @rdname bsSpike
#' @export
bsSpike_setup <- function(formula, data,
                          prior = NULL,
                          init = NULL,
                          monitors = NULL, niter = 10000, nburnin=1000,
                          thin = 1, nchain = 1, setSeed = FALSE){
  return(
    bsSpike.default(formula = formula, data  = data,
                    prior = prior,
                    init = init,
                    sampling = FALSE,
                    monitors = monitors, niter = niter, nburnin = nburnin,
                    thin = thin, nchain = nchain, setSeed = setSeed)
  )
}

bsSpike.default <- function(formula, data,
                    prior = NULL,
                    init = NULL,
                    sampling = TRUE,
                    monitors = NULL, niter = 10000, nburnin=1000,
                    thin = 1, nchain = 1, setSeed = FALSE
){
  start1 <- Sys.time()
  index_raw <- 0; nu <- 0

  # check sampling, prior, init parameters for independent execution
  checkOutput <- validate_and_finalize_args(
    sampling, niter, nburnin, thin,  nchain,
    prior, init, "spike", "bspline"
  )
  prior <- checkOutput$priorlist_final
  init <- checkOutput$initlist_final

  # environment
  envobj <- ls(envir=.GlobalEnv)
  on.exit(rm(list=ls(envir=.GlobalEnv)[which(!ls(envir=.GlobalEnv)%in%envobj)],envir=.GlobalEnv))

  .fns <- c(
    # a_common
    "quickSortOrderIndexOnly", "nimOrder", "nimSort",
    "sampleQuantile_nim", "quantile_nimble",

    # aa_bspline_ver3
    "SplineState", "any_duplicated", "mat_wo_col1", "update_spline_df", "update_x_index",
    "update_knot_sequence", "get_basis_simple", "simplify_knots", "get_inside_x",
    "gen_default_internal_knots", "SplineBase1", "SplineBase2", "basis", "bsNimble", "bsBasis",

    # bsFisher
    "postll_bspline_fisher","nimNorm","rW","besselI_nimble","Stheta",
    "estBeta_fisher","gvcCV","transX_fisher",
    "pred_bsplineFisher","indexSampler_bspline_fisher","betaSampler_bspline_fisher",


    "pred_fitted"
  )

  pkg <- "BayesSIM"
  ns <- asNamespace(pkg)
  list2env(mget(.fns, envir = ns, inherits = FALSE), envir = globalenv())


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

  # data dimension
  N <- length(Y)
  p <- ncol(X)


  # Model code
  Rmodel <- nimbleCode({
    # thetastar: spike and slab
    pi ~ dbeta(a0, b0)
    for (j in 1:p) {
      nu[j] ~ dbern(pi)
      index_raw[j] ~ dnorm(0, sd = sigma_theta)
      index_temp[j] <- index_raw[j] * nu[j] + 1e-10
    }
    index[1:p] <- index_temp[1:p]/sqrt(sum(index_temp[1:p]^2))

    # Linear predictor
    for (i in 1:N){
      Xlin[i] <- sum(X[i, 1:p] * index[1:p])
    }

    # Design matrix - b spline basis
    Xmat[1:N, 1:df] <- transX_fisher(Xlin[1:N], df = df, degree = degree, delta = delta)

    # likelihoods
    sigma2 ~ dinvgamma(a_sig, b_sig)
    beta[1:df] ~ dmnorm(mubeta[1:df], cov = covbeta[1:df, 1:df])
    linkFunction[1:N, 1] <- Xmat[1:N, 1:df] %*% matrix(beta[1:df], ncol = 1)
    for (i in 1:N){
      Y[i, 1] ~ dnorm(linkFunction[i, 1], var = sigma2)
    }

  })

  # Prior parameters
  ## check data dimension and save
  ## index: theta
  if (is.null(prior$index$nu$r1)||
      length(prior$index$nu$r1) >= 2 ||
      prior$index$nu$r1 < 0){
    stop("Prior index (r1) has incorrect value.")
  } else{
    a0 <- prior$index$nu$r1
  }

  if (is.null(prior$index$nu$r2)||
      length(prior$index$nu$r2) >= 2 ||
      prior$index$nu$r2 < 0){
    stop("Prior index (r2) has incorrect value.")
  } else{
    b0 <- prior$index$nu$r2
  }

  if (is.null(prior$index$sigma_theta)||
      length(prior$index$sigma_theta) >= 2 ||
      prior$index$sigma_theta < 0){
    stop("Prior index (sigma_theta) has incorrect value.")
  } else{
    sigma_theta <- prior$index$sigma_theta
  }

  ## link- basis
  if (is.null(prior$link$basis$df)||length(prior$link$basis$df) >= 2 || prior$link$basis$df< 0){
    stop("Prior bspline (df) has incorrect value.")
  } else{
    df <- prior$link$basis$df
  }

  if (is.null(prior$link$basis$degree)||length(prior$link$basis$degree) >= 2 || prior$link$basis$degree< 0){
    stop("Error: Prior bspline (degree) has incorrect value.")
  } else{
    degree <- prior$link$basis$degree
  }

  if (is.null(prior$link$basis$delta)||length(prior$link$basis$delta) >= 2 || prior$link$basis$delta< 0){
    stop("Prior bspline (delta) has incorrect value.")
  } else{
    delta <- prior$link$basis$delta
  }

  ## link - beta
  if (is.null(prior$link$beta$mu)){
    mubeta <- rep(0, df)
  }

  if (!is.vector(prior$link$beta$mu) & !is.null(prior$link$beta$mu)){
    stop("Prior beta (mu) should be vector.")
  }
  if (!is.null(prior$link$beta$mu) & length(prior$link$beta$mu) != (df)){
    stop("Incorrect dimension on prior beta (mu).")
  }

  if (is.null(prior$link$beta$cov)){
    covbeta <- diag(df)
  }
  if (!is.null(prior$link$beta$cov) & !is.matrix(prior$link$beta$cov)){
    stop("Prior beta (cov) should be matrix.")
  }

  ## sigma2
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

  # Initialize
  ## sigma2
  init_sigma2 <- init$sigma2
  if (length(init$sigma2) >= 2 || !is.numeric(init$sigma2)){
    stop("Initial value of sigma2 should be scalar.")
  }
  if (init$sigma2 < 0){
    stop("Initial value of sigma2 should be positive.")
  }

  ## beta
  if (!is.null(init$link$beta) & length(init$link$beta) != df){
    stop("Incorrect dimention on initial value of beta.")
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




  inits_list <- lapply(seq_len(nchain), function(j) initFunction_bSpike(X = X, Y = Y,
                                                                        pi = init$index$pi, nu = init$index$nu,
                                                                        index = init$index$index,
                                                                        sigma2 = init$sigma2,
                                                                        beta = init$link$beta,
                                                                        df = df, degree = degree, delta = delta,
                                                                        setSeed = seedNum[j]))
  firstInit <- inits_list[[1]]

  # Build model
  message("Build Model")
  suppressMessages(simpleModel <- nimbleModel(Rmodel,
                             data = list(X = X, Y = Y),
                             constants = list(p = p, N = N,
                                              a0 = a0, b0 = b0, sigma_theta = sigma_theta,
                                              a_sig = a_sig, b_sig = b_sig,
                                              df = df, degree = degree,
                                              delta = delta, mubeta = mubeta, covbeta = covbeta),
                             inits = firstInit))

  # Assign samplers
  message("Assign samplers")
  monitorsList <- c("nu", "index", "sigma2", "linkFunction", "Xlin", "beta", "index_raw", "pi")
  suppressMessages(mcmcConf <- configureMCMC(simpleModel,
                                             monitors = monitorsList,
                                             print = FALSE))

  mcmcConf$removeSamplers(c("beta"))
  mcmcConf$addSampler(target = c("beta"),
                      type   = "betaSampler_bspline_fisher")

  # mcmcConf$setSamplerExecutionOrder(c(1, 2, 3, 5, 6, 7, 8, 4))
  # c(1, 2:(p+1),(p+2), (p+3):(2p+2), (2p+3))
  mcmcConf$setSamplerExecutionOrder(c(1, (p+3):(2*p+2), 2:(p+1), (2*p+3), (p+2)))

  mcmc1 <- buildMCMC(mcmcConf)


  if (!sampling){
    mcmc.out <- NULL
    sampMCMC <- NULL
    samples <- NULL

  } else{
    start2 <- Sys.time()
    # Compile
    message("Compile Model")
    suppressMessages(CsimpleModel <- compileNimble(simpleModel))
    message("Compile MCMC")
    suppressMessages(Cmcmc <- compileNimble(mcmc1, project = simpleModel,
                           resetFunctions = TRUE))
    end2 <- Sys.time()

    # Sampling
    message("Run MCMC")
    mcmc.out <- NULL
    if (setSeed == FALSE){
      seedNum <- setSeed
    }
    mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                        thin = thin,
                        nchains = nchain, setSeed = seedNum,
                        inits = inits_list,
                        summary = FALSE, samplesAsCodaMCMC = TRUE)

    # output
    ## combine all chains
    sampMCMC <- mcmc.out
    samples <- sampleBind(sampMCMC, nchain)

    # if (fitted){ # posterior fitted value output (mean, median, sd)
    #   message("Compute posterior fitted value")
    #   # namesBeta <- paste0("index", 1:p)
    #   namesLink <- paste0("linkFunction[", 1:N, ", 1]")
    #   namesSigma <- "sigma2"
    #   LinkFunction_samples <- samples[, namesLink]
    #   sigma2_samples <- samples[, namesSigma]
    #   n <- nrow(LinkFunction_samples)
    #   p <- ncol(LinkFunction_samples)
    #
    #   message("Compile function..")
    #   suppressMessages(cpred_fitted <- compileNimble(pred_fitted))
    #   message("Computing predicted value..")
    #   fittedValue <- cpred_fitted(LinkFunction_samples,
    #                                sigma2_samples)
    #   fittedResult <- fittedValue
    #
    # } else{
    #   fittedResult <- NULL
    # }

  }
  end1 <- Sys.time()

  ## Input options 정리
  if (!sampling){
    time <- NULL
  } else{
    # total time - compile time
    time_wo_compile <- difftime(end1, start1, units = "secs") -
      difftime(end2, start2, units = "secs")

    time <- list(time_wo_compile = time_wo_compile)

  }

  # Results ------------------------------------------------------------------
  #Model point estimation
  modelESTLIST <- bsimFit_pointest(samples, X, Y)
  inputOptions <- list(origdata = list(x = X, y = Y), formula = formula,
                       prior = list(index = list(nu = list(r1 = a0, r2 = b0),sigma_theta = sigma_theta),
                                    link = list(basis = list(df = df, degree = degree, delta = delta),
                                                beta = list(mu = mubeta, cov = covbeta)),
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
                samples = sampMCMC, # return of runMCMC
                input = inputOptions,
                defModel = simpleModel, defSampler = mcmc1,
                modelName = "bsSpike")

    class(out) = "bsim"


  } else{
    out <- list(input = inputOptions,
                defModel = simpleModel,
                defSampler = mcmc1,
                modelName = "bsSpike")

    class(out) = "bsimSetup"

  }
  return(out)
}
