#' Bayesian single-index regression with B-spline link and von Mises-Fisher prior
#' @importFrom mvtnorm rmvnorm
#' @description Fits a single-index model \eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n}
#' where the link \eqn{f(\cdot)} is represented by B-spline and the
#' index vector \eqn{\theta} has von Mises–Fisher prior.
#'
#' @param x Numeric data.frame/matrix of predictors. Each row is an observation.
#' @param y Numeric response vector/matrix.
#' @param prior Optional named list of prior settings with sublists:
#' \describe{
#'   \item{\code{index}}{von Mises--Fisher prior for the projection vector \eqn{\theta}.
#'     \code{direction} is a unit direction vector of the von Mises--Fisher distribution.
#'     If \code{direction} is \code{NULL}, the estimated vector from projection pursuit regression is assigned.
#'     \code{dispersion} is the concentration parameter \eqn{c_{\mathrm{prior}} > 0}. (default \code{150})
#'   }
#'   \item{\code{link}}{B-spline basis and coefficient of B-spline setup.
#'          \enumerate{
#'          \item{\code{basis} For the basis of B-spline, \code{df} is the number of basis functions (default \code{21}), \code{degree} is the spline degree (default \code{2}) and
#'          \code{delta} is a small jitter for boundary-knot spacing control (default \code{0.001}).}
#'          \item{\code{beta} For the coefficient of B-spline, multivariate normal prior is assigned with mean \code{mu}, and covariance \code{cov}.
#'           By default, \eqn{\mathcal{N}_p\!\big(0, \mathrm{I}_p\big)}}
#' }}
#'   \item{\code{sigma2}}{Error-variance prior hyperparameters. An Inverse-Gamma prior is assigned to \eqn{\sigma^2}
#'         where \code{shape} is shape parameter and \code{rate} is rate parameter of inverse gamma distribution.
#'         (default \code{shape = 0.001, rate = 100})
#'   }
#' }
#' @param init Optional named list of initial values. If the values are not assigned, they are randomly sampled from prior.
#' \describe{
#'   \item{\code{index}}{Initial unit index vector \eqn{\theta}. By default, the vector is ranomly sampled from the von Mises--Fisher prior.}
#'    \item{\code{link}}{Initial spline coefficients(\code{beta}). By default,
#'          \eqn{\big(X_{\theta}^\top X_{\theta} + \rho\, \mathrm{I}\big)^{-1} X_{\theta}^\top Y} is computed,
#'     where \eqn{X_{\theta}} is the B-spline basis design matrix.
#'   }
#'   \item{\code{sigma2}}{Initial scalar error variance (default \code{0.01}).}
#'
#' }
#' @param sampling Logical. If \code{TRUE} (default), run MCMC; otherwise return prepared nimble model objects without sampling.
#' @param fitted Logical. If \code{TRUE} (default), fitted values drawn from posterior distribution are included in the output and \code{c("Xlin", "linkFunction", "beta")} is monitored for prediction.
#' @param monitors2 Optional character vector of additional monitor nodes. To check the names of the nodes, set \code{fit <- bsFisher(x, y, sampling = FALSE)} and then inspect the variable names stored in the model object using \code{fit$model$getVarNames()}.
#' @param niter Integer. Total MCMC iterations (default \code{10000}).
#' @param nburnin Integer. Burn-in iterations (default \code{1000}).
#' @param thin Integer. Thinning for monitors1 (default \code{1}).
#' @param thin2 Integer. Optional thinning for \code{monitors2} (default \code{1}).
#' @param nchain Integer. Number of MCMC chains (default \code{1}).
#' @param setSeed Logical or numeric argument.  Further details are provided in \link[nimble]{runMCMC}.
#'
#' @details
#' \strong{Model} The single–index model uses a \eqn{m}-order polynomial spline with \eqn{k} interior knots as follows:
#' \eqn{f(t) = \sum_{j=1}^{1+m+k} B_j(t)\,\beta_j} on \eqn{[a, b]} with \eqn{t_i = x_i' \theta, i = 1,\cdots, n}
#' and \eqn{\|\theta\|_2 = 1}. \eqn{\{\beta_j\}_{j=1}^{m+k}} are spline coefficient and \eqn{a_\theta} and \eqn{ b_\theta} are boundary knots where \eqn{a = min(t_i, i = 1, \cdots, n) - \delta},
#' and \eqn{b = max(t_i, i = 1,\cdots, n) + \delta}.
#'
#' \strong{Priors}
#' \itemize{
#' \item von Mises–Fisher prior on the index \eqn{\theta}: direction \code{prior$index$direction}, concentration \code{prior$index$dispersion}.
#' \item Inverse-Gamma prior on \eqn{\sigma^2}: \eqn{\sigma^2 \sim \mathrm{IG}(a_\sigma, b_\sigma)}.
#' \item Conditional on \eqn{\theta} and \eqn{\sigma^2}, the link coefficients follow
#'       \eqn{\beta = (\beta_1,\ldots,\beta_K)^\top \sim \mathcal{N}\!\big(\hat{\beta}_{\theta},\, \sigma^2 (X_{\theta}^\top X_{\theta})^{-1}\big)}.
#' }
#'
#' \strong{Sampling}
#' Random walk metropolis algorithm is used for index vector \eqn{\theta}. Given \eqn{\theta}, \eqn{\sigma^2} and \eqn{\beta} are sampled from posterior distribution. Further sampling method is described in Antoniadis et al.(2004).
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
#' # One tool version
#' fit <- bsFisher(X, y)
#'
#' # Split version
#' models <- bsFisher(X, y, sampling = FALSE)
#' Ccompile <- compileModelAndMCMC(models)
#' mcmc.out <- runMCMC(Ccompile$mcmc, niter = 5000, nburnin = 1000, thin = 1,
#'                    nchains = 1, setSeed = TRUE, inits = models$input$init,
#'                    summary = TRUE, samplesAsCodaMCMC = TRUE)
#' }
#'
#' @references
#' Antoniadis, A., Grégoire, G., & McKeague, I. W. (2004).
#' Bayesian estimation in single-index models. \emph{Statistica Sinica}, 1147-1164.
#'
#' Hornik, K., & Grün, B. (2014). movMF: an R package for fitting mixtures of von Mises-Fisher distributions.
#' \emph{Journal of Statistical Software}, 58, 1-31.
#'
#' @export


bsFisher <- function(x, y,
                     prior = list(
                       index = list(direction = NULL, dispersion = 150),
                       link = list(basis = list(df = 21, degree = 2, delta = 0.001),
                                   beta = list(mu = NULL, cov = NULL)),
                       sigma2  = list(shape = 0.001, rate = 100)),
                     init = list(index = NULL, link = list(beta = NULL), sigma2 = 0.01),
                     sampling = TRUE, fitted = TRUE,
                     monitors2 = NULL, niter = 10000, nburnin=1000,
                     thin = 1, thin2 = NULL, nchain = 1, setSeed = FALSE
){
  start1 <- Sys.time()
  index <- 0

  # check sampling, prior, init parameters for independent execution
  checkOutput <- validate_and_finalize_args(
    sampling, fitted, niter, nburnin, thin, thin2, nchain,
    prior, init, "fisher", "bspline"
  )
  prior <- checkOutput$priorlist_final
  init <- checkOutput$initlist_final

  envobj <- ls(envir=.GlobalEnv)
  on.exit(rm(list=ls(envir=.GlobalEnv)[which(!ls(envir=.GlobalEnv)%in%envobj)],envir=.GlobalEnv))

  # assign functions global environment
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
    "estBeta_fisher","gvcCV",
    "pred_bsplineFisher","indexSampler_bspline_fisher","betaSampler_bspline_fisher")

  pkg <- "BayesSIM"
  ns <- asNamespace(pkg)
  list2env(mget(.fns, envir = ns, inherits = FALSE), envir = .GlobalEnv)


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
  Rmodel <- nimble::nimbleCode({
    # index
    index[1:p] ~ dvMFnim(theta = index_prior[1:p])

    # Linear predictor
    for (i in 1:N){
      Xlin[i] <- sum(X[i, 1:p] * index[1:p])
    }

    # Design matrix - b spline basis
    Xmat[1:N, 1:df] <- transX_fisher(Xlin[1:N], df = df, degree = degree, delta = delta)

    # likelihood
    sigma2 ~ dinvgamma(a_sig, b_sig)
    beta[1:df] ~ dmnorm(mubeta[1:df], cov = covbeta[1:df, 1:df])
    linkFunction[1:N, 1] <- Xmat[1:N, 1:df] %*% matrix(beta[1:df], ncol = 1)
    for (i in 1:N){
      Y[i, 1] ~ dnorm(linkFunction[i, 1], var = sigma2)
    }
  })

  # Prior parameters
  ## check data dimension and save
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
  # index_direction <- index_direct_input

  ## link - basis
  if (is.null(prior$link$basis$df)||length(prior$link$basis$df) >= 2 || prior$link$basis$df< 0){
    stop("Prior basis (df) has incorrect value.")
  } else{
    df <- prior$link$basis$df
  }

  if (is.null(prior$link$basis$degree)||length(prior$link$basis$degree) >= 2 || prior$link$basis$degree< 0){
    stop("Prior basis (degree) has incorrect value.")
  } else{
    degree <- prior$link$basis$degree
  }

  if (is.null(prior$link$basis$delta)||length(prior$link$basis$delta) >= 2 || prior$link$basis$delta< 0){
    stop("Prior basis (delta) has incorrect value.")
  } else{
    delta <- prior$link$basis$delta
  }

  ## beta
  # mu
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
  # if (dim(prior$beta$cov)[1] != (df) || dim(prior$beta$cov)[2] != (df)){
  #   stop("Error: Incorrect dimention on prior beta (cov).")
  # }

  ## sigma
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
  ## index
  if (!is.null(init$index) & length(init$index) != p){
    stop("Incorrect dimension on initial value of index.")
  }

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

  inits_list <- lapply(seq_len(nchain), function(j) initFunction_bF(X = X, Y = Y,
                                                                    index = init$index,
                                                                    sigma2 = init$sigma2,
                                                                    beta = init$link$beta,
                                                                    df = df, degree = degree, delta = delta,
                                                                    index_direction = index_direction,
                                                                    setSeed = seedNum[j]))
  firstInit <- inits_list[[1]]



  # covariance


  # Build model
  message("Build Model")
  suppressMessages(
    simpleModel <- nimbleModel(Rmodel,
                               data = list(X = X, Y = Y),
                               constants = list(p = p, N = N, index_prior = index_direction,
                                                a_sig = a_sig, b_sig = b_sig,
                                                df = df, degree = degree,
                                                delta = delta,mubeta = mubeta, covbeta = covbeta),
                               inits = firstInit)
  )


  # Assign samplers
  message("Assign samplers")
  monitorsList <- c("index", "sigma2")
  if (fitted){
    monitorsList <- c(monitorsList, "linkFunction", "Xlin", "beta")
  }
  if (is.null(monitors2)){
    suppressMessages(
      mcmcConf <- configureMCMC(simpleModel,
                                monitors = monitorsList,
                                print = FALSE)
    )

  } else{
    suppressMessages(
      mcmcConf <- configureMCMC(simpleModel,
                                monitors = monitorsList,
                                monitors2 = monitors2,
                                print = FALSE)
    )

  }


  # mcmcConf1$printSamplers(executionOrder= TRUE)
  mcmcConf$removeSamplers(c("index"))
  mcmcConf$addSampler(target = c("index"),
                      type   = indexSampler_bspline_fisher)

  mcmcConf$removeSamplers(c("beta"))
  mcmcConf$addSampler(target = c("beta"),
                      type   = betaSampler_bspline_fisher)

  mcmcConf$setSamplerExecutionOrder(c(2, 3, 1))


  mcmc1 <- buildMCMC(mcmcConf)
  end1 <- Sys.time()

  if (!sampling){ # compile해서 output
    mcmc.out<- NULL
    fittedResult <- NULL
    sampMCMC <- NULL

  } else{
    # Compile
    message("Compile Model")
    suppressMessages(CsimpleModel <- compileNimble(simpleModel))
    message("Compile MCMC")
    suppressMessages(Cmcmc <- compileNimble(mcmc1, project = simpleModel, resetFunctions = TRUE))

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
                          nchains = nchain, setSeed = seedNum,
                          inits = inits_list,
                          summary = FALSE, samplesAsCodaMCMC = TRUE)
    } else{
      mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin,
                          thin = thin, thin2 = thin2,
                          nchains = nchain, setSeed = seedNum,
                          inits = inits_list,
                          summary = FALSE, samplesAsCodaMCMC = TRUE)
    }
    end2 <- Sys.time()

    # output

    start3 <- Sys.time()
    # if (enableWAIC){
    #   sampMCMC <- mcmc.out$samples
    # } else{
    #   sampMCMC <- mcmc.out
    # }
    sampMCMC <- mcmc.out
    samples <- NULL
    if (nchain > 1){
      for (i in 1:nchain){
        samples <- rbind(samples, sampMCMC[[i]])
      }
    } else if (nchain == 1){
      samples <- sampMCMC
    }

    if (fitted){ # posterior fitted value output (mean, median, sd)
      message("Compute posterior fitted value")
      # namesBeta <- paste0("index", 1:p)
      namesLink <- paste0("linkFunction[", 1:N, ", 1]")
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
                       prior = list(index = list(direction = index_direct_input, dispersion = dispersion_prior),
                                    link = list(basis = list(df = df, degree = degree, delta = delta),
                                                beta = list(mu = mubeta, cov = covbeta)),
                                    sigma = list(shape = a_sig, rate = b_sig)),
                       init = inits_list,
                       time = time)

  out <- list(model = simpleModel,
              sampler = mcmc1,
              sampling = sampMCMC,
              fitted = fittedResult,
              input = inputOptions,
              modelName = "bsFisher")
  class(out) = "bsimSpline"


  return(out)
}
