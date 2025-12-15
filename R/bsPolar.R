#' Bayesian Single-Index Regression with B-Spline Link and One-to-One Polar Transformation
#'
#' @description Fits a single-index model \eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n}
#' where the link \eqn{f(\cdot)} is represented by B-spline link function and the
#' index vector \eqn{\theta} is parameterized on the unit sphere via a one-to-one polar transformation.
#'
#' @inheritParams bsFisher
#' @details
#' \strong{Model} The single–index model uses a \eqn{m}-order polynomial spline with \eqn{k} interior knots as follows:
#' \deqn{f(t) = \sum_{j=1}^{m+k} B_j(t)\,\beta_j} on \eqn{[a, b]} with \eqn{t_i = X_i' \theta, i = 1,\cdots, n}
#' and \eqn{\|\theta\|_2 = 1}. \eqn{\{\beta_j\}_{j=1}^{m+k}} are spline coefficient and \eqn{a_\theta}, \eqn{ b_\theta} are boundary knots where \eqn{a_\theta = min(t_i, i = 1, \cdots, n) - \delta},
#' and \eqn{b_\theta = max(t_i, i = 1,\cdots, n) + \delta}. \eqn{\theta} lies on the unit sphere (\eqn{\|\theta\|_2=1})
#' to ensure identifiability and is parameterized via a one-to-one polar transformation with angle \eqn{\psi_1,...,\psi_{p-1}}, where \eqn{p} is the dimension of predictor.
#' Sampling is  performed on the angular parameters \eqn{\\psi} defining
#' the index vector.
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
#' \strong{Priors}
#' \itemize{
#' \item \eqn{\psi} is \eqn{p-1} dimension of polar angle of index vector and re-scaled Beta distribution on \eqn{[0, \pi]} is assigned.
#' \item Conditioned on \eqn{\theta} and \eqn{\sigma^2}, the link coefficients \eqn{\beta = (\beta_1,\ldots,\beta_{m+k})^\top} follow
#'        normal distribution with estimated mean vector \eqn{\hat{\beta}_{\theta} = (X_{\theta}'X_{\theta})^{-1}X_{\theta}'Y} and
#'        covariance \eqn{\sigma^2 (X_{\theta}^\top X_{\theta})^{-1}}, where \eqn{X_{\theta}} is the B-spline basis design matrix.
#' \item Inverse gamma prior on \eqn{\sigma^2} with shape parameter \eqn{a_\sigma} and rate parameter \eqn{b_\sigma}.
#' }
#'
#' \strong{Sampling}
#' Samplers are automatically assigned by nimble.
#'
#'
#' \strong{Prior hyper-parameters}
#' These are the prior hyper-parameters set in the function. You can define new values for each parameter in \link{prior_param}.
#' \enumerate{
#' \item Index vector:
#'     Only shape parameter \code{index_psi_alpha} of \eqn{p-1} dimension vector is needed since rate parameters is computed to satisfy \eqn{\theta_{j, \text{MAP}}}.
#'     By default, the shape parameter for each element of the polar vector is set to \code{5000}.
#'     \item Link function: B-spline basis and coefficient of B-spline setup.
#'          \itemize{
#'          \item{\code{basis}: For the basis of B-spline, \code{link_basis_df} is the number of basis functions (default \code{21}), \code{link_basis_degree} is the spline degree (default \code{2}) and
#'          \code{link_basis_delta} is a small jitter for boundary-knot spacing control (default \code{0.001}).}
#'          \item{\code{beta}: For the coefficient of B-spline, multivariate normal prior is assigned with mean \code{link_beta_mu}, and covariance \code{link_beta_cov}.
#'           By default, \eqn{\mathcal{N}_p\!\big(0, \mathrm{I}_p\big)}} is assigned.
#'           }
#'
#'   \item Error variance (\code{sigma2}): An Inverse gamma prior is assigned to \eqn{\sigma^2}
#'         where \code{sigma2_shape} is shape parameter and \code{sigma2_rate} is rate parameter of inverse gamma distribution.
#'         (default \code{sigma2_shape = 0.001, sigma2_rate = 100})
#'}
#'
#' \strong{Initial values}
#' These are the initial values set in the function. You can define new values for each initial value in \link{init_param}.
#' \enumerate{
#' \item Index vector: Initial vector of polar angle \code{index_psi} (\eqn{p-1} dimension). Each element of angle is between 0 and \eqn{\pi}.
#' By default, it is randomly draw from uniform distribution.
#'      \item Link function: Initial spline coefficients(\code{link_beta}). By default,
#'          \eqn{\big(X_{\theta}^\top X_{\theta} + \rho\, \mathrm{I}\big)^{-1} X_{\theta}^\top Y} is computed,
#'     where \eqn{X_{\theta}} is the B-spline basis design matrix.
#'   \item Error variance (\code{sigma2}): Initial scalar error variance (default \code{0.01}).
#' }
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
#' fit1 <- bsPolar(y ~ ., data = simdata,
#'                 niter = 5000, nburnin = 1000, nchain = 1)
#'
#' # Split version
#' models <- bsPolar_setup(y ~ ., data = simdata)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                    nchains = 1, setSeed = TRUE, inits = initList,
#'                    samplesAsCodaMCMC = TRUE)
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
#' Dhara, K., Lipsitz, S., Pati, D., & Sinha, D. (2019). A new Bayesian single index model with or without covariates missing at random.
#' \emph{Bayesian analysis}, 15(3), 759.
#'
#' @name bsPolar
#' @export
bsPolar <- function(formula, data,
                    prior = NULL,
                    init = NULL,
                    monitors = NULL, niter = 10000, nburnin=1000,
                    thin = 1, nchain = 1, setSeed = FALSE
){
  return(
    bsPolar.default(formula = formula, data  = data,
                    prior = prior,
                    init = init,
                    sampling = TRUE,
                    monitors = monitors, niter = niter, nburnin = nburnin,
                    thin = thin, nchain = nchain, setSeed = setSeed)
  )

}

#' @rdname bsPolar
#' @export
bsPolar_setup <- function(formula, data,
                          prior = NULL,
                          init = NULL,
                          monitors = NULL, niter = 10000, nburnin=1000,
                          thin = 1, nchain = 1, setSeed = FALSE
){
  return(
    bsPolar.default(formula = formula, data  = data,
                    prior = prior,
                    init = init,
                    sampling = FALSE,
                    monitors = monitors, niter = niter, nburnin = nburnin,
                    thin = thin, nchain = nchain, setSeed = setSeed)
  )
}

bsPolar.default <- function(formula, data,
                    prior = NULL,
                    init = NULL,
                    sampling = TRUE,
                    monitors = NULL, niter = 10000, nburnin=1000,
                    thin = 1, nchain = 1, setSeed = FALSE
){
  start1 <- Sys.time()
  psi <- 0

  # check sampling, prior, init parameters for independent execution
  checkOutput <- validate_and_finalize_args(
    sampling, niter, nburnin, thin,  nchain,
    prior, init,"polar", "bspline"
  )
  prior <- checkOutput$priorlist_final
  init <- checkOutput$initlist_final

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

    "pred_fitted", "Xlinear", "alphaTheta", "transX_fisher"
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
    # c is constant
    for(j in 1:(p-1)) {
      d[j] ~ dunif(1e-10, 1e10)
      psi[j] ~ dbeta(c[j], d[j])
    }

    index[1:p] <- alphaTheta(psi[1:(p-1)]*pi)

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
  ## psi
  if (!is.null(prior$index$psi$alpha) & length(prior$index$psi$alpha) != (p-1))
  {stop("Prior psi has incorrect dimension")}
  if (is.null(prior$index$psi$alpha)){
    psi_c <- rep(5000,(p-1))
  } else{
    psi_c <- prior$index$psi$alpha
  }

  ## link-basis
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
  init_psi <- init$index$psi
  if (!is.null(init_psi) & length(init_psi) != (p-1)){
    stop("Initial psi has incorrect dimension")
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

  inits_list <- lapply(seq_len(nchain), function(j) initFunction_bP(X = X, Y = Y,
                                                                    psi = init_psi,
                                                                    sigma2 = init_sigma2,
                                                                    beta = init$link$beta,
                                                                    df = df, degree = degree, delta = delta,
                                                                    setSeed = seedNum[j]))

  firstInit <- inits_list[[1]]

  # Build model
  message("Build Model")
  suppressMessages(simpleModel <- nimbleModel(Rmodel,
                             data = list(X = X, Y = Y),
                             constants = list(p = p, N = N, c = psi_c,
                                              a_sig = a_sig, b_sig = b_sig,
                                              df = df, degree = degree, pi = pi,
                                              delta = delta,mubeta = mubeta, covbeta = covbeta),
                             inits = firstInit))

  # Assign samplers
  message("Assign samplers")
  monitorsList <- c("index", "sigma2", "linkFunction", "Xlin", "beta", "d", "psi")
  suppressMessages(mcmcConf <- configureMCMC(simpleModel,
                                             monitors = monitorsList,
                                             print = FALSE))


  mcmcConf$removeSamplers(c("beta"))
  mcmcConf$addSampler(target = c("beta"),
                      type   = "betaSampler_bspline_fisher")

  # mcmcConf$setSamplerExecutionOrder(c(1, 2, 3, 5, 6, 7, 8, 4))

  mcmc1 <- buildMCMC(mcmcConf)


  mcmc.out <- NULL
  if (!sampling){
    sampMCMC <- NULL
    samples <- NULL
  } else{
    # Compile
    start2 <- Sys.time()
    message("Compile Model")
    suppressMessages(CsimpleModel <- compileNimble(simpleModel))
    message("Compile MCMC")
    suppressMessages(Cmcmc <- compileNimble(mcmc1, project = simpleModel, resetFunctions = TRUE))
    end2 <- Sys.time()

    # Sampling
    message("Run MCMC")

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
  }
  end1 <- Sys.time()

  ## Input options
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
                       prior = list(index = list(psi = list(alpha = psi_c)),
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
                modelName = "bsPolar")

    class(out) = "bsim"


  } else{
    out <- list(input = inputOptions,
                defModel = simpleModel,
                defSampler = mcmc1,
                modelName = "bsPolar")

    class(out) = "bsimSetup"

  }
  return(out)

}
