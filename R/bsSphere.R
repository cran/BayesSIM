#' Bayesian Single-Index Regression with B-Spline Link and Half-Unit Sphere Prior
#'
#' @description
#' Fits a single-index model \eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n}
#' where the link \eqn{f(\cdot)} is represented by B-spline link and the
#' index vector \eqn{\theta} is on half-unit sphere.
#'
#' @inheritParams bsFisher
#'
#' @details
#' \strong{Model} The single–index model uses a \eqn{m}-order polynomial spline with \eqn{k} interior knots and intercept as follows:
#' \deqn{f(t) = \sum_{j=1}^{1+m+k} B_j(t)\,\beta_j} on \eqn{[a, b]} with \eqn{t_i = X_i' \theta, i = 1,\cdots, n}
#' and \eqn{\|\theta\|_2 = 1}. \eqn{\{\beta_j\}_{j=1}^{m+k+1}} are spline coefficients and \eqn{a_\theta}, \eqn{ b_\theta} are boundary knots where \eqn{a_{\theta} = min(t_i, i = 1, \cdots, n)},
#' and \eqn{b_{\theta} = max(t_i, i = 1,\cdots, n)}. Variable selection is encoded by a binary vector \eqn{\nu}, equivalently
#' setting components of \eqn{\theta} to zero when \eqn{\nu_i = 0}.
#'
#' \strong{Priors}
#' \itemize{
#' \item The variable selection indicator \eqn{\nu} has Beta–Bernoulli hierarchy prior. Beta hyper-prior on the Bernoulli–inclusion probability \eqn{w},
#'     inducing \eqn{p(\nu) \propto \mathrm{Beta}(r_1 + n_\nu, r_2 + p - n_\nu)} where  \eqn{n_\nu = \Sigma_{i=1}^{p}I(\nu_i = 1)}.
#'     \eqn{r_1, r_2} are shape and rate parameter of beta distribution.
#'   \item Free‑knot prior: the number of knots \eqn{k} with mean \eqn{\lambda_k}. The knot locations \eqn{\xi_i, i = 1,...,k} have a Dirichlet prior on the scaled interval \eqn{[0, 1]}.
#'   \item Index vector prior is uniform on the half‑sphere of dimension \eqn{n_\nu} with first nonzero positive.
#'   \item Conjugate normal–inverse-gamma on \eqn{(\beta, \sigma^2)} enables analytic integration for RJMCMC with covariance \eqn{\tau \Sigma_0}.
#' }
#'
#' \strong{Sampling} Posterior exploration follows a hybrid RJMCMC with six move types:
#' add/remove predictor \eqn{\nu}, update \eqn{\theta}, add/delete/relocate a knot. The \eqn{\theta} update is a random‑walk
#' Metropolis via local rotations in a two‑coordinate subspace. Knot changes use simple proposals with tractable acceptance ratios.
#' Further sampling method is described in Wang (2009).
#'
#' \strong{Prior hyper-parameters}
#' These are the prior hyper-parameters set in the function. You can define new values for each parameter in \link{prior_param}.
#'  \enumerate{
#'  \item Index vector: \code{index_nu_r1, index_nu_r2} gives the shape and rate parameter of beta-binomial prior, respectively. (default \code{index_nu_r1 = 1, index_nu_r2 = 1}).
#' \item Link function: B-spline knots, basis and coefficient setup.
#' \itemize{
#'   \item{knots: Free-knot prior for the spline. \code{link_knots_lambda_k} is the Poisson mean for the number of
#'     interior knot \eqn{k} (default \code{5}). \code{link_knots_maxknots} is the maximum number of interior knots.
#'     If \code{link_knots_maxknots} is \code{NULL}, the number of interior knots is randomly drawn from a Poisson distribution.}
#'   \item{basis: For the basis of B-spline,  \code{link_basis_degree} is the spline
#'     degree (default \code{2}).}
#'   \item{beta: For the coefficient of B-spline,
#'    By default, \code{link_beta_mu} is a zero vector, \code{link_beta_tau} is set to the sample size,
#'     and \code{link_beta_Sigma0} is the identity matrix of dimension \eqn{1 + k + m},
#'     where \eqn{k} is the number of interior knots and \eqn{m} is the spline order.}}
#' \item Error variance (\code{sigma2}): Inverse gamma prior is assigned to \eqn{\sigma^2}
#'  where \code{sigma2_shape} is shape parameter and \code{sigma2_rate} is rate parameter of inverse gamma distribution.
#'     Small values for shape and rate parameters yield a weakly-informative prior on \eqn{\sigma^{2}}. (default \code{sigma2_shape = 0.001, sigma2_rate = 0.001})
#'  }
#'
#' \strong{Initial values}
#' These are the initial values set in the function. You can define new values for each initial value in \link{init_param}.
#'  \enumerate{
#'  \item Index vector: \code{index_nu} is binary vector indicating active predictors for the index.
#' \code{index} is initial unit-norm index vector \eqn{\theta} (automatically normalized if necessary, with the first nonzero element made positive for identifiability).
#'     The vector length must match the number of predictor.
#'     Ideally, positions where \code{index_nu} has a value of 1 should correspond to nonzero elements in \eqn{\theta}; elements corresponding to \code{index_nu} = 0 will be set to zero.
#'  \item Link function: \code{link_k} is initial number of interior knots. \code{link_knots} is initial vector of interior knot positions in \eqn{[0, 1]}, automatically scaled to the true boundary.
#'     Length of this vector should be equal to the initial value of \code{k}.
#'     \code{link_beta} is initial vector of spline coefficients. Length should be equal to the initial number of basis functions with intercept (\eqn{1 + k + m}).
#'   \item Error variance (\code{sigma2}): Initial scalar error variance. (default \code{0.01})
#'  }
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
#' fit1 <- bsSphere(y ~ ., data = simdata,
#'                  niter = 5000, nburnin = 1000, nchain = 1)
#'
#' # Split version
#' models <- bsSphere_setup(y ~ ., data = simdata)
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
#' Wang, H.-B. (2009). Bayesian estimation and variable selection for single index models.
#' \emph{Computational Statistics & Data Analysis}, 53, 2617–2627.
#'
#' Hornik, K., & Grün, B. (2014). movMF: an R package for fitting mixtures of von Mises-Fisher distributions.
#' \emph{Journal of Statistical Software}, 58, 1-31.
#'
#' @name bsSphere
#' @export
bsSphere <- function(formula, data,
                     prior = NULL,
                     init = NULL,
                     monitors = NULL, niter = 10000, nburnin=1000,
                     thin = 1, nchain = 1, setSeed = FALSE){

  return(bsSphere.default(formula = formula, data  = data,
                          prior = prior,
                          init = init,
                          sampling = TRUE,
                          monitors = monitors, niter = niter, nburnin = nburnin,
                          thin = thin, nchain = nchain, setSeed = setSeed))
}

#' @rdname bsSphere
#' @export
bsSphere_setup <- function(formula, data,
                           prior = NULL,
                           init = NULL,
                           monitors = NULL, niter = 10000, nburnin=1000,
                           thin = 1, nchain = 1, setSeed = FALSE){
  return(bsSphere.default(formula = formula, data  = data,
                          prior = prior,
                          init = init,
                          sampling = FALSE,
                          monitors = monitors, niter = niter, nburnin = nburnin,
                          thin = thin, nchain = nchain, setSeed = setSeed))
}

bsSphere.default <- function(formula, data,
                     prior = NULL,
                     init = NULL,
                     sampling = TRUE,
                     monitors = NULL, niter = 10000, nburnin=1000,
                     thin = 1, nchain = 1, setSeed = FALSE
){

  start1 <- Sys.time()
  index <- 0; k <- 0; knots <- 0; sigma2 <- 0

  # check sampling, prior, init parameters for independent execution --------
  checkOutput <- validate_and_finalize_args(
    sampling, niter, nburnin, thin, nchain,
    prior, init, "sphere", "bspline"
  )
  prior <- checkOutput$priorlist_final
  init <- checkOutput$initlist_final

  # global environment

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

    # bsSphere
    "transX_sp","estBetaInit_sp","computeSig","ComputeS","logdet_nim",
    "postll_bspline_sphere","postll_knots","betaFunction","computeA1_nu1",
    "computeA1_nu0","computeA2_add","computeA2_delete","rnormTrun","newKnots",
    "prop_add","prop_delete","rtruncnorm","pred_bsplineSphere",
    "nuSampler_bspline_sphere","indexSampler_bspline_sphere","knotsSampler_bspline_sphere",
    "betaSampler_bspline_sphere","sigma2Sampler_bspline_sphere", "transX_fisher",

    # distribution
    "dunitSphere", "runitSphere","dKnotsSimple", "rKnotsSimple",

    # utils
    "pred_fitted"
  )

  pkg <- "BayesSIM"
  ns <- asNamespace(pkg)
  list2env(mget(.fns, envir = ns, inherits = FALSE), envir = globalenv())

  suppressMessages(
  nimble::registerDistributions(list(
    dKnotsSimple = list(
      BUGSdist = "dKnotsSimple(a, b, k, alpha)",
      types = c("value = double(1)", "a = double(0)", "b = double(0)", "k = double(0)",
                "alpha = double(1)"),
      discrete = FALSE
    )
  ), verbose = FALSE))

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

  # data dimension
  N <- length(Y)
  p <- ncol(X)


  # Model code
  Rmodelcode <- nimbleCode({
    for (i in 1:p) {
      nu[i] ~ dbern(r1/(r1+r2))
    }

    # index
    index[1:p] ~ dunitSphere(p) ## half unit sphere

    # Linear predictor
    for (i in 1:N){
      Xlin[i] <- sum(X[i, 1:p] * index[1:p])
    }

    # boundary knots
    a_alpha <- min(Xlin[1:N])
    b_alpha <- max(Xlin[1:N])

    # knots
    k ~ dpois(lambda_k)
    knots[1:maxknots] ~ dKnotsSimple(a = a_alpha, b = b_alpha, k = k, alpha = alpha[1:maxknots])

    # Design matrix - b spline basis
    numBasis <- degree + k + 1 # real dimension
    Xmat[1:N, 1:maxBasis] <- transX_sp(Xlin[1:N], degree, knots[1:maxknots], k, maxBasis, a_alpha, b_alpha)

    # likelihood
    sigma2 ~ dinvgamma(a_sig, b_sig)
    covBeta[1:maxBasis, 1:maxBasis] <- (sigma2 * tau) * Sigma0[1:maxBasis, 1:maxBasis]
    beta[1:maxBasis, 1] ~ dmnorm(mubeta[1:maxBasis],
                                 cov = covBeta[1:maxBasis, 1:maxBasis])

    linkFunction[1:N, 1] <- Xmat[1:N, 1:maxBasis] %*% beta[1:maxBasis, 1]
    covY[1:N, 1:N] <- sigma2 * Sigma[1:N, 1:N]
    Y[1:N, 1] ~ dmnorm(linkFunction[1:N, 1], cov = covY[1:N, 1:N])

  })


  # Prior parameters
  ## check data dimension and save
  ## index - nu
  if (is.null(prior$index$nu$r1) || is.null(prior$index$nu$r2)){
    stop("nu prior is not available.")
  }
  if (prior$index$nu$r1 <= 0 || prior$index$nu$r2 <= 0){
    stop("nu prior should be positive.")
  }

  r1 <- prior$index$nu$r1
  r2 <- prior$index$nu$r2

  ## link - knots
  if (is.null(prior$link$knots$lambda_k) || is.null(prior$link$knots$lambda_k)){
    stop("knots prior is not available.")
  }
  if (prior$link$knots$lambda_k <= 0){
    stop("knots prior should be positive.")
  }

  lambda_k <- prior$link$knots$lambda_k

  maxknots <- ifelse(is.null(prior$link$knots$maxknots),
                     qpois(1-0.01, lambda = lambda_k), prior$link$knots$maxknots)

  if (length(maxknots) >= 2 || maxknots < 0){
    stop("Prior link (maxknots) has incorrect value.")
  }

  ## link - basis
 if (is.null(prior$link$basis$degree)||length(prior$link$basis$degree) >= 2 || prior$link$basis$degree< 0){
    stop("Prior basis (degree) has incorrect value.")
  } else{
    degree <- prior$link$basis$degree
  }

  alpha <- rep(1/maxknots, maxknots)
  maxBasis <- degree + 1 + maxknots


  ## beta
  ## tau
  tau <- ifelse(is.null(prior$link$beta$tau), N, prior$link$beta$tau)
  if (length(tau) >= 2 || tau < 0){
    stop("Prior beta (tau) has incorrect value.")
  }

  ## mu
  if (!is.vector(prior$link$beta$mu) & !is.null(prior$link$beta$mu)){
    stop("Prior beta (mu) should be vector.")
  }
  if (!is.null(prior$link$beta$mu) & length(prior$link$beta$mu) < maxBasis ){
    stop("Incorrect dimention on prior beta (mu).")
  }

  if (is.null(prior$link$beta$mu)){
    mubeta <- rep(0, maxBasis)
  } else{
    mubeta <- c(prior$link$beta$mu, rep(0, maxBasis - length(prior$link$beta$mu)))
  }



  ## Sigma0

  if (is.null(prior$link$beta$Sigma0)){
    Sigma0 <- diag(1, maxBasis)
  } else{
    Sigma0 <- prior$link$beta$Sigma0
  }
  if (!is.null(prior$link$beta$Sigma0) & !is.matrix(prior$link$beta$Sigma0)){
    stop("Prior beta (Sigma0) should be matrix.")
  }

  Sigma <- diag(1, N)


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
  ## nu
  if (!is.null(init$index$nu) & (length(init$index$nu) != p)){
    stop("Incorrect dimension on initial value of nu.")
  }

  if (!is.null(init$index$nu) & (!all(init$index$nu %in% c(0, 1)))){
    stop("nu should consist of 0 or 1.")
  }

  ## index
  if (!is.null(init$index$index) & (length(init$index$index) != p)){
    stop("Incorrect dimension on initial value of index")
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
                       function(j) initFunction_bS(nu = init$index$nu, index = init$index$index,
                                                   k = init$link$k, beta = init$link$beta,
                                                   sigma2 = init$sigma2, knots = init$link$knots,
                                                   p = p, X = X, Y = Y, lambda_k = lambda_k,
                                                   maxknots = maxknots, degree = degree,
                                                   maxBasis = maxBasis, tau = tau,
                                                   Sigma0 = Sigma0,setSeed = seedNum[j]))
  firstInit <- inits_list[[1]]


  # Build model
  message("Build Model")
  suppressMessages(simpleModel <- nimbleModel(Rmodelcode,
                                              data = list(X = X, Y = Y),
                                              constants = list(p = p, N = N, r1 = r1, r2 = r2,
                                                               lambda_k = lambda_k,maxknots = maxknots,
                                                               a_sig = a_sig, b_sig = b_sig, alpha= alpha,
                                                               tau = tau, degree = degree, maxBasis = maxBasis,
                                                               mubeta = mubeta, Sigma0 = Sigma0,Sigma = Sigma),
                                              inits = firstInit))

  # Assign samplers
  message("Assign samplers")
  # List that is needed for fitting the model
  monitorsList <- c("index", "nu", "sigma2", "linkFunction", "beta",
                    "k", "knots", "numBasis", "a_alpha", "b_alpha", "Xlin")

  suppressMessages(mcmcConf <- configureMCMC(simpleModel,
                                             monitors = monitorsList,
                                             print = FALSE))


  # mcmcConf1$printSamplers(executionOrder= TRUE)
  mcmcConf$removeSamplers(c("nu"))
  mcmcConf$addSampler(target = c("nu"),
                       type   = nuSampler_bspline_sphere)

  mcmcConf$removeSamplers(c("index"))
  mcmcConf$addSampler(target = c("index"),
                       type   = indexSampler_bspline_sphere)

  mcmcConf$removeSamplers(c("k", "knots"))
  mcmcConf$addSampler(target = c("k","knots"),
                       type   = knotsSampler_bspline_sphere)

  mcmcConf$removeSamplers(c("beta"))
  mcmcConf$addSampler(target = c("beta"),
                       type   = betaSampler_bspline_sphere)

  mcmcConf$removeSamplers(c("sigma2"))
  mcmcConf$addSampler(target = c("sigma2"),
                       type   = sigma2Sampler_bspline_sphere)


  mcmc1 <- buildMCMC(mcmcConf)



  if (!sampling){ # not mcmc sampling
    samples <- NULL

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

  }

  end1 <- Sys.time()

  # inputOptions <- NULL
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

  # input options
  inputOptions <- list(origdata = list(x = X, y = Y), formula = formula,
                       prior = list(index = list(nu = list(r1 = r1, r2 = r2)),
                                    link = list(knots = list(lambda_k = lambda_k, maxknots = maxknots),
                                                basis = list(degree = degree),
                                                beta = list(mu = mubeta, tau = tau, sigma0 = Sigma0)),
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
                modelName = "bsSphere")

    class(out) = "bsim"


  } else{
    out <- list(input = inputOptions,
                defModel = simpleModel,
                defSampler = mcmc1,
                modelName = "bsSphere")

    class(out) = "bsimSetup"

  }
  return(out)

}
