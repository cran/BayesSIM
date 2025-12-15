#' Bayesian Single-Index Regression with B-Spline Link and von Mises-Fisher Prior
#' @importFrom mvtnorm rmvnorm
#' @description Fits a single-index model \eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n}
#' where the link \eqn{f(\cdot)} is represented by B-spline and the
#' index vector \eqn{\theta} has von Mises–Fisher prior.
#'
#' @param formula an object of class \link{formula}. Interaction term is not allowed for single-index model.
#' @param data an data frame.
#' @param prior Optional named list of prior settings. Further descriptions are in \strong{Details} section.
#' @param init Optional named list of initial values. If the values are not assigned, they are randomly sampled from prior or designated value.
#' Further descriptions are in \strong{Details} section.
#' @param monitors Optional character vector of additional monitor nodes. To check the names of the nodes, fit the \code{model_setup} function and
#' then inspect the variable names stored in the model object using \link{getVarMonitor}.
#' @param niter Integer. Total MCMC iterations (default \code{10000}).
#' @param nburnin Integer. Burn-in iterations (default \code{1000}).
#' @param thin Integer. Thinning for monitors (default \code{1}).
#' @param nchain Integer. Number of MCMC chains (default \code{1}).
#' @param setSeed Logical or numeric argument.  Further details are provided in \link[nimble]{runMCMC} \code{setSeed} argument.
#'
#' @details
#' \strong{Model} The single–index model uses a \eqn{m}-order polynomial spline with \eqn{k} interior knots as follows:
#' \deqn{f(t) = \sum_{j=1}^{m+k} B_j(t)\,\beta_j} on \eqn{[a, b]} with \eqn{t_i = X_i' \theta, i = 1,\cdots, n}
#' and \eqn{\|\theta\|_2 = 1}. \eqn{\{\beta_j\}_{j=1}^{m+k}} are spline coefficients and \eqn{a_\theta}, \eqn{b_\theta} are boundary knots where \eqn{a_{\theta} = min(t_i, i = 1, \cdots, n) - \delta},
#' and \eqn{b_{\theta} = max(t_i, i = 1,\cdots, n) + \delta}.
#'
#' \strong{Priors}
#' \itemize{
#' \item von Mises–Fisher prior on the index \eqn{\theta} with direction and concentration.
#' \item Conditioned on \eqn{\theta} and \eqn{\sigma^2}, the link coefficients \eqn{\beta = (\beta_1,\ldots,\beta_{m+k})^\top} follow
#'        normal distribution with estimated mean vector \eqn{\hat{\beta}_{\theta} = (X_{\theta}'X_{\theta})^{-1}X_{\theta}'Y} and
#'        covariance \eqn{\sigma^2 (X_{\theta}^\top X_{\theta})^{-1}}, where \eqn{X_{\theta}} is the B-spline basis design matrix.
#' \item Inverse gamma prior on \eqn{\sigma^2} with shape parameter \eqn{a_\sigma} and rate parameter \eqn{b_\sigma}.
#' }
#'
#' \strong{Sampling}
#' Random walk metropolis algorithm is used for index vector \eqn{\theta}. Given \eqn{\theta}, \eqn{\sigma^2} and \eqn{\beta} are sampled from posterior distribution. Further sampling method is described in Antoniadis et al(2004).
#'
#' \strong{Prior hyper-parameters}
#' These are the prior hyper-parameters set in the function. You can define new values for each parameter in \link{prior_param}.
#' \enumerate{
#'   \item Index vector: von Mises--Fisher prior for the projection vector \eqn{\theta} (\code{index}).
#'         \code{index_direction} is a unit direction vector of the von Mises--Fisher distribution.
#'         By default, the estimated vector from projection pursuit regression is assigned.
#'         \code{index_dispersion} is the positive concentration parameter. By default, \code{150} is assigned.
#'
#'   \item Link function: B-spline basis and coefficient of B-spline setup.
#'         \itemize{
#'           \item \code{basis}: For the basis of B-spline, \code{link_basis_df} is the number of basis functions (default \code{21}),
#'                 \code{link_basis_degree} is the spline degree (default \code{2}) and
#'                 \code{link_basis_delta} is a small jitter for boundary knots spacing control (default \code{0.001}).
#'           \item \code{beta}: For the coefficient of B-spline, multivariate normal prior is assigned with mean \code{link_beta_mu}, and covariance \code{link_beta_cov}.
#'                 By default, \eqn{\mathcal{N}_p\!\big(0, \mathrm{I}_p\big)} is assigned.
#'         }
#'
#'   \item Error variance (\code{sigma2}): An Inverse gamma prior is assigned to \eqn{\sigma^2}
#'         where \code{sigma2_shape} is shape parameter and \code{sigma2_rate} is rate parameter of inverse gamma distribution.
#'         (default \code{sigma2_shape = 0.001, sigma2_rate = 100})
#' }
#'
#'
#' \strong{Initial values}
#' These are the initial values set in the function. You can define new values for each initial value in \link{init_param}.
#' \enumerate{
#'  \item Index vector: Initial unit index vector \eqn{\theta}. By default, the vector is randomly sampled from the von Mises--Fisher prior.
#'    \item Link function: Initial spline coefficients (\code{link_beta}). By default,
#'          \eqn{\big(X_{\theta}^\top X_{\theta} + \rho\, \mathrm{I}\big)^{-1} X_{\theta}^\top Y} is computed,
#'     where \eqn{X_{\theta}} is the B-spline basis design matrix.
#'
#'   \item Error variance (\code{sigma2}): Initial scalar error variance (default \code{0.01}).
#'   }
#'
#'
#' @return A \code{list} typically containing:
#' \describe{
#' \item{\code{coefficients}}{Mean estimates of index vector. Return list of \code{model_setup} does not include it.}
#' \item{\code{ses}}{Mean standard error of index vector. Return list of \code{model_setup} does not include it.}
#' \item{\code{residuals}}{Residuals with mean estimates of fitted values. Return list of \code{model_setup} does not include it.}
#' \item{\code{fitted.values}}{Mean estimates of fitted values. Return list of \code{model_setup} does not include it.}
#' \item{\code{linear.predictors}}{Mean estimates of single-index values. Return list of \code{model_setup} does not include it.}
#' \item{\code{gof}}{Goodness-of-fit. Return list of \code{model_setup} function does not include it.}
#' \item{\code{samples}}{Posterior draws of variables for computing fitted values of the model, including \eqn{\theta}, \eqn{\sigma^2}.
#' Return list of \code{model_setup} does not include it.}
#' \item{\code{input}}{List of data used in modeling, formula, input values for prior and initial values, and computation time without compiling.}
#' \item{\code{defModel}}{Nimble model object.}
#' \item{\code{defSampler}}{Nimble sampler object.}
#' \item{\code{modelName}}{Name of the model.}
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
#' simdata <- data.frame(x = X, y = y)
#' colnames(simdata) <- c(paste0("X", 1:4), "y")
#'
#' # One tool version
#' fit1 <- bsFisher(y ~ ., data = simdata,
#'                  niter = 5000, nburnin = 1000, nchain = 1)
#'
#' # Split version
#' models <- bsFisher_setup(y ~ ., data = simdata)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                    nchains = 1, setSeed = TRUE, inits = initList,
#'                    summary = TRUE, samplesAsCodaMCMC = TRUE)
#' fit2 <- as_bsim(models, mcmc.out)
#' summary(fit2)
#'
#'
#' }
#'
#' @references
#' Antoniadis, A., Grégoire, G., & McKeague, I. W. (2004).
#' Bayesian estimation in single-index models. \emph{Statistica Sinica}, 1147-1164.
#'
#' Hornik, K., & Grün, B. (2014). movMF: an R package for fitting mixtures of von Mises-Fisher distributions.
#' \emph{Journal of Statistical Software}, 58, 1-31.
#' @name bsFisher
#' @export
bsFisher <- function(formula, data,
                     prior = NULL,
                     init = NULL,
                     monitors = NULL, niter = 10000, nburnin = 1000,
                     thin = 1, nchain = 1, setSeed = FALSE){

  return(bsFisher.default(formula = formula, data  = data,
                          prior = prior,
                          init = init,
                          sampling = TRUE,
                          monitors = monitors, niter = niter, nburnin = nburnin,
                          thin = thin, nchain = nchain, setSeed = setSeed))
}



#' @name bsFisher
#' @export
bsFisher_setup <- function(formula, data,
                           prior = NULL,
                           init = NULL,
                           monitors = NULL, niter = 10000, nburnin = 1000,
                           thin = 1, nchain = 1, setSeed = FALSE){


  return(bsFisher.default(formula = formula, data  = data,
                          prior = prior,
                          init = init,
                          sampling = FALSE,
                          monitors = monitors, niter = niter, nburnin = nburnin,
                          thin = thin, nchain = nchain, setSeed = setSeed))
}

bsFisher.default <- function(formula, data,
                             prior = NULL,
                             init = NULL,
                             sampling = TRUE, monitors = NULL,
                             niter = 10000, nburnin = 1000,
                             thin = 1, nchain = 1, setSeed = FALSE
){
  start1 <- Sys.time()
  index <- 0;

  # check sampling, prior, init parameters for independent execution
  checkOutput <- validate_and_finalize_args(
    sampling, niter, nburnin, thin, nchain,
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

    # distribution
    "dvMFnim", "rvMFnim",

    # aa_bspline_ver3
    "SplineState", "any_duplicated", "mat_wo_col1", "update_spline_df", "update_x_index",
    "update_knot_sequence", "get_basis_simple", "simplify_knots", "get_inside_x",
    "gen_default_internal_knots", "SplineBase1", "SplineBase2", "basis", "bsNimble", "bsBasis",

    # bsFisher
    "postll_bspline_fisher","nimNorm","rW","besselI_nimble","Stheta",
    "estBeta_fisher","gvcCV", "transX_fisher",
    "pred_bsplineFisher","indexSampler_bspline_fisher","betaSampler_bspline_fisher")



  pkg <- "BayesSIM"
  ns <- asNamespace(pkg)
  list2env(mget(.fns, envir = ns, inherits = FALSE), envir = .GlobalEnv)
  suppressMessages(
    nimble::registerDistributions(list(
      dvMFnim = list(
        BUGSdist = "dvMFnim(theta)",
        types    = c("value = double(1)", "theta = double(1)"),
        discrete = FALSE
      )
    ), verbose = FALSE)
  )

  # --------------------------------------------------------------------------
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
  } else if(!sum(duplicated(c(colnames(data),tot.name))[-c(1:ncol(data))])==
            length(tot.name)){
    stop(paste(paste(tot.name[duplicated(c(colnames(data),
                                           tot.name))[-c(1:ncol(data))]],collapse=","),
               " is/are not in your data"))
  } else{
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
  N <- nrow(X)
  p <- ncol(X)

  # Model code
  Rmodel <- nimble::nimbleCode({
    # index
    index[1:p] ~ dvMFnim(theta = index_prior[1:p])

    # Linear predictor
    for (i in 1:N){
      Xlin[i] <- sum(X[i, 1:p] * index[1:p])
    }

    # Design matrix - B-spline basis
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
    index_direct_input <- as.vector(ppr(X, as.matrix(Y), nterms = 1)$alpha)
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
  monitorsList_essential <- c("index", "sigma2")
  # Add parameters for fitting
  monitorsList <- c(monitorsList_essential, "linkFunction", "Xlin", "beta")
  suppressMessages(
    mcmcConf <- configureMCMC(simpleModel,
                              monitors = monitorsList,
                              print = FALSE)
  )


  # mcmcConf1$printSamplers(executionOrder= TRUE)
  mcmcConf$removeSamplers(c("index"))
  mcmcConf$addSampler(target = c("index"),
                      type   = indexSampler_bspline_fisher)

  mcmcConf$removeSamplers(c("beta"))
  mcmcConf$addSampler(target = c("beta"),
                      type   = betaSampler_bspline_fisher)

  mcmcConf$setSamplerExecutionOrder(c(2, 3, 1))


  mcmc1 <- buildMCMC(mcmcConf)

  samples <- NULL
  if (!sampling){} else{
    # Compile
    start2 <- Sys.time()
    message("Compile Model")
    suppressMessages(CsimpleModel <- compileNimble(simpleModel))
    message("Compile MCMC")
    suppressMessages(Cmcmc <- compileNimble(mcmc1, project = simpleModel, resetFunctions = TRUE))
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
    sampMCMC <- mcmc.out
    samples <- sampleBind(sampMCMC, nchain)
  }
  end1 <- Sys.time()


  # Input options
  ## Calculate execution time
  if (!sampling){
    time <- NULL
  } else{
    # total time - compile time
    time <- difftime(end1, start1, units = "secs") -
      difftime(end2, start2, units = "secs")
  }

  # Results ------------------------------------------------------------------
  #Model point estimation
  modelESTLIST <- bsimFit_pointest(samples, X, Y)

  # input options
  inputOptions <- list(origdata = list(x = X, y = Y), formula = formula,
                       prior = list(index = list(direction = index_direct_input, dispersion = dispersion_prior),
                                    link = list(basis = list(df = df, degree = degree, delta = delta),
                                                beta = list(mu = mubeta, cov = covbeta)),
                                    sigma = list(shape = a_sig, rate = b_sig)),
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
                modelName = "bsFisher")

    class(out) = "bsim"


  } else{
    out <- list(input = inputOptions,
                defModel = simpleModel,
                defSampler = mcmc1,
                modelName = "bsFisher")

    class(out) = "bsimSetup"

  }
  return(out)


}

