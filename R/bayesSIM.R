#' Integrated Function for Bayesian Single-Index Regression
#' @importFrom utils modifyList
#' @description
#' Fitting a single–index model \eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n} in single integrated function.
#' @inheritParams bsFisher
#' @param indexprior Index vector prior among \code{"fisher"} (default), \code{"sphere"}, \code{"polar"}, \code{"spike"}.
#' @param link Link function among \code{"bspline"} (default), \code{"gp"}
#' @param prior Optional named list of prior settings. Further descriptions are in every specific model's \strong{Details} section.
#' @param init Optional named list of initial values. If the values are not assigned, they are randomly sampled from prior or designated value.
#' Further descriptions are in every specific model's \strong{Details} section.
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
#' @param x A fitted `BayesSIM` object.
#' @param digits Number of digits to display.
#' @param ... Additional arguments.
#'
#'
#' @inherit bsFisher return
#'
#' @details
#' Integrated function for Bayesian single-index model. Default model is von-Mises Fisher distribution for index vector with B-spline link function.
#'
#' @seealso [bsFisher()], [bsSphere()], [bsPolar()], [bsSpike()],
#' [gpFisher()], [gpSphere()], [gpPolar()], [gpPolarHigh()], [gpSpike()]
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
#' # One tool version - bsFisher
#' fit1 <- BayesSIM(y ~ ., data = simdata,
#'                  niter = 5000, nburnin = 1000,
#'                  nchain = 1)
#'
#' # Split version- bsFisher
#' models <- BayesSIM_setup(y ~ ., data = simdata)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                    nchains = 1, setSeed = TRUE, inits = initList,
#'                    summary = TRUE, samplesAsCodaMCMC = TRUE)
#' fit2 <- as_bsim(models, mcmc.out)
#' summary(fit2)
#' }
#' @name BayesSIM
#' @export
BayesSIM <- function(formula, data,
                     indexprior = "fisher", link = "bspline",
                     prior = NULL, init = NULL,
                     method = "FB",
                     lowerB = NULL, upperB = NULL,
                     monitors = NULL, niter = 10000, nburnin=1000,
                     thin = 1, nchain = 1, setSeed = FALSE){
  # prior_list <- 0; init_list <- 0



  if (link == "bspline"){
    if (indexprior == "fisher"){

      fit <- bsFisher(formula = formula, data = data,
                      prior = prior, init = init,
                      monitors = monitors, niter = niter, nburnin = nburnin,
                      thin = thin, nchain = nchain, setSeed = setSeed)

    } else if (indexprior == "sphere"){
      fit <- bsSphere(formula = formula, data = data,
                      prior = prior, init = init,
                      monitors = monitors, niter = niter, nburnin = nburnin,
                      thin = thin, nchain = nchain, setSeed = setSeed)

    } else if (indexprior == "polar"){
      fit <- bsPolar(formula = formula, data = data,
                     prior = prior, init = init,
                     monitors = monitors, niter = niter, nburnin = nburnin,
                     thin = thin, nchain = nchain, setSeed = setSeed)

    } else if (indexprior == "spike"){
      fit <- bsSpike(formula = formula, data = data,
                     prior = prior, init = init,
                     monitors = monitors, niter = niter, nburnin = nburnin,
                     thin = thin, nchain = nchain, setSeed = setSeed)

    } else{
      stop("Wrong index prior name!")
    }

  } else if (link == "gp"){
    if (indexprior == "fisher"){
      fit <- gpFisher(formula = formula, data = data,
                      prior = prior, init = init,
                      monitors = monitors, niter = niter, nburnin = nburnin,
                      thin = thin, nchain = nchain, setSeed = setSeed)

    } else if (indexprior == "sphere"){
      fit <- gpSphere(formula = formula, data = data,
                      prior = prior, init = init,
                      method = method, lowerB = lowerB, upperB = upperB,
                      monitors = monitors, niter = niter, nburnin = nburnin,
                      thin = thin, nchain = nchain, setSeed = setSeed)

    } else if (indexprior == "polar"){
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
      if (ncol(X) >= 5){
        fit <- gpPolarHigh(formula = formula, data = data,
                           prior = prior, init = init,
                           monitors = monitors, niter = niter, nburnin = nburnin,
                           thin = thin, nchain = nchain, setSeed = setSeed)
      } else{
        fit <- gpPolar(formula = formula, data = data,
                       prior = prior, init = init,
                       monitors = monitors, niter = niter, nburnin = nburnin,
                       thin = thin, nchain = nchain, setSeed = setSeed)
      }
    } else if (indexprior == "spike"){
      fit <- gpSpike(formula = formula, data = data,
                     prior = prior, init = init,
                     monitors = monitors, niter = niter, nburnin = nburnin,
                     thin = thin, nchain = nchain, setSeed = setSeed)

    } else{
      stop("Wrong index prior name!")
    }

  } else{
    stop("Wrong link function name!")
  }

  return(fit)

}
#' @rdname BayesSIM
#' @export
BayesSIM_setup <- function(formula, data,
                     indexprior = "fisher", link = "bspline",
                     prior = NULL, init = NULL,
                     method = "FB",
                     lowerB = NULL, upperB = NULL,
                     monitors = NULL, niter = 10000, nburnin=1000,
                     thin = 1, nchain = 1, setSeed = FALSE){
  # modifyList <- 0
  # prior_list <- 0; init_list <- 0

  if (link == "bspline"){
    if (indexprior == "fisher"){

      fit <- bsFisher_setup(formula = formula, data = data,
                            prior = prior, init = init,
                            monitors, niter, nburnin,
                            thin, nchain, setSeed)

    } else if (indexprior == "sphere"){
      fit <- bsSphere_setup(formula = formula, data = data,
                            prior = prior, init = init,
                            monitors, niter, nburnin,
                            thin, nchain, setSeed)

    } else if (indexprior == "polar"){
      fit <- bsPolar_setup(formula = formula, data = data,
                           prior = prior, init = init,
                           monitors, niter, nburnin,
                           thin, nchain, setSeed)

    } else if (indexprior == "spike"){
      fit <- bsSpike_setup(formula = formula, data = data,
                           prior = prior, init = init,
                           monitors, niter, nburnin,
                           thin, nchain, setSeed)

    } else{
      stop("Wrong index prior name!")
    }

  } else if (link == "gp"){
    if (indexprior == "fisher"){
      fit <- gpFisher_setup(formula = formula, data = data,
                            prior = prior, init = init,
                            monitors, niter, nburnin,
                            thin, nchain, setSeed)

    } else if (indexprior == "sphere"){
      fit <- gpSphere_setup(formula = formula, data = data,
                            prior = prior, init = init,
                            method, lowerB, upperB,
                            monitors, niter, nburnin,
                            thin, nchain, setSeed)

    } else if (indexprior == "polar"){
      if (ncol(x) >= 5){
        fit <- gpPolarHigh_setup(formula = formula, data = data,
                                 prior = prior, init = init,
                                 monitors, niter, nburnin,
                                 thin, nchain, setSeed)
      } else{
        fit <- gpPolar_setup(formula = formula, data = data,
                             prior = prior, init = init,
                             monitors, niter, nburnin,
                             thin, nchain, setSeed)
      }
    } else if (indexprior == "spike"){
      fit <- gpSpike_setup(formula = formula, data = data,
                           prior = prior, init = init,
                           monitors, niter, nburnin,
                           thin, nchain, setSeed)

    } else{
      stop("Wrong index prior name!")
    }

  } else{
    stop("Wrong link function name!")
  }

  return(fit)

}
