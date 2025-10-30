#' Integrated function for Bayesian single-index regression
#' @importFrom utils modifyList
#' @description
#' Fits a single–index model \eqn{Y_i \sim \mathcal{N}(f(X_i'\theta), \sigma^2), i = 1,\cdots,n} in one function.
#'
#' @param x Numeric data.frame/matrix of predictors. Each row is an observation.
#' @param y Numeric response numeric vector/matrix. Other types  are not available.
#' @param indexprior Index vector prior among \code{"fisher"} (default), \code{"sphere"}, \code{"polar"}, \code{"spike"}.
#' @param link Link function among \code{"bspline"} (default), \code{"gp"}
#' @param prior List of prior hyperparameters of index, link function, and sigma2. For further details,
#' refer to \code{help()} of designated function.
#' @param init List of initial values of index, link function, and sigma2. For further details,
#' refer to \code{help()} of designated function.
#' @param sampling Logical. If \code{TRUE} (default), run MCMC; otherwise return prepared nimble model objects without sampling.
#' @param fitted Logical. If \code{TRUE} (default), fitted values drawn from posterior distribution are included in the output and \code{c("Xlin", "linkFunction", "beta")} is monitored for prediction.
#' @param method Character, gpSphere model has 3 different types of sampling method, fully Bayesian method (\code{"FB"}), empirical Bayes approach (\code{"EB"}), and empirical Gibbs sampler (\code{"EG"}).
#' Assign one sampler method. Empirical sampling approach is recommended in high-dimensional data. By default, fully Bayesian approach is assigned.
#' @param lowerB This parameter is only for gpSphere model. Numeric vector of element-wise lower bounds for the \code{"L-BFGS-B"} method.
#' When the empirical Bayes or Gibbs sampler method is used, the marginal likelihood is optimized via \code{optim(method = "L-BFGS-B")}.
#' The vector must be ordered as \code{c(index vector, lengthscale, amp, sigma2)}; note that \code{sigma2} is included only for the empirical Bayes method (omit it for Gibbs).
#' By default, the lower bounds are \code{-1} for the index vector and \code{-1e2} for logarithm of \code{lengthscale}, \code{amp}, and (if present) \code{sigma2}.
#' @param upperB This parameter is only for gpSphere model. Numeric vector of element-wise upper bounds for the \code{"L-BFGS-B"} method.
#' When the empirical Bayes or Gibbs sampler method is used, the marginal likelihood is optimized via \code{optim(method = "L-BFGS-B")}.
#' The vector must be ordered as \code{c(index vector, lengthscale, amp, sigma2)}; note that \code{sigma2} is included only for the empirical Bayes method (omit it for Gibbs).
#' By default, the upper bounds are \code{1} for the index vector and \code{1e2} for logarithm of \code{lengthscale}, \code{amp}, and (if present) \code{sigma2}.
#' @param monitors2 Optional character vector of additional monitor nodes. Available: \code{c("Xlin", "k", "knots", "beta")}.
#' @param niter Integer. Total MCMC iterations (default \code{10000}).
#' @param nburnin Integer. Burn-in iterations (default \code{1000}).
#' @param thin Integer. Thinning for monitors1 (default \code{1}).
#' @param thin2 Integer. Optional thinning for \code{monitors2} (default \code{1}).
#' @param nchain Integer. Number of MCMC chains (default \code{1}). If >1, different initial values are assigned for each chain.
#' @param setSeed Logical or numeric argument.  Further details are provided in \link[nimble]{runMCMC}.
#'
#' @return A \code{list} typically containing:
#' \describe{
#'   \item{\code{model}}{Nimble model}
#'   \item{\code{sampler}}{Nimble sampler}
#'   \item{\code{sampling}}{Posterior draws of samples with \code{coda mcmc} object. \eqn{\nu}(spike-and slab prior), \eqn{\theta}, \eqn{\sigma^2}, \code{monitors2} variables are included.}
#'   \item{\code{fitted}}{If \code{fitted = TRUE}, in-sample fitted values is given.}
#'   \item{\code{input}}{List of data,input values for prior and initial values, and computation time without compiling.}
#' }
#'
#' @details
#' Integrated function for Bayesian single-index model. Default model is von-Mises Fisher distribution for index vector with B-spline link function.
#'
#' @seealso [bsFisher()], [bsSphere()], [bsPolar()], [bsSpike()],
#' [gpFisher()], [gpSphere()], [gpPolar()], [gpSpike()]
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 100; d <- 4
#' theta <- c(2, 1, 1, 1); theta <- theta / sqrt(sum(theta^2))
#' f <- function(u) u^2 * exp(u)
#' sigma <- 0.5
#' X <- matrix(runif(n * d, -1, 1), nrow = n)
#' index_vals <- as.vector(X %*% theta)
#' y <- f(index_vals) + rnorm(n, 0, sigma)
#'
#' # One-call version
#' fit <- BayesSIM(X, y, indexprior = "sphere", link = "bspline")
#'
#' # Split version
#' models <- BayesSIM(X, y, indexprior = "sphere", link = "bspline", sampling = FALSE)
#' Ccompile <- compileModelAndMCMC(models)
#' mcmc.out <- runMCMC(Ccompile$mcmc, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, init = models$input$init,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#' }
#' @export
BayesSIM <- function(x, y, indexprior = "fisher", link = "bspline",
                     prior = NULL, init = NULL,
                     sampling = TRUE, fitted = TRUE, method = "FB",
                     lowerB = NULL, upperB = NULL,
                     monitors2 = NULL, niter = 10000, nburnin=1000,
                     thin = 1, thin2 = NULL, nchain = 1, setSeed = FALSE){
  # modifyList <- 0

  # check sampling, prior, init parameters for independent execution
  # checkOutput <- validate_and_finalize_args(
  #   sampling, fitted, niter, nburnin, thin, thin2, nchain,
  #   prior, init, indexprior, link
  # )
  # priorlist_final <- checkOutput$priorlist_final
  # initlist_final <- checkOutput$initlist_final



  if (link == "bspline"){
    if (indexprior == "fisher"){

      fit <- bsFisher(x, y,
                      prior = prior, init = init,
                      sampling, fitted,monitors2, niter, nburnin,
                      thin, thin2, nchain, setSeed)

    } else if (indexprior == "sphere"){
      fit <- bsSphere(x, y,
                      prior = prior, init = init,
                      sampling, fitted,
                      monitors2, niter, nburnin,
                      thin, thin2, nchain, setSeed)

    } else if (indexprior == "polar"){
      fit <- bsPolar(x, y,
                     prior = prior, init = init,
                     sampling, fitted,
                     monitors2, niter, nburnin,
                     thin, thin2, nchain, setSeed)

    } else if (indexprior == "spike"){
      fit <- bsSpike(x, y,
                     prior = prior, init = init,
                     sampling, fitted,
                     monitors2, niter, nburnin,
                     thin, thin2, nchain, setSeed)

    } else{
      stop("Wrong index prior name!")
    }

  } else if (link == "gp"){
    if (indexprior == "fisher"){
      fit <- gpFisher(x, y,
                      prior = prior, init = init,
                      sampling, fitted,
                      monitors2, niter, nburnin,
                      thin, thin2, nchain, setSeed)

    } else if (indexprior == "sphere"){
      fit <- gpSphere(x, y,
                      prior = prior, init = init,
                      sampling, fitted, method, lowerB, upperB,
                      monitors2, niter, nburnin,
                      thin, thin2, nchain, setSeed)

    } else if (indexprior == "polar"){
      if (ncol(x) >= 5){
        fit <- gpPolarHigh(x, y,
                           prior = prior, init = init,
                           sampling, fitted,
                           monitors2, niter, nburnin,
                           thin, thin2, nchain, setSeed)
      } else{
        fit <- gpPolar(x, y,
                       prior = prior, init = init,
                       sampling, fitted,
                       monitors2, niter, nburnin,
                       thin, thin2, nchain, setSeed)
      }
    } else if (indexprior == "spike"){
      fit <- gpSpike(x, y,
                     prior = prior, init = init,
                     sampling, fitted,
                     monitors2, niter, nburnin,
                     thin, thin2, nchain, setSeed)

    } else{
      stop("Wrong index prior name!")
    }

  } else{
    stop("Wrong link function name!")
  }

  return(fit)

}
