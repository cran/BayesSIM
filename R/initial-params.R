#' Build an Initial Value List for BayesSIM Models
#'
#' @description
#' `init_param` is a convenience helper that constructs a nested initial value list
#' for a given combination of index vector and link function.
#' It starts from the model-specific default prior, and then overwrites only those components for
#' which the user supplies non-`NULL` arguments.
#'
#' This allows users to modify selected hyper-parameters without having to know
#' or manually reconstruct the full nested prior list structure.
#'
#' @param indexprior Character scalar indicating the prior for the index.
#'   Typically one of `"fisher"`, `"sphere"`, `"polar"`, or `"spike"`.
#'   The valid options mirror those used in the corresponding model functions.
#' @param link Character scalar indicating the link function family.
#'   Typically `"bspline"` for B-spline link functions or `"gp"` for Gaussian
#'   process link functions. The valid options mirror those used in the
#'   corresponding model functions.
#'
#' @param index,index_nu,index_psi,index_pi
#'   Optional initial values for index and related parameter values.
#'
#' @param link_beta,link_k,link_knots,link_lengthscale,link_amp,link_kappa,link_inv_lambda
#'   Optional initial values for components under link functions.
#'
#' @param sigma2 Optional numeric scalar giving the initial value of \eqn{\sigma^2}.
#'
#' @details
#' \code{init_param(indexprior, link)} can be used to obtain the random initial values
#' list for the requested combination of index prior and link function.
#' For any argument that is not `NULL`, the corresponding field in the nested prior list is overwritten.
#'
#' The detailed meaning and recommended choices for each initial values depend
#' on the specific model, index vector and link function.
#' For those details, please refer to the documentation of the corresponding
#' model-fitting functions.
#'
#' @seealso [bsFisher()], [bsSphere()], [bsPolar()], [bsSpike()],
#' [gpFisher()], [gpSphere()], [gpPolar()], [gpPolarHigh()], [gpSpike()]
#'
#' @return A nested list with components \code{index}, \code{link}, and
#'   \code{sigma2}.
#' @examples
#' ## Default initial values for Fisher index + B-spline link:
#' i0 <- init_param("fisher", "bspline")
#'
#' ## Modify only a few initial values:
#' i1 <- init_param(
#'   indexprior = "fisher",
#'   link       = "bspline",
#'   index      = c(1, 0, 0),      # initial direction of the index
#'   link_beta  = rep(0, 21),      # initial values for spline coefficients
#'   sigma2     = 0.1              # initial value for sigma^2
#' )
#'
#' ## Example with GP link:
#' i2 <- init_param(
#'   indexprior        = "sphere",
#'   link              = "gp",
#'   link_lengthscale  = 0.2,      # initial GP length-scale
#'   link_amp          = 1.5,      # initial GP amplitude
#'   sigma2            = 1         # initial variance
#' )
#'
#' @export
init_param <- function(
    indexprior, link, index = NULL, index_nu = NULL, index_psi = NULL,
    index_pi  = NULL, link_beta = NULL,link_k  = NULL,
    link_knots  = NULL, link_lengthscale = NULL, link_amp = NULL,
    link_kappa  = NULL, link_inv_lambda  = NULL, sigma2 = NULL
) {

  ## 1. default init list
  init <- init_param_default(indexprior = indexprior,
                             link       = link)

  ## 2. mapping flat args → nested components
  arg_map <- c(
    index        = "index",
    index_nu     = "index.nu",
    index_psi    = "index.psi",
    index_pi     = "index.pi",

    link_beta        = "link.beta",
    link_k           = "link.k",
    link_knots       = "link.knots",
    link_lengthscale = "link.lengthscale",
    link_amp         = "link.amp",
    link_kappa       = "link.kappa",
    link_inv_lambda  = "link.inv_lambda",

    sigma2       = "sigma2"
  )

  if (indexprior == "spike"){ #bsSpike, gpSpike
    arg_map$index <- "index.index"
  }

  ## 3. apply overrides
  env <- environment()

  for (arg_name in names(arg_map)) {
    val <- get(arg_name, envir = env, inherits = FALSE)

    if (!is.null(val)) {
      path <- strsplit(arg_map[[arg_name]], ".", fixed = TRUE)[[1L]]
      init <- set_nested(init, path, val, strict = TRUE)
    }
  }

  init
}

init_param_default <- function(indexprior, link){
  initList <- list(index = NULL, link = NULL, sigma2 = NULL)
  if (link == "bspline"){
    if (indexprior == "fisher"){
      initList <- list(index = NULL, link = list(beta = NULL), sigma2 = 0.01)

    } else if (indexprior == "sphere"){
      initList <-list(index = list(nu = NULL, index = NULL),
                      link = list(k = NULL, knots = NULL, beta = NULL),
                      sigma2 = 0.01)

    } else if (indexprior == "polar"){
      initList <- list(index = list(psi = NULL), link = list(beta = NULL), sigma2 = 0.01)

    } else if (indexprior == "spike"){
      initList <- list(index = list(pi = 0.5, nu = NULL, index = NULL),
                       link = list(beta = NULL),
                       sigma2 = 0.01)
    } else{
      stop("Error: Wrong index prior name!")
    }

  } else if (link == "gp"){
    if (indexprior == "fisher"){
      initList <- list(index = NULL,
                       link = list(lengthscale = 0.1, amp = 1),
                       sigma2 = 1)

    } else if (indexprior == "sphere"){
      initList <- list(index = NULL,
                       link = list(lengthscale = 0.1, amp = 1),
                       sigma2 = 1)

    } else if (indexprior == "polar"){
      initList <- list(index = list(psi = NULL),
                       link = list(kappa = 2),
                       sigma2 = 0.01)

    } else if (indexprior == "spike"){
      initList <- list(index = list(pi = 0.5, nu = NULL, index = NULL),
                       link = list(inv_lambda = NULL),
                       sigma2 = NULL)

    } else{
      stop("Error: Wrong index prior name!")
    }

  } else{
    stop("Error: Wrong link function name!")
  }

  return(initList)

}
