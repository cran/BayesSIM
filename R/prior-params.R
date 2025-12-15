#' Build a Prior List for BayesSIM Models
#'
#' @description
#' `prior_param` is a convenience helper that constructs a nested prior list
#' for a given combination of index prior and link function.
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
#' @param index_direction,index_dispersion,index_nu_r1,index_nu_r2,index_psi_alpha,index_sigma_theta,index_r1,index_r2
#'   Optional overrides for hyper-parameters related to the index prior.
#'
#' @param link_basis_df,link_basis_degree,link_basis_delta
#'   Optional overrides for the B-spline basis hyper-parameters, such as the effective degrees of freedom,
#'   spline degree, and penalty parameter.
#' @param link_knots_lambda_k,link_knots_maxknots
#'   Optional overrides for the B-spline knot-selection hyper-parameters in, used for models with adaptive knot placement.
#'
#' @param link_beta_mu,link_beta_cov,link_beta_tau,link_beta_Sigma0
#'   Optional overrides for the prior on spline coefficients. The detailed interpretation of these
#'   hyper-parameters depends on the specific model and is described in the
#'   documentation of each model-fitting function.
#'
#' @param link_lengthscale_shape,link_lengthscale_rate
#'   Optional overrides for the hyper-parameters of the GP length-scale prior.
#'
#' @param link_amp_a,link_amp_b
#'   Optional overrides for the hyper-parameters of the GP amplitude prior.
#'
#' @param link_kappa_min,link_kappa_max,link_kappa_grid_width
#'   Optional overrides for the hyper-parameters in used in models with polar index and GP-type link,
#'   to control the grid or support for the concentration parameter \eqn{\kappa} in `gpPolar`.
#'
#' @param link_inv_lambda_shape,link_inv_lambda_rate
#'   Optional overrides for spike-and-slab–type GP link priors.
#'
#' @param sigma2_shape,sigma2_rate
#'   Optional overrides for the inverse-gamma prior on the observation
#'   variance \eqn{\sigma^2}.
#'
#' @details
#'  \code{prior_param(indexprior, link)} can be used to obtain the default prior
#' list for the requested combination of index prior and link function.
#' For any argument that is not `NULL`, the corresponding field in the nested prior list is overwritten.
#'
#'
#' The detailed meaning and recommended choices for each hyper-parameter depend
#' on the specific model, prior of index vector and link function.
#' For those details, please refer to the documentation of the corresponding
#' model-fitting functions.
#'
#' @return
#' A nested list with top-level elements \code{index}, \code{link}, and
#' \code{sigma2}, suitable for passing to the \code{prior} argument of the
#' various BayesSIM model fitting functions.
#'
#' @seealso [bsFisher()], [bsSphere()], [bsPolar()], [bsSpike()],
#' [gpFisher()], [gpSphere()], [gpPolar()], [gpPolarHigh()], [gpSpike()]
#'
#' @examples
#' ## Default prior for Fisher index + B-spline link:
#' p0 <- prior_param("fisher", "bspline")
#'
#' ## Modify only a few hyper-parameters:
#' p1 <- prior_param(
#'   indexprior          = "fisher",
#'   link                = "bspline",
#'   sigma2_shape        = 0.5,
#'   link_basis_df       = 30,
#'   index_direction     = c(1, 0, 0)
#' )
#' @export
prior_param <- function(indexprior, link,
                        index_direction = NULL, index_dispersion = NULL,
                        index_nu_r1 = NULL, index_nu_r2 = NULL,
                        index_psi_alpha = NULL, index_sigma_theta = NULL,
                        index_r1 = NULL, index_r2 = NULL,
                        link_basis_df = NULL, link_basis_degree = NULL,
                        link_basis_delta = NULL, link_knots_lambda_k = NULL,
                        link_knots_maxknots = NULL, link_beta_mu = NULL,
                        link_beta_cov = NULL, link_beta_tau = NULL,
                        link_beta_Sigma0 = NULL, link_lengthscale_shape = NULL,
                        link_lengthscale_rate = NULL, link_amp_a = NULL,
                        link_amp_b = NULL, link_kappa_min = NULL,
                        link_kappa_max = NULL, link_kappa_grid_width = NULL,
                        link_inv_lambda_shape = NULL, link_inv_lambda_rate = NULL,
                        sigma2_shape = NULL, sigma2_rate = NULL) {

  ## 1. default prior
  prior <- prior_param_default(indexprior = indexprior,
                               link       = link)

  ## 2. mapping table
  arg_map <- c(
    index_direction        = "index.direction",
    index_dispersion       = "index.dispersion",
    index_nu_r1            = "index.nu.r1",
    index_nu_r2            = "index.nu.r2",
    index_psi_alpha        = "index.psi.alpha",
    index_r1               = "index.r1",
    index_r2               = "index.r2",
    index_sigma_theta      = "index.sigma_theta",

    link_basis_df          = "link.basis.df",
    link_basis_degree      = "link.basis.degree",
    link_basis_delta       = "link.basis.delta",

    link_knots_lambda_k    = "link.knots.lambda_k",
    link_knots_maxknots    = "link.knots.maxknots",

    link_beta_mu           = "link.beta.mu",
    link_beta_cov          = "link.beta.cov",
    link_beta_tau          = "link.beta.tau",
    link_beta_Sigma0       = "link.beta.Sigma0",

    link_lengthscale_shape = "link.lengthscale.shape",
    link_lengthscale_rate  = "link.lengthscale.rate",
    link_amp_a         = "link.amp.a_amp",
    link_amp_b         = "link.amp.b_amp",

    link_kappa_min   = "link.kappa.min",
    link_kappa_max   = "link.kappa.max",
    link_kappa_grid_width  = "link.kappa.grid_width",

    link_inv_lambda_shape  = "link.inv_lambda_shape",
    link_inv_lambda_rate   = "link.inv_lambda_rate",

    sigma2_shape           = "sigma2.shape",
    sigma2_rate            = "sigma2.rate"
  )

  env <- environment()


  for (arg_name in names(arg_map)) {
    val <- get(arg_name, envir = env, inherits = FALSE)
    if (!is.null(val)) {
      path <- strsplit(arg_map[[arg_name]], ".", fixed = TRUE)[[1L]]
      prior <- set_nested(
        lst   = prior,
        path  = path,
        value = val,
        strict = TRUE
      )
    }
  }


  return(prior)


}



set_nested <- function(lst, path, value, parent_path = NULL, strict = TRUE) {
  nm <- path[1L]
  full_path <- c(parent_path, nm)

  # Last stage: assignment
  if (length(path) == 1L) {
    if (strict && !is.null(lst) && !nm %in% names(lst)) {
      stop("Unknown prior field: ",
           paste(full_path, collapse = " -> "),
           call. = FALSE)
    }
    if (is.null(lst)) lst <- list()
    lst[[nm]] <- value
    return(lst)
  }

  # middle state: check the path
  if (is.null(lst[[nm]])) {
    if (strict) {
      stop("Unknown prior path: ",
           paste(full_path, collapse = " -> "),
           call. = FALSE)
    } else {
      lst[[nm]] <- list()
    }
  } else if (!is.list(lst[[nm]])) {
    stop("Cannot descend into non-list prior field: ",
         paste(full_path, collapse = " -> "),
         call. = FALSE)
  }

  lst[[nm]] <- set_nested(
    lst[[nm]],
    path        = path[-1L],
    value       = value,
    parent_path = full_path,
    strict      = strict
  )
  lst
}

prior_param_default <- function(indexprior, link){
  priorList <- list(index = NULL, link = NULL, sigma2 = NULL)
  if (link == "bspline"){
    if (indexprior == "fisher"){
      priorList <- list(index = list(direction = NULL, dispersion = 150),
                        link = list(basis = list(df = 21, degree = 2, delta = 0.001),
                                    beta = list(mu = NULL, cov = NULL)),
                        sigma2 = list(shape = 0.001, rate = 100))

    } else if (indexprior == "sphere"){
      priorList <- list(index = list(nu = list(r1 = 1, r2 = 1)),
                        link = list(knots = list(lambda_k = 5, maxknots = NULL),
                                    basis = list(degree = 2),
                                    beta = list(mu = NULL, tau = NULL,Sigma0 = NULL)),
                        sigma2 = list(shape = 0.001, rate = 0.001))

    } else if (indexprior == "polar"){
      priorList <- list(index = list(psi = list(alpha = NULL)),
                        link = list(basis = list(df = 21,degree = 2, delta = 0.001),
                                    beta = list(mu = NULL, cov = NULL)),
                        sigma2 = list(shape = 0.001, rate = 100))

    } else if (indexprior == "spike"){
      priorList <- list(index = list(nu = list(r1 = 1, r2 = 1), sigma_theta = 0.25),
                        link = list(basis = list(df = 21, degree = 2, delta = 0.01),
                                    beta = list(mu = NULL, cov = NULL)),
                        sigma2 = list(shape = 0.001, rate = 100))


    } else{
      stop("Error: Wrong index prior name!")
    }

  } else if (link == "gp"){
    if (indexprior == "fisher"){
      priorList <- list(index = list(direction = NULL, dispersion = 150),
                        link = list(lengthscale = list(shape = 1/8, rate = 1/8),
                                    amp = list(a_amp = -1, b_amp = 1)),
                        sigma2 = list(shape = 1, rate = 1))

    } else if (indexprior == "sphere"){
      priorList <- list(index = NULL,
                        link = list(lengthscale = list(shape = 1/8, rate = 1/8),
                                    amp = list(a_amp = -1, b_amp = 1)),
                        sigma2 = list(shape = 1, rate = 1))

    } else if (indexprior == "polar"){
      priorList <- list(index = list(psi = list(alpha = NULL)),
                        link = list(kappa = list(min= 0.5, max = 4, grid_width = 0.1)),
                        sigma2 = list(shape = 2, rate = 0.01))

    } else if (indexprior == "spike"){
      priorList <- list(index = list(r1 = 1, r2 = 1, sigma_theta = 0.25),
                        link = list(inv_lambda_shape = 1, inv_lambda_rate = 0.1),
                        sigma2 = list(shape = 0.001, rate = 0.001))

    } else{
      stop("Error: Wrong index prior name!")
    }

  } else{
    stop("Error: Wrong link function name!")
  }

  return(priorList)

}


# prior1 <- prior_param(
#   indexprior = "fisher",
#   link       = "bspline",
#   sigma2_shape    = 0.5,               # sigma2$shape
#   sigma2_rate     = 10,                # sigma2$rate
#   link_basis_df   = 30,                # link$basis$df
#   index_direction = c(1, 0, 0, 0)      # index$direction
# )
#
# str(prior1)
