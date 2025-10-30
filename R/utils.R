initFunction_bP <- function(X, Y,
                            psi, sigma2, beta,
                            df, degree, delta,
                            index_direction,
                            setSeed){

  if (setSeed != FALSE){
    set.seed(setSeed)
  }

  n <- dim(X)[1]
  p <- dim(X)[2]

  # hyper-parameter d
  init_d <- rep(6000,p-1)


  # index
  if (is.null(psi)){
    init_psi <- runif(p-1,0.01,1)
    init_index <- alphaTheta(init_psi)
  } else{
    init_psi <- psi
    init_index <- alphaTheta(init_psi * pi)
  }
  if (sum((init_index)^2) != 1){
    init_index <- init_index/sqrt(sum(init_index^2))
  }


  # sigma2
  init_sigma2 <- sigma2
  init_Xlin <- as.vector(X %*% matrix(init_index, ncol = 1))
  init_Xmat <- transX_fisher(init_Xlin, df = df, degree = degree, delta = delta)
  init_rho <- gvcCV(init_Xmat, Y)

  # beta
  if (is.null(beta)){
    init_beta <- estBeta_fisher(init_Xmat, Y, init_rho)
  } else{
    init_beta <- beta
  }
  init_linkFunction <- init_Xmat %*%  init_beta

  return(list(d = init_d, psi = init_psi,
              index = init_index, Xlin = init_Xlin,
              Xmat = init_Xmat,
              sigma2 = init_sigma2, beta = init_beta,
              linkFunction = init_linkFunction))

}


initFunction_bSpike <- function(X, Y,
                                pi, nu, index,
                                sigma2, beta,
                                df, degree, delta,
                                setSeed){

  if (setSeed != FALSE){
    set.seed(setSeed)
  }
  n <- dim(X)[1]
  p <- dim(X)[2]

  ## pi, nu, theta
  init_pi <- pi
  if (length(pi) >= 2 || !is.numeric(pi)){
    stop("Error: Initial value of pi should be scalar.")
  }
  if (pi < 0){
    stop("Error: Initial value of pi should be positive.")
  }

  ## nu, theta
  if (!is.null(nu) & (length(nu) != p)){
    stop("Error: Incorrect dimension on initial value of nu.")
  }

  if (!is.null(nu) & (!all(nu %in% c(0, 1)))){
    stop("Error: Delta should consist of 0 or 1.")
  }

  if (is.null(nu) & is.null(index)){
    init_nu <- c(1, rbinom(p-1, 1 ,prob = init_pi))
  } else if (is.null(nu) & !is.null(index)){
    init_nu <- as.numeric(index!=0)
  } else{
    init_nu <- nu
  }


  if (!is.null(index) & (length(index) != p)){
    stop("Error: Incorrect dimension on initial value of index")
  }

  if (is.null(index)){
    nonzero <- sum(init_nu)
    nonzeroidx <- which(init_nu == 1)
    temp <- rtruncnorm(n = sum(nonzeroidx), 0, 1, lower = 0, upper = Inf)
    init_index <- rep(0, p)
    init_index[nonzeroidx] <-temp[nonzeroidx]
  } else{
    bad_zero  <- which(init_nu == 0 & index != 0)
    bad_nonzero <- which(init_nu == 1 & index == 0)

    if (length(bad_zero) > 0) {
      warning(sprintf(
        "Theta is nonzero at positions where nu = 0: %s",
        paste(bad_zero, collapse = ", ")
      ))
    }

    if (length(bad_nonzero) > 0) {
      warning(sprintf(
        "Theta is zero at positions where nu = 1: %s",
        paste(bad_nonzero, collapse = ", ")
      ))
    }

    init_index <- index * init_nu
  }

  if (sum(init_index) < 0){
    init_index <- (-1) * init_index
  }

  init_Xlin <- as.vector(X %*% matrix(init_index, nrow = p))


  # sigma2
  init_sigma2 <- sigma2
  init_Xlin <- as.vector(X %*% matrix(init_index, ncol = 1))
  init_Xmat <- transX_fisher(init_Xlin, df = df, degree = degree, delta = delta)
  init_rho <- gvcCV(init_Xmat, Y)

  # beta
  if (is.null(beta)){
    init_beta <- estBeta_fisher(init_Xmat, Y, init_rho)
  } else{
    init_beta <- beta
  }
  init_linkFunction <- init_Xmat %*%  init_beta

  return(list(pi = init_pi, nu = init_nu,
              index_raw = init_index,
              index_temp = init_index,
              index = init_index,
              Xlin = init_Xlin,
              Xmat = init_Xmat,
              sigma2 = init_sigma2, beta = init_beta,
              linkFunction = init_linkFunction))

}


initFunction_gpFisher <- function(X, Y,
                                  index, lengthscale,
                                  amp, sigma2,
                                  setSeed){
  if (setSeed != FALSE){
    set.seed(setSeed)
  }

  N <- nrow(X)
  p <- ncol(X)

  # index
  if (is.null(index)){
    # init_index <- rvMFnim(1, index_direction)
    init_index <- as.vector(ppr(X, Y, nterms = p)$beta)
  } else{
    init_index <- index
  }
  if (sum((init_index)^2) != 1){
    init_index <- init_index/sqrt(sum(init_index^2))
  }

  init_Xlin <- as.vector(X %*% matrix(init_index, nrow = p))

  # lengthscale
  init_lengthscale <- lengthscale
  if (length(init_lengthscale) >= 2 || !is.numeric(init_lengthscale)){
    stop("Error: Initial value of lengthscale should be scalar.")
  }
  if (init_lengthscale < 0){
    stop("Error: Initial value of lengthscale should be positive.")
  }

  # amp (eta)
  init_amp <- amp
  if (length(init_amp) >= 2 || !is.numeric(init_amp)){
    stop("Error: Initial value of amp should be scalar.")
  }
  if (init_amp < 0){
    stop("Error: Initial value of amp should be positive.")
  }

  init_cov <- expcov_gpSphere(init_Xlin, init_lengthscale, init_amp)
  init_linkFunction <- mvtnorm::rmvnorm(1, rep(0, N), sigma = init_cov)[1,]

  # sigma2
  init_sigma2 <- sigma2
  if (length(init_sigma2) >= 2 || !is.numeric(init_sigma2)){
    stop("Error: Initial value of sigma2 should be scalar.")
  }
  if (init_sigma2 < 0){
    stop("Error: Initial value of sigma2 should be positive.")
  }
  init_Simga <- init_sigma2 * diag(1, N)

  return(list(index = init_index, index_temp = init_index,
              lengthscale = init_lengthscale, amp = init_amp,
              Xlin = init_Xlin,cov = init_cov,
              sigma2 = init_sigma2,
              linkFunction = init_linkFunction, Sigma = init_Simga))

}


prior.param.default <- function(indexprior, link){
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
      priorList <- list(index = list(nu = list(r1 = 1, r2 = 1), index = list(sigma_theta = 0.25)),
                        link = list(basis = list(df = 21, degree = 2, delta = 0.001),
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
                        link = list(kappa = list(min_kappa= 0.5, max_kappa = 4, grid.width = 0.1)),
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

}



prior.param.default <- function(indexprior, link){
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
      priorList <- list(index = list(nu = list(r1 = 1, r2 = 1), index = list(sigma_theta = 0.25)),
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
                        link = list(kappa = list(min_kappa= 0.5, max_kappa = 4, grid.width = 0.1)),
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


init.param.default <- function(indexprior, link){
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
      initList <- list(index = list(index = NULL),
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


param.check <- function(user, template){
  old <- options(error = NULL)
  on.exit(options(old))

  wrong <- setdiff(names(user), names(template))
  if (length(wrong) > 0) {
    stop("Invalid field(s): ", paste(wrong, collapse = ", "), call. = FALSE)
  }

  for (nm in names(user)) {
    if (is.list(user[[nm]]) && is.list(template[[nm]])) {
      param.check(user[[nm]], template[[nm]])
    }
  }

  return(TRUE)
}

# check parameters - function
validate_and_finalize_args <- function(
    sampling, fitted, niter, nburnin, thin, thin2, nchain,
    prior, init, indexprior, link
) {
  if (!is.logical(sampling)){
    stop("'sampling' argument should be logical.")
  }
  if (is.logical(sampling) & (length(sampling) > 1)){
    stop("'sampling' argument should be scalar.")
  }
  if (!is.logical(fitted)){
    stop("'fitted' argument should be logical.")
  }
  if (is.logical(fitted) & (length(fitted) > 1)){
    stop("'fitted' argument should be scalar.")
  }
  if (!is.numeric(niter) || length(niter) != 1 || is.na(niter)){
    stop("'niter' argument should be numeric scalar.")
  }
  if (!is.numeric(nburnin) || length(nburnin) != 1 || is.na(nburnin)){
    stop("'nburnin' argument should be numeric scalar.")
  }
  if (!is.numeric(thin) || length(thin) != 1 || is.na(thin)){
    stop("'thin' argument should be numeric scalar.")
  }
  if (!is.null(thin2)){
    if (!is.numeric(thin2) || length(thin2) != 1 || is.na(thin2)){
      stop("'thin2' argument should be numeric scalar.")
    }
  }
  if (!is.numeric(nchain) || length(nchain) != 1 || is.na(nchain)){
    stop("'nchain' argument should be numeric scalar.")
  }

  # Prior & init
  priorlist <- prior.param.default(indexprior, link)
  if (!is.null(prior) & param.check(prior, priorlist)){
    priorlist_final <- modifyList(priorlist, prior, keep.null = TRUE)
  } else if (is.null(prior)){
    priorlist_final <- priorlist
  }

  # Initial value
  initlist <- init.param.default(indexprior, link)
  if (!is.null(init) & param.check(init, initlist)){
    initlist_final <- modifyList(initlist, init, keep.null = TRUE)
  } else if (is.null(init)){
    initlist_final <- initlist
  }

  list(
    priorlist_final = priorlist_final,
    initlist_final  = initlist_final
  )
}


# prediction
pred_fitted <- nimbleFunction(
  run = function(linkFunction_samples = double(2),
                 sigma2_samples = double(1)){
    returnType(double(2))

    nsamp <- nimDim(linkFunction_samples)[1]
    N <- nimDim(linkFunction_samples)[2]
    pred <- nimMatrix(0, nrow = nsamp, ncol = N)

    for (i in 1:nsamp){
      linkPred <- linkFunction_samples[i, ]
      predsigma <- chol(diag(rep(sigma2_samples[i], N)))
      pred[i,] <- rmnorm_chol(1, mean = linkPred,
                              cholesky = predsigma, prec_param  = FALSE)
    }
    return(pred)
  }
)

