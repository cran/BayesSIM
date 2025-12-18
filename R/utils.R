#' @importFrom ggplot2 ggplot aes geom_point geom_line labs theme_minimal theme geom_ribbon
# Initial values setting for bsPolar
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

# Initial values setting for bsSpike
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

# Initial values setting for gpFisher
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
    sampling, niter, nburnin, thin, nchain,
    prior, init, indexprior, link
) {
  if (!is.logical(sampling)){
    stop("'sampling' argument should be logical.")
  }
  if (is.logical(sampling) & (length(sampling) > 1)){
    stop("'sampling' argument should be scalar.")
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

  if (!is.numeric(nchain) || length(nchain) != 1 || is.na(nchain)){
    stop("'nchain' argument should be numeric scalar.")
  }

  # Prior & init
  priorlist <- prior_param_default(indexprior, link)
  if (!is.null(prior) & param.check(prior, priorlist)){
    priorlist_final <- modifyList(priorlist, prior, keep.null = TRUE)
  } else if (is.null(prior)){
    priorlist_final <- priorlist
  }

  # Initial value
  initlist <- init_param_default(indexprior, link)
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

# Make samples into one matrix
sampleBind <- function(mcmc.out, nchain){
  samples <- NULL
  if (nchain == 1 & !is.list(mcmc.out)){
    samples <- mcmc.out
  } else{
    for (i in 1:nchain){
      samples <- rbind(samples, mcmc.out[[i]])
    }
  }
  return(samples)
}

# Get monitorable variables - used in getVarMonitor()
expand_monitor <- function(vars) {
  out <- character(0)
  baseName <- character(0)

  for (v in vars) {

    # If there is no [], add
    if (!grepl("\\[", v)) {
      if (v == "sigma2") next
      out <- c(out, v)
      baseName <- c(baseName, v)
      next
    }

    # extract variable name
    base <- sub("\\[.*$", "", v)        # "index[1:4]" -> "index"
    if (grepl("index", base)) next
    if (base == "nu") next
    baseName <- c(baseName, base)

    # extract inside []
    inside <- sub("^.*\\[", "", v)      # "index[1:4]"     -> "1:4]"
    inside <- sub("\\]$", "", inside)   # "1:4]"           -> "1:4"

    parts <- strsplit(inside, ",")[[1]]
    parts <- trimws(parts)

    first <- parts[1]

    if (!grepl(":", first)) {
      out <- c(out, v)
      next
    }

    # a:b parsing
    ab <- as.integer(strsplit(first, ":")[[1]])
    a  <- ab[1]
    b  <- ab[2]

    rest <- ""
    if (length(parts) > 1) {
      rest <- paste0(", ", paste(parts[-1], collapse = ","))
    }

    for (k in seq(a, b)) {
      out <- c(out, sprintf("%s[%d%s]", base, k, rest))
    }
  }

  return(list(baseName = unique(baseName), out = out))
}

# Inside all the modeling functions (Extract coef, fitted, residuals etc.)
bsimFit_pointest <- function(samples, dataX, datay){
  # samples: combined samples of chain
  # relevant information
  p <- ncol(dataX)
  N <- nrow(dataX)
  namesIndex <- paste0("index[", 1:p, "]")
  namesXlin <- paste0("Xlin[", 1:N, "]")
  # link function
  colsN <- colnames(samples)
  namesLink <- grep("linkFunction", colsN, value = TRUE)
  # has_comma <- grepl(",", lf_cols)
  # namesLink <- ifelse(
  #   has_comma,
  #   paste0(lf_cols, paste0("linkFunction[", 1:N, ", 1]")),
  #   paste0(lf_cols, paste0("linkFunction[", 1:N, "]"))
  # )
  namesSigma <- "sigma2"
  namesVar <- colnames(dataX)
  indexName <- colnames(dataX)
  if (is.null(samples)){
    coefficients <- ses_coef <- se <- linear.predictors <- fitted.values <- residuals <- gof <- NULL
  } else{
    ALLsamp <- samples


    # Compute coefficients, sd, fitted.values, linear.predictors - mean
    coefficients <- apply(ALLsamp[ ,namesIndex], 2, mean)
    attr(coefficients, "names") <- indexName
    ses_coef <- apply(ALLsamp[ ,namesIndex], 2, sd)
    attr(ses_coef, "names") <- indexName
    se <- mean(ALLsamp[ ,namesSigma])
    linear.predictors <- apply(ALLsamp[ ,namesXlin], 2, mean)
    fitted.values <- apply(ALLsamp[ ,namesLink], 2, mean)
    residuals <- datay[,1] - as.vector(fitted.values)
    gof <- mean(residuals^2) # residual sum of squares
  }


  outList <- list(coefficients = coefficients,
                  ses_coef = ses_coef,
                  se = se,
                  residuals = residuals,
                  linear.predictors = as.vector(linear.predictors),
                  fitted.values = as.vector(fitted.values),
                  gof = gof)

  return(outList)

}

# Fitted plot algorithm
plot_bsim_fitted <- function(result, interval){ # result: data.frame
  idxValue <- 0; pred <- 0; LB <- 0; UB <- 0; truey <- 0

  level <- result$level
  if (is.null(result$truey)){

    if (!interval){
      if (is.vector(result$fitted)){
        pred <- result$fitted
      } else{
        pred <- result$fitted[ ,1]
      }
      results <- data.frame(pred = pred, idxValue = result$idxValue)
      df_pred_line <- results[order(results$idxValue), ]
      A2 <- ggplot(df_pred_line, aes(x = idxValue)) +
        geom_point(aes(y = pred), color = "grey20", alpha = 0.55, size = 1.6) +
        geom_line(aes(y = pred), color = "firebrick", linewidth = 1) +
        labs(x = "Index value", y = "Predicted",
             title = paste0("Fitted curve"),
             subtitle = "Points: observed, Line: posterior mean") +
        theme_minimal(base_size = 12) +
        theme(
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(face = "bold")
        )

    } else{
      results2 <- data.frame(pred = result$fitted[,1],
                             idxValue = result$idxValue, LB = result$fitted$LB,
                             UB = result$fitted$UB)
      df_pred_line2 <- results2[order(results2$idxValue), ]
      A2 <- ggplot(df_pred_line2, aes(x = idxValue)) +
        geom_ribbon(aes(ymin = LB, ymax = UB),
                    fill = "firebrick", alpha = 0.12, colour = NA) +
        geom_point(aes(y = pred), color = "grey20", alpha = 0.55, size = 1.6) +
        geom_line(aes(y = pred), color = "firebrick", linewidth = 1) +
        labs(x = "Index value", y = "Predicted",
             title = paste0("Fitted curve with ", level, "% interval"),
             subtitle = "Points: predicted, Line: posterior mean") +
        theme_minimal(base_size = 12) +
        theme(
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(face = "bold")
        )
    }
    return(A2)
  } else{
    if (!interval){
      if (is.vector(result$fitted)){
        pred <- result$fitted
      } else{ # se.fit = TRUE
        pred <- result$fitted[ ,1]
      }
      results <- data.frame(truey = result$truey, pred = pred,
                            idxValue = result$idxValue)
    } else{
      results <- data.frame(truey = result$truey, pred = result$fitted[,1],
                            idxValue = result$idxValue)
    }


    # 1) fitted vs. observed
    if (!is.null(result$truey)){
      rmse <- mean((results$truey-results$pred)^2)

      A1 <- ggplot(results, aes(x = truey, y = pred)) +
        geom_point(color = "black", alpha = 0.6) +
        geom_abline(slope = 1, intercept = 0, color = "firebrick", linewidth = 0.8) +
        theme_minimal(base_size = 14) +
        labs(
          title = "Fitted Plot",
          x = "True",
          y = "Predicted",
          subtitle = paste("RMSE:", round(rmse, 4))
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(face = "bold"),
          panel.grid.minor = element_blank()
        )
    }


    # 2) index vs. y
    if (!interval){
      df_pred_line <- results[order(results$idxValue), ]
      A2 <- ggplot(df_pred_line, aes(x = idxValue)) +
        geom_point(aes(y = truey), color = "grey20", alpha = 0.55, size = 1.6) +
        geom_line(aes(y = pred), color = "firebrick", linewidth = 1) +
        labs(x = "Index value", y = "Response",
             title = paste0("Fitted curve"),
             subtitle = "Points: observed, Line: posterior mean") +
        theme_minimal(base_size = 12) +
        theme(
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(face = "bold")
        )

    } else{
      results2 <- data.frame(truey = result$truey, pred = result$fitted[,1],
                             idxValue = result$idxValue, LB = result$fitted$LB,
                             UB = result$fitted$UB)
      df_pred_line2 <- results2[order(results2$idxValue), ]
      A2 <- ggplot(df_pred_line2, aes(x = idxValue)) +
        geom_ribbon(aes(ymin = LB, ymax = UB),
                    fill = "firebrick", alpha = 0.12, colour = NA) +
        geom_point(aes(y = truey), color = "grey20", alpha = 0.55, size = 1.6) +
        geom_line(aes(y = pred), color = "firebrick", linewidth = 1) +
        labs(x = "Index value", y = "Response",
             title = paste0("Fitted curve with ", level, "% interval"),
             subtitle = "Points: observed, Line: posterior mean") +
        theme_minimal(base_size = 12) +
        theme(
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(face = "bold")
        )
    }
    return(A1 + A2)
  }

}



