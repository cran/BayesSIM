#' @import methods
#' @import patchwork
#' @importFrom magrittr %>%
#' @import ggplot2
#' @importFrom tidyr pivot_longer all_of
#' @importFrom dplyr select
#' @importFrom stats mad

# Common functions and S3 class for the BayesSIM
## bsim: Total results for fitted model
## bsimSetup: Model setup without sampling
## bsimPred: object for prediction -> used in fitted plots
bsim <- function(coefficients, ses, residuals,
                 fitted.values, linear.predictors, gof,
                 input, samples, model, sampler, modelName) {
  structure(
    list(
      coefficients = coefficients,
      ses = ses, residuals = residuals,
      fitted.values = fitted.values,
      linear.predictors = linear.predictors,
      gof = gof,
      input = input,
      samples = samples,
      defModel = model, defSampler = sampler,
      modelName = modelName
    ),
    class = "bsim"
  )
}

bsimSetup <- function(coefficients, ses, residuals,
                      fitted.values, linear.predictors, gof,
                      input, samples, model, sampler, modelName) {
  structure(
    list(
      input = input,
      defModel = model,
      defSampler = sampler,
      modelName = modelName
    ),
    class = "bsimSetup"
  )
}

bsim.predict <- function(fitted, truey, idxValue, level){
  structure(
    list(
      fitted = fitted,
      truey = truey,
      idxValue = idxValue,
      level = level
    ),
    class = "bsimPred"
  )
}

#' Construct a Fitted Model Object from Model Setup and MCMC Output
#'
#' Create a fitted \code{bsim} object by combining a `BayesSIM`
#' setup object with MCMC samples returned by \code{runMCMC()}.
#'
#' @param setup A `BayesSIM` setup object, typically the output of a
#'   \code{_setup} function.
#' @param mcmc.out MCMC output corresponding to the result of a call to \code{runMCMC()}.
#'
#' @return
#' An object of class \code{"bsim"} containing posterior samples,
#' point estimates, fitted values, and related model information.
#'
#' @details
#' This function is mainly intended for workflows where the model structure
#' and the MCMC sampling are performed separately. It collects the MCMC draws across chains, and returns an object of class \code{"bsim"}
#' that can be used with generic functions such as \code{summary()}, \code{plot()}, and \code{predict()}.
#'
#'
#' @examples
#' \donttest{
#' simdata2 <- data.frame(DATA1$X, y = DATA1$y)
#' models <- BayesSIM_setup(y ~ ., data = simdata2)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                    nchains = 1, setSeed = TRUE, inits = initList,
#'                    summary = TRUE, samplesAsCodaMCMC = TRUE)
#' fit2 <- as_bsim(models, mcmc.out)
#' summary(fit2)
#' }
#'
#' @export
as_bsim <- function(setup, mcmc.out){
  nchain <- setup$input$samplingOptions$nchain
  # samples <- sampleBind(mcmc.out, nchain)
  if ("samples" %in% names(mcmc.out)){
    temp <- mcmc.out[["samples"]]
  } else{
    temp <- mcmc.out
  }

  # gpSphere - EB
  if (setup$modelName == "gpSphere"){
    if (setup$input$samplingOptions$method == "EB"){
      map_index <- setup$defModel$getConstants()$index
      map_amp <- setup$defModel$getConstants()$amp
      map_lengthscale <- setup$defModel$getConstants()$lengthscale
      map_sigma2 <- setup$defModel$getConstants()$sigma2
      p <- length(map_index)

      sampMCMC <- as.matrix(temp)
      nmcmcsamp <- nrow(sampMCMC)
      index_mat <- matrix(rep(map_index, nmcmcsamp), ncol = p, byrow = TRUE)
      colnames(index_mat) <- paste0("index[",1:p, "]")
      hype_mat <- matrix(rep(c(map_lengthscale, map_amp, map_sigma2), nmcmcsamp),
                         ncol = 3, byrow = TRUE)
      colnames(hype_mat) <- c("lengthscale", "amp", "sigma2")
      temp <- cbind(sampMCMC, index_mat, hype_mat)
    }
  }

  # gpSpike
  samples <- sampleBind(temp, nchain)


  X <- setup$input$origdata$x
  Y <- setup$input$origdata$y
  modelESTLIST <- bsimFit_pointest(samples, X, Y)
  out <- list(
    coefficients = modelESTLIST$coefficients,
    ses_coef = modelESTLIST$ses_coef, se = modelESTLIST$se,
    residuals = modelESTLIST$residuals,
    fitted.values = modelESTLIST$fitted.values,
    linear.predictors = modelESTLIST$linear.predictors,
    gof = modelESTLIST$gof,
    input = setup$input,
    samples = temp,
    defModel = setup$defModel, defSampler = setup$defSampler,
    modelName = setup$modelName
  )

  class(out) = "bsim"
  return(out)

}

#' @rdname BayesSIM
#' @method print bsim
#' @export
print.bsim <- function(x, digits = 3, ...){
  # Model definition if possible
  # 1) Header
  cat("BayesSIM model\n")

  # 2) Model info
  cat("formula:     ", x$input$formula[c(2, 1, 3)], "\n")
  cat("observations:", nrow(x$input$origdata$x), "\n")
  cat("predictors:  ", ncol(x$input$origdata$x), "\n")
  cat("------\n")
  if(is.null(x$samples)){
    cat("Try compiling the object and draw posterior samples from `runMCMC()`")
  } else{
    # estimated coefficient
    cat("Single-index estimates: \n")
    # cat(x$coefficients)
    A <- data.frame(mean = x$coefficients,
                    empirical.se = x$ses_coef)

    rownames(A) <- names(x$coefficients)

    print(round(A, digits))

    cat("\n")

    # Goodness of fit: Residual sum of square
    cat("Goodness of fit (Mean residual sum of square): ", round(x$gof, digits))
    cat("\n")
  }

}

#' Extract Index Vector Coefficients from `BayesSIM`
#'
#' Computes posterior summaries of the single-index model index vector
#' from a fitted `BayesSIM`. Users may choose either the
#' posterior mean or median as the point estimate.
#'
#' @param object A fitted object of \code{BayesSIM} or individual model.
#' @param method Character string indicating the summary statistic to compute.
#'   Options are \code{"mean"} or \code{"median"}. Default is \code{"mean"}.
#' @param se Logical value whether computing standard error for index estimates.
#' If `method` is \code{"mean"}, standard deviation of index vector MCMC samples is gained.
#' If  `method` is \code{"median"}, median absolute deviation of index vector MCMC samples is gained.
#' `FALSE` is default.
#' @param ... Additional arguments passed to other methods.
#'
#' @return A numeric vector or data.frame of estimated coefficient and standard error of the index vector.
#'
#' @examples
#' \donttest{
#' simdata2 <- data.frame(DATA1$X, y = DATA1$y)
#'
#' # 1. One tool version
#' fit_one <- BayesSIM(y ~ ., data = simdata2,
#'                     niter = 5000, nburnin = 1000, nchain = 1)
#'
#' # Check median index vector estimates with standard errors
#' coef(fit_one, method = "median", se = TRUE)
#'
#' # Fitted index values of median prediction
#' fitted(fit_one, type = "linpred", method = "median")
#'
#' # Residuals of median prediction
#' residuals(fit_one, method = "median")
#'
#' # Summary of the model
#' summary(fit_one)
#'
#' # Convergence diagnostics
#' nimTraceplot(fit_one)
#'
#' # Goodness of fit
#' GOF(fit_one)
#'
#' # Fitted plot
#' plot(fit_one)
#'
#' # Prediction with 95% credible interval at new data
#' newx <- data.frame(X1 = rnorm(10), X2 = rnorm(10), X3 = rnorm(10), X4 = rnorm(10))
#' pred <- predict(fit_one, newdata = newx, interval = "credible", level = 0.95)
#' plot(pred)
#'
#'
#' # 2. Split version
#' models <- BayesSIM_setup(y ~ ., data = simdata2)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, inits = initList,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#'
#' # "fit_split" becomes exactly the same as the class of "fit_one" object and apply generic functions.
#' fit_split <- as_bsim(models, mcmc.out)
#'
#' }
#'
#' @method coef bsim
#' @export
coef.bsim <- function(object, method = c("mean", "median"), se = FALSE, ...){
  method <- match.arg(method)

  if (sum(method == c("mean", "median")) == 2){
    method <- NULL
  }

  coefVector <- NULL

  if (is.null(method) || method == "mean"){
    coefVector <- object$coefficients

  } else if (method == "median"){
    dataX <- object$input$origdata$x
    p <- ncol(dataX)
    namesIndex <- paste0("index[", 1:p, "]")
    indexName <- colnames(dataX)
    nchain <- object$input$samplingOptions$nchain
    ALLsamp <- sampleBind(object$samples, nchain)

    coefficients_median <- apply(ALLsamp[ ,namesIndex], 2, median)
    attr(coefficients_median, "names") <- indexName
    coefVector <- coefficients_median

  } else{
    stop("method could be either mean or median.")
  }

  if (se){
    if (is.null(method) || method == "mean"){
      se.est <- object$ses_coef
    } else{
      se.est <- apply(ALLsamp[ ,namesIndex], 2, function(x) mad(x))
    }
    matResult <- rbind(coefVector, se.est)
    rownames(matResult) <- c("est", "Std.error")
    return(matResult)
  } else{
    return(coefVector)
  }

}


#' @rdname gof
#' @export
GOF <- function(object){
  UseMethod("GOF")
}

#' Goodness of Fit for `BayesSIM`
#'
#' Generic function applied to \code{BayesSIM}.
#' It extracts goodness of fit of the `BayesSIM`.
#'
#' @param object A fitted object of \code{BayesSIM} or individual model.
#' @param ... Additional arguments passed to other methods.
#'
#' @return Mean squared error of model with mean of MCMC sample is applied.
#' @method GOF bsim
#' @examples
#' \donttest{
#' simdata2 <- data.frame(DATA1$X, y = DATA1$y)
#'
#' # 1. One tool version
#' fit_one <- BayesSIM(y ~ ., data = simdata2,
#'                     niter = 5000, nburnin = 1000, nchain = 1)
#'
#' # Check median index vector estimates with standard errors
#' coef(fit_one, method = "median", se = TRUE)
#'
#' # Fitted index values of median prediction
#' fitted(fit_one, type = "linpred", method = "median")
#'
#' # Residuals of median prediction
#' residuals(fit_one, method = "median")
#'
#' # Summary of the model
#' summary(fit_one)
#'
#' # Convergence diagnostics
#' nimTraceplot(fit_one)
#'
#' # Goodness of fit
#' GOF(fit_one)
#'
#' # Fitted plot
#' plot(fit_one)
#'
#' # Prediction with 95% credible interval at new data
#' newx <- data.frame(X1 = rnorm(10), X2 = rnorm(10), X3 = rnorm(10), X4 = rnorm(10))
#' pred <- predict(fit_one, newdata = newx, interval = "credible", level = 0.95)
#' plot(pred)
#'
#'
#' # 2. Split version
#' models <- BayesSIM_setup(y ~ ., data = simdata2)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, inits = initList,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#'
#' # "fit_split" becomes exactly the same as the class of "fit_one" object and apply generic functions.
#' fit_split <- as_bsim(models, mcmc.out)
#'
#' }
#' @name gof
#' @export
GOF.bsim <- function(object, ...){
  object$gof
}

#' Extract Residuals from `BayesSIM`
#'
#' Returns the model residuals based on the posterior fitted values of a
#' `BayesSIM`. Residuals can be computed using either the posterior mean or median fitted values.
#'
#' @param object A fitted object of \code{BayesSIM} or individual model.
#' @param method Character string specifying the summary statistic used to
#'   compute the fitted values. Options are \code{"mean"} or \code{"median"}.
#'   Default is \code{"mean"}.
#' @param ... Additional arguments passed to other methods.
#'
#' @return A numeric vector of residuals. (\eqn{\mathbf{r} = \mathbf{Y} - \hat{\mathbf{Y}}})
#' \eqn{\hat{\mathbf{Y}}} can be mean or median of MCMC samples.
#'
#' @examples
#' \donttest{
#' simdata2 <- data.frame(DATA1$X, y = DATA1$y)
#'
#' # 1. One tool version
#' fit_one <- BayesSIM(y ~ ., data = simdata2,
#'                     niter = 5000, nburnin = 1000, nchain = 1)
#'
#' # Check median index vector estimates with standard errors
#' coef(fit_one, method = "median", se = TRUE)
#'
#' # Fitted index values of median prediction
#' fitted(fit_one, type = "linpred", method = "median")
#'
#' # Residuals of median prediction
#' residuals(fit_one, method = "median")
#'
#' # Summary of the model
#' summary(fit_one)
#'
#' # Convergence diagnostics
#' nimTraceplot(fit_one)
#'
#' # Goodness of fit
#' GOF(fit_one)
#'
#' # Fitted plot
#' plot(fit_one)
#'
#' # Prediction with 95% credible interval at new data
#' newx <- data.frame(X1 = rnorm(10), X2 = rnorm(10), X3 = rnorm(10), X4 = rnorm(10))
#' pred <- predict(fit_one, newdata = newx, interval = "credible", level = 0.95)
#' plot(pred)
#'
#'
#' # 2. Split version
#' models <- BayesSIM_setup(y ~ ., data = simdata2)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, inits = initList,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#'
#' # "fit_split" becomes exactly the same as the class of "fit_one" object and apply generic functions.
#' fit_split <- as_bsim(models, mcmc.out)
#'
#' }
#' @name genBasic
#' @method residuals bsim
#' @export
residuals.bsim <- function(object, method = c("mean", "median"), ...){
  method <- match.arg(method)
  # if (class(object) != "bsim"){
  #   stop("object should be 'bsim'.")
  # }

  if (sum(method == c("mean", "median")) == 2){
    method <- NULL
  }


  if (is.null(method) || method == "mean"){
    return(object$residuals)

  } else if (method == "median"){
    N <- nrow(object$input$origdata$x)
    nchain <- object$input$samplingOptions$nchain
    ALLsamp <- sampleBind(object$samples, nchain)
    colsN <- colnames(ALLsamp)
    namesLink <- grep("linkFunction", colsN, value = TRUE)
    datay <- object$input$origdata$y


    fitted.median <- apply(ALLsamp[ ,namesLink], 2, median)
    residuals.median <- datay[,1] - as.vector(fitted.median)
    return(residuals.median)

  } else{
    stop("method could be either mean or median.")
  }

}

#' Extract Fitted Values from `BayesSIM`
#'
#' Computes fitted values from a `BayesSIM`, using either the
#' posterior mean or median of the estimated link function with index values.
#' Fitted values can be returned on the latent scale or on the linear
#' predictor scale.
#'
#' @inheritParams genBasic
#' @param type Character string indicating the scale on which fitted values
#'   are returned. Default is \code{"latent"}.
#'   \itemize{
#'     \item \code{"latent"}: fitted response values \eqn{\hat{y} = E(\mathbf{Y}|\mathbf{X})}.
#'     \item \code{"linpred"}: linear predictor \eqn{X'\theta}.
#'   }
#' @return A numeric vector of fitted values.
#'
#' @examples
#' \donttest{
#' simdata2 <- data.frame(DATA1$X, y = DATA1$y)
#'
#' # 1. One tool version
#' fit_one <- BayesSIM(y ~ ., data = simdata2,
#'                     niter = 5000, nburnin = 1000, nchain = 1)
#'
#' # Check median index vector estimates with standard errors
#' coef(fit_one, method = "median", se = TRUE)
#'
#' # Fitted index values of median prediction
#' fitted(fit_one, type = "linpred", method = "median")
#'
#' # Residuals of median prediction
#' residuals(fit_one, method = "median")
#'
#' # Summary of the model
#' summary(fit_one)
#'
#' # Convergence diagnostics
#' nimTraceplot(fit_one)
#'
#' # Goodness of fit
#' GOF(fit_one)
#'
#' # Fitted plot
#' plot(fit_one)
#'
#' # Prediction with 95% credible interval at new data
#' newx <- data.frame(X1 = rnorm(10), X2 = rnorm(10), X3 = rnorm(10), X4 = rnorm(10))
#' pred <- predict(fit_one, newdata = newx, interval = "credible", level = 0.95)
#' plot(pred)
#'
#'
#' # 2. Split version
#' models <- BayesSIM_setup(y ~ ., data = simdata2)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, inits = initList,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#'
#' # "fit_split" becomes exactly the same as the class of "fit_one" object and apply generic functions.
#' fit_split <- as_bsim(models, mcmc.out)
#'
#' }
#'
#' @method fitted bsim
#' @export
fitted.bsim <- function(object,
                        type = c("latent", "linpred"),
                        method = c("mean", "median"),...){ # point estimate
  type <- match.arg(type)
  method <- match.arg(method)

  # if (class(object) != "bsim"){
  #   stop("object should be 'bsim'.")
  # }
  if (sum(type == c("latent", "linpred")) == 2){
    type <- NULL
  }

  if (sum(method == c("mean", "median")) == 2){
    method <- NULL
  }


  if (is.null(type) || type == "latent"){
    if (is.null(method) || method == "mean"){
      return(object$fitted.values)

    } else if (method == "median"){
      N <- nrow(object$input$origdata$x)
      nchain <- object$input$samplingOptions$nchain
      ALLsamp <- sampleBind(object$samples, nchain)
      colsN <- colnames(ALLsamp)
      namesLink <- grep("linkFunction", colsN, value = TRUE)

      fitted.median <- apply(ALLsamp[ ,namesLink], 2, median)
      return(as.vector(fitted.median))

    } else{ # method error
      stop("method could be either mean or median.")
    }
  } else if (type == "linpred"){
    N <- nrow(object$input$origdata$x)
    namesXlin <- paste0("Xlin[", 1:N, "]")
    nchain <- object$input$samplingOptions$nchain
    ALLsamp <- sampleBind(object$samples, nchain)

    if (is.null(method) || method == "mean"){
      linpred.mean <- apply(ALLsamp[ ,namesXlin], 2, mean)
      return(as.vector(linpred.mean))

    } else if (method == "median"){
      linpred.median <- apply(ALLsamp[ ,namesXlin], 2, median)
      return(as.vector(linpred.median))

    } else{ # method error
      stop("method could be either mean or median.")
    }

  } else{ # type error
    stop("type could be either latent or linpred.")
  }
}

#' Summarize `BayesSIM`
#' @name summary.bsim
#' @description
#' Provides a \code{summary} for `BayesSIM`.
#'
#' @param object A fitted object of \code{BayesSIM} or individual model.
#' @param x A summary output of \code{BayesSIM} or individual model.
#' @param digits The minimum number of significant digits to be printed.
#' @param ... Further arguments passed.
#'
#' @details
#' A \code{list} of summary statistics for MCMC samples, including \code{data.frame} table for the results.
#' Each row corresponds to a model parameter, and columns report the statistics.
#'
#' @return
#' The function summarizes posterior MCMC samples by reporting key statistics, including:
#' \itemize{
#'   \item Posterior mean and median
#'   \item Empirical standard deviation
#'   \item 95% credible interval (lower and upper quantiles)
#'   \item Potential scale reduction factor (\code{gelman}) for multiple chains
#'   \item Effective sample size (\code{ESS})
#' }
#'
#' By default, the index vector and error variance are only included in the summary.
#' If variable selection methods are used, such as uniform sphere and spike-and-slab prior,
#' the indicator vector (\code{nu}) is also included.
#' Note that the potential scale reduction factor for \code{nu} can be reported as
#' \code{NaN} or \code{Inf}, since the indicator rarely changes during the MCMC run.
#'
#' If the model is fitted with single chain, both \code{all.chain} and \code{chain} have identical information.
#'
#' @method summary bsim
#' @examples
#' \donttest{
#' simdata2 <- data.frame(DATA1$X, y = DATA1$y)
#'
#' # 1. One tool version
#' fit_one <- BayesSIM(y ~ ., data = simdata2,
#'                     niter = 5000, nburnin = 1000, nchain = 1)
#'
#' # Check median index vector estimates with standard errors
#' coef(fit_one, method = "median", se = TRUE)
#'
#' # Fitted index values of median prediction
#' fitted(fit_one, type = "linpred", method = "median")
#'
#' # Residuals of median prediction
#' residuals(fit_one, method = "median")
#'
#' # Summary of the model
#' summary(fit_one)
#'
#' # Convergence diagnostics
#' nimTraceplot(fit_one)
#'
#' # Goodness of fit
#' GOF(fit_one)
#'
#' # Fitted plot
#' plot(fit_one)
#'
#' # Prediction with 95% credible interval at new data
#' newx <- data.frame(X1 = rnorm(10), X2 = rnorm(10), X3 = rnorm(10), X4 = rnorm(10))
#' pred <- predict(fit_one, newdata = newx, interval = "credible", level = 0.95)
#' plot(pred)
#'
#'
#' # 2. Split version
#' models <- BayesSIM_setup(y ~ ., data = simdata2)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, inits = initList,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#'
#' # "fit_split" becomes exactly the same as the class of "fit_one" object and apply generic functions.
#' fit_split <- as_bsim(models, mcmc.out)
#'
#' }
#'
#' @seealso \code{\link[coda]{gelman.diag}}, \code{\link[coda]{effectiveSize}}
#' @export
#'
summary.bsim <- function(object, ...){
  p <- ncol(object$input$origdata$x)
  nchain <- object$input$samplingOptions$nchain
  namesVar <- colnames(object$input$origdata$x)
  printIndex <- paste0("index_", namesVar)
  namesIndex <- paste0("index[", 1:p, "]")
  namesSigma <- "sigma2"
  namesExtra <- character(0)
  temp <- object$input$samplingOptions$monitors
  if (!is.null(temp)){
    varList <- getVarMonitor(object, type = "list")
    for (v in temp){
      if (!grepl("\\[", v)){ # do not include "["
        # list up all name relevant to monitors variable
        for (j in varList){
          base <- sub("\\[.*$", "", j)
          if (v == base){
            namesExtra <- c(namesExtra, j)
          }
        }
      } else{
        namesExtra <- c(namesExtra, v)
      }
    }
  }


  # Names
  if (object$modelName %in% c("bsSphere", "bsSpike", "gpSpike")){ # Variable selection methods
    namesIndicator <- paste0("nu[", 1:p, "]")
    namesIndex_Var <- paste0("nu_",namesVar)

    if (!is.null(temp)){
      # param_names <- c(namesSigma, namesIndicator, namesIndex, namesExtra)
      param_names_order <- c(namesIndicator,  namesIndex, namesSigma, namesExtra)
      param_names_list <- c(namesIndex_Var, printIndex, namesSigma, namesExtra)
    } else{
      # param_names <- c(namesSigma, namesIndicator, namesIndex)
      param_names_order <- c(namesIndicator,  namesIndex, namesSigma)
      param_names_list <- c(namesIndex_Var, printIndex, namesSigma)
    }
  } else{
    if (!is.null(temp)){
      param_names_order <- c(namesIndex, namesSigma, namesExtra)
      param_names_list <- c(printIndex, namesSigma, namesExtra)
    } else{
      param_names_order <- c(namesIndex, namesSigma)
      param_names_list <- c(printIndex, namesSigma)
    }
  }

  calc_stat <- function(nchain, sampObj, namesIndex,
                       param_names_order, param_names_list){
    summarize_one_chain <- function(sampMat) {
      # sampMat: iter × param (matrix)
      # index standardized
      realTheta <- t(apply(sampMat[, namesIndex, drop = FALSE], 1, function(x) {
        x / sqrt(sum(x^2))
      }))
      colnames(realTheta) <- namesIndex

      # Rest of the parameters
      restNames <- setdiff(param_names_order, namesIndex)
      subTemp   <- cbind(sampMat[, restNames, drop = FALSE], realTheta)
      subTemp   <- subTemp[, param_names_order, drop = FALSE]

      # Compute mean, median, sd, CI
      stat_mat <- apply(subTemp, 2, function(x) {
        c(
          mean           = mean(x, na.rm = TRUE),
          median         = median(x, na.rm = TRUE),
          Standard.error = sd(x, na.rm = TRUE),
          LB             = quantile(x, probs = 0.025, na.rm = TRUE),
          UB             = quantile(x, probs = 0.975, na.rm = TRUE)
        )
      })

      stat_df <- as.data.frame(t(stat_mat))
      stat_df$ess <- coda::effectiveSize(subTemp)
      rownames(stat_df) <- param_names_list

      return(stat_df)
    }

    finalSummaryList <- list(
      chain     = vector("list", nchain),
      all.chain = NULL
    )

    if (nchain == 1){
      tmpResult <- summarize_one_chain(sampObj)
      finalSummaryList$chain <- tmpResult
      finalSummaryList$all.chain <- tmpResult
    } else{
      for (i in seq_len(nchain)) {
        sampTemp <- object$samples[[i]]
        tmpResult <- summarize_one_chain(sampTemp)
        finalSummaryList$chain[[i]] <- tmpResult
      }
      all_sub <- sampleBind(sampObj, nchain) # sampObj is MCMC
      tempALL <- summarize_one_chain(all_sub)
      gel     <- as.vector(coda::gelman.diag(sampObj[,param_names_order],
                                             multivariate = FALSE)$psrf[, 1])
      tempALL$gelman <- gel
      finalSummaryList$all.chain <- tempALL
    }
    return(finalSummaryList)

  }

  totalTable <- calc_stat(nchain, object$samples, namesIndex,
                          param_names_order, param_names_list)

  class(totalTable) <- "summary.bsim"
  return(totalTable)

}

#' @rdname summary.bsim
#' @export
print.summary.bsim <- function(x, digits = 3, ...){
  cat("Summary statistics of all chains\n")
  tempTable <- x$all.chain
  if ("gelman" %in% colnames(tempTable)){
    tempTable$gelman[is.nan(tempTable$gelman)] <- NA
    tempTable$gelman <- round(tempTable$gelman, digits = digits)
    tempTable$gelman <- as.character(tempTable$gelman)
    tempTable$gelman[is.na(tempTable$gelman)] <- "-"
  }
  print(tempTable, digits = digits)
}

#' Traceplot for `BayesSIM`
#'
#' @description
#' Provides traceplot for objects of `BayesSIM`.
#'
#' @param x A fitted object of \code{BayesSIM} or individual model.
#' @param ... Further arguments passed to \code{\link{plot}}.
#' @return
#' Traceplots for MCMC samples are displayed.
#' By default, the index vector and error variance are only included in the summary.
#' Extra monitor variables are included, if specified.
#' @examples
#' \donttest{
#' simdata2 <- data.frame(DATA1$X, y = DATA1$y)
#'
#' # 1. One tool version
#' fit_one <- BayesSIM(y ~ ., data = simdata2,
#'                     niter = 5000, nburnin = 1000, nchain = 1)
#'
#' # Check median index vector estimates with standard errors
#' coef(fit_one, method = "median", se = TRUE)
#'
#' # Fitted index values of median prediction
#' fitted(fit_one, type = "linpred", method = "median")
#'
#' # Residuals of median prediction
#' residuals(fit_one, method = "median")
#'
#' # Summary of the model
#' summary(fit_one)
#'
#' # Convergence diagnostics
#' nimTraceplot(fit_one)
#'
#' # Goodness of fit
#' GOF(fit_one)
#'
#' # Fitted plot
#' plot(fit_one)
#'
#' # Prediction with 95% credible interval at new data
#' newx <- data.frame(X1 = rnorm(10), X2 = rnorm(10), X3 = rnorm(10), X4 = rnorm(10))
#' pred <- predict(fit_one, newdata = newx, interval = "credible", level = 0.95)
#' plot(pred)
#'
#'
#' # 2. Split version
#' models <- BayesSIM_setup(y ~ ., data = simdata2)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, inits = initList,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#'
#' # "fit_split" becomes exactly the same as the class of "fit_one" object and apply generic functions.
#' fit_split <- as_bsim(models, mcmc.out)
#'
#' }
#' @export
nimTraceplot <- function(x, ...){
  Iteration <- 0; Value <- 0; num.chain <- 0; n <- 0
  p <- ncol(x$input$origdata$x)

  nchain <- x$input$samplingOptions$nchain
  varnames <- colnames(x$input$origdata$x)
  index_var <- paste0("index_", varnames)

  if (nchain > 1){
    samples <- NULL
    nchain <- length(x$samples)
    for (i in 1:nchain){
      n <- nrow(x$samples[[i]])
      temp <- cbind(chain.num = rep(i, n), x$samples[[i]])
      samples <- rbind(samples, temp)
      num.chain = samples[,1]
    }
  } else{
    n <- nrow(x$samples)
    samples <- x$samples
    num.chain = rep(1,n)
  }
  ##
  # index
  namesTheta <- paste0("index[", 1:p, "]")
  theta_samples <- samples[, namesTheta]

  # sigma2
  namesSigma <- "sigma2"
  sigma2_samples <- samples[, namesSigma]
  sigma2_df <- data.frame(sigma2 = sigma2_samples)

  # monitor variables
  temp <- x$input$samplingOptions$monitors
  if (!is.null(temp)){
    namesExtra <- character(0)
    temp <- x$input$samplingOptions$monitors
    varList <- getVarMonitor(x, type = "list")

     for (v in temp){
       if (!grepl("\\[", v)){ # do not include "["
         # list up all name relevant to monitors variable
         for (j in varList){
           base <- sub("\\[.*$", "", j)
           if (v == base){
             namesExtra <- c(namesExtra, j)
           }
           }
         } else{
           namesExtra <- c(namesExtra, v)
           }
      }

    extra_samples <- samples[, namesExtra]
  }

  if (x$modelName %in% c("bsSphere", "bsSpike", "gpSpike")){
    namesNu <- paste0("nu[", 1:p, "]")
    varNu <- paste0("nu_", varnames)
    nu_samples <- samples[, namesNu]
    if (!is.null(temp)){
      df <- data.frame(Iteration = rep(1:n, nchain),
                       num.chain =  num.chain,
                       nu_samples,
                       theta_samples,
                       sigma2_df, extra_samples)
      totalName <- c(varNu, index_var, namesSigma, namesExtra)
    } else{
      df <- data.frame(Iteration = rep(1:n, nchain),
                       num.chain =  num.chain,
                       nu_samples,
                       theta_samples,
                       sigma2_df)
      totalName <- c(varNu, index_var, namesSigma)
    }
  } else{
    if (!is.null(temp)){
      df <- data.frame(Iteration = rep(1:n, nchain),
                       num.chain =  num.chain,
                       theta_samples,
                       sigma2_df, extra_samples)
      totalName <- c(index_var, namesSigma, namesExtra)
    } else{
      df <- data.frame(Iteration = rep(1:n, nchain),
                       num.chain =  num.chain,
                       theta_samples,
                       sigma2_df)
      totalName <- c(index_var, namesSigma)
    }
  }

  colnames(df)[-c(1,2)] <- totalName

  if (nchain == 1){

    df_long <- df %>% dplyr::select(-num.chain) %>%
      pivot_longer(
        cols = -Iteration,
        names_to = "Parameter",
        values_to = "Value"
      )
    df_long$Parameter <- factor(df_long$Parameter,
                                levels = totalName)
    ggplot2::ggplot(df_long, aes(x = Iteration, y = Value)) +
      geom_line() +
      facet_wrap(~Parameter, scales = "free_y") +
      theme_bw() +
      labs(title = "Traceplots for parameters",
           x = "Iteration", y = "Value")

  } else{
    df_long <- df %>%
      pivot_longer(
        cols = all_of(totalName),
        names_to = "Parameter",
        values_to = "Value"
      )
    df_long$Parameter <- factor(df_long$Parameter,
                                levels = totalName)
    df_long$num.chain <- factor(df_long$num.chain,
                                levels = 1:nchain)
    # traceplot
    ggplot2::ggplot(df_long, aes(x = Iteration, y = Value, col = num.chain)) +
      geom_line() +
      facet_wrap(~Parameter, scales = "free_y") +
      theme_bw() +
      labs(title = "Traceplots for parameters",
           x = "Iteration", y = "Value")
  }

}

#' Plot Method for `BayesSIM`
#'
#' Produce diagnostic plots for a fitted Bayesian single-index model.
#'
#' @param x A fitted object of \code{BayesSIM} or individual model.
#' @param method Character string specifying the summary used for the
#'   posterior fitted values. Options are `"mean"` or `"median"`.
#'   Default is `"mean"`.
#' @param interval A logical value indicating whether a credible interval is included in the fitted plot. Default is `TRUE`.
#' @param alpha Numeric value between 0 and 1 specifying the credible level.
#'   By default, `alpha = 0.95` produces a 95% credible interval.
#' @param ... Additional arguments passed to underlying plotting functions.
#'
#'
#' @details
#' The function internally calls [predict()] on the fitted model object
#' to obtain posterior summaries of \eqn{\hat{y}}. Predicted value of \eqn{y} is \eqn{\hat{f}(X'\hat{\theta})}.
#'
#' - If `interval = TRUE`, the function requests posterior credible intervals
#'   and overlays them on the fitted curve.
#' - If `interval = FALSE`, only the posterior mean or median curve is drawn.
#' @return The output consists of two plots:
#' \enumerate{
#'   \item **Observed vs Predicted plot**: a diagnostic scatter plot comparing
#'   actual outcomes with posterior fitted values to visually assess model fit.
#'   \item **Fitted curve plot**: posterior \eqn{\hat{y}} as a function of the
#'   estimated single index, optionally with \eqn{100 \times \alpha} % credible intervals.
#' }
#'
#' @seealso [predict.bsim()], [summary.bsim()]
#'
#' @examples
#' \donttest{
#' simdata2 <- data.frame(DATA1$X, y = DATA1$y)
#'
#' # 1. One tool version
#' fit_one <- BayesSIM(y ~ ., data = simdata2,
#'                     niter = 5000, nburnin = 1000, nchain = 1)
#'
#' # Check median index vector estimates with standard errors
#' coef(fit_one, method = "median", se = TRUE)
#'
#' # Fitted index values of median prediction
#' fitted(fit_one, type = "linpred", method = "median")
#'
#' # Residuals of median prediction
#' residuals(fit_one, method = "median")
#'
#' # Summary of the model
#' summary(fit_one)
#'
#' # Convergence diagnostics
#' nimTraceplot(fit_one)
#'
#' # Goodness of fit
#' GOF(fit_one)
#'
#' # Fitted plot
#' plot(fit_one)
#'
#' # Prediction with 95% credible interval at new data
#' newx <- data.frame(X1 = rnorm(10), X2 = rnorm(10), X3 = rnorm(10), X4 = rnorm(10))
#' pred <- predict(fit_one, newdata = newx, interval = "credible", level = 0.95)
#' plot(pred)
#'
#'
#' # 2. Split version
#' models <- BayesSIM_setup(y ~ ., data = simdata2)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, inits = initList,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#'
#' # "fit_split" becomes exactly the same as the class of "fit_one" object and apply generic functions.
#' fit_split <- as_bsim(models, mcmc.out)
#'
#' }
#' @method plot bsim
#' @name plot
#' @export
plot.bsim <- function(x, method = c("mean", "median"),
                            interval = TRUE, alpha = 0.95, ...){

  method <- match.arg(method)

  if (!is.logical(interval)){
    stop("interval should be logical value.")
  }

  if (interval == TRUE){
    pred <- predict.bsim(x, method = method,
                    interval = "credible",
                    alpha = alpha)
  } else{
    pred <- predict.bsim(x, method = method,
                    alpha = alpha)
  }

  # summaryResult <- list(fitted = NULL, truey = NULL,
  #                       idxValue = NULL, level = level)
  # summaryResult$truey <- as.vector(x$input$origdata$y)
  # summaryResult$idxValue <- fitted(x, type = "linpred", method = method)
  # pred <- fitted(x, method = method)
  #
  # if (interval == TRUE){
  #   fitted <- data.frame(pred = pred, )
  # }

  # object: fitting object

  plot_bsim_fitted(pred, interval)
}


#' Prediction Method for `BayesSIM`
#'
#' Generate predictions from a fitted Bayesian single-index model.
#'
#' @param object A fitted object of \code{BayesSIM} or individual model.
#' @param newdata Optional data frame for which predictions should be made.
#'   If `NULL`, predictions are returned for the training data.
#' @param se.fit A logical value indicating whether standard errors are required.
#'   Default is `FALSE`.
#' @param type Character string specifying the type of prediction. `"response"`
#'   is the default.
#'   \describe{
#'     \item{`"response"`}{
#'       Posterior predictive summaries of the response \eqn{Y}.
#'       This corresponds to draws from the posterior predictive distribution
#'       \eqn{Y^{(m)} \sim N(f(X'\theta^{(m)}), \sigma^{2(m)})}
#'       and therefore incorporates both the uncertainty in the link function
#'       and the variability of the error term for each \eqn{m^{th}} MCMC sample.}
#'     \item{`"latent"`}{
#'       Posterior summaries of the latent mean structure
#'       \eqn{E(Y \mid X) = f^{(m)}(t^{(m)})}, where \eqn{t^{(m)} = X'\theta^{(m)}}.
#'       Unlike `"response"`, it excludes the noise term and calculated by
#'       \eqn{f^{(m)}(X'\theta^{(m)})} for each \eqn{m^{th}} MCMC sample of
#'       \eqn{\theta}.}
#'     \item{`"index"`}{
#'       Posterior summaries of the single index \eqn{t^{(m)} = X'\theta^{(m)}}.}
#'   }
#'
#' @param method Character string determining the posterior summary used for
#'   point predictions. Options are `"mean"` or `"median"`. Default is `"mean"`.
#'
#' @param interval Character string indicating whether a credible interval should be returned. Default is `"none"`.
#'   \describe{
#'     \item{`"none"`}{Return only point predictions.}
#'     \item{`"credible"`}{Return a \eqn{100 \times \text{level}\%} credible interval.}
#'   }
#'
#' @param level Numeric value between 0 and 1 specifying the credible level.
#'   `level = 0.95` yields a 95% credible interval. Default is `0.95`.
#'
#' @param ... Additional arguments.
#'
#' @details
#' This method extracts MCMC posterior samples stored in a BayesSIM object
#' and computes posterior summaries of:
#' \itemize{
#'   \item the **posterior predictive response** \eqn{Y \mid X} (type `"response"`),
#'   \item the **latent link function** evaluation \eqn{E(Y \mid X) = f(X'\theta)} (type `"latent"`), or
#'   \item the **single index** \eqn{X'\theta} (type `"index"`).
#' }
#'
#' The key distinction is that `"response"` incorporates the posterior
#' variability of the error term \eqn{\epsilon}, whereas `"latent"` represents
#' the noiseless conditional expectation \eqn{E(Y \mid X)} computed directly
#' from the link function and the posterior draws of \eqn{\theta}.
#'
#' When `interval = "credible"`, the returned object includes lower and upper
#' credible bounds computed via posterior quantiles for the chosen prediction
#' scale.
#'
#' If `newdata` is supplied, predictions are evaluated at the new covariate
#' values by computing the corresponding posterior index \eqn{t = X'\theta}
#' and applying the link function.
#'
#' @return
#' A list containing at least the following components:
#' \describe{
#'   \item{`fitted`}{Posterior mean or median predictions. If \code{se.fit = TRUE} or
#'   \code{interval = "credible"}, standard error and lower, upper bounds of the credible interval is also included.}
#'   \item{`truey`}{True response value of test data. When true response value is not available, \code{NULL} is saved.}
#'   \item{`idxValue`}{Linear index value is saved where \eqn{\theta} is estimated by \code{method}.}
#'   \item{`level`}{Credible level.}
#' }
#'
#' @examples
#' \donttest{
#' simdata2 <- data.frame(DATA1$X, y = DATA1$y)
#'
#' # 1. One tool version
#' fit_one <- BayesSIM(y ~ ., data = simdata2,
#'                     niter = 5000, nburnin = 1000, nchain = 1)
#'
#' # Check median index vector estimates with standard errors
#' coef(fit_one, method = "median", se = TRUE)
#'
#' # Fitted index values of median prediction
#' fitted(fit_one, type = "linpred", method = "median")
#'
#' # Residuals of median prediction
#' residuals(fit_one, method = "median")
#'
#' # Summary of the model
#' summary(fit_one)
#'
#' # Convergence diagnostics
#' nimTraceplot(fit_one)
#'
#' # Goodness of fit
#' GOF(fit_one)
#'
#' # Fitted plot
#' plot(fit_one)
#'
#' # Prediction with 95% credible interval at new data
#' newx <- data.frame(X1 = rnorm(10), X2 = rnorm(10), X3 = rnorm(10), X4 = rnorm(10))
#' pred <- predict(fit_one, newdata = newx, interval = "credible", level = 0.95)
#' plot(pred)
#'
#'
#' # 2. Split version
#' models <- BayesSIM_setup(y ~ ., data = simdata2)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, inits = initList,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#'
#' # "fit_split" becomes exactly the same as the class of "fit_one" object and apply generic functions.
#' fit_split <- as_bsim(models, mcmc.out)
#'
#' }
#' @method predict bsim
#' @export
predict.bsim <- function(object, newdata = NULL,
                         se.fit = FALSE,
                         type = c("response", "latent", "index"),
                         method = c("mean", "median"),
                         interval = c("none", "credible"),
                         level = 0.95, ...){
  start1 <- Sys.time()
  truey <- 0; pred <- 0; LB <- 0; UB <- 0

  # environment
  envobj <- ls(envir=.GlobalEnv)
  on.exit(rm(list=ls(envir=.GlobalEnv)[which(!ls(envir=.GlobalEnv)%in%envobj)],envir=.GlobalEnv))

  # aa_bspline_ver3
  .fns <- c(
    # a_common
    "quickSortOrderIndexOnly", "nimOrder", "nimSort",
    "sampleQuantile_nim", "quantile_nimble",

    # aa_bspline_ver3
    "SplineState", "any_duplicated", "mat_wo_col1", "update_spline_df", "update_x_index",
    "update_knot_sequence", "get_basis_simple", "simplify_knots", "get_inside_x",
    "gen_default_internal_knots", "SplineBase1", "SplineBase2", "basis", "bsNimble", "bsBasis",

    # bsFisher
    "pred_bsplineFisher",

    # bsSphere
    "pred_bsplineSphere","expcovTest_gpSphere",
    "expcovTest_gpPolar","expcovnn_gpSpike",

    # extra
    "expcov_gpPolar", "expcov_gpSphere", "expcov_gpSpike",
    "transX_fisher", "transX_sp", "Xlinear", "alphaTheta",

    # distribution
    "dvMFnim", "dvMFnim", "dunitSphere", "runitSphere"

  )

  pkg <- "BayesSIM"
  ns <- asNamespace(pkg)
  list2env(mget(.fns, envir = ns, inherits = FALSE), envir = globalenv())

  suppressMessages(
    nimble::registerDistributions(list(
      dvMFnim = list(
        BUGSdist = "dvMFnim(theta)",
        types    = c("value = double(1)", "theta = double(1)"),
        discrete = FALSE
      )
    ), verbose = FALSE)
  )

  suppressMessages(
    nimble::registerDistributions(list(
      dKnotsSimple = list(
        BUGSdist = "dKnotsSimple(a, b, k, alpha)",
        types = c("value = double(1)", "a = double(0)", "b = double(0)", "k = double(0)",
                  "alpha = double(1)"),
        discrete = FALSE
      )
    ), verbose = FALSE)
  )

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


  type <- match.arg(type)
  method <- match.arg(method)
  interval <- match.arg(interval)

  # check - modify
  # if (is.null(object$samples)){
  #   stop("No MCMC samples available. Make sure 'sample = TRUE'")
  # }

  if (!is.character(method) || length(method) != 1){
    stop("'method' argument should be correctly specified.")
  }

  if (!is.logical(se.fit)){
    stop("'se.fit' argument should be correctly specified.")
  }

  if (!is.numeric(level) || length(level) != 1){
    stop("'level' argument should be correctly specified.")
  }

  if (!(type %in% c("index", "latent", "response"))){
    stop("type arugment is incorrect.")
  }

  # ------------------------------------------------------------------
  # Samples ready
  nchain <- object$input$samplingOptions$nchain
  samples <- sampleBind(object$samples, nchain)

  # Data formation for new data set
  if (!is.null(newdata)){ # test data
    # test new data
    if (!is.data.frame(newdata)){
      stop("newdata should be data.frame.")
    }

    formula <- as.character(object$input$formula)
    class.name <- formula[2]
    data.name <- strsplit(formula[3]," \\+ ")[[1]]
    # int.flag <- any(strsplit(formula[3]," \\* ")[[1]] == formula[3])

    if(data.name[1] == "."){
      allX <- colnames(newdata)[colnames(newdata) != class.name]
      tot.name <- c(class.name, allX)
    } else{
      tot.name <- c(class.name, data.name)
    }
    varset <- setdiff(tot.name, class.name)
    newdataX <- newdata[ ,varset]

    if (!(class.name %in% colnames(newdata))) { # new Y is not included
      predIncluded <- FALSE
      newdataY <- NULL
    } else{ # new Y is included
      predIncluded <- TRUE
      newdataY <- as.vector(newdata[,class.name])
    }

  } else{
    predIncluded <- TRUE
    newdataX <- object$input$origdata$x
    newdataY <- object$input$origdata$y
  }

  train_p <- ncol(newdataX)
  # N <- nrow(newdataX)
  N <- nrow(object$input$origdata$x)
  trainY <- as.vector(object$input$origdata$y)


  ## newdataY, newdataX exists (prediction starts) ----------------------------
  nsamp <- nrow(samples)
  namesIndex <- paste0("index[", 1:train_p, "]")
  namesXlin <- paste0("Xlin[", 1:N, "]") ##

  indexSample <- samples[, namesIndex]
  XlinSample <- samples[, namesXlin]
  sigma2_samples <- samples[, "sigma2"]
  newdataMat <- as.matrix(newdataX)

  indexValueMatrix <- indexSample %*% t(newdataX) # N(MCMCsamp) * nrow(X)

  # What - posterior samples
  if (type == "index"){ # newdataX %*% indexSample
    testPred <- indexValueMatrix # N(MCMCsamp) * nrow(X)
  }

  if (type %in% c("response", "latent")){
    # B-spline
    if (object$modelName == "bsSphere"){
      maxknots <- object$input$prior$link$knots$maxknots
      k <- samples[ ,"k"]
      knotsSample <- samples[ ,paste0("knots[", 1:maxknots, "]")]
      dfSample <- samples[ ,"numBasis"]
      degree <- object$input$prior$link$basis$degree
      namesBeta <- paste0("beta[", 1:maxknots, ", 1]")
      betaSample <- samples[,namesBeta]
      mina <- samples[,"a_alpha"]
      maxb <- samples[,"b_alpha"]


      message("Compile function..")
      suppressMessages(
        cpred_bsplineSphere <- compileNimble(pred_bsplineSphere)
      )

      if (type == "latent"){

        start2 <- Sys.time()
        message("Computing predicted value..")
        testPred <- cpred_bsplineSphere(newdataMat, indexSample,
                                        knotsSample, XlinSample,
                                        betaSample, sigma2_samples, k,
                                        dfSample,
                                        mina, maxb, nsamp,degree,
                                        prediction = 1)
        end2 <- Sys.time()


      } else if (type == "response"){
        start2 <- Sys.time()
        message("Computing predicted value..")
        testPred <- cpred_bsplineSphere(newdataMat, indexSample,
                                        knotsSample, XlinSample,
                                        betaSample, sigma2_samples, k,
                                        dfSample,
                                        mina, maxb, nsamp,degree,
                                        prediction = 2)
        end2 <- Sys.time()

      }

    } else if (object$modelName %in% c("bsFisher", "bsPolar", "bsSpike")){ # c("bsFisher", "bsPolar", "bsSpike")
      df <- object$input$prior$link$basis$df
      degree <- object$input$prior$link$basis$degree
      delta <- object$input$prior$link$basis$delta
      n_internal_knots <- df - degree
      prob_vec <- seq(0, 1, length.out = n_internal_knots+2)[2:(n_internal_knots+1)]
      knotsSample <- t(apply(XlinSample, 1, function(x) quantile(x, prob_vec))) # n_internal_knots * nsamp
      namesBeta <- paste0("beta[", 1:df, "]")
      betaSample <- samples[ ,namesBeta]
      end1 <- Sys.time()


      message("Compile function..")
      suppressMessages(
        cpred_bsplineFisher <- compileNimble(pred_bsplineFisher)
      )

      if (type == "latent"){
        start2 <- Sys.time()
        message("Computing predicted value..")
        testPred <- cpred_bsplineFisher(newdataMat, indexSample,
                                        knotsSample, XlinSample,
                                        betaSample, sigma2_samples,
                                        nsamp, df, degree, delta,
                                        prediction = 1)
        end2 <- Sys.time()
      } else if (type == "response"){
        start2 <- Sys.time()
        message("Computing predicted value..")
        testPred <- cpred_bsplineFisher(newdataMat, indexSample,
                                        knotsSample, XlinSample,
                                        betaSample, sigma2_samples,
                                        nsamp, df, degree, delta,
                                        prediction = 2)
        end2 <- Sys.time()
      }
    }

    # Gp
    if (object$modelName %in% c("gpSphere", "gpFisher")){
      # Hyper-parameters of expcov
      lengthSample <- samples[ ,"lengthscale"]
      ampSample <- samples[ ,"amp"]
      # newdataMat <- as.matrix(newdataX)

      message("Compile function..")
      suppressMessages(
        cpred_gpSphere <- compileNimble(pred_gpSphere)
      )

      if (type == "latent"){ # latent
        start2 <- Sys.time()
        message("Computing predicted value..")
        testPred <- cpred_gpSphere(newdataMat, nsamp, trainY, indexSample,
                                   XlinSample, sigma2_samples,
                                   lengthSample, ampSample,
                                   prediction = 1)
        end2 <- Sys.time()
      } else{# response
        start2 <- Sys.time()
        message("Computing predicted value..")
        testPred <- cpred_gpSphere(newdataMat, nsamp, trainY, indexSample,
                                   XlinSample, sigma2_samples,
                                   lengthSample, ampSample,
                                   prediction = 2)
        end2 <- Sys.time()
      }

    } else if (object$modelName == "gpPolar"){
      # Hyper-parameters of expcov
      kappaSample <- samples[ , "kappa"]
      newdataMat <- as.matrix(newdataX)

      message("Compile function..")
      suppressMessages(
        cpred_gpPolar <- compileNimble(pred_gpPolar)
      )

      if (type == "latent"){
        start2 <- Sys.time()
        message("Computing predicted value..")
        testPred <- cpred_gpPolar(newdataMat, nsamp, trainY, indexSample,
                                  XlinSample, sigma2_samples,
                                  kappaSample, prediction = 1)
        end2 <- Sys.time()
      } else{ # response
        start2 <- Sys.time()
        message("Computing predicted value..")
        testPred <- cpred_gpPolar(newdataMat, nsamp, trainY, indexSample,
                                  XlinSample, sigma2_samples,
                                  kappaSample, prediction = 2)
        end2 <- Sys.time()
      }

    } else if (object$modelName == "gpSpike"){
      namesIndex <- paste0("index[", 1:train_p, "]")
      indexstarSample <- samples[, namesIndex]
      invlambdaSample <- samples[,"invlambda"]
      newdataMat <- as.matrix(scale(newdataX))

      message("Compile function..")
      suppressMessages(
        cpred_gpSpike <- compileNimble(pred_gpSpike)
      )

      if (type == "latent"){
        # latent
        start2 <- Sys.time()
        message("Computing predicted value..")
        testPred <- cpred_gpSpike(newdataMat, nsamp, trainY, indexstarSample,
                                  XlinSample, sigma2_samples, invlambdaSample,
                                  prediction = 1)
        end2 <- Sys.time()
      } else{ # response
        start2 <- Sys.time()
        message("Computing predicted value..")
        testPred <- cpred_gpSpike(newdataMat, nsamp, trainY, indexstarSample,
                                  XlinSample, sigma2_samples, invlambdaSample,
                                  prediction = 2)
        end2 <- Sys.time()
      }

    }


  }
  # testpred: samples

  # Compute results
  alpha <- 1-level
  fittedResult <- data.frame(mean = apply(testPred, 2, mean),
                             median = apply(testPred, 2, median),
                             se = apply(testPred, 2, sd),
                             LB = apply(testPred, 2, quantile, probs = alpha/2, na.rm = TRUE),
                             UB = apply(testPred, 2, quantile, probs = 1-alpha/2, na.rm = TRUE))

  # interval, se
  Name <- c("mean", "median", "se", "LB", "UB")
  if (se.fit){
    Name <- c(method, "se")
  } else{
    Name <- method
  }

  if (interval == "credible"){
    Name <- c(Name, "LB", "UB")
  }
  end1 <- Sys.time()

  summaryResult <- list(fitted = NULL, truey = NULL,
                        idxValue = NULL, level = level)


  summaryResult$fitted <- fittedResult[ ,Name]
  summaryResult$truey <- newdataY # if predIncluded = FALSE, truey = NULL
  if (method == "mean"){
    summaryResult$idxValue <- apply(indexValueMatrix, 2, mean)
  } else if (method == "median"){
    summaryResult$idxValue <- apply(indexValueMatrix, 2, median)
  }

  class(summaryResult) <- "bsimPred"
  invisible(summaryResult)

}

# predict structure
#' @rdname plot
#' @export
plot.bsimPred <- function(x, ...){ # predict object

  # object: predict object
  if (sum(c("LB", "UB") %in% colnames(x$fitted)) == 2){
    interval = TRUE
  } else{
    interval = FALSE
  }

  plot_bsim_fitted(x, interval)

}



#' Retrieve Monitorable Variables
#'
#' @description
#' Functions for retrieving the variables that can be monitored.
#'
#' @param object A fitted object of \code{BayesSIM}, \code{BayesSIM_setup} or individual model.
#' @param type Options for variables. By default, `type = "name"` is used that it only prints the name of the node.
#' If you put name of the nodes, the MCMC outputs gave you all elements of the variable, in case of the vector.
#' If `type = "list"`, the dimension of the nodes are printed. If you put name and dimension of the nodes, only specific location of vector or matrix can be seen in `summary` or `nimTraceplot`.
#'
#' @details
#' The function returns a list of variables that can be used in \code{monitors2} in the \code{bayesSIM} function.
#' You can also refer to \code{\link{getModelDef}} to understand the model structure and designate necessary variables.
#' Stochastic nodes of the model are recommended to be monitored.
#'
#' @return A vector of variables that can be monitored in the model.
#' @seealso
#' \code{\link{getModelDef}}
#'
#' @examples
#' simdata2 <- data.frame(DATA1$X, y = DATA1$y)
#' models <- BayesSIM_setup(y ~ ., data = simdata2)
#'
#' # Get monitorable variables
#' getVarMonitor(models)
#' # Get the list of variables with dimension
#' getVarMonitor(models, type = "list")
#'
#' @export
getVarMonitor <- function(object, type = c("name", "list")){
  type <- match.arg(type)
  possibleList <- object$defModel$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  tempList <- expand_monitor(possibleList)
  if (type == "name"){
    return(tempList$baseName)
  } else if (type == "list"){
    return(tempList$out)
  } else{
    stop("Argument `type` is wrong.")
  }


}

#' Get Definition of the Model
#'
#' @description
#' Functions for identifying definition of the \pkg{nimble} model.
#'
#' @param object A fitted object of \code{BayesSIM} or individual model.
#'
#' @details
#' The function that contain Bayes SIM model structure in \pkg{nimble}.
#' This function is for advanced users.
#' There are several functions used in the model definition.
#' \itemize{
#' \item `transX_fisher`, `transX_sp`: Making B-spline basis.
#' \item `dvMFnim`: Distribution of von Mises-Fisher.
#' \item `dKnotsSimple`: Distribution of the free knots for `bsSphere`.
#' \item `dunitSphere`: Distribution of unit sphere.
#' \item `alphaTheta`: One-to-one polar transformation. Making index vector from individual angular vector `psi`.
#' \item `expcov_gpSphere`, `expcov_gpPolar`, `expcov_gpSpike`: Covariance kernel of each model.
#' `expcov_gpSphere` uses squared-exponential kernel, `expcov_gpPolar` uses OU process kernel, and `expcov_gpSpike` uses squared-exponential including its own parameter, \eqn{\lambda^{-1}}.
#' \item `Xlinear`: Making linear combination with index vector.
#' }
#'
#' @return BUGS code of the model definition.
#'
#'
#' @seealso
#' \code{\link{getVarMonitor}}
#'
#' @examples
#' simdata2 <- data.frame(DATA1$X, y = DATA1$y)
#' models <- bsFisher_setup(y ~ ., data = simdata2)
#' # Get model definition
#' getModelDef(models)
#'
#' @export
getModelDef <- function(object){
  return(object$defModel$modelDef$BUGScode)
}


#' Get Initial Value of the Model
#'
#' @description
#' Functions for getting list of initial values of the \pkg{nimble} model.
#'
#' @param object A fitted object of \code{BayesSIM}, `BayesSIM_setup` or individual model.
#'
#' @details
#' The list of initial values are returned. This can be helpful to use when you use `BayesSIM_setup`.
#' You should be aware of that if you want to get more than 1 chain of MCMC samples, you should change `nchain` argument in `BayesSIM_setup`.
#' The output of initial values would be different, depending on the number of chain.
#'
#' You can apply `BayesSIM` object when you want to check the list of initial values.
#'
#' @return BUGS code of the model definition.
#'
#'
#' @examples
#' \donttest{
#' simdata2 <- data.frame(DATA1$X, y = DATA1$y)
#'
#'
#' # Split version
#' models <- BayesSIM_setup(y ~ ., data = simdata2)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, inits = initList,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#'
#' # "fit_split" becomes exactly the same as the class of "fit_one" object and apply generic functions.
#' fit_split <- as_bsim(models, mcmc.out)
#'
#' }
#'
#'
#' @export
getInit <- function(object){

  return(object$input$init)
}


#' Compile a 'nimble' Model and Its MCMC
#'
#' @description
#' Compiles a nimble model object and a corresponding (uncompiled) MCMC
#' algorithm and returns the compiled pair.
#'
#' @param object Class `BayesSIM_setup` object
#'
#' @details
#' The function first compiles `nimble` model object, then compiles `nimble` sampler.
#' Both compiled model and compiled MCMC samplers are returned as a list.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{model}}{the compiled NIMBLE model (external pointer object).}
#'   \item{\code{mcmc}}{the compiled MCMC function/algorithm bound to the model.}
#' }
#'
#' @seealso
#' \code{\link[nimble]{nimbleModel}},
#' \code{\link[nimble]{configureMCMC}},
#' \code{\link[nimble]{buildMCMC}},
#' \code{\link[nimble]{compileNimble}},
#' \code{\link[nimble]{runMCMC}}
#'
#' @examples
#' \donttest{
#' simdata2 <- data.frame(DATA1$X, y = DATA1$y)
#'
#' # 1. One tool version
#' fit_one <- BayesSIM(y ~ ., data = simdata2,
#'                     niter = 5000, nburnin = 1000, nchain = 1)
#'
#'
#' # 2. Split version
#' models <- BayesSIM_setup(y ~ ., data = simdata2)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, inits = initList,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#'
#' # "fit_split" becomes exactly the same as the class of "fit_one" object and apply generic functions.
#' fit_split <- as_bsim(models, mcmc.out)
#'
#' }
#'
#' @export
compileModelAndMCMC <- function(object) {
  model <- object$defModel; mcmc <- object$defSampler
  # global assign
  envobj <- ls(envir=.GlobalEnv)
  on.exit(rm(list=ls(envir=.GlobalEnv)[which(!ls(envir=.GlobalEnv)%in%envobj)],envir=.GlobalEnv))


  # aa_bspline_ver3
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

    # bsSphere
    "transX_sp","estBetaInit_sp","computeSig","ComputeS","logdet_nim",
    "postll_bspline_sphere","postll_knots","betaFunction","computeA1_nu1",
    "computeA1_nu0","computeA2_add","computeA2_delete","rnormTrun","newKnots",
    "prop_add","prop_delete","rtruncnorm","pred_bsplineSphere",
    "nuSampler_bspline_sphere","indexSampler_bspline_sphere","knotsSampler_bspline_sphere",
    "betaSampler_bspline_sphere","sigma2Sampler_bspline_sphere",

    # gpPolar
    "alphaTheta","Xlinear","invcov","expcov_gpPolar","expcovTest_gpPolar",
    "obj_btt_theta","thetaPrior","pred_gpPolar","gibbsSampler_sigma2",
    "gibbsSampler_kappa","MH_thetaeta",

    # gpspike
    "expcov_gpSpike","expcovnn_gpSpike","computeA1","computeA2",
    "llFunLambda","llFunThetaV1","llFunThetaV2","transitionTheta",
    "multiplyMatrixByConstant","pred_gpSpike",
    "SamplingLambda_gp_spike","SamplingThetaV_gp_spike",
    "gibbsSigma2_gp_spike","gibbsGam_gp_spike",

    # gpSphere
    "expcov_gpSphere","expcovTest_gpSphere","conBeta","obj_btt",
    "obj_btt_EB","pred_gpSphere","indexSampler_gpSphere","sigma2Sampler_gpSphere",
    "optSampler",

    # utils
    "pred_fitted",
    # distribution
    "dvMFnim", "rvMFnim", "dunitSphere", "runitSphere","dKnotsSimple", "rKnotsSimple"

  )

  pkg <- "BayesSIM"
  ns <- asNamespace(pkg)
  list2env(mget(.fns, envir = ns, inherits = FALSE), envir = globalenv())

  suppressMessages(
    nimble::registerDistributions(list(
      dvMFnim = list(
        BUGSdist = "dvMFnim(theta)",
        types    = c("value = double(1)", "theta = double(1)"),
        discrete = FALSE
      )
    ), verbose = FALSE)
  )

  suppressMessages(
    nimble::registerDistributions(list(
      dKnotsSimple = list(
        BUGSdist = "dKnotsSimple(a, b, k, alpha)",
        types = c("value = double(1)", "a = double(0)", "b = double(0)", "k = double(0)",
                  "alpha = double(1)"),
        discrete = FALSE
      )
    ), verbose = FALSE)
  )

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

  message("Compile model")
  suppressMessages(Cmodel <- nimble::compileNimble(model))
  message("Compile samplers")
  suppressMessages(Cmcmc  <- nimble::compileNimble(mcmc, project = model, resetFunctions = TRUE))
  return(list(model = Cmodel, sampler = Cmcmc))
}


#' Get nimbleModel and nimbleSampler Object from the Result of `compileModelAndMCMC`
#' @description
#' Return compiled nimble model object and a corresponding MCMC samplers.
#'
#'
#' @param object The result of `compileModelAndMCMC` function.
#'
#' @return
#' `get_model` returns compiled Nimble model object.
#' `get_sampler` returns compiled Nimble sampler object, directly using in `runMCMC` function.
#'
#' @seealso
#' \code{\link[nimble]{nimbleModel}},
#' \code{\link[nimble]{configureMCMC}},
#' \code{\link[nimble]{buildMCMC}},
#' \code{\link[nimble]{compileNimble}},
#' \code{\link[nimble]{runMCMC}}
#'
#' @examples
#' \donttest{
#' simdata2 <- data.frame(DATA1$X, y = DATA1$y)
#'
#' # 1. One tool version
#' fit_one <- BayesSIM(y ~ ., data = simdata2,
#'                     niter = 5000, nburnin = 1000, nchain = 1)
#'
#'
#' # 2. Split version
#' models <- BayesSIM_setup(y ~ ., data = simdata2)
#' Ccompile <- compileModelAndMCMC(models)
#' nimSampler <- get_sampler(Ccompile)
#' initList <- getInit(models)
#' mcmc.out <- runMCMC(nimSampler, niter = 5000, nburnin = 1000, thin = 1,
#'                     nchains = 1, setSeed = TRUE, inits = initList,
#'                     summary = TRUE, samplesAsCodaMCMC = TRUE)
#'
#' # "fit_split" becomes exactly the same as the class of "fit_one" object and apply generic functions.
#' fit_split <- as_bsim(models, mcmc.out)
#'
#' }
#' @name getFunction
#' @export
get_model <- function(object){
  return(object$model)
}

#' @rdname getFunction
#' @export
get_sampler <- function(object){
  return(object$sampler)
}
