#' @import methods
bsimSpline <- function(model, sampler, sampling, fitted, input, modelName) {
  structure(
    list(
      model = model,
      sampler = sampler,
      sampling = sampling,
      fitted = fitted,
      input = input,
      modelName = modelName
    ),
    class = "bsimSpline"
  )
}

#' Summarize bsimSpline Results
#'
#' @description
#' Provides a \code{summary} method for objects of class \code{bsimSpline}, corresponding to a
#' single-index model with a Gaussian process link function.
#'
#' @param object An object of class \code{bsimSpline}.
#' @param verbose Logical. If \code{TRUE}, the summary values for all chain is printed.
#' @param ... Further arguments passed to \code{\link{summary}}.
#'
#' @details
#' A \code{data.frame} of summary statistics for MCMC samples. Each row corresponds to
#' a model parameter, and columns report the statistics.
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
#' By default, only the index vector and error variance are included in the summary.
#' If a uniform sphere prior is used, the indicator vector (\code{nu}) is also summarized.
#' Note that the potential scale reduction factor for \code{nu} can be reported as
#' \code{NaN} or \code{Inf}, since the indicator rarely changes during the MCMC run.
#'
#' @method summary bsimSpline
#'
#' @seealso \code{\link[coda]{gelman.diag}}, \code{\link[coda]{effectiveSize}}
#' @export
summary.bsimSpline <- function(object, verbose = TRUE, ...){
  if (verbose){
    message("Summary statistics of all chains\n")
  }
  p <- ncol(object$input$data$x)
  finalSummaryList <- list(chain = NULL, all.chain = NULL)
  namesIndex <- paste0("index[", 1:p, "]")

  if (object$modelName %in% c("bsSphere", "bsSpike")){
    namesSigma <- "sigma2"
    namesIndicator <- paste0("nu[", 1:p, "]")
    param_names <- c(namesSigma, namesIndicator,  namesIndex)
    param_names_order <- c(namesIndicator,  namesIndex, namesSigma)

    if (is.list(object$sampling)){ # Multiple chain
      nchain <- length(object$sampling)
      summaryList <- vector("list", nchain)
      finalSummaryList$chain <- summaryList
      samplesFULL <- NULL
      for (i in 1:nchain){ # Each chain
        sampTemp <- object$sampling[[i]]
        thetastar_samp <- sampTemp[, namesIndex, drop = FALSE]
        realTheta <- t(apply(thetastar_samp, 1, function(x) x/sqrt(sum(x^2))))
        theta_summary <- apply(realTheta, 2, function(x) {
          c(
            mean = mean(x),
            median = median(x),
            sd   = sd(x),
            LB = quantile(x, probs = 0.025, na.rm = TRUE),
            UB = quantile(x, probs = 1-0.05/2, na.rm = TRUE)
          )
        })
        restNames <- setdiff(param_names, namesIndex)
        result <- apply(sampTemp[, restNames, drop = FALSE], 2, function(x) {
          c(
            mean = mean(x),
            median = median(x),
            sd   = sd(x),
            LB = quantile(x, probs = 0.025, na.rm = TRUE),
            UB = quantile(x, probs = 1-0.05/2, na.rm = TRUE)
          )
        })
        result <- as.data.frame(t(result))
        result <- rbind(result, t(theta_summary))
        result$ess <- coda::effectiveSize(sampTemp[,param_names])
        # summaryList[[i]] <- result[param_names_order, ]
        finalSummaryList$chain[[i]] <- result[param_names_order, ]
        samplesFULL <- rbind(samplesFULL, object$sampling[[i]])
      }

      # Total
      thetastar_samp <- samplesFULL[, namesIndex, drop = FALSE]
      realTheta <- t(apply(thetastar_samp, 1, function(x) x/sqrt(sum(x^2))))
      theta_summary <- apply(realTheta, 2, function(x) {
        c(
          mean = mean(x),
          median = median(x),
          sd   = sd(x),
          LB = quantile(x, probs = 0.025, na.rm = TRUE),
          UB = quantile(x, probs = 1-0.05/2, na.rm = TRUE)
        )
      })
      restNames <- setdiff(param_names, namesIndex)

      resultTotal <- apply(samplesFULL[, restNames, drop = FALSE], 2, function(x) {
        c(
          mean = mean(x),
          median = median(x),
          sd   = sd(x),
          LB = quantile(x, probs = 0.025, na.rm = TRUE),
          UB = quantile(x, probs = 1-0.05/2, na.rm = TRUE)
        )
      })
      resultTotal <- as.data.frame(t(resultTotal))
      resultTotal <- rbind(resultTotal, t(theta_summary))
      resultTotal$gelman <-  coda::gelman.diag(object$sampling[,param_names],
                                               multivariate = FALSE)$psrf[,1]
      resultTotal$ESS <- coda::effectiveSize(object$sampling[,param_names])
      # summaryList[[nchain+1]] <- resultTotal[param_names_order, ]
      finalSummaryList$all.chain <- resultTotal[param_names_order, ]
      names(finalSummaryList$chain) <- paste0("chain", 1:nchain)

      if (verbose){
        print(finalSummaryList$all.chain, digits = 2)
      }

      invisible(finalSummaryList)

    } else{
      samples <- object$sampling
      thetastar_samp <- samples[, namesIndex, drop = FALSE]
      realTheta <- t(apply(thetastar_samp, 1, function(x) x/sqrt(sum(x^2))))
      theta_summary <- apply(realTheta, 2, function(x) {
        c(
          mean = mean(x, na.rm = TRUE),
          median = median(x, na.rm = TRUE),
          sd   = sd(x, na.rm = TRUE),
          LB = quantile(x, probs = 0.025, na.rm = TRUE),
          UB = quantile(x, probs = 1-0.05/2, na.rm = TRUE)
        )
      })
      restNames <- setdiff(param_names, namesIndex)
      summary_df <- apply(samples[, restNames, drop = FALSE], 2, function(x) {
        c(
          mean = mean(x, na.rm = TRUE),
          median = median(x, na.rm = TRUE),
          sd   = sd(x, na.rm = TRUE),
          LB = quantile(x, probs = 0.025, na.rm = TRUE),
          UB = quantile(x, probs = 1-0.05/2, na.rm = TRUE)
        )
      })

      result <- as.data.frame(t(summary_df))
      result <- rbind(result, t(theta_summary))
      result$ESS <- coda::effectiveSize(object$sampling[,param_names])
      finalSummaryList$all.chain <- result[param_names_order,]
      if (verbose){
        print(finalSummaryList$all.chain, digits = 2)
      }
      invisible(finalSummaryList)
    }


  } else{ # others
    namesIndex <- paste0("index[", 1:p, "]")
    namesSigma <- "sigma2"
    param_names <- c(namesIndex, namesSigma)

    if (is.list(object$sampling)){ # Multiple chain
      nchain <- length(object$sampling)
      summaryList <- vector("list", nchain)
      finalSummaryList$chain <- summaryList
      samplesFULL <- NULL
      for (i in 1:nchain){ # Single chain
        sampTemp <- object$sampling[[i]]
        result <- apply(sampTemp[, param_names, drop = FALSE], 2, function(x) {
          c(
            mean = mean(x),
            median = median(x),
            sd   = sd(x),
            LB = quantile(x, probs = 0.025, na.rm = TRUE),
            UB = quantile(x, probs = 1-0.05/2, na.rm = TRUE)
          )
        })
        result <- as.data.frame(t(result))
        # result <- rbind(result, restNames)
        result$ESS <- coda::effectiveSize(sampTemp[,param_names])
        # summaryList[[i]] <- result
        finalSummaryList$chain[[i]] <- result
        samplesFULL <- rbind(samplesFULL, object$sampling[[i]])
      }

      # Total
      resultTotal <- apply(samplesFULL[, param_names, drop = FALSE], 2, function(x) {
        c(
          mean = mean(x),
          median = median(x),
          sd   = sd(x),
          LB = quantile(x, probs = 0.025, na.rm = TRUE),
          UB = quantile(x, probs = 1-0.05/2, na.rm = TRUE)
        )
      })
      resultTotal <- as.data.frame(t(resultTotal))
      resultTotal$gelman <-  coda::gelman.diag(object$sampling[,param_names],
                                               multivariate = FALSE)$psrf[,1]
      resultTotal$ESS <- coda::effectiveSize(object$sampling[,param_names])
      # summaryList[[nchain+1]] <- resultTotal
      finalSummaryList$all.chain <- resultTotal

      names(finalSummaryList$chain) <- paste0("chain", 1:nchain)
      if (verbose){
        print(finalSummaryList$all.chain, digits = 2)
      }
      invisible(finalSummaryList)

    } else{
      samples <- object$sampling
      summary_df <- apply(samples[, param_names, drop = FALSE], 2, function(x) {
        c(
          mean = mean(x),
          median = median(x),
          sd   = sd(x),
          LB = quantile(x, probs = 0.025, na.rm = TRUE),
          UB = quantile(x, probs = 1-0.05/2, na.rm = TRUE)
        )
      })
      result <- as.data.frame(t(summary_df))
      result$ESS <- coda::effectiveSize(object$sampling[,param_names])
      finalSummaryList$all.chain <- result
      if (verbose){
        print(finalSummaryList$all.chain, digits = 2)
      }
      invisible(finalSummaryList)
    }
  }



}



#' Traceplot for bsimSpline Results
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer all_of
#' @import ggplot2
#' @description
#' Provides traceplot for objects of class \code{bsimSpline}, corresponding to a
#' single-index model with a Gaussian process link function.
#'
#' @param x An object of class \code{bsimSpline}.
#' @param ... Further arguments passed to \code{\link{plot}}.
#' @return
#' Traceplots for MCMC samples are displayed.
#' By default, only the index vector and error variance are included in the summary.
#'
#' @method plot bsimSpline
#' @export
plot.bsimSpline <- function(x, ...){
  Iteration <- 0; Value <- 0; num.chain <- 0
  p <- ncol(x$input$data$x)
  if (is.list(x$sampling)){
    samples <- NULL
    nchain <- length(x$sampling)
    for (i in 1:nchain){
      n <- nrow(x$sampling[[i]])
      temp <- cbind(chain.num = rep(i, n), x$sampling[[i]])
      samples <- rbind(samples, temp)
    }
    ##
    namesTheta <- paste0("index[", 1:p, "]")
    theta_samples <- samples[, namesTheta]

    namesSigma <- "sigma2"
    sigma2_samples <- samples[, namesSigma]
    sigma2_df <- data.frame(sigma2 = sigma2_samples)

    df <- data.frame(Iteration = rep(1:n, nchain),
                     num.chain =  samples[,1],
                     theta_samples,
                     sigma2_df)
    colnames(df)[-c(1,2)] <- c(namesTheta, namesSigma)

    df_long <- df %>%
      pivot_longer(
        cols = all_of(c(namesTheta, namesSigma)),
        names_to = "Parameter",
        values_to = "Value"
      )

    df_long$Parameter <- factor(df_long$Parameter, levels = c(namesTheta, namesSigma))
    df_long$num.chain <- factor(df_long$num.chain, levels = 1:nchain)

    # traceplot
    ggplot2::ggplot(df_long, aes(x = Iteration, y = Value, col = num.chain)) +
      geom_line() +
      facet_wrap(~ Parameter, scales = "free_y") +
      theme_bw() +
      labs(title = "Traceplots for parameters",
           x = "Iteration", y = "Value")

  } else{ # one chain only
    samples <- x$sampling
    namesTheta <- paste0("index[", 1:p, "]")
    theta_samples <- samples[, namesTheta]
    namesSigma <- "sigma2"

    sigma2_samples <- samples[, namesSigma]
    sigma2_df <- data.frame(sigma2 = sigma2_samples)

    df <- data.frame(Iteration = 1:nrow(theta_samples),
                     theta_samples,
                     sigma2_df)
    colnames(df)[-1] <- c(namesTheta, namesSigma)

    df_long <- df %>%
      pivot_longer(
        cols = -Iteration,
        names_to = "Parameter",
        values_to = "Value"
      )

    df_long$Parameter <- factor(df_long$Parameter, levels = c(namesTheta, namesSigma))

    # traceplot
    ggplot2::ggplot(df_long, aes(x = Iteration, y = Value)) +
      geom_line() +
      facet_wrap(~ Parameter, scales = "free_y") +
      theme_bw() +
      labs(title = "Traceplots for parameters",
           x = "Iteration", y = "Value")

  }

}

#########################################################

#' Predict method for \code{bsimSpline} objects
#' @import patchwork
#' @importFrom magrittr %>%
#' @import ggplot2
#'
#' @description
#' Computes posterior predictive summaries from a fitted single–index B-spline model
#' (class \code{bsimSpline}). Optionally returns a diagnostic plot and RMSE when
#' observed responses are provided.
#'
#' @param object An object of class \code{bsimSpline} created with MCMC sampling
#'   enabled.
#' @param newdata A \code{data.frame} of predictors for prediction. The columns
#'   must match the columns of train data. If a column named \code{"y"} is present, it is treated as the observed response
#'   and is used to compute RMSE and to draw fitted plots. If \code{NULL},
#'   fitted values with train data are summarized.
#' @param method A character scalar selecting the point estimator returned in the
#'   summary: \code{"mean"} (posterior mean) or \code{"median"} (posterior median).
#' @param se.fit Logical; if \code{TRUE}, include posterior standard deviation and
#'   interval bounds in the returned summary and shaded area for the interval
#'   in the fitted curve plot.
#' @param level Numeric scalar in \eqn{[0, 1]}; nominal coverage for intervals. By default, \code{level = 0.95}.
#' @param ... Further arguments passed to \code{\link{predict}}.
#'
#' @details
#' The function first gathers posterior draws and computes predictive draws for each row of
#' \code{newdata} (or for the training data when \code{newdata = NULL}). From
#' these draws it forms.:
#' \itemize{
#'   \item posterior mean and median,
#'   \item empirical standard deviation (\code{sd}),
#'   \item lower/upper quantiles for a \code{level} interval (\code{LB}, \code{UB}).
#' }
#' If \code{newdata} supplies a response column \code{y}, the function also
#' computes the mean squared prediction error (reported as \code{rmse}) and
#' creates two ggplot objects: (i) True vs Predicted with a 45° line and
#' (ii) Index value vs Response showing the fitted curve. If
#' \code{se.fit = TRUE}, a shaded area visualizes the \code{level} interval.
#'
#' The index used in the fitted-curve plot is formed by projecting predictors
#' onto the estimated index vector extracted from \code{summary(object)}.
#'
#' @return A \code{list} with components:
#' \itemize{
#'   \item \code{summary}: a \code{data.frame} with columns matching
#'         \code{method}, either \code{mean} or \code{median}; if \code{se.fit = TRUE}
#'         the columns \code{sd}, \code{LB}, \code{UB} are included.
#'   \item \code{plot}: a combined \code{ggplot} object with the fitted plots when an observed
#'         \code{y} is available, otherwise \code{NULL}.
#'   \item \code{rmse}: mean squared prediction error when \code{y} is available, otherwise \code{NULL}.
#' }
#'
#'
#' @method predict bsimSpline
#' @seealso \code{\link{summary.bsimSpline}}
#' @export
predict.bsimSpline <- function(object, newdata = NULL,
                               method = "mean", se.fit = FALSE,
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
    "pred_bsplineSphere"

  )

  pkg <- "BayesSIM"
  ns <- asNamespace(pkg)
  list2env(mget(.fns, envir = ns, inherits = FALSE), envir = globalenv())



  # check
  if (is.null(object$sampling)){
    stop("No MCMC samples available. Make sure 'sample = TRUE'")
  }

  if (!is.character(method) || length(method) != 1){
    stop("'method' argument should be correctly specified.")
  }

  if (method != "mean" & method != "median"){
    stop("'method' argument should be mean or median.")
  }

  if (!is.logical(se.fit)){
    stop("'se.fit' argument should be correctly specified.")
  }

  if (!is.numeric(level) || length(level) != 1){
    stop("'level' argument should be correctly specified.")
  }

  N <- nrow(object$input$data$x)


  # no newdata
  if (is.null(newdata)){
    testPred <- object$fitted
    predIncluded <- TRUE
  } else{
    # newdata
    train_p <- ncol(object$input$data$x)
    y <- as.vector(object$input$data$y)


    # test new data
    if (!is.data.frame(newdata)){
      stop("newdata should be data.frame.")
    }


    if (!("y" %in% colnames(newdata))) {
      predIncluded <- FALSE
      newdataX <- newdata
    } else{
      predIncluded <- TRUE
      newdataX <- newdata[,setdiff(names(newdata), "y")]
      newdataY <- newdata$y
    }
    if (ncol(newdataX) != train_p){
      stop("Dimension of newdata is incorrect with train data.")
    }

    ## 0. samples ready
    if (is.list(object$sampling)){
      samples <- NULL
      nchain <- length(object$sampling)
      for (i in 1:nchain){
        samples <- rbind(samples, object$sampling[[i]])
      }
    } else{
      samples <- object$sampling
    }

    nsamp <- nrow(samples)
    namesIndex <- paste0("index[", 1:train_p, "]")
    namesXlin <- paste0("Xlin[", 1:N, "]")
    indexSample <- samples[, namesIndex]
    XlinSample <- samples[, namesXlin]
    sigma2_samples <- samples[, "sigma2"]
    newdataMat <- as.matrix(newdataX)
    # newdataMat <- rbind(object$input$data$x, temp)
    testPred <- matrix(0, nrow = nsamp, ncol = nrow(newdata))

    if (object$modelName %in% c("bsFisher", "bsPolar", "bsSpike")){
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


      start2 <- Sys.time()
      message("Computing predicted value..")
      testPred <- cpred_bsplineFisher(newdataMat, indexSample,
                                      knotsSample, XlinSample,
                                      betaSample, sigma2_samples,
                                      nsamp, df, degree, delta)
      end2 <- Sys.time()


    } else if (object$modelName == "bsSphere"){
      maxknots <- object$input$prior$link$knots$maxknots
      k <- samples[ ,"k"]
      knotsSample <- samples[ ,paste0("knots[", 1:maxknots, "]")]
      dfSample <- samples[ ,"numBasis"]
      degree <- object$input$prior$link$basis$degree
      namesBeta <- paste0("beta[", 1:maxknots, ", 1]")
      betaSample <- samples[,namesBeta]
      mina <- samples[,"a_alpha"]
      maxb <- samples[,"b_alpha"]
      end1 <- Sys.time()


      message("Compile function..")
      suppressMessages(
        cpred_bsplineSphere <- compileNimble(pred_bsplineSphere)
      )


      start2 <- Sys.time()
      message("Computing predicted value..")
      testPred <- cpred_bsplineSphere(newdataMat, indexSample,
                                      knotsSample, XlinSample,
                                      betaSample, sigma2_samples, k,
                                      dfSample,
                                      mina, maxb, nsamp,degree)
      end2 <- Sys.time()
    }

  }


  alpha <- 1-level
  # testPred <- testPred[-c(1:N),]
  fittedResult <- data.frame(mean = apply(testPred, 2, mean),
                             median = apply(testPred, 2, median),
                             sd = apply(testPred, 2, sd),
                             LB = apply(testPred, 2, quantile, probs = alpha/2, na.rm = TRUE),
                             UB = apply(testPred, 2, quantile, probs = 1-alpha/2, na.rm = TRUE))

  if (se.fit){
    summaryOuput <- fittedResult[,c(method, "sd", "LB", "UB")]
  } else{
    summaryOuput <- fittedResult[,c(method)]
  }

  if (!predIncluded){
    plots <- NULL
    rmse <- NULL
  } else{
    if (method == "mean"){
      predY <- fittedResult$mean
    } else if (method == "median"){
      predY <- fittedResult$median
    }

    if (is.list(object$sampling) & length(object$sampling) > 1){ # multiple chain
      TEMP <- summary(object, verbose = FALSE)$all.chain
      estidx <- TEMP[grepl("index", rownames(TEMP)),method]
    } else {
      TEMP <- summary(object, verbose = FALSE)$all.chain
      estidx <- TEMP[grepl("index", rownames(TEMP)),method]
    }

    ## projection
    if (is.null(newdata)){
      matX <- object$input$data$x
      idxValue <- matX %*% matrix(estidx, ncol = 1)
      TRUEY <- object$input$data$y
    } else if (predIncluded){
      idxValue <- newdataMat %*% matrix(estidx, ncol = 1)
      TRUEY <- newdataY
    }

    results <- data.frame(truey = TRUEY, pred = predY,
                          idxValue = idxValue[,1])

    rmse <- sqrt(mean((results$truey - results$pred)^2))

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

    # index vs. y
    if (!se.fit){
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
      results2 <- data.frame(truey = TRUEY, pred = predY,
                             idxValue = idxValue[,1], LB = fittedResult$LB,
                             UB = fittedResult$UB)
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



    plots <- A1 + A2

  }

  if (is.null(newdata)){
    time <- NULL
  } else{
    time <- difftime(end2, start2, units = "secs") +
      difftime(end1, start1, units = "secs")
  }

  print(summaryOuput)
  invisible(list(summary = summaryOuput,
                 plot = plots,
                 rmse = rmse,
                 time = time))


}

