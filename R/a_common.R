#' @importFrom stats dbeta dnorm median knots pnorm ppr qnorm qpois quantile
#' @importFrom stats rbeta rbinom rgamma rnorm rpois runif sd
#' @importFrom MASS mvrnorm

# Sorting algorithm: quick sort
quickSortOrderIndexOnly <- nimble::nimbleFunction(
  run = function(x = double(1), decreasing = logical(0, default = FALSE)) {
    n <- length(x)
    idx <- integer(n)
    for(i in 1:n) idx[i] <- i  # initial index

    # stack
    stackLeft <- integer(n)
    stackRight <- integer(n)
    top <- 1
    stackLeft[top] <- 1
    stackRight[top] <- n

    while(top > 0) {
      left <- stackLeft[top]
      right <- stackRight[top]
      top <- top - 1

      if(left < right) {
        pivotIdx <- floor((left + right) / 2)
        pivotVal <- x[idx[pivotIdx]]
        i <- left
        j <- right

        if(decreasing) {
          while(i <= j) {
            while(x[idx[i]] > pivotVal) i <- i + 1
            while(x[idx[j]] < pivotVal) j <- j - 1
            if(i <= j) {
              tmp <- idx[i]; idx[i] <- idx[j]; idx[j] <- tmp
              i <- i + 1; j <- j - 1
            }
          }
        } else {
          while(i <= j) {
            while(x[idx[i]] < pivotVal) i <- i + 1
            while(x[idx[j]] > pivotVal) j <- j - 1
            if(i <= j) {
              tmp <- idx[i]; idx[i] <- idx[j]; idx[j] <- tmp
              i <- i + 1; j <- j - 1
            }
          }
        }

        if(left < j) {
          top <- top + 1
          stackLeft[top] <- left
          stackRight[top] <- j
        }
        if(i < right) {
          top <- top + 1
          stackLeft[top] <- i
          stackRight[top] <- right
        }
      }
    }

    returnType(integer(1))
    return(idx)
  }
)

nimOrder <- nimble::nimbleFunction(
  run = function(x = double(1), decreasing = logical(0, default = FALSE)) {
    returnType(integer(1))
    return(quickSortOrderIndexOnly(x, decreasing))
  }
)


nimSort <- nimble::nimbleFunction(
  run = function(x = double(1), decreasing = logical(0, default = FALSE)) {
    idx <- nimOrder(x, decreasing)
    out <- numeric(length(x))
    for(i in 1:length(x)) {
      out[i] <- x[idx[i]]
    }
    returnType(double(1))
    return(out)
  }
)


# Quantile function
sampleQuantile_nim <- nimble::nimbleFunction(
  run = function(x   = double(1),
                 prob= double(1)) {
    returnType(double(1))
    n   <- nimDim(x)[1]
    sortx <- nimSort(x)
    np  <- nimDim(prob)[1]
    out <- numeric(np)
    for(p in 1:np) {
      k <- ceiling(prob[p] * (n + 1))
      if(k < 1)        k <- 1
      if(k > n)        k <- n
      out[p] <- sortx[k]
    }
    return(out)
  }
)



quantile_nimble <- nimble::nimbleFunction(
  run = function(x   = double(1),
                 prob= double(1)) {
    returnType(double(1))

    n   <- nimDim(x)[1]
    sortx <- nimSort(x)
    np  <- nimDim(prob)[1]
    out <- numeric(np)

    for (j in 1:np){

      h  <- (n - 1) * prob[j]
      i  <- floor(h)          # 0-based
      a  <- h - i             # alpha

      x_left  <- sortx[i + 1]
      x_right <- sortx[i + 2]

      if (prob[j] == 0){
        out[j] <- sortx[1]}
      else if (prob[j] == 1) {
        out[j] <- sortx[n]}
      else{
        out[j] <- (1 - a) * x_left + a * x_right
      }
    }

    return(out)
  }
)



#' Compile a NIMBLE model and its MCMC
#'
#' Compiles a NIMBLE model object and a corresponding (uncompiled) MCMC
#' algorithm and returns the compiled pair.
#'
#' @param fullmodel Class "bsimSpline" and "bsimGp" object with sampling = FALSE
#'
#' @details
#' The function first compiles \code{fullmodel$model} with
#' \code{nimble::compileNimble()}, then compiles \code{fullmodel$sampler} with
#' \code{nimble::compileNimble(project = model)}, where \code{model} is the
#' uncompiled model used to build the sampler. The compiled model and
#' compiled MCMC are returned as a list.
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
#' # Split version
#' models <- bsplineFisher(DATA1$X, DATA1$y, sampling = FALSE)
#' Ccompile <- compileModelAndMCMC(models)
#' mcmc.out <- runMCMC(Ccompile$mcmc, niter = 5000, nburnin = 1000, thin = 1,
#'                    nchains = 1, setSeed = TRUE, inits = models$input$initial,
#'                    summary = TRUE, samplesAsCodaMCMC = TRUE)
#'
#' }
#'
#' @export
compileModelAndMCMC <- function(fullmodel) {
  model <- fullmodel$model;mcmc <- fullmodel$sampler
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
    "pred_fitted"
  )

  pkg <- "BayesSIM"
  ns <- asNamespace(pkg)
  list2env(mget(.fns, envir = ns, inherits = FALSE), envir = globalenv())




  message("Compile model")
  suppressMessages(Cmodel <- nimble::compileNimble(model))
  message("Compile samplers")
  suppressMessages(Cmcmc  <- nimble::compileNimble(mcmc, project = model))
  return(list(model = Cmodel, mcmc = Cmcmc))
}
