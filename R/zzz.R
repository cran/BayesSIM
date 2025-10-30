# R/zzz.R
#' @keywords internal
"_PACKAGE"
.onLoad <- function(lib, pkg) {
  requireNamespace("nimble", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
}
.onAttach <- function(libname, pkgname) {
  # ns <- asNamespace(pkgname)
  # if (!identical(parent.env(.nimble_env), ns)) {
  #   parent.env(.nimble_env) <- ns
  # }
  packageStartupMessage("Loading BayesSIM Registering multiple variants of the following distributions:\n ",
                        "dvMFnim", ", KnotsSimple", ", unitSphere.\n")

  # Register the distributions explicitly for two reasons:
  # 1. Avoid message to user about automatic registrations upon first use in a nimbleModel
  # 2. Establish default len = 0 via reparameterization mechanism.
  suppressMessages({
    nimble::registerDistributions(list(
      dvMFnim = list(
        BUGSdist = "dvMFnim(theta)",
        types    = c("value = double(1)", "theta = double(1)"),
        discrete = FALSE
      )
    ), verbose = FALSE)

    nimble::registerDistributions(list(
      dKnotsSimple = list(
        BUGSdist = "dKnotsSimple(a, b, k, alpha)",
        types = c("value = double(1)", "a = double(0)", "b = double(0)", "k = double(0)",
                  "alpha = double(1)"),
        discrete = FALSE
      )
    ), verbose = FALSE)

    nimble::registerDistributions(list(
      dunitSphere = list(
        BUGSdist     = "dunitSphere(dim)",
        types        = c("value = double(1)",
                         "dim   = double(0)"),
        discrete     = FALSE
      )
    ), verbose = FALSE)

    # # aa_bspline_ver3
    # .fns <- c(
    #   # a_common
    #   "quickSortOrderIndexOnly", "nimOrder", "nimSort",
    #   "sampleQuantile_nim", "quantile_nimble",
    #
    #   # aa_bspline_ver3
    #   "SplineState", "any_duplicated", "mat_wo_col1", "update_spline_df", "update_x_index",
    #   "update_knot_sequence", "get_basis_simple", "simplify_knots", "get_inside_x",
    #   "gen_default_internal_knots", "SplineBase1", "SplineBase2", "basis", "bsNimble", "bsBasis",
    #
    #   # bsFisher
    #   "postll_bspline_fisher","nimNorm","rW","besselI_nimble","Stheta",
    #   "estBeta_fisher","gvcCV","transX_fisher",
    #   "pred_bsplineFisher","indexSampler_bspline_fisher","betaSampler_bspline_fisher",
    #
    #   # bsSphere
    #   "transX_sp","estBetaInit_sp","computeSig","ComputeS","logdet_nim",
    #   "postll_bspline_sphere","postll_knots","betaFunction","computeA1_nu1",
    #   "computeA1_nu0","computeA2_add","computeA2_delete","rnormTrun","newKnots",
    #   "prop_add","prop_delete","rtruncnorm","pred_bsplineSphere",
    #   "nuSampler_bspline_sphere","indexSampler_bspline_sphere","knotsSampler_bspline_sphere",
    #   "betaSampler_bspline_sphere","sigma2Sampler_bspline_sphere",
    #
    #   # gpPolar
    # "alphaTheta","Xlinear","invcov","expcov_gpPolar","expcovTest_gpPolar",
    # "obj_btt_theta","thetaPrior","pred_gpPolar","gibbsSampler_sigma2",
    # "gibbsSampler_kappa","MH_thetaeta",
    #
    # # gpspike
    # "expcov_gpSpike","expcovnn_gpSpike","computeA1","computeA2",
    # "llFunLambda","llFunThetaV1","llFunThetaV2","transitionTheta",
    # "multiplyMatrixByConstant","pred_gpSpike",
    # "SamplingLambda_gp_spike","SamplingThetaV_gp_spike",
    # "gibbsSigma2_gp_spike","gibbsGam_gp_spike",
    #
    # # gpSphere
    # "expcov_gpSphere","expcovTest_gpSphere","conBeta","obj_btt",
    # "obj_btt_EB","pred_gpSphere","indexSampler_gpSphere","sigma2Sampler_gpSphere",
    # "optSampler",
    #
    # # utils
    # "pred_fitted"
    #   )
    #
    # pkg <- "BayesSIM"
    # ns <- asNamespace(pkg)
    # list2env(mget(.fns, envir = ns, inherits = FALSE), envir = globalenv())
    #

    })}


