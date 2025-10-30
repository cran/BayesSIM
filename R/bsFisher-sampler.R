# Sampler
## 1. theta random walk metropolis algorithm
indexSampler_bspline_fisher <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control){
    calcNodes <- model$getDependencies(target)
    # vMFnimsample <- vMFnim()
    Rpostll <- postll_bspline_fisher(model)
    propDens <- 1000
    Y <- model$Y
    rho <- gvcCV(model$Xmat, Y)
  },
  run = function(){
    currTheta <- model$index
    beforell <- Rpostll$run(rho)
    proposal <- rvMFnim(1, theta = currTheta * propDens)
    model$index <<- nimC(proposal)
    model$calculate(calcNodes)
    afterll <- Rpostll$run(rho)

    L <- (afterll - beforell)
    ratio <- min(c(1, exp(L)))
    u <- runif(1)

    # accept/reject
    if(u < ratio) {
      nimCopy(from = model, to = mvSaved, row = 1, nodes  = calcNodes, logProb = TRUE)
    }
    if(u > ratio){
      nimCopy(from = mvSaved, to = model, row = 1, nodes  = calcNodes, logProb = TRUE)
    }
  },
  methods = list(reset = function(){})
)

## 2. sampling beta - gibbs
betaSampler_bspline_fisher <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control){
    calcNodes <- model$getDependencies(target)
    Y <- model$Y
    N <- model$getConstants()$N
    rho <- gvcCV(model$Xmat, Y)
  },
  run = function(){
    Xmat <- model$Xmat
    temp <- nimDim(Xmat)[2]
    sigma2 <- model$sigma2
    cov <- chol(inverse(t(Xmat) %*% Xmat) *
                  nimMatrix(sigma2, ncol = temp, nrow = temp))
    betasamp <- rmnorm_chol(n = 1, mean = estBeta_fisher(Xmat, Y, rho),
                           cholesky = cov, prec_param = FALSE)
    model$beta <<- betasamp

    model$calculate(calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1,
            nodes  = calcNodes, logProb = TRUE)

  },
  methods = list(reset = function(){})
)
