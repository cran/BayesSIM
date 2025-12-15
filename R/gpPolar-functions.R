
# Additional functions for gpPolar
initfunction_gpPolar <- function(x, y, kappa_init, sigma2_init, psi_init = NULL, grid.with = 0.1,
                         sig_a = 1, sig_b = 0.01, setSeed){

  if (setSeed != FALSE){
    set.seed(setSeed)
  }

  n <- dim(x)[1]
  p <- dim(x)[2]
  theta <- matrix(0,p-1,1)
  theta[,1] <- ifelse(is.null(psi_init), runif(p-1, 0, pi), psi_init)
  thetat <- theta
  parabeta <- betaalpha1(thetat)

  #Converting theta to alpha
  alphahat<-matrix(0,p,1) # p*1 matrix
  alphahat[,1] <- alphaTheta(theta)
  Xlin <-x %*% alphahat
  storeorder <- order(Xlin)

  #Ordering x and y for fitting OU process
  yobs <- y[order(Xlin)]
  xobs <- x[order(Xlin)]

  #Finding inverse of C
  Cinv <- invcov(Xlin, kappa_init)

  #generate from f
  sigma2 <-(sigma2_init^2) #Initial sigma^2
  Sigma <- solve(Cinv+diag(n)/sigma2)
  mu <- (Sigma %*% yobs)/sigma2
  f <- mvrnorm(1,mu=mu,Sigma=Sigma,tol=1e-8)

  #Generate from sigma2
  a <- sig_a;
  b <- sig_b;
  err <- sum((yobs-f)^2)/2
  sigma2 <-1/rgamma(1,shape=(a+n/2),rate=(b+err))

  thetaNum <- seq(0.01, pi-0.01, grid.with)
  thetas  <- do.call(expand.grid,
                     replicate(p-1, thetaNum, simplify = FALSE))
  thetas<-as.matrix(thetas)
  num.calc<-dim(thetas)[1]
  g <- f[order(storeorder)]

  PP1 <- apply(thetas,1,function(theta) optimized(g,x,kappa_init,theta,p))
  theta <- thetas[which.max(PP1),]

  thetat <- theta
  parabeta <- betaalpha1(thetat)

  #Converting theta to alpha
  alphahat<-matrix(0,p,1)
  alphahat[,1] <- alphaTheta(theta)


  Xlin <- x %*% alphahat
  storeorder <- order(Xlin)
  yobs<-y[order(Xlin)]
  xobs<-x[order(Xlin)]

  #Finding inverse of C
  Cinv <- invcov(Xlin, kappa_init)


  #generate from f
  sigma2<-(sigma2_init^2)
  Sigma<-solve(Cinv+diag(n)/sigma2)
  mu<-(Sigma %*% yobs)/sigma2
  f <- mvrnorm(1,mu=mu,Sigma=Sigma,tol=1e-8)

  #Generate from sigma2
  # a=7;
  # b=2;
  err<-sum((yobs-f)^2)/2
  sigma2<-1/rgamma(1,shape=(a+n/2),rate=(b+err))

  TEMPlist <- list(psi_init = theta, index_init = alphahat[,1],
            kappa_init = kappa_init, Xlin_init = Xlin, d_init = parabeta[,2],
            sigma2_init = sigma2, linkFunction_init = f)

  init_psi <- TEMPlist$psi_init/pi
  init_index <- TEMPlist$index_init
  init_kappa <- TEMPlist$kappa_init
  init_Xlin <- TEMPlist$Xlin_init[,1]
  init_cov <- expcov_gpPolar(init_Xlin, init_kappa)
  init_sigma2 <- TEMPlist$sigma2_init
  init_linkFunction <- TEMPlist$linkFunction_init
  init_d <- TEMPlist$d_init

  return(list(psi = init_psi, index = init_index,
              kappa = init_kappa, d = init_d,
              Xlin = init_Xlin, cov = init_cov,
              sigma2 = init_sigma2, linkFunction = init_linkFunction))

}


optimized<-function(f,xobs,kappa,thetahat,p)
{
  n=length(f)
  thetahat<-as.numeric(as.matrix(thetahat))
  alphahat<-matrix(0,p,1)
  alphahat[,1] <- alphaTheta(thetahat)
  #alphahat<-as.matrix(alphahat)
  x<-xobs %*% alphahat
  f<-f[order(x)]
  x<-x[order(x)]

  r=matrix(0,n+1,1)
  r[2:n]=exp(-kappa*(x[2:n]-x[1:(n-1)]))
  r[1]=0
  r[n+1]=0

  e<-matrix(0,n+1,1)
  e<-r/(1-(r^2))
  d<-matrix(0,n,1)
  d[1:n]<-1+r[1:n]*e[1:n]+r[2:(n+1)]*e[2:(n+1)]

  Cinv<-matrix(0,n,n)
  diag(Cinv)=d

  for(i in 1:(n-1))
  {
    Cinv[i,i+1]=-e[i+1]
    Cinv[i+1,i]=-e[i+1]
  }


  return(0.5*log((det(Cinv)))+(-0.5*f%*%Cinv%*%f))


}

betaalpha1<-function(theta)
{
  p <- length(theta)
  t <- matrix(0,p,1)
  t[,1] <- theta/pi
  abeta <- matrix(0,p,1)
  abeta[,1] <- 5000
  bbeta <- matrix(0,p,1)
  for(i in 1:p){
    bbeta[i]<-(abeta[i]-1-abeta[i]*t[i]+2*t[i])/t[i]}
  return(cbind(abeta,bbeta))
}


alphaTheta <- nimbleFunction(
  run = function(theta = double(1)){
    returnType(double(1))
    p <- length(theta) + 1
    alpha <- numeric(p)
    for (i in 1:p){
      if (i == 1){
        alpha[i] <- sin(theta[i])
      } else if (i == p){
        alpha[i] <- prod(cos(theta[1:(i-1)]))
      } else {
        alpha[i] <- prod(cos(theta[1:(i-1)]))*sin(theta[i])
      }
    }

    alpha <- alpha/sqrt(sum(alpha^2))

    if (alpha[1] < 0) {
      alpha <- -alpha
    }

    return(alpha)
  }
)

############################ Linear predictor ####################################
Xlinear <- nimbleFunction(
  run = function(betavec = double(1), dataX = double(2)){
    returnType(double(1))
    n <- nimDim(dataX)[1]
    p <- nimDim(dataX)[2]

    Xlin <- numeric(n)
    for (i in 1:n) {
      Xlin[i] <- 0
      for (j in 1:p) {
        Xlin[i] <- Xlin[i] + dataX[i, j] * betavec[j]
      }
    }

    return(Xlin)

  }
)


invcov <- nimbleFunction(
  run = function(xlin = double(1), kappa = double(0)){
    returnType(double(2))
    # sort Xlin
    sortedXlin <- nimSort(xlin) ## needs to be coded
    n <- length(xlin)

    r <- nimNumeric(n+1, init = TRUE)
    r[2:n] <- exp(-kappa*(sortedXlin[2:n]-sortedXlin[1:(n-1)]))


    q <- nimNumeric(n+1, init = TRUE)
    denom <- 1 - r^2
    for (i in 1:(n+1)) {
      if (abs(denom[i]) < 1e-8) denom[i] <- 1e-8
    }
    q <- r / denom

    s <- nimNumeric(n, init = TRUE)
    s <- 1 + r[1:n]*q[1:n] + r[2:(n+1)]*q[2:(n+1)]

    # tridiagonal
    Cinv <- nimMatrix(0,nrow = n, ncol = n)
    for(i in 1:(n-1)){
      Cinv[i, i] <- s[i]
      Cinv[i,i+1] <- -q[i+1]
      Cinv[i+1,i] <- -q[i+1]
    }
    Cinv[n, n] <- s[n]


    return(Cinv)

  }
)


expcov_gpPolar <- nimbleFunction(
  run = function(vec = double(1), kappa = double(0)){
    returnType(double(2))
    n = length(vec)
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- exp(-kappa * abs(vec[i] - vec[j]))
      }
    }
    result = (result + t(result))/2 + diag(rep(1e-4, n))

    return(result)

  }
)


expcovTest_gpPolar <- nimbleFunction(
  run = function(vec1 = double(1), vec2 = double(1), kappa = double(0)){
    returnType(double(2))
    n1 = length(vec1); n2 = length(vec2)
    result <- matrix(nrow = n1, ncol = n2, init = FALSE)
    for(i in 1:n1){
      for(j in 1:n2){
        result[i, j] <- exp(-kappa * abs(vec1[i] - vec2[j]))
      }
    }
    return(result)

  }
)


obj_btt_theta <- nimbleFunction(
  setup = function(model, p) {
    p_local <- p
  },
  run = function(par = double(1)) {
    returnType(double(0))
    theta <- par[1:p_local]

    beta <- alphaTheta(theta*pi)
    n = nimDim(model$x)[1]
    Xlin <- Xlinear(beta, model$x)

    kappa <- model$kappa[1]
    inv_cov <- invcov(Xlin, kappa)
    orderf <- nimOrder(Xlin)
    f <- model$linkFunction[orderf]
    matF <- nimMatrix(f, ncol = 1)
    lp <- dmnorm_chol(f, mean = rep(0, n),
                      cholesky = chol(inv_cov), prec_param = TRUE, log=TRUE)
    return(lp)
  }
)


thetaPrior <- nimbleFunction(
  run = function(theta = double(1)) {
    returnType(double(0))
    logP <- 0.0
    p_1    <- nimDim(theta)[1]
    p    <- p_1 + 1
    for(j in 1:p_1) {
      logP <- logP + (p - j) * log(abs(cos(theta[j])))
    }
    return(logP)
  }
)

initfunction_gpPolar2 <- function(X, kappa, sigma2, psi , setSeed){

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
  init_kappa <- kappa
  init_Xlin <- as.vector(X %*% matrix(init_index, ncol = 1))
  init_cov <-  expcov_gpPolar(init_Xlin, init_kappa)
  init_linkFunction <- mvtnorm::rmvnorm(1, rep(0, n), sigma = init_cov)[1,]

  return(list(psi = init_psi, index = init_index,
              kappa = init_kappa, d = init_d,
              Xlin = init_Xlin, cov = init_cov,
              sigma2 = init_sigma2, linkFunction = init_linkFunction))

}


pred_gpPolar <- nimbleFunction(
  run = function(newdata = double(2), nsamp = integer(0), y = double(1),
                 indexSample = double(2),
                 XlinSample = double(2), sigma2_samples = double(1),
                 kappaSample = double(1), prediction = integer(0)){
    returnType(double(2))

    new_ncol <- nimDim(newdata)[1]
    orig_ncol <- length(y)

    testPred <- matrix(0, nrow = nsamp, ncol = new_ncol)
    for (i in 1:nsamp){
      currKappa <- kappaSample[i]
      currSigma <- sigma2_samples[i]
      # 1. compute index
      sampleZ <- newdata %*% matrix(indexSample[i, ], ncol = 1)

      # 2. compute covariance
      ## 1) C(ori, ori)
      cov_ori_ori <- expcov_gpPolar(XlinSample[i,], currKappa)
      ## 2) C(new, ori)
      cov_new_ori <- expcovTest_gpPolar(sampleZ[,1], XlinSample[i,], currKappa)
      ## 3) c(new, new)
      cov_new_new <- expcov_gpPolar(sampleZ[,1], currKappa)

      # 3. compute mu, Sigma
      midMatrix <- inverse(cov_ori_ori + diag(rep(currSigma, orig_ncol)))
      mu <- cov_new_ori %*% midMatrix %*% matrix(y, ncol = 1)
      Sigcov <- cov_new_new - cov_new_ori %*% midMatrix %*%  t(cov_new_ori)

      # 4. Sampling
      if (prediction == 1){ # latent
        cholpredcov <- chol(Sigcov)
        testPred[i,] <- t(rmnorm_chol(1, mean = mu[,1], cholesky = cholpredcov,
                                      prec_param  = FALSE))
      } else{ # response
        predcov <- Sigcov + diag(rep(currSigma, new_ncol))
        cholpredcov <- chol(predcov)
        testPred[i,] <- t(rmnorm_chol(1, mean = mu[,1], cholesky = cholpredcov,
                                      prec_param  = FALSE))
      }

    }
    return(testPred)
  }
)
