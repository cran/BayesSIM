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
