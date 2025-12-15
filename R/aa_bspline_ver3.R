#' Additional Functions for Spline
SplineState <- nimble::nimbleList(
  x                       = double(1),
  internal_knots          = double(1),
  boundary_knots          = double(1),
  degree                  = integer(0, default = 3),
  order                   = integer(0, default = 4),
  spline_df               = integer(0, default = 4),

  knot_sequence           = double(1),
  has_internal_multiplicity = logical(0, default = 0),
  is_knot_sequence_latest   = logical(0, default = 0),
  is_extended_knot_sequence = logical(0, default = 0),

  x_index                 = integer(1),
  is_x_index_latest       = logical(0),

  basis = double(2)
)


any_duplicated <- nimble::nimbleFunction(
  run = function(x = double(1)) {
    returnType(logical(0))
    n <- length(x)
    if(n <= 1)
      return(FALSE)

    for(i in 1:n) {
      if(i > 1) {
        for(j in 1:(i-1)) {
          if(x[i] == x[j])
            return(TRUE)            # at least one duplicated element
        }
      }
    }
    return(FALSE)                   # no duplicated
  }
)






mat_wo_col1 <- nimble::nimbleFunction(
  run = function(mat = double(2)){
    returnType(double(2))
    x_ncol <- nimDim(mat)[2]
    if (x_ncol > 1){
      idxCol <- 2:x_ncol
      return(mat[, idxCol])
    } else{
      nimStop("No column left in the matrix.")
    }
  }
)


update_spline_df = nimble::nimbleFunction(
  run = function(temp = SplineState()){
    temp$spline_df <- nimDim(temp$internal_knots)[1]+ temp$order
  }
)


update_x_index = nimble::nimbleFunction(
  run = function(temp = SplineState()){
    if (temp$is_x_index_latest & nimDim(temp$x_index)[1] > 0){
      return()
    }

    n <- nimDim(temp$x)[1]
    x_index <- numeric(n)
    m <- nimDim(temp$internal_knots)[1]

    if (m == 0) {
      x_index <- nimRep(1, nimDim(temp$x)[1])
    } else{
      idx <- nimOrder(temp$x)
      sortedX <- nimSort(temp$x)

      j <- 0
      for(i in 1:n) {
        xi <- sortedX[i]
        while((j < m) & (xi >= temp$internal_knots[j+1]))
          j <- j + 1L
        x_index[i] <- (j)
      }

      x_index <- x_index[nimOrder(idx*1.0)]
    }
    # print(x_index)
    temp$x_index <- x_index
    temp$is_x_index_latest <- TRUE

  })


update_knot_sequence = nimble::nimbleFunction(
  run = function(temp = SplineState()){
    n <- temp$order
    out <- nimC(nimRep(temp$boundary_knots[1], n), temp$internal_knots,
                nimRep(temp$boundary_knots[2], n))
    temp$knot_sequence <- out
    temp$is_knot_sequence_latest <- TRUE
  })



get_basis_simple <- nimble::nimbleFunction(
  run = function(temp = SplineState()){
    returnType(double(2))

    update_spline_df(temp)
    update_x_index(temp)

    bmat <- nimMatrix(0, nimDim(temp$x)[1], temp$spline_df)

    for(i in 1:nimDim(temp$x)[1]){
      bmat[i, (temp$x_index[i]+1)] <- 1
    }

    if (temp$degree > 0){
      update_knot_sequence(temp)
    }

    for (k in 1:temp$degree){
      k_offset <- temp$degree - k
      for (i in 1:nimDim(temp$x)[1]){
        saved <- 0.0
        for (j in 1:k){
          j_index <- (temp$x_index[i]) + j
          i1 <- j_index + k_offset + 1
          i2 <- j_index + temp$order
          den <- temp$knot_sequence[i2] - temp$knot_sequence[i1]
          term <- bmat[i, j_index]/den
          bmat[i, j_index] = saved + (temp$knot_sequence[i2] - temp$x[i])*term
          saved <- (temp$x[i] - temp$knot_sequence[i1])*term
        }
        bmat[i, (temp$x_index[i]+k+1)] <- saved
      }
    }

    return(bmat)

  }
)



simplify_knots <- nimble::nimbleFunction(
  run = function(temp = SplineState(),
                 internal_knots = double(1),
                 boundary_knots = double(1)){

    # boundary knots
    if (nimDim(boundary_knots)[1] == 0){
      if(nimDim(boundary_knots)[1] != 2 & nimDim(temp$x)[1] > 0){
        left <- min(temp$x)
        right <- max(temp$x)
        if (left == right){
          nimStop("Cannot set boundary knots from x.")
        }
        temp$boundary_knots <- nimC(left, right)

      }
    }
    else{
      if (sum(is.na(boundary_knots)) > 0){
        nimStop("Cannot set boundary knots from x.")
      }

      left <- boundary_knots[1]
      right <- boundary_knots[2]

      if (left == right){
        nimStop("Need two distinct boundary knots.")
      }
      temp$boundary_knots <- nimC(left, right)
    }

    # Inner knots
    if (sum(is.na(internal_knots)) > 0){
      nimStop("Internal knots cannot contain NA.")
    }

    if (nimDim(internal_knots)[1] > 0){
      # sort knots
      sortedKnots <- nimSort(internal_knots)
      minKnots <- sortedKnots[1]
      maxKnots <- sortedKnots[nimDim(internal_knots)[1]]

      if (nimDim(temp$boundary_knots)[1] == 2 &
          (temp$boundary_knots[1] >= minKnots | temp$boundary_knots[2] <= maxKnots)){
        nimStop("Internal knots must be set inside boundary.")
      }
      allKnots <- nimC(internal_knots, temp$boundary_knots)
      temp$has_internal_multiplicity <- any_duplicated(allKnots)
      temp$internal_knots <- sortedKnots
    }
    else{
      temp$has_internal_multiplicity <- FALSE
      temp$internal_knots <- numeric(0)
    }
  }
)


get_inside_x = nimble::nimbleFunction(
  run = function(x = double(1), boundary_knots = double(1)){
    returnType(double(1))
    minKnots <- min(boundary_knots)
    maxKnots <- max(boundary_knots)
    elements <- (x >= minKnots & x <= maxKnots)
    res <- x[which(elements)]

    return(res)
  })


gen_default_internal_knots = nimble::nimbleFunction(
  run = function(inside_x = double(1),
                 boundary_knots = double(1),
                 n_internal_knots = integer(0)){
    returnType(double(1))

    prob_vec <- nimSeq(0, 1, length.out = n_internal_knots+2)[2:(n_internal_knots+1)]
    internal_knots <- quantile_nimble(inside_x, prob_vec)
    any_dup  <- any_duplicated(internal_knots)

    if (any_dup){
      nimPrint("Set equidistant internal knots \n")
      nimPrint("(found duplicated knots from quantiles).")
      return(nimSeq(boundary_knots[1], boundary_knots[2],
                    length.out = n_internal_knots+2)[2:(n_internal_knots+1)])
    }

    min_int_knots <- internal_knots[1]
    max_int_knots <- internal_knots[n_internal_knots]
    if (boundary_knots[1] >= min_int_knots |
        boundary_knots[2] <= max_int_knots) {
      return(nimSeq(boundary_knots[1], boundary_knots[2],
                    length.out = n_internal_knots+2)[2:(n_internal_knots+1)])
    }
    return(internal_knots)
  }
)




SplineBase1 = nimble::nimbleFunction(
  run = function(x = double(1), spline_df = integer(0), degree = integer(0, default = 3),
                 intercept = integer(0), Boundary.knots =  double(1, default = numeric(0))){
    returnType(SplineState())

    temp <- SplineState$new()

    temp$x <- x
    temp$degree <- degree
    temp$order <- degree + 1
    temp$spline_df <- spline_df

    # determine internal knots
    num_knots <- spline_df - temp$order
    simplify_knots(temp, numeric(0), Boundary.knots)

    # Internal knots
    if(num_knots > 0){
      inside_x <- get_inside_x(temp$x, temp$boundary_knots)
      temp$internal_knots <- gen_default_internal_knots(inside_x, temp$boundary_knots,num_knots)
    }

    return(temp)

  }
)


SplineBase2 = nimble::nimbleFunction(
  run = function(x = double(1), degree = integer(0, default = 3),
                 internal_knots = double(1),
                 Boundary.knots =  double(1, default = numeric(0))){
    returnType(SplineState())

    temp <- SplineState$new()
    temp$x <- x
    temp$degree <- degree
    simplify_knots(temp, internal_knots, Boundary.knots)
    temp$order <- degree + 1

    return(temp)

  }
)


basis <- nimble::nimbleFunction(
  run = function(temp = SplineState(), complete_basis = logical(0)){
    returnType(SplineState())

    b_mat <- get_basis_simple(temp)

    if (complete_basis){
      temp$basis <- b_mat
    } else{
      b_mat_fix <- mat_wo_col1(b_mat)
      temp$basis <- b_mat_fix
    }

    return(temp)

  }
)


bsNimble <- nimble::nimbleFunction(
  run = function(x = double(1), df = integer(0, default = 0),
                 knots = double(1, default = numeric(0)),
                 degree = integer(0, default = 3),
                 intercept = logical(0, default = FALSE),
                 boundary_knots = double(1, default = numeric(0))){
    # returnType(double(2))
    returnType(SplineState())

    # Decide knots
    if (df > 0 & length(knots) == 0){
      wo_intercept <- !intercept
      spline_df <- df + wo_intercept

      B <- SplineBase1(x              = x,
                       spline_df      = spline_df,
                       degree         = degree,
                       intercept      = wo_intercept,
                       Boundary.knots = boundary_knots)
    } else{
      B <- SplineBase2(x              = x,
                            internal_knots = knots,
                            degree         = degree,
                            Boundary.knots = boundary_knots)
    }

    # Calculate basis
    basisT <- basis(B, complete_basis = intercept)


    return(basisT)

  }
)


bsBasis <- nimble::nimbleFunction(
  run = function(x = double(1), df = integer(0, default = 0),
                 knots = double(1),
                 degree = integer(0, default = 3),
                 intercept = logical(0, default = FALSE),
                 boundary_knots = double(1)){
    returnType(double(2))
    A <- bsNimble(x = x, df = df,
                  knots = knots, degree = degree,
                  intercept = intercept,
                  boundary_knots = boundary_knots)
    return(A$basis)

  }
)
