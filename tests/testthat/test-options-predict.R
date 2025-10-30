test_that("predict: input validation errors are thrown early", {
  # predict.bsimGp
  obj0 <- structure(list(sampling = NULL), class = "bsimGp")
  expect_error(BayesSIM:::predict.bsimGp(obj0), "No MCMC samples")
  obj <- structure(
    list(
      modelName = "gpFisher",
      input = list(data = list(x = matrix(rnorm(12), 4, 3), y = rnorm(4))),
      sampling = matrix(0, 1, 1),
      fitted   = matrix(rnorm(20), 5, 4)
    ),
    class = "bsimGp"
  )

  expect_error(BayesSIM:::predict.bsimGp(obj, method = 1L), "method.*specified")
  expect_error(BayesSIM:::predict.bsimGp(obj, method = "mode"), "mean or median")
  expect_error(BayesSIM:::predict.bsimGp(obj, se.fit = "yes"), "se.fit.*specified")
  expect_error(BayesSIM:::predict.bsimGp(obj, level = c(0.9, 0.95)), "level.*specified")
  expect_error(BayesSIM:::predict.bsimGp(obj, newdata = 1:3), "newdata should be data.frame")
  bad_df <- data.frame(a = rnorm(4), b = rnorm(4))
  expect_error(BayesSIM:::predict.bsimGp(obj, newdata = bad_df), "Dimension of newdata")

  # predict.bsimSpline
  obj1 <- structure(list(sampling = NULL), class = "bsimSpline")
  expect_error(BayesSIM:::predict.bsimSpline(obj1), "No MCMC samples")
  obj <- structure(
    list(
      modelName = "bsFisher",
      input = list(data = list(x = matrix(rnorm(12), 4, 3), y = rnorm(4))),
      sampling = matrix(0, 1, 1),
      fitted   = matrix(rnorm(20), 5, 4)
    ),
    class = "bsimSpline"
  )

  expect_error(BayesSIM:::predict.bsimSpline(obj, method = 1L), "method.*specified")
  expect_error(BayesSIM:::predict.bsimSpline(obj, method = "mode"), "mean or median")
  expect_error(BayesSIM:::predict.bsimSpline(obj, se.fit = "yes"), "se.fit.*specified")
  expect_error(BayesSIM:::predict.bsimSpline(obj, level = c(0.9, 0.95)), "level.*specified")
  expect_error(BayesSIM:::predict.bsimSpline(obj, newdata = 1:3), "newdata should be data.frame")
  bad_df <- data.frame(a = rnorm(4), b = rnorm(4))
  expect_error(BayesSIM:::predict.bsimSpline(obj, newdata = bad_df), "Dimension of newdata")

})

test_that("predict: happy path (newdata = NULL) returns summary/plot/rmse", {
  skip_if_not_installed("ggplot2")

  set.seed(1)
  N <- 6; p <- 3; nsamp <- 10
  obj <- structure(
    list(
      modelName = "gpFisher",
      input = list(data = list(
        x = matrix(rnorm(N * p), N, p),
        y = rnorm(N)
      )),
      sampling = matrix(0, 1, 1),
      fitted   = matrix(rnorm(nsamp * N), nsamp, N)
    ),
    class = "bsimGp"
  )

  fake_summary <- function(object, verbose = FALSE, ...) {
    ac <- data.frame(
      mean   = seq_len(p),
      median = seq_len(p) / 2,
      row.names = paste0("index[", 1:p, "]")
    )
    list(all.chain = ac)
  }

  res <- with_mocked_bindings(
    {
      out <- expect_invisible(BayesSIM:::predict.bsimGp(
        obj, newdata = NULL, method = "mean", se.fit = TRUE, level = 0.9
      ))
      expect_type(out, "list")
      expect_s3_class(out$plot, "ggplot")
      expect_true(is.data.frame(out$summary))
      expect_true(all(c("mean", "sd", "LB", "UB") %in% colnames(out$summary)))
      expect_true(is.numeric(out$rmse))
    },
    summary.bsimGp = fake_summary,
    .package = "BayesSIM"
  )

  obj2 <- structure(
    list(
      modelName = "bsFisher",
      input = list(data = list(
        x = matrix(rnorm(N * p), N, p),
        y = rnorm(N)
      )),
      sampling = matrix(0, 1, 1),
      fitted   = matrix(rnorm(nsamp * N), nsamp, N)
    ),
    class = "bsimSpline"
  )

  fake_summary <- function(object, verbose = FALSE, ...) {
    ac <- data.frame(
      mean   = seq_len(p),
      median = seq_len(p) / 2,
      row.names = paste0("index[", 1:p, "]")
    )
    list(all.chain = ac)
  }

  res <- with_mocked_bindings(
    {
      out <- expect_invisible(BayesSIM:::predict.bsimSpline(
        obj2, newdata = NULL, method = "mean", se.fit = TRUE, level = 0.9
      ))
      expect_type(out, "list")
      expect_s3_class(out$plot, "ggplot")
      expect_true(is.data.frame(out$summary))
      expect_true(all(c("mean", "sd", "LB", "UB") %in% colnames(out$summary)))
      expect_true(is.numeric(out$rmse))
    },
    summary.bsimSpline = fake_summary,
    .package = "BayesSIM"
  )


})
