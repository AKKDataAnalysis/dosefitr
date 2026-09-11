make_rout_model_selection_data <- function() {
  logc <- seq(-9, -5, length.out = 12L)
  response <- 20 + (100 - 20) / (1 + 10^((logc + 7) * 1))

  data.frame(
    "log(inhibitor).[M]" = logc,
    "Compound A" = response,
    "Compound A.2" = response + 0.2,
    check.names = FALSE
  )
}

test_that("explicit ROUT model modes do not switch silently", {
  dat <- make_rout_model_selection_data()

  fit3 <- suppressWarnings(rout_outliers(
    dat, model = "3pl", direction = "inhibition", verbose = FALSE
  ))
  fit4 <- suppressWarnings(rout_outliers(
    dat, model = "4pl", direction = "inhibition", verbose = FALSE
  ))

  expect_identical(fit3$params$model, "3pl")
  expect_identical(fit4$params$model, "4pl")
  expect_true(nrow(fit3$results) > 0L)
  expect_true(nrow(fit4$results) > 0L)
  expect_true(all(fit3$results$model_used == "3PL"))
  expect_true(all(fit4$results$model_used == "4PL"))
})

test_that("auto is the default and is explicit in metadata", {
  dat <- make_rout_model_selection_data()

  auto <- suppressWarnings(rout_outliers(
    dat, model = "auto", direction = "inhibition", verbose = FALSE
  ))
  default_fit <- suppressWarnings(rout_outliers(
    dat, direction = "inhibition", verbose = FALSE
  ))

  expect_identical(auto$params$model, "auto")
  expect_identical(default_fit$params$model, "auto")
  expect_false("n_param" %in% names(formals(rout_outliers)))
  expect_false("n_param" %in% names(auto$params))
})

test_that("batch ROUT propagates the requested fixed model", {
  dat <- make_rout_model_selection_data()
  batch <- list(
    plate_01 = list(
      result = list(modified_ratio_table = dat)
    )
  )

  fit <- suppressWarnings(rout_outliers_batch(
    batch,
    model = "3pl",
    direction = "inhibition",
    verbose = FALSE
  ))

  expect_identical(fit$params$model, "3pl")
  expect_false("n_param" %in% names(formals(rout_outliers_batch)))
  expect_false("n_param" %in% names(fit$params))
  expect_true(all(
    fit$plate_01$result$rout_results$results$model_used == "3PL"
  ))
})

test_that("invalid ROUT model modes are rejected", {
  dat <- make_rout_model_selection_data()

  expect_error(
    rout_outliers(dat, model = "5pl", verbose = FALSE),
    'model must be one of "3pl", "4pl", or "auto"'
  )
})
