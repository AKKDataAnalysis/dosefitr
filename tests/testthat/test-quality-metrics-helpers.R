test_that("shared Z'-factor helper calculates from replicated controls", {
  out <- .dosefitr_z_prime(
    control_0_values   = c(10, 12),
    control_100_values = c(100, 102)
  )

  expected <- 1 - 3 * (stats::sd(c(10, 12)) + stats::sd(c(100, 102))) /
    (mean(c(100, 102)) - mean(c(10, 12)))

  expect_equal(out$value, expected)
  expect_identical(out$comment, "high (>0.7)")
})

test_that("shared Z'-factor helper returns NA rather than error for one control", {
  one_zero <- .dosefitr_z_prime(
    control_0_values   = 10,
    control_100_values = c(100, 102)
  )
  one_hundred <- .dosefitr_z_prime(
    control_0_values   = c(10, 12),
    control_100_values = 100
  )

  expect_true(is.na(one_zero$value))
  expect_true(is.na(one_zero$comment))
  expect_true(is.na(one_hundred$value))
  expect_true(is.na(one_hundred$comment))
})

test_that("fixed 0% control requires a supplied SD", {
  without_sd <- .dosefitr_z_prime(
    control_100_values   = c(100, 102),
    fixed_control_0_mean = 10
  )
  with_sd <- .dosefitr_z_prime(
    control_100_values   = c(100, 102),
    fixed_control_0_mean = 10,
    fixed_control_0_sd   = 2
  )

  expect_true(is.na(without_sd$value))
  expect_true(is.finite(with_sd$value))
})
