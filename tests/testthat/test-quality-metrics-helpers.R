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

test_that("luciferase signal thresholds preserve defaults and boundaries", {
  defaults <- c(1000, 10000, 100000)

  expect_identical(.classify_luciferase_signal(NA_real_, defaults),
                   "insufficient luciferase signal")
  expect_identical(.classify_luciferase_signal(1000, defaults),
                   "insufficient (<=1000)")
  expect_identical(.classify_luciferase_signal(1000.1, defaults),
                   "low (1000<x<=10000)")
  expect_identical(.classify_luciferase_signal(10000, defaults),
                   "low (1000<x<=10000)")
  expect_identical(.classify_luciferase_signal(10000.1, defaults),
                   "medium (10000<x<=100000)")
  expect_identical(.classify_luciferase_signal(100000, defaults),
                   "medium (10000<x<=100000)")
  expect_identical(.classify_luciferase_signal(100000.1, defaults),
                   "high (>100000)")
})

test_that("luciferase signal thresholds can be customised in one argument", {
  custom <- c(insufficient = 2000, low = 20000, medium = 150000)

  expect_identical(.classify_luciferase_signal(1500, custom),
                   "insufficient (<=2000)")
  expect_identical(.classify_luciferase_signal(15000, custom),
                   "low (2000<x<=20000)")
  expect_identical(.classify_luciferase_signal(120000, custom),
                   "medium (20000<x<=150000)")
  expect_identical(.classify_luciferase_signal(200000, custom),
                   "high (>150000)")
})

test_that("invalid luciferase signal thresholds fail clearly", {
  expect_error(.validate_luciferase_signal_thresholds(c(1000, 10000)),
               "three finite")
  expect_error(.validate_luciferase_signal_thresholds(c(1000, 1000, 100000)),
               "strictly increasing")
  expect_error(.validate_luciferase_signal_thresholds(c(-1, 10000, 100000)),
               "non-negative")
  expect_error(.validate_luciferase_signal_thresholds(c(1000, 10000, Inf)),
               "finite")
})

test_that("Z-prime thresholds preserve defaults and exact boundaries", {
  defaults <- c(0.25, 0.5, 0.7)

  expect_identical(.classify_three_level_quality(0.25, defaults),
                   "insufficient (<=0.25)")
  expect_identical(.classify_three_level_quality(0.3, defaults),
                   "low (0.25<x<=0.5)")
  expect_identical(.classify_three_level_quality(0.5, defaults),
                   "low (0.25<x<=0.5)")
  expect_identical(.classify_three_level_quality(0.6, defaults),
                   "medium (0.5<x<=0.7)")
  expect_identical(.classify_three_level_quality(0.7, defaults),
                   "medium (0.5<x<=0.7)")
  expect_identical(.classify_three_level_quality(0.8, defaults),
                   "high (>0.7)")
})

test_that("Z-prime helper applies custom thresholds", {
  controls_0 <- c(10, 12)
  controls_100 <- c(100, 102)
  custom <- c(-0.5, 0, 0.95)

  out <- .dosefitr_z_prime(
    control_0_values = controls_0,
    control_100_values = controls_100,
    thresholds = custom
  )

  expect_true(out$value > 0 && out$value <= 0.95)
  expect_identical(out$comment, "medium (0<x<=0.95)")
})

test_that("assay-window thresholds preserve defaults and can be customised", {
  defaults <- c(1.5, 2, 3)

  expect_identical(.classify_assay_window(1.5, defaults),
                   "insufficient (<=1.5)")
  expect_identical(.classify_assay_window(2, defaults),
                   "low (1.5<x<=2)")
  expect_identical(.classify_assay_window(3, defaults),
                   "medium (2<x<=3)")
  expect_identical(.classify_assay_window(3.1, defaults),
                   "high (>3)")

  custom <- c(insufficient = 1.2, low = 2.5, medium = 4)
  expect_identical(.classify_assay_window(2, custom),
                   "low (1.2<x<=2.5)")
  expect_identical(.classify_assay_window(3, custom),
                   "medium (2.5<x<=4)")
})

test_that("invalid Z-prime and assay-window thresholds fail clearly", {
  expect_error(.validate_z_prime_thresholds(c(0.25, 0.7, 0.5)),
               "strictly increasing")
  expect_silent(.validate_z_prime_thresholds(c(-1, 0, 0.7)))
  expect_error(.validate_assay_window_thresholds(c(-1, 2, 3)),
               "non-negative")
  expect_error(.validate_assay_window_thresholds(c(1, 2, Inf)),
               "finite")
})
