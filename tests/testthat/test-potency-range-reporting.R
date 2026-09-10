# -----------------------------------------------------------------------------
# Out-of-range potency reporting
# -----------------------------------------------------------------------------

test_that("potency range classification uses the actual tested log10 interval", {
  tested <- c(NA, -9, -8, -7, -6, -5, NA)

  below <- dosefitr:::.classify_potency_range(-10, tested)
  within <- dosefitr:::.classify_potency_range(-7, tested)
  at_lower <- dosefitr:::.classify_potency_range(-9, tested)
  at_upper <- dosefitr:::.classify_potency_range(-5, tested)
  above <- dosefitr:::.classify_potency_range(-4, tested)

  expect_identical(below$status, "Below tested range")
  expect_identical(below$relation, "<")
  expect_true(below$is_outside)
  expect_identical(within$status, "Within tested range")
  expect_false(within$is_outside)
  expect_false(at_lower$is_outside)
  expect_false(at_upper$is_outside)
  expect_identical(above$status, "Above tested range")
  expect_identical(above$relation, ">")
  expect_equal(below$lowest_log10_M, -9)
  expect_equal(above$highest_log10_M, -5)
})

test_that("estimate, censor, and na policies format both sides correctly", {
  tested <- -9:-5

  estimate <- dosefitr:::.potency_reporting(-10, tested, "estimate")
  below <- dosefitr:::.potency_reporting(-10, tested, "censor")
  above <- dosefitr:::.potency_reporting(-4, tested, "censor")
  missing <- dosefitr:::.potency_reporting(-10, tested, "na")

  expect_identical(estimate$reported_log10_M, "-10")
  expect_identical(below$reported_log10_M, "<-9")
  expect_identical(below$reported_p, ">9")
  expect_true(startsWith(below$reported_M, "<"))
  expect_identical(above$reported_log10_M, ">-5")
  expect_identical(above$reported_p, "<5")
  expect_true(startsWith(above$reported_M, ">"))
  expect_true(is.na(missing$reported_log10_M))
  expect_true(is.na(missing$reported_M))
  expect_true(is.na(missing$reported_p))
})

test_that("batch_drc_analysis exposes a backward-compatible reporting policy", {
  expect_true("outside_range" %in% names(formals(batch_drc_analysis)))
  expect_identical(
    eval(formals(batch_drc_analysis)$outside_range),
    c("legacy", "estimate", "censor", "na")
  )
})
