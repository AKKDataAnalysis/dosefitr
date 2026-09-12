test_that("Hill warning ranges are assay- and model-specific", {
  expect_equal(.hill_warning_limits("3pl", "nanobret"), c(0.5, 1.5))
  expect_equal(.hill_warning_limits("4pl", "nanobret"), c(0.5, 2.0))
  expect_equal(.hill_warning_limits("3pl", "viability"), c(0.3, 3.0))
  expect_equal(.hill_warning_limits("4pl", "viability"), c(0.3, 3.0))
})

test_that("Hill warning messages use signed direction-specific ranges", {
  expect_equal(
    .hill_warning_message(-1.8, "inhibition", "3pl", "nanobret"),
    "Hill Slope (expected -1.5 to -0.5): -1.800"
  )
  expect_equal(
    .hill_warning_message(-1.8, "inhibition", "4pl", "nanobret"),
    ""
  )
  expect_equal(
    .hill_warning_message(2.5, "activation", "4pl", "nanobret"),
    "Hill Slope (expected 0.5 to 2): 2.500"
  )
  expect_equal(
    .hill_warning_message(-2.5, "inhibition", "4pl", "viability"),
    ""
  )
  expect_equal(
    .hill_warning_message(-3.5, "inhibition", "4pl", "viability"),
    "Hill Slope (expected -3 to -0.3): -3.500"
  )
})

test_that("Hill warning limits reject unsupported inputs", {
  expect_error(.hill_warning_limits("5pl", "nanobret"), "model")
  expect_error(.hill_warning_limits("4pl", "other"), "assay_type")
})
