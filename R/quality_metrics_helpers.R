# Internal quality-metric helpers shared by the NanoBRET and viability
# pipelines. These functions are intentionally not exported.

.dosefitr_z_prime <- function(control_0_values = NULL,
                              control_100_values = NULL,
                              fixed_control_0_mean = NULL,
                              fixed_control_0_sd = NULL,
                              min_replicates = 2L) {
  min_replicates <- as.integer(min_replicates)
  if (length(min_replicates) != 1L || is.na(min_replicates) ||
      min_replicates < 2L) {
    min_replicates <- 2L
  }

  finite_values <- function(x) {
    x <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE)))
    x[is.finite(x)]
  }

  control_100_values <- finite_values(control_100_values)
  if (length(control_100_values) < min_replicates) {
    return(list(value = NA_real_, comment = NA_character_))
  }

  if (!is.null(fixed_control_0_mean)) {
    mean_0 <- suppressWarnings(as.numeric(fixed_control_0_mean)[1L])
    sd_0   <- suppressWarnings(as.numeric(fixed_control_0_sd)[1L])
    if (!is.finite(mean_0) || length(sd_0) == 0L || !is.finite(sd_0) ||
        sd_0 < 0) {
      return(list(value = NA_real_, comment = NA_character_))
    }
  } else {
    control_0_values <- finite_values(control_0_values)
    if (length(control_0_values) < min_replicates) {
      return(list(value = NA_real_, comment = NA_character_))
    }
    mean_0 <- mean(control_0_values)
    sd_0   <- stats::sd(control_0_values)
  }

  mean_100 <- mean(control_100_values)
  sd_100   <- stats::sd(control_100_values)
  difference <- mean_100 - mean_0

  if (!all(is.finite(c(mean_0, sd_0, mean_100, sd_100, difference))) ||
      difference == 0) {
    return(list(value = NA_real_, comment = NA_character_))
  }

  # Preserve the package's established control orientation: the 100% control
  # is expected to have the larger signal. This intentionally matches the
  # historical NanoBRET calculation.
  value <- 1 - (3 * (sd_100 + sd_0) / difference)

  comment <- if (!is.finite(value)) {
    NA_character_
  } else if (value > 0.7) {
    "high (>0.7)"
  } else if (value > 0.5) {
    "medium (0.5<x<0.7)"
  } else if (value > 0.25) {
    "low (0.25<x<0.5)"
  } else {
    "insufficient (<=0.25)"
  }

  list(value = value, comment = comment)
}
