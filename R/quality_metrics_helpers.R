# Internal quality-metric helpers shared by the NanoBRET and viability
# pipelines. These functions are intentionally not exported.

.validate_three_quality_thresholds <- function(x, argument,
                                               require_non_negative = FALSE) {
  invalid <- !is.numeric(x) || length(x) != 3L || any(!is.finite(x)) ||
    any(diff(x) <= 0)
  if (require_non_negative && is.numeric(x) && any(x < 0, na.rm = TRUE)) {
    invalid <- TRUE
  }
  if (invalid) {
    sign_requirement <- if (require_non_negative) "non-negative, " else ""
    stop(
      paste0(
        argument, " must contain three finite, ", sign_requirement,
        "strictly increasing numeric values: ",
        "c(insufficient_max, low_max, medium_max)."
      ),
      call. = FALSE
    )
  }
  unname(as.numeric(x))
}

.format_quality_threshold <- function(x) {
  format(x, trim = TRUE, scientific = FALSE)
}

.classify_three_level_quality <- function(value, thresholds,
                                          na_label = "insufficient") {
  if (length(value) != 1L || is.na(value) || !is.finite(value)) {
    return(na_label)
  }

  insufficient_max <- thresholds[1L]
  low_max          <- thresholds[2L]
  medium_max       <- thresholds[3L]

  if (value > medium_max) {
    sprintf("high (>%s)", .format_quality_threshold(medium_max))
  } else if (value > low_max) {
    sprintf(
      "medium (%s<x<=%s)",
      .format_quality_threshold(low_max),
      .format_quality_threshold(medium_max)
    )
  } else if (value > insufficient_max) {
    sprintf(
      "low (%s<x<=%s)",
      .format_quality_threshold(insufficient_max),
      .format_quality_threshold(low_max)
    )
  } else {
    sprintf(
      "insufficient (<=%s)",
      .format_quality_threshold(insufficient_max)
    )
  }
}

.validate_z_prime_thresholds <- function(x) {
  .validate_three_quality_thresholds(
    x, "z_prime_thresholds", require_non_negative = FALSE
  )
}

.validate_assay_window_thresholds <- function(x) {
  .validate_three_quality_thresholds(
    x, "assay_window_thresholds", require_non_negative = TRUE
  )
}

.classify_assay_window <- function(value, thresholds = c(1.5, 2, 3)) {
  thresholds <- .validate_assay_window_thresholds(thresholds)
  .classify_three_level_quality(value, thresholds)
}

.dosefitr_z_prime <- function(control_0_values = NULL,
                              control_100_values = NULL,
                              fixed_control_0_mean = NULL,
                              fixed_control_0_sd = NULL,
                              min_replicates = 2L,
                              thresholds = c(0.25, 0.5, 0.7)) {
  thresholds <- .validate_z_prime_thresholds(thresholds)
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

  comment <- if (!is.finite(value)) NA_character_ else
    .classify_three_level_quality(value, thresholds)

  list(value = value, comment = comment)
}


# Validate and classify the mean donor/luciferase signal. The three values are
# the inclusive upper bounds for insufficient, low, and medium, respectively;
# values above the third bound are high.
.validate_luciferase_signal_thresholds <- function(x) {
  .validate_three_quality_thresholds(
    x, "luciferase_signal_thresholds", require_non_negative = TRUE
  )
}

.classify_luciferase_signal <- function(value,
                                        thresholds = c(1000, 10000, 100000)) {
  thresholds <- .validate_luciferase_signal_thresholds(thresholds)
  .classify_three_level_quality(
    value, thresholds, na_label = "insufficient luciferase signal"
  )
}
