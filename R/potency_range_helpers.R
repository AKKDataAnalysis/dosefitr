# Internal helper used by batch_drc_analysis() to compare an estimated potency
# with the concentration interval that was actually tested. Concentrations are
# expected on the log10 molar scale.
.classify_potency_range <- function(log_potency, tested_log_concentrations) {
  tested <- suppressWarnings(as.numeric(tested_log_concentrations))
  tested <- tested[is.finite(tested)]

  empty <- list(
    status = "Not estimable",
    relation = NA_character_,
    is_outside = NA,
    lowest_log10_M = NA_real_,
    highest_log10_M = NA_real_,
    lowest_M = NA_real_,
    highest_M = NA_real_
  )

  if (length(tested) == 0L) {
    empty$status <- "Tested range unavailable"
    return(empty)
  }

  lower <- min(tested)
  upper <- max(tested)
  empty$lowest_log10_M <- lower
  empty$highest_log10_M <- upper
  empty$lowest_M <- 10^lower
  empty$highest_M <- 10^upper

  potency <- suppressWarnings(as.numeric(log_potency)[1L])
  if (!is.finite(potency)) {
    return(empty)
  }

  if (potency < lower) {
    empty$status <- "Below tested range"
    empty$relation <- "<"
    empty$is_outside <- TRUE
  } else if (potency > upper) {
    empty$status <- "Above tested range"
    empty$relation <- ">"
    empty$is_outside <- TRUE
  } else {
    empty$status <- "Within tested range"
    empty$relation <- "="
    empty$is_outside <- FALSE
  }

  empty
}

# Internal formatter for the three explicit outside-range reporting modes.
# It deliberately does not mutate fitted parameters or confidence intervals.
.potency_reporting <- function(log_potency, tested_logs,
                               mode = c("estimate", "censor", "na")) {
  mode <- match.arg(tolower(mode), c("estimate", "censor", "na"))
  range_info <- .classify_potency_range(log_potency, tested_logs)
  log_value <- suppressWarnings(as.numeric(log_potency)[1L])

  reported_log <- if (is.finite(log_value)) {
    as.character(round(log_value, 3))
  } else {
    NA_character_
  }
  reported_M <- if (is.finite(log_value)) {
    fmt_adaptive(10^log_value)
  } else {
    NA_character_
  }
  reported_p <- if (is.finite(log_value)) {
    as.character(round(-log_value, 3))
  } else {
    NA_character_
  }

  if (isTRUE(range_info$is_outside)) {
    if (mode == "censor") {
      boundary_log <- if (range_info$relation == "<") {
        range_info$lowest_log10_M
      } else {
        range_info$highest_log10_M
      }
      boundary_M <- if (range_info$relation == "<") {
        range_info$lowest_M
      } else {
        range_info$highest_M
      }
      reported_log <- paste0(range_info$relation,
                             as.character(round(boundary_log, 3)))
      reported_M <- paste0(range_info$relation, fmt_adaptive(boundary_M))
      # pIC50/pEC50 reverses the inequality because pIC50 = -logIC50.
      p_relation <- if (range_info$relation == "<") ">" else "<"
      reported_p <- paste0(p_relation,
                           as.character(round(-boundary_log, 3)))
    } else if (mode == "na") {
      reported_log <- NA_character_
      reported_M <- NA_character_
      reported_p <- NA_character_
    }
  }

  c(range_info, list(
    fitted_log10_M = log_value,
    fitted_M = if (is.finite(log_value)) 10^log_value else NA_real_,
    reported_log10_M = reported_log,
    reported_M = reported_M,
    reported_p = reported_p,
    mode = mode
  ))
}
