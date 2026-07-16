#' Analytic dose-response curve from a corrected parameter table
#'
#' Internal helper: evaluate a logistic dose-response curve from the
#' \code{parameters} data frame stored in a fit result.  Used by the plot
#' functions when
#' \code{result$biological_plausibility_check$needs_correction == TRUE}
#' so the rendered curve reflects the corrected parameters shown in the
#' batch report rather than the raw drm coefficients used by
#' \code{predict(model, ...)}.
#'
#' Supports both 3PL (rows: Bottom / Top / LogIC50 / IC50 / Span) and
#' 4PL (rows: Bottom / Top / LogIC50 / HillSlope / IC50 / Span) tables.
#' For 3PL, HillSlope defaults to \code{hill_default} (caller-supplied:
#' -1 for inhibition, +1 for activation) when not present.
#'
#' @param x Numeric vector of log10 inhibitor concentrations.
#' @param parameters Data frame with columns Parameter and Value.
#' @param hill_default Numeric, used when parameters lacks a HillSlope row.
#'
#' @return Numeric vector of predicted response values, one per \code{x}, or
#'   \code{NULL} if the parameter table lacks Bottom / Top / LogIC50.
#' @noRd
analytic_dose_response <- function(x, parameters, hill_default = -1) {
  if (is.null(parameters) || !is.data.frame(parameters)) return(NULL)
  need <- c("Bottom", "Top", "LogIC50")
  if (!all(need %in% parameters$Parameter)) return(NULL)
  b <- parameters$Value[parameters$Parameter == "Bottom"][1]
  t <- parameters$Value[parameters$Parameter == "Top"][1]
  li <- parameters$Value[parameters$Parameter == "LogIC50"][1]
  if (!is.finite(b) || !is.finite(t) || !is.finite(li)) return(NULL)
  hs <- if ("HillSlope" %in% parameters$Parameter) {
    v <- parameters$Value[parameters$Parameter == "HillSlope"][1]
    if (is.finite(v) && v != 0) v else hill_default
  } else hill_default
  b + (t - b) / (1 + 10^((x - li) * hs))
}
