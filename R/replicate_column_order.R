# Internal helper: group combined replicate columns by their base identifier
# and sort each group by replicate number (base, .2, .3, ...). The order of
# distinct identifiers follows their first appearance in the plate.
.order_combined_replicate_columns <- function(data) {
  if (is.null(dim(data)) || ncol(data) < 2L) {
    return(data)
  }

  column_names <- colnames(data)
  base_names <- sub("\\.\\d+$", "", column_names)
  base_order <- match(base_names, unique(base_names))

  replicate_number <- suppressWarnings(
    as.integer(sub("^.*\\.", "", column_names))
  )
  replicate_number[column_names == base_names] <- 1L

  # Non-standard names without a numeric suffix are treated as the first
  # replicate, while original order resolves any remaining ties.
  replicate_number[is.na(replicate_number)] <- 1L
  ordered_columns <- order(
    base_order,
    replicate_number,
    seq_along(column_names)
  )

  data[, ordered_columns, drop = FALSE]
}
