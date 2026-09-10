# Internal helper shared by panel and legend ordering. It returns positions
# rather than modifying the supplied plots or their source result objects.
.panel_order_index <- function(labels,
                               panel_order = c("original", "alphabetical")) {
  panel_order <- match.arg(tolower(panel_order),
                           c("original", "alphabetical"))

  n_items <- length(labels)
  if (n_items <= 1L || panel_order == "original") {
    return(seq_len(n_items))
  }

  labels <- as.character(labels)
  labels[is.na(labels)] <- ""
  labels_utf8 <- enc2utf8(labels)

  # The final integer key makes ties stable, preserving their input order.
  order(tolower(labels_utf8), labels_utf8, seq_along(labels_utf8),
        method = "radix")
}
