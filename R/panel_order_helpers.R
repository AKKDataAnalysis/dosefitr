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
  labels_utf8 <- tolower(enc2utf8(labels))

  # Build natural-sort keys without converting digit runs to numeric values.
  # Padding every run to the longest run in the selected labels makes, for
  # example, "Compound 3" sort before "Compound 10" while also supporting
  # arbitrarily long identifiers and preserving leading-zero ties.
  digit_matches <- gregexpr("[0-9]+", labels_utf8, perl = TRUE)
  digit_runs <- regmatches(labels_utf8, digit_matches)
  all_digit_runs <- unlist(digit_runs, use.names = FALSE)
  pad_width <- max(c(1L, nchar(all_digit_runs)))
  padded_runs <- lapply(digit_runs, function(runs) {
    if (length(runs) == 0L) return(character())
    paste0(vapply(nchar(runs), function(n) {
      strrep("0", pad_width - n)
    }, character(1L)), runs)
  })
  natural_keys <- labels_utf8
  regmatches(natural_keys, digit_matches) <- padded_runs

  # The integer key makes case-insensitive/numeric ties stable.
  order(natural_keys, seq_along(natural_keys), method = "radix")
}
