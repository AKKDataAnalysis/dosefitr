# -----------------------------------------------------------------------------
# test-plot_multiple_compounds.R -- coverage for the `legend_label` arg.
# -----------------------------------------------------------------------------
#
# `legend_label` controls the text shown in the side-legend of
# `plot_multiple_compounds()`.  Four modes:
#
#   * "auto"      (default) -- smart shortening based on whether construct
#                              or compound is constant across the selection.
#   * "compound"            -- always show only the compound part.
#   * "construct"           -- always show only the construct (target) part.
#   * "full"               -- always show the full "construct:compound" label.
#
# The other side-legend logic (colours, line-wrapping, etc.) is not exercised
# here.  These tests only assert what the four modes place in the legend's
# `labels` vector -- read straight from the ggplot object.

# --- Setup helpers -----------------------------------------------------------

# Build a DRC batch on the bundled NanoBRET fixture and return plate_01's
# `drc_result` slot (the object shape that plot_multiple_compounds accepts
# via `results = <drc_result>`).  Uses the shared stage_nanobret_dir()
# helper from helper-fixtures.R.
build_plate01_drc <- function() {
  work_dir <- stage_nanobret_dir("plmc_nb_")
  ratio_res <- batch_ratio_analysis(
    directory        = work_dir,
    info_file        = "nanobret_info.xlsx",
    data_pattern     = "nanobret_plate_\\d+\\.xlsx$",
    control_0perc    = "1",
    control_100perc  = "24",
    selected_columns = 1:24,
    generate_reports = FALSE,
    output_dir       = tempfile("plmc_nb_out_"),
    verbose          = FALSE
  )
  drc_res <- batch_drc_analysis(ratio_res, model = "3pl",
                                generate_reports = FALSE, verbose = FALSE)
  list(work_dir = work_dir,
       drc_result = drc_res$drc_results$plate_01$drc_result)
}

# Extract the colour-scale `labels` vector from a ggplot returned by
# plot_multiple_compounds().  Named character; names are the raw compound
# keys (e.g. "KinaseA:Cpd1"), values are what the user sees in the legend.
color_labels <- function(gg) {
  scale <- gg$scales$get_scales("colour")
  if (is.null(scale) || is.null(scale$labels))
    stop("no colour scale on this ggplot object")
  scale$labels
}

# Draw the ggplot with a specific legend_label mode.  Returns the ggplot
# so callers can peek at its scales.  save_plot must be a real path per
# plot_multiple_compounds' contract (FALSE is rejected).
render_with_mode <- function(drc_result, mode, out_dir,
                             indices = c(1, 4, 7, 10, 13)) {
  args <- list(
    results          = drc_result,
    compound_indices = indices,
    save_plot        = file.path(out_dir, paste0("mode_", mode, ".png")),
    verbose          = FALSE
  )
  if (!is.null(mode)) args$legend_label <- mode
  do.call(plot_multiple_compounds, args)
}

# --- Tests -------------------------------------------------------------------

test_that("legend_label = 'auto' (default) preserves the smart shortening", {
  bundle <- build_plate01_drc()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  out_dir <- tempfile("plmc_auto_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # Indices 1,4,7,10,13 span >= 2 constructs (KinaseA/KinaseB) AND all
  # compounds differ, so auto's rule "both differ -> show full label"
  # applies.  Expect all labels contain a colon.
  gg   <- render_with_mode(bundle$drc_result, mode = NULL, out_dir = out_dir)
  labs <- color_labels(gg)
  expect_true(all(grepl(":", labs, fixed = TRUE)),
              info = paste("labels:", paste(labs, collapse = ", ")))
  # Values should equal names (raw keys)
  expect_equal(unname(labs), unname(names(labs)))
})

test_that("legend_label defaults to 'auto' when the argument is omitted", {
  bundle <- build_plate01_drc()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  out_dir <- tempfile("plmc_def_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  gg_omitted <- render_with_mode(bundle$drc_result, mode = NULL,   out_dir = out_dir)
  gg_auto    <- render_with_mode(bundle$drc_result, mode = "auto", out_dir = out_dir)

  expect_equal(color_labels(gg_omitted), color_labels(gg_auto))
})

test_that("legend_label = 'compound' shows only the compound part", {
  bundle <- build_plate01_drc()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  out_dir <- tempfile("plmc_cmp_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  gg   <- render_with_mode(bundle$drc_result, mode = "compound", out_dir = out_dir)
  labs <- color_labels(gg)

  # No colon anywhere in the label VALUES (names still carry the raw key)
  expect_false(any(grepl(":", unname(labs), fixed = TRUE)),
               info = paste("labels:", paste(labs, collapse = ", ")))

  # Each value should equal the compound part of its raw key
  expected <- sub("^[^:]+:", "", names(labs))
  expect_equal(unname(labs), expected)
})

test_that("legend_label = 'construct' shows only the construct part", {
  bundle <- build_plate01_drc()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  out_dir <- tempfile("plmc_con_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  gg   <- render_with_mode(bundle$drc_result, mode = "construct", out_dir = out_dir)
  labs <- color_labels(gg)

  expect_false(any(grepl(":", unname(labs), fixed = TRUE)),
               info = paste("labels:", paste(labs, collapse = ", ")))

  # Each value should equal the target part of its raw key
  expected <- sub(":.*$", "", names(labs))
  expect_equal(unname(labs), expected)
})

test_that("legend_label = 'full' always shows the full construct:compound label", {
  bundle <- build_plate01_drc()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  out_dir <- tempfile("plmc_full_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # This scenario picks 5 compounds spanning 2 constructs, so `auto`
  # ALSO shows the full label.  To distinguish `full` from `auto`, we
  # also test on compounds 1:3 (all KinaseA), where auto would drop
  # the construct part but full must keep it.
  gg_multi <- render_with_mode(bundle$drc_result, mode = "full", out_dir = out_dir)
  labs_multi <- color_labels(gg_multi)
  expect_equal(unname(labs_multi), unname(names(labs_multi)))
  expect_true(all(grepl(":", labs_multi, fixed = TRUE)))

  # Same-construct scenario: KinaseA:Cpd1..Cpd3
  gg_same <- render_with_mode(bundle$drc_result, mode = "full",
                              out_dir = out_dir, indices = 1:3)
  labs_same <- color_labels(gg_same)
  expect_true(all(grepl(":", labs_same, fixed = TRUE)),
              info = paste("labels:", paste(labs_same, collapse = ", ")))
  # And `auto` on the same input DOES drop the construct: sanity check
  gg_same_auto <- render_with_mode(bundle$drc_result, mode = "auto",
                                   out_dir = out_dir, indices = 1:3)
  labs_same_auto <- color_labels(gg_same_auto)
  expect_false(any(grepl(":", labs_same_auto, fixed = TRUE)),
               info = paste("auto labels:", paste(labs_same_auto, collapse = ", ")))
  # Therefore auto and full DIFFER on this input, proving the override worked
  expect_false(identical(unname(labs_same), unname(labs_same_auto)))
})

test_that("legend_label rejects invalid modes with a clear error", {
  bundle <- build_plate01_drc()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  out_dir <- tempfile("plmc_bad_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # Random string
  expect_error(
    render_with_mode(bundle$drc_result, mode = "elephant", out_dir = out_dir),
    regexp = "legend_label"
  )
  # Wrong type
  expect_error(
    plot_multiple_compounds(
      results          = bundle$drc_result,
      compound_indices = 1:2,
      save_plot        = file.path(out_dir, "typebad.png"),
      legend_label     = 3L,
      verbose          = FALSE
    ),
    regexp = "legend_label"
  )
  # Length 2
  expect_error(
    plot_multiple_compounds(
      results          = bundle$drc_result,
      compound_indices = 1:2,
      save_plot        = file.path(out_dir, "lenbad.png"),
      legend_label     = c("auto", "compound"),
      verbose          = FALSE
    ),
    regexp = "legend_label"
  )
})

# =============================================================================
# `format` arg coverage
# =============================================================================
#
# `format` is the file-format hint used by save_plot.  Mirrors the same-named
# arg on batch_save_all_drc_plots(): passed straight through as the file
# extension, no whitelist, and ggsave() infers the graphics device from it.
#
# Precedence rules exercised here:
#   * save_plot = TRUE       -> auto-name uses <base>.<format>
#   * save_plot = "path"     -> if path has no extension, append .<format>
#                            -> if path has an extension, path wins and
#                               format is ignored (verbose message on mismatch)

test_that("format = default 'png' produces .png on auto-name", {
  bundle <- build_plate01_drc()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  out_dir <- tempfile("plmc_fmt_A_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  old_wd <- getwd(); setwd(out_dir); on.exit(setwd(old_wd), add = TRUE)

  # save_plot = TRUE lands in getwd(); default format = "png"
  plot_multiple_compounds(bundle$drc_result, compound_indices = 1:2,
                          save_plot = TRUE, verbose = FALSE)
  pngs <- list.files(out_dir, pattern = "^multiple_compounds_.*\\.png$")
  expect_length(pngs, 1L)
})

test_that("format = 'svg' produces .svg on auto-name", {
  skip_if_not_installed("svglite")
  bundle <- build_plate01_drc()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  out_dir <- tempfile("plmc_fmt_B_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  old_wd <- getwd(); setwd(out_dir); on.exit(setwd(old_wd), add = TRUE)

  plot_multiple_compounds(bundle$drc_result, compound_indices = 1:2,
                          save_plot = TRUE, format = "svg", verbose = FALSE)
  svgs <- list.files(out_dir, pattern = "^multiple_compounds_.*\\.svg$")
  expect_length(svgs, 1L)
  # And no PNG should have been produced
  pngs <- list.files(out_dir, pattern = "^multiple_compounds_.*\\.png$")
  expect_length(pngs, 0L)
})

test_that("format applies when save_plot path has no extension", {
  bundle <- build_plate01_drc()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  out_dir <- tempfile("plmc_fmt_C_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  target <- file.path(out_dir, "explicit_no_ext")
  plot_multiple_compounds(bundle$drc_result, compound_indices = 1:2,
                          save_plot = target, format = "pdf", verbose = FALSE)
  expect_true(file.exists(paste0(target, ".pdf")))
  expect_false(file.exists(target))          # no extensionless file
})

test_that("path extension wins when both path and format specify a format", {
  skip_if_not_installed("svglite")
  bundle <- build_plate01_drc()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  out_dir <- tempfile("plmc_fmt_D_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  target <- file.path(out_dir, "explicit.svg")
  # Capture messages: verbose = TRUE should yield an override notice
  msgs <- capture_messages(
    plot_multiple_compounds(bundle$drc_result, compound_indices = 1:2,
                            save_plot = target, format = "pdf", verbose = TRUE)
  )
  # File landed at the SVG path, not the PDF path
  expect_true(file.exists(target))
  expect_false(file.exists(sub("\\.svg$", ".pdf", target)))
  # Override message emitted (may be one of several messages)
  expect_true(any(grepl("ignoring `format", msgs, fixed = TRUE)),
              info = paste("messages:", paste(msgs, collapse = " || ")))
})

test_that("matching path/format combo does NOT emit override message", {
  bundle <- build_plate01_drc()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  out_dir <- tempfile("plmc_fmt_E_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  target <- file.path(out_dir, "explicit.png")
  msgs <- capture_messages(
    plot_multiple_compounds(bundle$drc_result, compound_indices = 1:2,
                            save_plot = target, format = "png", verbose = TRUE)
  )
  expect_true(file.exists(target))
  expect_false(any(grepl("ignoring `format", msgs, fixed = TRUE)),
               info = paste("messages:", paste(msgs, collapse = " || ")))
})

test_that("format = '.svg' (leading dot) is tolerated", {
  skip_if_not_installed("svglite")
  bundle <- build_plate01_drc()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  out_dir <- tempfile("plmc_fmt_G_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  old_wd <- getwd(); setwd(out_dir); on.exit(setwd(old_wd), add = TRUE)

  # Leading dot and whitespace should be stripped before use
  plot_multiple_compounds(bundle$drc_result, compound_indices = 1:2,
                          save_plot = TRUE, format = ".svg", verbose = FALSE)
  svgs <- list.files(out_dir, pattern = "^multiple_compounds_.*\\.svg$")
  expect_length(svgs, 1L)
  # No literal "..svg" extension anywhere
  double_dot <- list.files(out_dir, pattern = "\\.\\.svg$")
  expect_length(double_dot, 0L)
})

test_that("format rejects empty / NA / wrong type with a clear error", {
  bundle <- build_plate01_drc()
  on.exit(unlink(bundle$work_dir, recursive = TRUE), add = TRUE)
  out_dir <- tempfile("plmc_fmt_F_"); dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  target <- file.path(out_dir, "bad")
  # Empty string
  expect_error(
    plot_multiple_compounds(bundle$drc_result, compound_indices = 1:2,
                            save_plot = target, format = "", verbose = FALSE),
    regexp = "`format` must be"
  )
  # NA
  expect_error(
    plot_multiple_compounds(bundle$drc_result, compound_indices = 1:2,
                            save_plot = target, format = NA_character_,
                            verbose = FALSE),
    regexp = "`format` must be"
  )
  # Wrong type
  expect_error(
    plot_multiple_compounds(bundle$drc_result, compound_indices = 1:2,
                            save_plot = target, format = 42, verbose = FALSE),
    regexp = "`format` must be"
  )
  # Length > 1
  expect_error(
    plot_multiple_compounds(bundle$drc_result, compound_indices = 1:2,
                            save_plot = target, format = c("png", "svg"),
                            verbose = FALSE),
    regexp = "`format` must be"
  )
})
