# =============================================================================
# R/shared/utils.R  -  small, generic utilities used by more than one stage
# -----------------------------------------------------------------------------
# Only genuinely cross-stage helpers live here: tidy_org_code (used by the
# site-code build and the reference scripts) and check_input (a stop-early input
# guard, generic enough that any stage might want it). Stage-specific helpers
# stay in that stage's own _helpers.R.
#
# This file loads the one package its functions need, rather than relying on a
# stage having loaded it first, so it is safe to source on its own.
# =============================================================================

suppressPackageStartupMessages(library(stringr))

# Tidy an organisation code: upper case, drop anything that is not a letter or a
# digit, and treat an empty string as missing rather than as a value. On the
# 20260212 extracts every diagnosis_trust and site_code_of_diagnosis value is
# already a clean compact code, so this only guards against stray characters.
tidy_org_code <- function(x) {
  x <- str_to_upper(str_trim(as.character(x)))
  x <- str_replace_all(x, "[^A-Z0-9]", "")
  dplyr::na_if(x, "")
}

# Stop early and clearly if a file is not the one we think it is. A missing
# column noticed here saves a puzzling empty result later.
check_input <- function(df, needed, label, path, min_rows = 1000L) {
  absent <- setdiff(needed, names(df))
  if (length(absent))
    stop(label, " has no ", paste(absent, collapse = ", "), " column.",
         "\n  - check that ", basename(path), " is the expected file.",
         call. = FALSE)
  if (nrow(df) < min_rows)
    stop(label, " has only ", nrow(df), " rows, expected at least ", min_rows,
         ".\n  - this looks like a part-read or a test file rather than the real",
         " one: ", basename(path), call. = FALSE)
  invisible(TRUE)
}