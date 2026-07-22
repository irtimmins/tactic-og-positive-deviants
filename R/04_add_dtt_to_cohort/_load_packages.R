# =============================================================================
# _load_packages.R  -  packages the CWT merge and waiting-time derivation need
# -----------------------------------------------------------------------------
# Just the library() calls. Sourced by 01_define_parameters.R. The stage's own
# helpers (modality_group_of and the rest) still live in 01 for now - only the
# package loading and the shared paths/utilities have been split out so far.
# =============================================================================

suppressPackageStartupMessages({
  library(haven)     # read the .dta extracts
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(magrittr)  # %>%
})