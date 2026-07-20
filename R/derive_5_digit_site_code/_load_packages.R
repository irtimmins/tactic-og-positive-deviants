# =============================================================================
# _load_packages.R  -  packages the five-character site-code build needs
# -----------------------------------------------------------------------------
# Just the library() calls, so the environment this stage depends on can be read
# (and provisioned) in one glance, separately from the logic in _helpers.R and
# the settings in 01_define_parameters.R. Sourced by 01 before _helpers.R.
# =============================================================================

suppressPackageStartupMessages({
  library(haven)     # read the .dta extracts
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(magrittr)  # %>%
})