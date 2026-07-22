# =============================================================================
# _load_packages.R  -  packages the endoscopy-site-from-HES build needs
# -----------------------------------------------------------------------------
# Just the library() calls, so the environment this stage depends on can be read
# (and provisioned) in one glance, separately from the logic in _helpers.R and
# the settings in 01_define_parameters.R. Sourced by 01 before _helpers.R.
#
# tidyr is here and not in the COSD stage: this one reshapes the 24 operation
# slots long to find the endoscopy, which the COSD stage never has to do.
# =============================================================================

suppressPackageStartupMessages({
  library(haven)     # read the .dta extracts
  library(dplyr)
  library(tidyr)     # pivot_longer, to unfold the 24 operation slots
  library(stringr)
  library(tibble)
  library(purrr)
  library(magrittr)  # %>%
})