# =============================================================================
# _load_packages.R  -  packages the positive-deviance analysis needs
# -----------------------------------------------------------------------------
# The one package every script in this stage uses is loaded here (dplyr).
# The heavier, script-specific packages are loaded by the scripts that need
# them, so a stage that only builds the cohort does not pay to load the modelling
# stack:
#   02_estimation_weights     balancer, ggplot2    (ggplot2 guarded - figure only)
#   03_shrinkage              rstan
#   04_ranks_caterpillars     rstan, ggplot2       (ggplot2 guarded)
# The cohort-side scripts these settings build on live in stage 05 and need
# none of the above.
# Keeping those in their own scripts is what lets each be run - and the cohort
# build tested - without the full modelling stack installed.
# =============================================================================

suppressPackageStartupMessages(library(dplyr))
