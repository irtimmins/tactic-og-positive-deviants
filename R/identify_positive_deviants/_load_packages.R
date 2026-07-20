# =============================================================================
# _load_packages.R  -  packages the positive-deviance analysis needs
# -----------------------------------------------------------------------------
# The one package every script in this stage uses is loaded here (dplyr).
# The heavier, script-specific packages are loaded by the scripts that need
# them, so a stage that only builds the cohort does not pay to load the modelling
# stack:
#   02_build_cohort           lubridate
#   03_table1_characteristics flextable, officer   (guarded - Word output only)
#   04_estimation_weights     balancer, ggplot2    (ggplot2 guarded - figure only)
#   05/06 shrinkage & ranks   rstan, ggplot2       (ggplot2 guarded)
# Keeping those in their own scripts is what lets each be run - and the cohort
# build tested - without the full modelling stack installed.
# =============================================================================

suppressPackageStartupMessages(library(dplyr))
