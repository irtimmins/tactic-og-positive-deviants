# 00  master - run the OG positive-deviance analysis end to end
# -----------------------------------------------------------------------------
# Run from the project root (the folder holding R/). Each step writes to out_dir
# and the next reads it, so order matters. Settings live in 01_config.R; set
# base_dir (where og_cohort_cwt.rds is) and out_dir (where results go) there, or
# assign them in this session before sourcing.
#
# The analysis identifies diagnosing sites with consistently short times from
# diagnostic endoscopy to the decision to treat, among curatively-treated stage
# 1-3 oesophago-gastric patients, adjusting for age, comorbidity, season and
# calendar year, with Bayesian shrinkage of the per-site estimates.
#
# Prerequisites: dplyr, tidyr, lubridate, balancer, rstan; ggplot2 for the
# figures; flextable and officer for the Word Table 1. The Stan model file
# dp_normal_cont.stan sits alongside these scripts.

r_dir <- "R/identify_positive_deviants"
step <- function(file) {
  message("\n========== ", file, " ==========")
  source(file.path(r_dir, file), local = new.env())
}

step("02_build_cohort.R")            # inclusion + attrition + hospital volume floor
step("03_table1_characteristics.R")  # Table 1 and alliance representativeness
step("04_estimation_weights.R")      # balancing-weights direct standardisation
step("05_shrinkage.R")               # Bayesian shrinkage of the per-site estimates
step("06_ranks_caterpillars.R")      # ranks, caterpillar, candidates

message("\nPositive-deviance analysis complete. Outputs in the configured out_dir.")
