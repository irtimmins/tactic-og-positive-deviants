# =============================================================================
# 06  Identify the positive deviants  -  run this step
# -----------------------------------------------------------------------------
# WHERE: the analysis server. Reads the analysis cohort stage 05 wrote. Needs
# the modelling stack (balancer, rstan); if either is missing this stops with a
# clear message rather than failing obscurely partway through.
#
# Estimates each hospital's case-mix-standardised endoscopy-to-DTT time with
# balancing weights, stabilises the estimates with Bayesian shrinkage, then
# ranks hospitals and flags the sustained fast performers.
#
# Needs:    pd_cohort.rds (stage 05)
# Produces: fit_primary.rds, ranks, caterpillar plots, candidate list
# =============================================================================

source("R/config/directories.R")

run_this_step <- function() {
  for (pkg in c("balancer", "rstan")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("stage 06 needs the '", pkg, "' package, which is not installed. ",
           "Stages up to 05 run without it; install it for the modelling.",
           call. = FALSE)
  }
  dir_build <- "R/06_identify_positive_deviants"
  step <- function(file) {
    cat("\n========== ", file, " ==========\n", sep = "")
    source(file.path(dir_build, file), local = new.env())
  }
  step("02_estimation_weights.R")   # balancing weights, standardised means
  step("03_shrinkage.R")            # Bayesian stabilisation
  step("04_ranks_caterpillars.R")   # ranks, plots, sustained candidates
  cat("\nStage 06 complete.\n")
  invisible(NULL)
}
run_this_step()
