# =============================================================================
# 01  configuration for the positive-deviance modelling
# -----------------------------------------------------------------------------
# Sourced at the top of every script in this stage. It builds on stage 05:
# everything about the cohort - the paths, the funnel, the outcome, the volume
# floor - comes from that stage's config, sourced first, so the two stages can
# never disagree about what the cohort is. Only the modelling settings are here.
#
# What this stage does, in plain terms:
#   For each diagnosing hospital, estimate the average endoscopy-to-decision-to-
#   treat time adjusted for patient mix (case-mix standardisation by balancing
#   weights), stabilise small-hospital estimates with a Bayesian shrinkage step,
#   then rank hospitals and flag consistently fast ("positive deviant")
#   providers.
#
# It reads the analysis cohort (pd_cohort.rds) produced by stage 05, so run
# that stage first.
# ============================================================================

source("R/05_derive_analysis_cohort/01_config.R")
source("R/06_identify_positive_deviants/_load_packages.R")
source("R/06_identify_positive_deviants/_helpers.R")

# the shrinkage model file, alongside these scripts in the repo.
stan_dir  <- "R/06_identify_positive_deviants"
stan_file <- file.path(stan_dir, "dp_normal_cont.stan")

# patient-level intermediate (carries patient_pseudo_id): dir_out, as in 05.
fit_rds <- file.path(dir_out, "fit_primary.rds")   # holds the weighted patient frame

# site/hospital-level intermediates, no patient rows: dir_debug, in the repo.
site_rds <- file.path(dir_debug, "site_sustained.rds")
stan_rds <- file.path(dir_debug, "stan_sustained.rds")

# case-mix adjustment --------------------------------------------------------
# the model adjusts for age at diagnosis (linear), comorbidity (RCS Charlson,
# coded 0 / 1 / 2+), season (quarter of diagnosis) and calendar-year half. No
# other patient factors are adjusted for.
#
#   age_coding : "cont" linear age, or "band" (<50 / 50-59 / 60-69 / 70-79 / 80+)
#   cci_coding : "0_1_2plus" (0 / 1 / 2+), "0_1_2_3plus", "cont", or "none"
age_coding <- "cont"
cci_coding <- "0_1_2plus"

# season (quarter) and calendar-period dummy names, built in 02.
season_terms <- c("q2", "q3", "q4")   # quarter-of-diagnosis dummies (ref = Q1)
year_term    <- "yr_late"             # later calendar-year half of the window

# build the case-mix covariates for a given coding. Returns the augmented data
# and the continuous / binary covariate names for the weighting step. Reference
# categories are 60-69 for age bands and 0 conditions for comorbidity.
code_covariates <- function(d, age = age_coding, cci = cci_coding) {
  cont <- character(0); bin <- character(0)
  if (age == "cont") {
    cont <- c(cont, "agediag")
  } else if (age == "band") {
    b <- cut(d$agediag, c(-Inf, 50, 60, 70, 80, Inf),
             labels = c("u50", "50_59", "60_69", "70_79", "80p"))
    for (lv in c("u50", "50_59", "70_79", "80p")) {
      col <- paste0("age_", lv); d[[col]] <- as.integer(b == lv); bin <- c(bin, col)
    }
  }
  if (cci == "cont") {
    cont <- c(cont, "cci_n_conditions")
  } else if (cci == "0_1_2plus") {
    d$cci_1  <- as.integer(d$cci_n_conditions == 1)
    d$cci_2p <- as.integer(d$cci_n_conditions >= 2)
    bin <- c(bin, "cci_1", "cci_2p")
  } else if (cci == "0_1_2_3plus") {
    d$cci_1  <- as.integer(d$cci_n_conditions == 1)
    d$cci_2  <- as.integer(d$cci_n_conditions == 2)
    d$cci_3p <- as.integer(d$cci_n_conditions >= 3)
    bin <- c(bin, "cci_1", "cci_2", "cci_3p")
  }
  list(data = d, cont = cont, bin = bin)
}

# case-mix balance strictness (balancing weights) ---------------------------
# lambda controls how hard the weights push for an exact match to the overall
# patient mix. 0 = match as closely as possible; larger = a gentler match that
# keeps more effective sample size. lambda_main is the working value; lambda_grid
# is scanned to show the trade-off.
lambda_grid <- c(0, .01, .05, .1, .25, .5, 1, 1.5, 2, 2.5, 3)
lambda_main <- 0.01

# Bayesian shrinkage priors (in days) ----------------------------------------
# prior_mu_sd: how far the overall average wait could plausibly sit from the data
#   mean. prior_tau_scale: the typical size of genuine between-hospital
#   differences. Both weakly informative on the day scale.
prior_mu_sd     <- 50
prior_tau_scale <- 10

# effective sample size ------------------------------------------------------
# after weighting, a hospital with fewer than this many effective patients has
# little usable information; flagged in diagnostics.
ess_threshold <- 5

# candidate selection --------------------------------------------------------
# a sustained candidate has at least this posterior probability of sitting in the
# fastest 20% of hospitals.
prob_cut <- 0.80
# ----------------------------------------------------------------------------# ----------------------------------------------------------------------------
# Helper functions (standardisation, weighting, shrinkage, ranking) live in
# _helpers.R in this folder, sourced above.
