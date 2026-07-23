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
#   then rank hospitals and flag providers.
#
#   Two estimands, reported side by side:
#     SUSTAINED  the level over the whole window, standardised to the whole
#                window's patient mix - who is consistently fast.
#     IMPROVED   the change from period 1 to period 2, with BOTH periods
#                standardised to period 1's patient mix - who got faster. Fixing
#                the reference at period 1 is what makes this a change in speed
#                rather than a change in who was diagnosed.
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

# site/hospital-level intermediates, no patient rows.
#
# NOTE: site_rds and stan_rds are NOT read by any script in this stage - 02, 03
# and 04 all build their own paths from out_dir. They are left here because
# something outside this stage may use them, but they point at dir_debug while
# the scripts write to out_dir, and on a real run with dir_transfer set those
# are different folders. Treat them as stale rather than as the real location.
site_rds        <- file.path(dir_debug, "site_sustained.rds")
stan_rds        <- file.path(dir_debug, "stan_sustained.rds")

# the improvement pair, written only when enough sites have patients in both
# periods (see 02). These use out_dir, matching where the scripts actually read
# and write, so they can be used directly.
site_improve_rds <- file.path(out_dir, "site_improve.rds")
stan_improve_rds <- file.path(out_dir, "stan_improve.rds")

# ---------------------------------------------------------------------------
# the two periods for the improvement estimand
# ---------------------------------------------------------------------------
# Written out in full. Nothing is inferred from the data and nothing is derived
# from the cohort's own date range: the periods are whatever is set here.
#
# This matters because stage 05 also carries a `period` column, split at the
# MIDPOINT of the observed endoscopy dates. That midpoint moves whenever the
# extract changes, and it does not fall on a year boundary, so a patient
# diagnosed in late December can land in the "later" period while a table
# captioned by calendar year says otherwise. 02 re-derives period and yr_late
# from the dates below, so this stage never depends on where that midpoint fell.
#
# period 1 is the BASELINE: the earlier arm of the change, and the reference
# population both arms are standardised to. period 2 is the later arm. The
# sustained estimand ignores the split and uses the whole window.
#
# Set these to the two calendar years of the audit window before running.
period_1_start <- as.Date("2023-01-01")
period_1_end   <- as.Date("2023-12-31")
period_2_start <- as.Date("2024-01-01")
period_2_end   <- as.Date("2024-12-31")

# a site needs this many patients in EACH period to get a change estimate. Sites
# below it keep their sustained estimate and are dropped from the improvement
# analysis only - reported in 02 rather than silently missing.
min_per_period <- min_per_year

# --- check the periods before anything reads them --------------------------
# Ordering and overlap are always mistakes and stop the run. Periods that are
# not whole calendar years are allowed but warned about, because the captions in
# 04 are written in years.
local({
  dates <- list(period_1_start = period_1_start, period_1_end = period_1_end,
                period_2_start = period_2_start, period_2_end = period_2_end)
  for (nm in names(dates)) {
    d <- dates[[nm]]
    if (!inherits(d, "Date") || length(d) != 1 || is.na(d))
      stop(nm, " must be a single, non-missing Date - set it in this stage's ",
           '01_config.R, e.g. as.Date("2023-01-01").', call. = FALSE)
  }
  if (period_1_start > period_1_end || period_2_start > period_2_end)
    stop("each period must start before it ends - check the four period dates.",
         call. = FALSE)
  if (period_2_start <= period_1_end)
    stop("period 2 starts on or before period 1 ends (", period_1_end, " / ",
         period_2_start, "). The periods must not overlap: a patient would ",
         "otherwise sit in both arms of the change.", call. = FALSE)
  if (period_2_start != period_1_end + 1)
    warning("there is a gap between period 1 (ends ", period_1_end, ") and ",
            "period 2 (starts ", period_2_start, "). Patients diagnosed in the ",
            "gap contribute to the sustained estimand but to neither arm of the ",
            "improvement one.", call. = FALSE)
  whole_year <- function(a, b)
    format(a, "%m-%d") == "01-01" && format(b, "%m-%d") == "12-31" &&
    format(a, "%Y") == format(b, "%Y")
  if (!whole_year(period_1_start, period_1_end) ||
      !whole_year(period_2_start, period_2_end))
    warning("the periods are not whole calendar years. That is allowed, but the ",
            "captions in 04 are written in years and will need checking.",
            call. = FALSE)
})

# labels for the captions, derived from the dates so a caption cannot disagree
# with the analysis under it. A whole calendar year reads as "2023".
.period_label <- function(a, b) {
  if (format(a, "%m-%d") == "01-01" && format(b, "%m-%d") == "12-31" &&
      format(a, "%Y") == format(b, "%Y")) format(a, "%Y")
  else paste(format(a, "%b %Y"), "to", format(b, "%b %Y"))
}
period_1_label <- .period_label(period_1_start, period_1_end)
period_2_label <- .period_label(period_2_start, period_2_end)
window_label   <- if (grepl("^[0-9]{4}$", period_1_label) &&
                      grepl("^[0-9]{4}$", period_2_label)) {
  paste0(period_1_label, "-", period_2_label)
} else {
  paste(format(period_1_start, "%b %Y"), "to", format(period_2_end, "%b %Y"))
}

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