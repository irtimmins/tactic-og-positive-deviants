# 01  configuration and shared helpers for the positive-deviance analysis
# -----------------------------------------------------------------------------
# Sourced at the top of every script in this stage. Part A holds the settings you
# may want to change; Part B holds helper functions you should not need to edit.
#
# What this stage does, in plain terms:
#   Among curatively-treated oesophago-gastric patients, we estimate for each
#   diagnosing hospital the average time from diagnostic endoscopy to the
#   decision to treat, adjusted for differences in patient mix (case-mix
#   standardisation by balancing weights), then stabilise small-hospital
#   estimates with a Bayesian shrinkage step before ranking hospitals and
#   flagging consistently fast ("positive deviant") providers.
#
# It reads the CWT-merged cohort (og_cohort_cwt.rds) produced by the
# merge_cwt_to_get_dtt stage, so run that stage first.

# ============================================================================
# PART A  -  settings
# ============================================================================

# paths ----------------------------------------------------------------------
# Shared paths - the same file every stage in the pipeline sources. Defines
# dir_raw, dir_out (patient-level data, on the restricted drive), dir_ref
# (reference lookups), dir_debug (aggregate, non-patient intermediates, in the
# repo) and dir_transfer (the results-transfer area, off by default). See
# R/config/directories.R for the full explanation.
source("R/config/directories.R")
source("R/identify_positive_deviants/_load_packages.R")
source("R/identify_positive_deviants/_helpers.R")

# input: the CWT-merged patient cohort. f_cohort_cwt is defined in the shared
# paths file, so this stage never has to guess what the merge stage called it.
in_rds <- f_cohort_cwt

# reference lookups (site->trust map, valid-trust list) - dir_ref, from the
# shared paths file.
site_trust_map_csv <- file.path(dir_ref, "site_trust_map.csv")
valid_trusts_csv   <- file.path(dir_ref, "valid_diagnosing_trusts.csv")

# the shrinkage model file, alongside these scripts in the repo.
stan_dir  <- "R/identify_positive_deviants"
stan_file <- file.path(stan_dir, "dp_normal_cont.stan")

# --- where this stage's own outputs go, on the same three-way split as the rest
# of the pipeline -------------------------------------------------------------
# patient-level (carries patient_pseudo_id): dir_out, the restricted drive.
cohort_rds <- file.path(dir_out, "pd_cohort.rds")   # 02 -> 03/04, patient rows
fit_rds    <- file.path(dir_out, "fit_primary.rds") # holds the weighted patient frame

# site/hospital-level intermediates, no patient rows: dir_debug, in the repo.
site_rds <- file.path(dir_debug, "site_sustained.rds")
stan_rds <- file.path(dir_debug, "stan_sustained.rds")

# key shareable outputs (aggregate tables and figures): dir_transfer, if it has
# been set for a real run (see run_4_positive_deviance.R). If it has not - a
# simulated or exploratory run - fall back to a local folder under dir_debug, so
# nothing here ever requires the S: drive to exist just to test the code.
if (is.null(dir_transfer)) {
  warning("dir_transfer is not set - writing this stage's results to a local ",
          "folder instead of the S: transfer area. Set dir_transfer (see ",
          "run_4_positive_deviance.R) before sourcing for a real run.",
          call. = FALSE)
  out_dir <- file.path(dir_debug, "positive_deviance_test_output")
} else {
  out_dir <- dir_transfer
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

flow_csv       <- file.path(out_dir, "pd_flow.csv")            # patient attrition
hosp_excl_csv  <- file.path(out_dir, "pd_hospitals_excluded.csv")

# analysis window ------------------------------------------------------------
# the diagnosis-date window is NOT set here. sotn_cohort (applied in 02) already
# carries the SoTN audit period as part of its own definition, so the analysis
# window is whatever period sotn_cohort represents - no separate window filter
# is applied on top of it. A future, different reporting period (e.g. Jan 2023
# to Mar 2025) would replace the sotn_cohort restriction, not add to it; see the
# commented-out block in 02_build_cohort.R where that would go.

# inclusion: which treatment pathways count as curative --------------------
# The tx_pathway values (from the merge stage) that define the curative cohort.
# Resection, neoadjuvant-then-surgery, adjuvant, and the non-surgical definitive
# pathways are in; EMR/ESD only is included as a curative endoscopic pathway.
# "EMR/ESD then surgery" is held out for now (see pathways_flagged) pending
# clinical confirmation of what that pathway represents. Everything else
# (palliative, SACT only, no treatment) is excluded as non-curative.
pathways_include <- c(
  "Surgery only",
  "Surgery + neoadjuvant chemo",
  "Surgery + neoadjuvant RT",
  "Surgery + neoadjuvant chemoRT",
  "Surgery + adjuvant chemo",
  "Definitive chemoRT",
  "Curative RT only",
  "EMR/ESD only")

# flagged, and excluded until the pathway is clarified. Reported in the attrition
# so the count is visible rather than silently dropped.
pathways_flagged <- c("EMR/ESD then surgery")

# tumour site restriction ----------------------------------------------------
# the analysis is oesophageal only: keep C15, drop C16 (gastric). tumour_site is
# the 3-character ICD-10 site in the registry. Set include_sites to c("C15","C16")
# to analyse both, or c("C16") for gastric.
include_sites <- c("C15")

# hospital unit and volume floor ---------------------------------------------
# the unit of analysis is the site of diagnosis - the 5-character site code
# derived in the site-build stage. Its operating trust is resolved from the
# ODS-derived site->trust map (site_trust_map.csv, built in the reference stage):
# parent_trust is what ODS says operates the site, which differs from the site's
# own first three characters where ODS relocated it. A site must have at least
# min_per_year diagnoses in EVERY calendar year of the window to be analysed;
# smaller sites are excluded (and listed in pd_hospitals_excluded.csv).
hosp_var       <- "site_dx_code"
min_per_year   <- 5

# the ODS site->trust map from the reference stage. Keyed on site_code, giving
# parent_trust (the operating trust). If absent, the trust falls back to the
# site's first three characters, with a warning.
site_trust_map_csv <- file.path(dir_ref, "site_trust_map.csv")

# valid diagnosing trusts ----------------------------------------------------
# only sites whose operating trust is a recognised OG-cancer diagnoser are
# analysed. The list is the NOGCA State of the Nation data-quality table (English
# NHS trusts with recorded OG diagnoses), read from valid_trusts_csv - a file
# with a trust_code column, prepared from the SoTN workbook. Kept in dir_ref
# alongside the other reference lookups (it is a general lookup, not specific to
# a single analysis run). If the file is absent the restriction is skipped with
# a warning, so the pipeline still runs.
valid_trusts_csv <- file.path(dir_ref, "valid_diagnosing_trusts.csv")

# outcome --------------------------------------------------------------------
# time from diagnostic endoscopy to the decision to treat, in days. The merge
# stage anchors this strictly on the diagnostic endoscopy; patients with no
# endoscopy date have no endoscopy-anchored wait and are excluded (and counted)
# in 02. Only positive (> 0) waits are analysed, and waits beyond max_wait are
# set aside as implausible.
outcome_var    <- "wt_endo_to_dtt"
max_wait       <- 180
drop_zero_wait <- TRUE

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
# ----------------------------------------------------------------------------
# Helper functions (standardisation, weighting, shrinkage, ranking) live in
# _helpers.R, sourced at the top of this file. Kept separate so this file reads
# as configuration only.