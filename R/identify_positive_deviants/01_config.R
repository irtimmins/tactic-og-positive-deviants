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
# in_rds points at the CWT-merged cohort. out_dir is where every result from this
# stage is written; change it to send results elsewhere. stan_file is the
# normal-normal shrinkage model, kept alongside these scripts.
if (!exists("base_dir")) base_dir <- "Data/OG"
in_rds  <- file.path(base_dir, "og_cohort_cwt.rds")

# reference lookups built once in the reference stage (site->trust map, etc.)
if (!exists("dir_ref")) dir_ref <- "Data/reference"

# out_dir: all outputs (tables, figures, intermediate rds) land here. Set it
# before sourcing this file to redirect results without editing the script.
if (!exists("out_dir")) out_dir <- "Output/positive_deviance"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

stan_dir  <- "R/identify_positive_deviants"
stan_file <- file.path(stan_dir, "dp_normal_cont.stan")

# hand-off files between scripts
cohort_rds     <- file.path(out_dir, "pd_cohort.rds")          # 02 -> 03/04...
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
# stage anchors this on endoscopy where present and diagnosis otherwise; only
# positive (> 0) waits are analysed, and waits beyond max_wait are set aside as
# implausible.
outcome_var    <- "wt_anchor_to_dtt"
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

# ============================================================================
# PART B  -  helper functions (no need to edit)
# ============================================================================
# Carried across, unchanged in substance, from the colon positive-deviance
# analysis; each is written to be read straight through and checked by hand.

# --- standardising covariates -----------------------------------------------
# a continuous covariate, shifted to mean 0 and scaled to sd 1
z_std <- function(x) (x - mean(x)) / sd(x)

# a 0/1 covariate, centred on its proportion p and scaled by sqrt(p(1-p)). The
# 0.05 floor applies to the scale only, so a very rare category is not blown up
# enough to swamp the balance; the centre stays the true proportion.
bin_std <- function(x) {
  p     <- mean(x)
  scale <- sqrt(max(p, 0.05) * (1 - max(p, 0.05)))
  (x - p) / scale
}

# build the standardised covariate matrix
make_std_matrix <- function(df, cont, bin) {
  std <- df[c(cont, bin)]
  for (v in cont) std[[v]] <- z_std(df[[v]])
  for (v in bin)  std[[v]] <- bin_std(df[[v]])
  as.matrix(std)
}

# the mean and sd of each continuous covariate, and the proportion of each binary
# one, in a chosen reference population.
ref_moments <- function(df, cont, bin) {
  cont_mean <- numeric(0); cont_sd <- numeric(0); bin_p <- numeric(0)
  for (v in cont) {
    cont_mean[v] <- mean(df[[v]])
    cont_sd[v]   <- sd(df[[v]])
  }
  for (v in bin) bin_p[v] <- mean(df[[v]])
  list(cont_mean = cont_mean, cont_sd = cont_sd, bin_p = bin_p)
}

# standardised matrix centred and scaled to externally supplied reference moments
make_std_matrix_ref <- function(df, cont, bin, ref) {
  std <- df[c(cont, bin)]
  for (v in cont) {
    std[[v]] <- (df[[v]] - ref$cont_mean[v]) / ref$cont_sd[v]
  }
  for (v in bin) {
    p        <- ref$bin_p[v]
    scale    <- sqrt(max(p, 0.05) * (1 - max(p, 0.05)))
    std[[v]] <- (df[[v]] - p) / scale
  }
  as.matrix(std)
}

# pull the per-patient weight out of a balancer standardize() result
extract_weights <- function(std_out) apply(std_out$weights, 1, max)

# weighted quantile, robust in small or degenerate cells
w_quantile <- function(x, w, p = 0.5) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x  <- x[ok]; w <- w[ok]
  if (length(x) == 0) return(NA_real_)
  if (length(unique(x)) == 1) return(x[1])
  o  <- order(x); x <- x[o]; w <- w[o]
  cum_w <- cumsum(w) / sum(w)
  keep  <- !duplicated(cum_w, fromLast = TRUE)
  cum_w <- cum_w[keep]; x <- x[keep]
  if (length(x) < 2) return(x[length(x)])
  approx(cum_w, x, xout = p, rule = 2, ties = "ordered")$y
}

# per-hospital standardised summaries for a continuous outcome
site_summary <- function(df) {
  library(dplyr)
  s <- df %>%
    group_by(hosp) %>%
    summarise(
      n          = n(),
      n_eff      = sum(w)^2 / sum(w^2),
      raw_mean   = mean(y),
      raw_median = median(y),
      stand      = weighted.mean(y, w),
      stand_med  = w_quantile(y, w, 0.5),
      stand_adj  = weighted.mean(resid, w) + mean(canonical),
      sd_w       = sqrt(sum(w^2 * (y - weighted.mean(y, w))^2) / sum(w^2)),
      sd_adj     = sqrt(sum(w^2 * resid^2) / sum(w^2)),
      .groups = "drop"
    )
  sd_pool     <- sqrt(weighted.mean(s$sd_w^2,   s$n_eff))
  sd_pool_adj <- sqrt(weighted.mean(s$sd_adj^2, s$n_eff))
  s %>% mutate(
    se          = sd_w        / sqrt(n_eff),
    se_pool     = sd_pool     / sqrt(n_eff),
    se_adj      = sd_adj      / sqrt(n_eff),
    se_adj_pool = sd_pool_adj / sqrt(n_eff),
    sd_pool     = sd_pool,
    sd_pool_adj = sd_pool_adj
  )
}

# Directly standardise hospital waiting times with balancing weights. balancer
# must be loaded by the calling script.
run_standardise <- function(patient_data,
                            continuous_covariates,
                            binary_covariates,
                            lambda             = lambda_main,
                            reference_moments  = NULL,
                            target_population  = NULL,
                            weight_upper_limit = NULL) {
  
  data <- patient_data %>% arrange(hosp)
  
  if (is.null(reference_moments)) {
    reference_moments <- ref_moments(data, continuous_covariates, binary_covariates)
  }
  
  balance_matrix <- make_std_matrix_ref(data, continuous_covariates,
                                        binary_covariates, reference_moments)
  
  covariate_varies <- apply(balance_matrix, 2,
                            function(column) { s <- sd(column); is.finite(s) && s > 0 })
  balance_matrix  <- balance_matrix[, covariate_varies, drop = FALSE]
  covariates_used <- colnames(balance_matrix)
  
  standardize_args <- list(X = balance_matrix, target = rep(0, ncol(balance_matrix)),
                           Z = data$hosp, lambda = lambda, exact_global = FALSE)
  if (!is.null(weight_upper_limit)) standardize_args$uplim <- weight_upper_limit
  weight_solution <- do.call(standardize, standardize_args)
  data$w <- extract_weights(weight_solution)
  
  data$y     <- data$wait
  outcome_model <- lm(reformulate(covariates_used, "wait"), data = data)
  data$resid <- resid(outcome_model)
  data$canonical <- if (is.null(target_population)) {
    mean(fitted(outcome_model))
  } else {
    mean(predict(outcome_model, newdata = target_population))
  }
  
  site <- site_summary(data) %>%
    left_join(distinct(data, hosp, diag_hosp = diag_hosp_canon), by = "hosp")
  
  list(site = site, data = data, weights = weight_solution, model = outcome_model)
}

# --- posterior ranking ------------------------------------------------------
# posterior ranking metrics from a draws-by-hospital matrix of latent means.
# shorter wait is better, so rank 1 is the best performer.
rank_metrics <- function(draws, tops = c(.05, .10, .20, .25, .50)) {
  J <- ncol(draws)
  ranks <- t(apply(draws, 1, rank, ties.method = "average"))
  out <- data.frame(
    exp_rank = colMeans(ranks),
    rank_lo  = apply(ranks, 2, quantile, 0.025),
    rank_hi  = apply(ranks, 2, quantile, 0.975)
  )
  for (p in tops) {
    threshold <- ceiling(p * J)
    out[[paste0("p_top", p * 100)]] <- colMeans(ranks <= threshold)
  }
  out
}

# --- Bayesian shrinkage -----------------------------------------------------
# The normal-normal shrinkage model. Feed it per-hospital point estimates y and
# their standard errors se; it pulls each hospital towards the overall mean by an
# amount set by its precision. rstan must be loaded by the calling script.
fit_shrink <- function(y, se, mu_mean = mean(y),
                       mu_sd = prior_mu_sd, tau_scale = prior_tau_scale,
                       seed = 8675309, refresh = 2000, cores = 1,
                       adapt_delta = 0.95) {
  dat <- list(J = length(y), y_site_obs = y, sigma_site_obs = se,
              prior_mu_mean = mu_mean, prior_mu_sd = mu_sd,
              prior_tau_scale = tau_scale)
  rstan::stan(stan_file, data = dat, seed = seed,
              chains = 4, iter = 4000, warmup = 2000, refresh = refresh, cores = cores,
              control = list(adapt_delta = adapt_delta, max_treedepth = 12))
}

# --- hospital display names -------------------------------------------------
# Title Case a name with small joining words left lower-case (except the first
# word) and NHS kept upper-case.
name_small_words <- c("and", "of", "the", "for", "in", "on", "at", "to", "by", "an", "a", "or")
name_acronyms    <- c("nhs")
title_case <- function(x) {
  vapply(x, function(one) {
    words <- strsplit(tolower(one), " ")[[1]]
    for (i in seq_along(words)) {
      if (words[i] %in% name_acronyms) {
        words[i] <- toupper(words[i])
      } else if (!(words[i] %in% name_small_words) || i == 1) {
        substr(words[i], 1, 1) <- toupper(substr(words[i], 1, 1))
      }
    }
    paste(words, collapse = " ")
  }, character(1), USE.NAMES = FALSE)
}