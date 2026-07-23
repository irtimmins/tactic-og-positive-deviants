# =============================================================================
# _helpers.R  -  functions for the positive-deviance analysis
# -----------------------------------------------------------------------------
# The standardisation, weighting, shrinkage and ranking machinery, kept separate
# from the settings in 01_config.R so that file reads as configuration only.
# Carried across, unchanged in substance, from the colon positive-deviance
# analysis; each is written to be read straight through and checked by hand.
#
# Assumes _load_packages.R has been sourced first (these use dplyr; the Stan and
# balancer calls are made by the scripts that source this, which load those
# packages themselves).
# =============================================================================


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
# --- hospital display names --------------------------------------------------
# A named character vector mapping five-character site code -> hospital name, for
# the ranking tables. The names come from the ODS site-to-trust map built in
# stage 01; ODS returns them in capitals, so title_case (R/shared/utils.R) tidies
# them here once rather than at each place they are printed.
#
# Returns NULL when the map is absent or has no usable names, and the callers
# fall back to showing the site code alone - a missing lookup should cost a
# nicer label, not the whole table.
hospital_names <- function(path = f_site_trust_map) {
  if (!file.exists(path)) return(NULL)
  map <- try(read.csv(path, colClasses = "character"), silent = TRUE)
  if (inherits(map, "try-error") ||
      !all(c("site_code", "name") %in% names(map))) return(NULL)
  map <- map[!is.na(map$site_code) & nzchar(map$site_code) &
               !is.na(map$name) & nzchar(map$name), c("site_code", "name"), drop = FALSE]
  map <- map[!duplicated(map$site_code), , drop = FALSE]
  if (!nrow(map)) return(NULL)
  setNames(title_case(map$name), map$site_code)
}

# --- the improvement estimand ------------------------------------------------
# For each site, the difference between its directly-standardised mean waiting
# time in period 2 and in period 1. The periods are the explicit date ranges set
# in 01_config.R, and a patient's period is assigned in 02; nothing is derived
# here.
#
# Both arms are standardised to the SAME fixed population - period 1's case-mix
# - which is what makes the difference a change in speed rather than a change in
# who was diagnosed. If a site's patients got older between the periods, an
# unstandardised comparison would read that as getting slower.
#
# Returns one row per site with the change (delta; NEGATIVE means faster in
# period 2, i.e. improved), its standard error, and the two period counts. Sites
# without min_patients_per_period in both periods are dropped: a change cannot be
# estimated from one arm. A covariate that is constant within a period (the
# later-period indicator) carries no information there and is dropped by
# run_standardise, so passing season + year here amounts to season only.
standardise_change <- function(patient_data,
                               continuous_covariates,
                               binary_covariates,
                               lambda = lambda_main,
                               min_patients_per_period = min_per_period) {
  
  patients_per_period <- patient_data %>%
    count(hosp, period) %>%
    tidyr::pivot_wider(names_from = period, values_from = n, values_fill = 0)
  for (nm in c("first", "second"))
    if (!nm %in% names(patients_per_period)) patients_per_period[[nm]] <- 0L
  
  hospitals_kept <- patients_per_period %>%
    filter(first >= min_patients_per_period, second >= min_patients_per_period) %>%
    pull(hosp)
  if (!length(hospitals_kept)) return(NULL)
  
  period_1_data <- patient_data %>% filter(period == "first",  hosp %in% hospitals_kept)
  period_2_data <- patient_data %>% filter(period == "second", hosp %in% hospitals_kept)
  
  # the fixed target both arms are standardised to: period 1's case-mix
  baseline_moments <- ref_moments(period_1_data, continuous_covariates,
                                  binary_covariates)
  
  earlier <- run_standardise(patient_data          = period_1_data,
                             continuous_covariates = continuous_covariates,
                             binary_covariates     = binary_covariates,
                             lambda                = lambda,
                             reference_moments     = baseline_moments,
                             target_population     = period_1_data)$site
  later   <- run_standardise(patient_data          = period_2_data,
                             continuous_covariates = continuous_covariates,
                             binary_covariates     = binary_covariates,
                             lambda                = lambda,
                             reference_moments     = baseline_moments,
                             target_population     = period_1_data)$site
  
  earlier %>%
    select(hosp, diag_hosp, mean_earlier = stand_adj, se_earlier = se_adj_pool,
           n1 = n, ess1 = n_eff) %>%
    inner_join(later %>% select(hosp, mean_later = stand_adj,
                                se_later = se_adj_pool, n2 = n, ess2 = n_eff),
               by = "hosp") %>%
    mutate(delta    = mean_later - mean_earlier,
           se_delta = sqrt(se_earlier^2 + se_later^2))
}

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
# title_case() (and its word lists) now live in R/shared/utils.R, since it is
# a generic name-formatting helper with no connection to modelling - stage 05
# needs it too, for the hospital and trust names in pd_hospital_counts.csv.