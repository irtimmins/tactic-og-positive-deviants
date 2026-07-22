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