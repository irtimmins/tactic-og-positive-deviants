# 02  directly standardised hospital waiting times via balancing weights
# -----------------------------------------------------------------------------
# For each diagnosing site, reweight its patients to the overall patient mix
# (balancer::standardize) so the standardised mean endoscopy-to-DTT time is
# comparable across sites. A pooled prognostic model supplies residuals for the
# augmented estimate and the population-mean prediction. Adjusts for age,
# comorbidity, season and calendar period. Writes the per-site estimates that
# the shrinkage step reads, for both estimands:
#
#   SUSTAINED  the level over the whole window, standardised to the whole
#              window's patient mix.
#   IMPROVED   the change from period 1 to period 2, both arms standardised to
#              period 1's patient mix.
#
# Reads : pd_cohort.rds
# Writes: fit_primary.rds, site_sustained.rds, site_improve.rds,
#         lambda_tradeoff.csv

library(balancer)
library(dplyr)

source("R/06_identify_positive_deviants/01_config.R")
df <- readRDS(cohort_rds)

# the two periods ------------------------------------------------------------
# Re-derived here from the explicit dates in 01_config.R, overwriting the period
# stage 05 carried. Stage 05 splits at the midpoint of the observed endoscopy
# dates, which moves with the extract and does not fall on a year boundary;
# these dates do not move. yr_late goes with it, so the sustained model's
# calendar term and the improvement split always mean the same thing.
if (!"endoscopy_date" %in% names(df))
  stop("pd_cohort.rds has no endoscopy_date, so the periods cannot be assigned.",
       call. = FALSE)

df <- df %>%
  mutate(endoscopy_date = as.Date(endoscopy_date),
         period = factor(case_when(
           endoscopy_date >= period_1_start & endoscopy_date <= period_1_end ~ "first",
           endoscopy_date >= period_2_start & endoscopy_date <= period_2_end ~ "second",
           TRUE                                                              ~ NA_character_),
           levels = c("first", "second")),
         yr_late = as.integer(period == "second"))

n_outside <- sum(is.na(df$period))
cat(sprintf("periods: %s (%s to %s) and %s (%s to %s)\n",
            period_1_label, period_1_start, period_1_end,
            period_2_label, period_2_start, period_2_end))
print(table(period = df$period, useNA = "ifany"))
if (n_outside)
  cat(sprintf("  %d patient(s) sit outside both periods: they contribute to the\n",
              n_outside),
      "  sustained estimand but to neither arm of the improvement one.\n", sep = "")

# yr_late is a covariate in the sustained model, so it cannot be missing there.
# A patient outside both periods has no later/earlier status; treat them as not
# in the later period rather than dropping them from the sustained analysis.
df$yr_late[is.na(df$yr_late)] <- 0L

# balance vs effective-sample-size trade-off across lambda -------------------
balance_tradeoff <- function(d, cont, bin, grid = lambda_grid) {
  d  <- d %>% arrange(hosp)
  X  <- make_std_matrix(d, cont, bin)
  Z  <- d$hosp
  pm <- lm(as.formula(paste("y_std ~", paste(colnames(X), collapse = " + "))),
           data = data.frame(y_std = d$wait, X))
  beta <- coef(pm)[-1]
  beta[is.na(beta)] <- 0
  
  hosp_means <- rowsum(X, Z) / as.numeric(table(Z))
  unw <- as.numeric(abs((hosp_means %*% beta)))
  
  res <- t(sapply(grid, function(l) {
    so <- standardize(X, rep(0, ncol(X)), Z, lambda = l, exact_global = FALSE)
    w  <- extract_weights(so)
    wm <- (t(so$weights) %*% X)
    wt <- as.numeric(abs(wm %*% beta))
    ne <- tapply(w, Z, function(x) sum(x)^2 / sum(x^2))
    c(lambda = l,
      bias_removed = 1 - mean(wt) / mean(unw),
      mean_eff_n   = mean(ne),
      mean_deff    = mean(ne / as.numeric(table(Z))))
  }))
  as.data.frame(res)
}

# primary standardisation: age + cci + season + calendar year ----------------
cv          <- code_covariates(df)
primary_bin <- c(cv$bin, season_terms, year_term)
trade       <- balance_tradeoff(cv$data, cv$cont, primary_bin)
cat("balance vs effective n by lambda:\n"); print(round(trade, 3))
write.csv(trade, file.path(out_dir, "lambda_tradeoff.csv"), row.names = FALSE)

fit_main <- run_standardise(patient_data          = cv$data,
                            continuous_covariates = cv$cont,
                            binary_covariates     = primary_bin,
                            lambda                = lambda_main)
site_sustained <- fit_main$site
saveRDS(fit_main,       file.path(out_dir, "fit_primary.rds"))
saveRDS(site_sustained, file.path(out_dir, "site_sustained.rds"))

cat(sprintf("\n%d sites standardised. Mean effective n per site %.1f (raw %.1f).\n",
            nrow(site_sustained), mean(site_sustained$n_eff), mean(site_sustained$n)))

# lambda trade-off figure, if ggplot2 is available ---------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  trade_curve <- trade %>%
    mutate(pct_bias_removed = 100 * bias_removed,
           is_main = abs(lambda - lambda_main) < 1e-9)
  p_trade <- ggplot(trade_curve, aes(mean_eff_n, pct_bias_removed)) +
    geom_path(colour = "grey60", linewidth = 0.6) +
    geom_point(aes(colour = is_main, size = is_main)) +
    scale_colour_manual(values = c(`FALSE` = "darkblue", `TRUE` = "firebrick"),
                        guide = "none") +
    scale_size_manual(values = c(`FALSE` = 2.8, `TRUE` = 3.8), guide = "none") +
    labs(x = "Average effective sample size per hospital",
         y = "Average percentage bias reduction (%)") +
    theme_classic(base_size = 13)
  ggsave(file.path(out_dir, "lambda_tradeoff.png"), p_trade,
         width = 120, height = 95, units = "mm", dpi = 300, bg = "white")
  cat("lambda trade-off figure written to lambda_tradeoff.png\n")
} else {
  cat("ggplot2 not installed - trade-off figure skipped (CSV written).\n")
}

# improvement estimand -------------------------------------------------------
# The change from period 1 to period 2, both arms standardised to period 1's
# case-mix. Sites without min_per_period patients in BOTH periods cannot have a
# change estimated and are dropped from this estimand only - they keep their
# sustained estimate. That loss is reported here rather than showing up later as
# a shorter caterpillar with no explanation.
site_improve <- standardise_change(patient_data          = cv$data,
                                   continuous_covariates = cv$cont,
                                   binary_covariates     = primary_bin)

if (is.null(site_improve) || !nrow(site_improve)) {
  cat("\nNo site has at least", min_per_period,
      "patients in both periods, so the improvement estimand is not estimable.\n",
      "Only the sustained results will be produced. Check the period dates in\n",
      "01_config.R against the cohort's date range.\n")
  if (file.exists(site_improve_rds)) file.remove(site_improve_rds)
} else {
  saveRDS(site_improve, site_improve_rds)
  n_dropped <- nrow(site_sustained) - nrow(site_improve)
  cat(sprintf(paste0("\nimprovement estimand: %d of %d sites have >= %d patients ",
                     "in both periods\n"),
              nrow(site_improve), nrow(site_sustained), min_per_period))
  if (n_dropped > 0)
    cat(sprintf("  %d site(s) too small in one period - sustained only\n", n_dropped))
  cat(sprintf("  mean change %.1f days (negative = faster in %s)\n",
              mean(site_improve$delta), period_2_label))
  cat(sprintf("  %d site(s) faster in %s, %d slower\n",
              sum(site_improve$delta < 0), period_2_label,
              sum(site_improve$delta > 0)))
}

cat("02 complete. Next: 03 shrinks the per-site estimates.\n")