# 04  directly standardised hospital waiting times via balancing weights
# -----------------------------------------------------------------------------
# For each diagnosing site, reweight its patients to the overall patient mix
# (balancer::standardize) so the standardised mean endoscopy-to-DTT time is
# comparable across sites. A pooled prognostic model supplies residuals for the
# augmented estimate and the population-mean prediction. Adjusts for age,
# comorbidity, season and calendar-year half. Writes the per-site estimates that
# the shrinkage step reads.
#
# Reads : pd_cohort.rds
# Writes: fit_primary.rds, site_sustained.rds, lambda_tradeoff.csv

library(balancer)
library(dplyr)

source("R/identify_positive_deviants/01_config.R")
df <- readRDS(cohort_rds)

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

cat("04 complete. Next: 05 shrinks the per-site estimates.\n")
