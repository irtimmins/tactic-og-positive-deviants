# 03  Bayesian shrinkage of the per-site estimates
# -----------------------------------------------------------------------------
# Fit the normal-normal shrinkage model to the augmented weighted per-site means
# and their pooled standard errors. The model pulls each site towards the overall
# mean by an amount set by its precision, stabilising small-site estimates before
# ranking. Vague (weakly informative) priors, set in 01_config.R. The shrinkage
# routine itself (fit_shrink) is shared, in _helpers.R.
#
# One fit per estimand. They differ in where the prior is centred: the sustained
# means are centred on their own average, whereas the improvement is a change,
# for which the natural null is zero change - so mu_mean = 0 there.
#
# Reads : site_sustained.rds, site_improve.rds (if 02 produced one)
# Writes: stan_sustained.rds, stan_improve.rds

library(rstan)
library(dplyr)

source("R/06_identify_positive_deviants/01_config.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

site <- readRDS(file.path(out_dir, "site_sustained.rds")) %>% arrange(hosp)
fit_sustained <- fit_shrink(site$stand_adj, site$se_adj_pool)
print(fit_sustained, pars = c("mu_true", "sigma_true"))
saveRDS(list(site = site, fit = fit_sustained),
        file.path(out_dir, "stan_sustained.rds"))

# improvement ----------------------------------------------------------------
# mu_mean = 0: the prior sits on "no change", so a site is pulled towards having
# improved no more than average rather than towards the average level of waiting.
if (file.exists(site_improve_rds)) {
  imp <- readRDS(site_improve_rds) %>% arrange(hosp)
  fit_improve <- fit_shrink(imp$delta, imp$se_delta, mu_mean = 0)
  print(fit_improve, pars = c("mu_true", "sigma_true"))
  saveRDS(list(site = imp, fit = fit_improve), stan_improve_rds)
  cat(sprintf("improvement shrinkage fitted for %d sites\n", nrow(imp)))
} else {
  cat("no site_improve.rds - 02 found the improvement estimand not estimable,\n",
      "so only the sustained fit is produced.\n", sep = "")
  if (file.exists(stan_improve_rds)) file.remove(stan_improve_rds)
}

cat("03 complete. Next: 04 ranks the sites and draws the caterpillar plot.\n")