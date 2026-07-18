# 05  Bayesian shrinkage of the per-site estimates
# -----------------------------------------------------------------------------
# Fit the normal-normal shrinkage model to the augmented weighted per-site means
# and their pooled standard errors. The model pulls each site towards the overall
# mean by an amount set by its precision, stabilising small-site estimates before
# ranking. Vague (weakly informative) priors, set in 01_config.R. The shrinkage
# routine itself (fit_shrink) is shared, in 01_config.R.
#
# Reads : site_sustained.rds
# Writes: stan_sustained.rds

library(rstan)
library(dplyr)

source("R/identify_positive_deviants/01_config.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

site <- readRDS(file.path(out_dir, "site_sustained.rds")) %>% arrange(hosp)
fit_sustained <- fit_shrink(site$stand_adj, site$se_adj_pool)
print(fit_sustained, pars = c("mu_true", "sigma_true"))
saveRDS(list(site = site, fit = fit_sustained),
        file.path(out_dir, "stan_sustained.rds"))

cat("05 complete. Next: 06 ranks the sites and draws the caterpillar plot.\n")
