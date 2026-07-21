# =============================================================================
# 91  Check the positive-deviance logic
# -----------------------------------------------------------------------------
# Proves the analysis machinery against small worked examples where the answer is
# known by hand, and the cohort build against simulated data. Split in two:
#
#   A. pure-logic checks - covariate coding, the weighted quantile, the
#      standardised per-site summary, rank metrics, title case, and the cohort
#      build's funnel (endoscopy exclusion, the hospital table, the volume
#      floor). These need only dplyr and run anywhere.
#   B. modelling checks - the balancing weights (balancer) and the shrinkage
#      (rstan). These need those packages and the compiled Stan model, so they
#      are guarded behind requireNamespace() and skipped, with a note, where the
#      packages are not installed. That keeps the script useful on a machine
#      without the modelling stack (where A still runs) as well as on the server.
#
# Run from the repository root:
#   Rscript R/identify_positive_deviants/91_check_positive_deviance.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

r_dir <- "R/identify_positive_deviants"

# load the settings and helpers without running an analysis. 01_config.R sources
# the shared paths file, which needs dir_out to look absolute; point everything
# at temporary folders so nothing real is touched.
assign("dir_out", tempfile("pd_check_out_"), envir = globalenv())
assign("dir_debug", tempfile("pd_check_dbg_"), envir = globalenv())
assign("dir_ref", tempdir(), envir = globalenv())
assign("dir_transfer", tempfile("pd_check_tr_"), envir = globalenv())
suppressWarnings(suppressMessages(source(file.path(r_dir, "01_config.R"))))

.checks <- new.env(); .checks$rows <- list()
expect <- function(label, cond) {
  ok <- isTRUE(cond)
  .checks$rows[[length(.checks$rows) + 1]] <- list(label = label, ok = ok)
  cat(if (ok) "  pass  " else "  FAIL  ", label, "\n")
}
near <- function(a, b, tol = 1e-8) is.finite(a) && is.finite(b) && abs(a - b) < tol

# ===========================================================================
cat("\nA. pure-logic checks\n")
# ===========================================================================

cat("\nCovariate coding\n")
d0 <- tibble(agediag = c(45, 55, 65, 75, 85),
             cci_n_conditions = c(0, 1, 2, 3, 4))
cc <- code_covariates(d0, age = "cont", cci = "0_1_2plus")
expect("continuous age is carried through as agediag",
       "agediag" %in% cc$cont && all(cc$data$agediag == d0$agediag))
expect("comorbidity 0/1/2+ makes two dummies",
       all(c("cci_1", "cci_2p") %in% cc$bin))
expect("the 2+ dummy captures 2 and above, not 1",
       all(cc$data$cci_2p == c(0, 0, 1, 1, 1)) &&
         all(cc$data$cci_1 == c(0, 1, 0, 0, 0)))

ccb <- code_covariates(d0, age = "band", cci = "0_1_2_3plus")
expect("age bands make the expected non-reference dummies",
       all(c("age_u50", "age_50_59", "age_70_79", "age_80p") %in% ccb$bin))
expect("60-69 is the omitted reference band (no dummy is set for age 65)",
       ccb$data$age_u50[3] == 0 && ccb$data$age_50_59[3] == 0 &&
         ccb$data$age_70_79[3] == 0 && ccb$data$age_80p[3] == 0)

cat("\nThe weighted quantile\n")
# w_quantile places the CDF at the upper edge of each point's weight and
# interpolates, so for 1..5 with equal weights the p=0.5 point falls halfway
# between the 2nd and 3rd values (2.5), not on the 3rd. That is the estimator's
# convention; check it behaves consistently rather than assuming a different one.
expect("equal weights interpolate to the expected weighted median",
       near(w_quantile(c(1, 2, 3, 4, 5), rep(1, 5), 0.5), 2.5))
expect("the weighted median is monotone in the requested probability",
       w_quantile(1:5, rep(1, 5), 0.25) < w_quantile(1:5, rep(1, 5), 0.75))
expect("a single distinct value is returned whatever the weights",
       near(w_quantile(c(7, 7, 7), c(2, 5, 1), 0.5), 7))
expect("weight concentrated on one value pulls the median to it",
       near(w_quantile(c(1, 100), c(1000, 1), 0.5), 1, tol = 1))

cat("\nThe standardised per-site summary\n")
# two sites, uniform weights, so the standardised mean is the raw mean and the
# summary maths can be checked by hand
site_df <- tibble(
  hosp = c(1, 1, 1, 2, 2, 2),
  y    = c(10, 20, 30, 40, 50, 60),
  w    = rep(1, 6),
  resid = c(-5, 0, 5, -5, 0, 5),
  canonical = rep(35, 6),
  diag_hosp_canon = c("A", "A", "A", "B", "B", "B"))
# site_summary expects a diag_hosp_canon-less frame joined later; it groups by
# hosp and uses y, w, resid, canonical
ss <- site_summary(site_df %>% mutate(diag_hosp_canon = NULL))
expect("the standardised mean equals the raw mean under uniform weights",
       near(ss$stand[ss$hosp == 1], 20) && near(ss$stand[ss$hosp == 2], 50))
expect("the effective n equals the count under uniform weights",
       near(ss$n_eff[ss$hosp == 1], 3) && near(ss$n_eff[ss$hosp == 2], 3))
expect("the adjusted mean adds the mean residual to the canonical mean",
       near(ss$stand_adj[ss$hosp == 1], 35) && near(ss$stand_adj[ss$hosp == 2], 35))

cat("\nRank metrics\n")
# three sites with cleanly separated latent means: site 1 fastest (smallest),
# site 3 slowest. Shorter wait is better, so site 1 should rank 1.
set.seed(1)
draws <- cbind(rnorm(2000, 10, 1), rnorm(2000, 20, 1), rnorm(2000, 30, 1))
rm_out <- rank_metrics(draws)
expect("the fastest site has the lowest expected rank",
       rm_out$exp_rank[1] < rm_out$exp_rank[2] &&
         rm_out$exp_rank[2] < rm_out$exp_rank[3])
expect("the fastest site is almost certainly in the top third",
       rm_out$p_top50[1] > 0.99)
expect("the slowest site is almost never in the top third",
       rm_out$p_top50[3] < 0.01)

cat("\nTitle case for hospital names\n")
expect("small joining words stay lower case, except the first",
       title_case("THE ROYAL MARSDEN") == "The Royal Marsden")
expect("NHS is kept upper case",
       title_case("guy's and st thomas nhs trust") == "Guy's and St Thomas NHS Trust")

cat("\nThe cohort build funnel\n")
# a tiny simulated cohort with known properties, run through 02, to prove the
# endoscopy exclusion and the hospital table behave. Build the input by hand so
# the expected drops are known.
set.seed(7)
mk <- function(n, site, endo = TRUE, stage = 1L, path = "Surgery only",
               wait = 30L) tibble(
                 patient_pseudo_id = paste0(site, "_", sample(1e6, n)),
                 site_dx_code = site, tumour_site = "C15", tx_pathway = path,
                 stage_RR = stage, diagnosis_year = sample(c(2023L, 2024L), n, TRUE),
                 endoscopy_date = if (endo) as.Date("2023-06-01") + sample(0:300, n, TRUE) else as.Date(NA),
                 wt_endo_to_dtt = if (endo) wait else NA_integer_,
                 age = 70L, gender = 1L, rcs_ch_score = 0L, ethnicity_grp = 1L,
                 imd_2019_RR = 3L, diagnosis_alliance = "Alliance A", sotn_cohort = 1L)

sim <- bind_rows(
  mk(40, "RXX01"),                       # big site, kept
  mk(40, "RXX02"),                       # second site in same trust, kept
  mk(20, "RYY01"),                       # single-site trust, kept
  mk(6,  "RZZ01", endo = FALSE),         # no endoscopy - all dropped early
  mk(4,  "RWW01"))                       # too small - excluded on the floor

# run 02 against it. 02 sources 01_config, which resets in_rds to
# file.path(dir_out, "og_cohort_cwt.rds"), so write the fixture there rather than
# rely on an in_rds override surviving. Point the reference lookups at absent
# files so trust resolution falls back to the code prefix and the valid-trust
# step is skipped - what we want for a self-contained logic check.
sim_dir <- tempfile("pd_check_data_"); dir.create(sim_dir)
saveRDS(sim, file.path(sim_dir, "og_cohort_cwt.rds"))
assign("dir_out", sim_dir, envir = globalenv())
assign("out_dir", tempfile("pd_out_"), envir = globalenv()); dir.create(out_dir)
assign("site_trust_map_csv", tempfile(), envir = globalenv())   # absent
assign("valid_trusts_csv",   tempfile(), envir = globalenv())   # absent
flow_env <- new.env()
invisible(capture.output(suppressMessages(
  sys.source(file.path(r_dir, "02_build_cohort.R"), flow_env))))
flow <- get("flow", envir = flow_env)
hosp_tab <- get("hosp_table", envir = flow_env)

endo_step <- flow$dropped[flow$step == "Diagnostic endoscopy date recorded"]
expect("the no-endoscopy site is dropped at the endoscopy step (6 patients)",
       endo_step == 6)
expect("the hospital table lists every site present before the floor",
       setequal(hosp_tab$hospital, c("RXX01", "RXX02", "RYY01", "RWW01")))
expect("multi-hospital trusts are sorted first in the table",
       hosp_tab$hospitals_in_trust[1] == 2)
expect("the hospital table shows full trust and hospital names",
       all(c("trust_name", "hospital_name") %in% names(hosp_tab)))
included_col <- grep("^Included in analysis", names(hosp_tab), value = TRUE)
expect("the small site is marked not included, the big ones are marked included",
       length(included_col) == 1 &&
         hosp_tab[[included_col]][hosp_tab$hospital == "RWW01"] == "No" &&
         hosp_tab[[included_col]][hosp_tab$hospital == "RXX01"] == "Yes")
final_hosps <- get("df", envir = flow_env)$diag_hosp_canon
expect("the too-small site is excluded from the final cohort",
       !("RWW01" %in% final_hosps) && all(c("RXX01", "RXX02", "RYY01") %in% final_hosps))

# ===========================================================================
cat("\nB. modelling checks (balancer, rstan)\n")
# ===========================================================================

if (requireNamespace("balancer", quietly = TRUE)) {
  suppressPackageStartupMessages(library(balancer))
  # an imbalanced two-site cohort: site 2 is older on average. Standardising
  # should move each site's weighted age mean towards the overall mean, so the
  # gap between the sites' weighted means is smaller than between their raw means.
  set.seed(3)
  n1 <- 200; n2 <- 200
  pd <- tibble(
    hosp = c(rep(1L, n1), rep(2L, n2)),
    diag_hosp_canon = c(rep("A", n1), rep("B", n2)),
    agediag = c(rnorm(n1, 65, 8), rnorm(n2, 75, 8)),
    cci_n_conditions = rpois(n1 + n2, 0.6),
    wait = c(rnorm(n1, 40, 10), rnorm(n2, 55, 10)))
  cv <- code_covariates(pd, age = "cont", cci = "0_1_2plus")
  fit <- run_standardise(cv$data, cv$cont, cv$bin, lambda = 0)
  
  raw_gap  <- abs(mean(pd$agediag[pd$hosp == 1]) - mean(pd$agediag[pd$hosp == 2]))
  wm <- fit$data %>% group_by(hosp) %>%
    summarise(m = weighted.mean(agediag, w), .groups = "drop")
  wt_gap <- abs(wm$m[1] - wm$m[2])
  expect("balancing weights shrink the between-site age gap",
         wt_gap < raw_gap)
  expect("every patient keeps a positive weight",
         all(fit$data$w > 0))
} else {
  cat("  (skipped - balancer not installed)\n")
}

if (requireNamespace("rstan", quietly = TRUE) && file.exists(stan_file)) {
  # shrinkage should pull an extreme, imprecise site towards the overall mean:
  # its posterior mean sits strictly between its noisy estimate and the grand
  # mean. Use few iterations - this is a behaviour check, not an inference.
  y  <- c(30, 40, 50, 90)          # site 4 is an outlier
  se <- c(2, 2, 2, 20)             # and imprecise
  fit <- suppressWarnings(fit_shrink(y, se, refresh = 0))
  post <- rstan::extract(fit, pars = "y_site_true")$y_site_true
  pm4 <- mean(post[, 4])
  grand <- mean(y)
  expect("shrinkage pulls the noisy outlier site towards the mean",
         pm4 < y[4] && pm4 > grand)
} else {
  cat("  (skipped - rstan not installed or Stan model absent)\n")
}

# ===========================================================================
res <- bind_rows(lapply(.checks$rows, as_tibble))
n_fail <- sum(!res$ok)
cat("\n", nrow(res), "checks,", n_fail, "failed\n")
if (n_fail) {
  cat("\nfailed:\n"); cat(paste0("  ", res$label[!res$ok], collapse = "\n"), "\n")
  if (!interactive()) quit(status = 1, save = "no") else stop(n_fail, " check(s) failed - see the output above.", call. = FALSE)
}
cat("All checks passed.\n")