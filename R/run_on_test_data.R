# =============================================================================
# run_on_test_data.R  -  the whole pipeline, end to end, on one simulated cohort
# -----------------------------------------------------------------------------
# WHERE: anywhere - a laptop, the internet machine, the server. Touches no real
#        data and needs no internet.
#
# This is the realistic dry run: it makes ONE simulated cohort in the shape of
# the real extracts, then runs all three real stages against it by the same file
# handoff the real pipeline uses -
#
#   1. site of diagnosis     (derive_hospital_code_from_cosd)   -> og_cohort_site.rds
#   2. CWT merge / waiting    (merge_cwt_to_get_dtt)       -> og_cohort_cwt.rds
#   3. positive deviants      (identify_positive_deviants) -> pd_cohort.rds, ...
#
# It is the dry-run counterpart to run_master.R: same stages, same handoffs, but
# on made-up data, so it answers "does the whole thing actually run on this
# machine, and what do the outputs look like?" without waiting on the server or
# touching real data. run_to_test_logic.R is the other half - it proves the
# stage LOGIC with pass/fail checks; this one produces real, inspectable output
# files you can open, profile, or eyeball.
#
# Everything is written under Data/OG_test (well away from the real Data/OG), and
# dir_transfer is left NULL so made-up numbers never reach the real S: transfer
# area. Nothing here is a formal test; for the pass/fail counts run
# run_to_test_logic.R.
#
# Note on how the paths are set: the stage scripts source R/config/directories.R
# with a plain source(), which evaluates in the global environment - so the path
# overrides here are set globally (not in a local scope) before each stage is
# sourced, the same way run_master.R does it. They are reset before each stage so
# the whole run stays inside Data/sim and Data/OG_test.
#
# Run from the repository root:
#   Rscript run_on_test_data.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
})

sim_scale <- 1        # 1 = full size; set to 0.2 for a quick smaller run
sim_seed  <- 20260212

dir_sim   <- "Data/sim"        # the simulated raw extracts and reference maps
dir_test  <- "Data/OG_test"    # all stage outputs land here, not in Data/OG
dir.create(dir_sim,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_test, recursive = TRUE, showWarnings = FALSE)

source("R/shared/simulate_cohort_extras.R")

stage_banner <- function(txt) cat("\n\n==================== ", txt,
                                  " ====================\n", sep = "")

# point all the pipeline paths at the test folders, in the global environment,
# so R/config/directories.R (sourced plainly by each stage, hence evaluated
# globally) leaves them alone rather than filling in the real W:/S: defaults.
# dir_transfer stays NULL so nothing reaches the real transfer area.
use_test_paths <- function(out = dir_test) {
  # Data/OG_test is repo-relative on purpose (synthetic data, kept next to the
  # code so it is easy to inspect); opt into the guard in directories.R that
  # otherwise blocks a repo-relative dir_out for real patient data.
  assign("allow_repo_relative_out", TRUE, envir = globalenv())
  assign("dir_raw",      dir_sim, envir = globalenv())
  assign("dir_out",      out,     envir = globalenv())
  assign("dir_ref",      dir_sim, envir = globalenv())
  assign("dir_transfer", NULL,    envir = globalenv())
}

# ---------------------------------------------------------------------------
# 1. site of diagnosis
# ---------------------------------------------------------------------------
# The site-code simulator writes the raw RAPID and COSD extracts and the two
# reference maps into dir_sim; the real build then reads them and writes
# og_cohort_site.rds. That output has the registry and the derived site, but not
# the treatment columns the later stages need - those are added just below.
stage_banner("1. site of diagnosis")
source("R/derive_hospital_code_from_cosd/90_simulate_inputs.R")

use_test_paths()
source("R/derive_hospital_code_from_cosd/00_master.R")

# enrich og_cohort_site.rds with treatment dates, stage and covariates, so the
# CWT and positive-deviance stages have everything they read. This is what makes
# the single cohort flow all the way through rather than each stage needing its
# own separate fake data.
site_rds <- file.path(dir_test, "og_cohort_site.rds")
cohort <- add_cohort_extras(readRDS(site_rds), seed = sim_seed)

# shape the sotn_cohort into something the analysis funnel can actually yield a
# cohort from - purely a dry-run presentation choice, not part of any stage's
# logic. The raw simulator spreads diagnoses thinly over 2018-2025 across ~120
# trusts, so no site clears the ">= 5 in every year" floor and the analysis
# cohort would come out empty. Concentrate the audit cohort into a recent two
# calendar years (matching the real Jan 2023 - Dec 2024 shape) and into the
# busier sites, so the floor is clearable and there is real output to inspect.
set.seed(sim_seed)
busy_sites <- cohort %>%
  filter(!is.na(site_dx_code)) %>%
  count(site_dx_code, sort = TRUE) %>%
  slice_head(n = 40) %>% pull(site_dx_code)
in_window <- cohort$site_dx_code %in% busy_sites
cohort$sotn_cohort <- 0
cohort$sotn_cohort[in_window] <- as.numeric(runif(sum(in_window)) < 0.85)
# give the audit-cohort patients diagnosis (and so endoscopy) dates inside a
# clean two-year window. The treatment dates were set by add_cohort_extras as
# offsets from the ORIGINAL diagnosis date, so shift them by the same amount the
# diagnosis date moves - otherwise treatment would sit years off from the new
# endoscopy date and the endoscopy-to-DTT chain would be invalid.
sotn <- cohort$sotn_cohort == 1
n_sotn <- sum(sotn)
new_dx <- as.Date("2023-01-01") + sample(0:729, n_sotn, TRUE)
shift  <- as.integer(new_dx - as.Date(cohort$diagnosisdate[sotn]))
shift_date <- function(v) { v[sotn] <- as.Date(v[sotn]) + shift; v }
cohort$diagnosisdate   <- shift_date(cohort$diagnosisdate)
cohort$endodateHES     <- shift_date(cohort$endodateHES)
cohort$surgery_date    <- shift_date(cohort$surgery_date)
cohort$sact_first_date <- shift_date(cohort$sact_first_date)
cohort$rt_first_date   <- shift_date(cohort$rt_first_date)
cohort$EMR_ESDdateHES  <- shift_date(cohort$EMR_ESDdateHES)
cohort$diagnosis_year[sotn] <- as.integer(format(cohort$diagnosisdate[sotn], "%Y"))

saveRDS(cohort, site_rds)
cat(sprintf("enriched %s (treatment, stage, covariates); shaped a %d-patient audit cohort over 2023-2024 across %d sites\n",
            basename(site_rds), n_sotn, length(busy_sites)))

# ---------------------------------------------------------------------------
# 1b. endoscopy hospital from HES
# ---------------------------------------------------------------------------
# The endoscopy stage re-finds each patient's diagnostic endoscopy in HES-APC and
# reads its five-character site. As with CWT, the HES episodes have to belong to
# the SAME patients as the cohort, so they are built from it with
# make_hes_extract rather than from the stage's own simulator (which invents its
# own separate patients). It runs after the date shifting above, so the episodes
# sit against the cohort's final endoscopy dates.
#
# The stage normally reads the rapid dta and the cut-down HES extract from disk.
# Here it is handed the enriched cohort and the built extract directly through
# read_rapid / read_hes, which is the same injection the stage's own checks use -
# the raw simulated rapid dta has no endoscopy columns, they were added by
# add_cohort_extras.
stage_banner("1b. endoscopy hospital from HES")
hes_extract <- make_hes_extract(cohort, seed = sim_seed)
cat(sprintf("built a HES-APC extract of %d episodes for %d of the cohort's patients\n",
            nrow(hes_extract), length(unique(hes_extract$patient_pseudo_id))))

use_test_paths()
assign("read_rapid", function() cohort,      envir = globalenv())
assign("read_hes",   function() hes_extract, envir = globalenv())
source("R/derive_hospital_code_from_hes/02_add_endoscopy_site.R")
source("R/derive_hospital_code_from_hes/03_endoscopy_diagnostics.R")
# put the injections back, so the later stages read from disk as they normally do
rm(read_rapid, read_hes, envir = globalenv())

# ---------------------------------------------------------------------------
# 2. CWT merge and waiting times
# ---------------------------------------------------------------------------
# The CWT stage reads og_cohort_site.rds (the enriched cohort just written) and a
# CWT extract. The extract has to belong to the SAME patients as the cohort, or
# the merge matches no one - so build it from our cohort with make_cwt_extract,
# rather than from the CWT stage's own simulator (which invents its own separate
# patients). Patients with no CWT row simply get no decision-to-treat date, which
# is the real situation.
stage_banner("2. CWT merge and waiting times")
cwt_extract <- make_cwt_extract(cohort, seed = sim_seed)
haven::write_dta(cwt_extract,
                 file.path(dir_sim, "20260212_all_cwt_rapid_202601_OG.dta"))
cat(sprintf("built a CWT extract of %d rows for %d of the cohort's patients\n",
            nrow(cwt_extract), length(unique(cwt_extract$patient_pseudo_id))))

use_test_paths()
f_cwt <- file.path(dir_sim, "20260212_all_cwt_rapid_202601_OG.dta")
source("R/merge_cwt_to_get_dtt/02_derive_pathway.R")
source("R/merge_cwt_to_get_dtt/03_cwt_merge.R")

# ---------------------------------------------------------------------------
# 3. positive-deviance analysis
# ---------------------------------------------------------------------------
# Reads og_cohort_cwt.rds and builds the analysis cohort, table 1, the weights,
# the shrinkage and the rankings - the real scripts, on the simulated cohort.
# The modelling steps (04-06) need balancer and rstan; where those are not
# installed the build (02) and table 1 (03) still run, and the rest is skipped
# with a note, so this is useful with or without the modelling stack.
stage_banner("3. positive-deviance analysis")
use_test_paths()
source("R/identify_positive_deviants/02_build_cohort.R")
source("R/identify_positive_deviants/03_table1_characteristics.R")
if (requireNamespace("balancer", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE)) {
  source("R/identify_positive_deviants/04_estimation_weights.R")
  source("R/identify_positive_deviants/05_shrinkage.R")
  source("R/identify_positive_deviants/06_ranks_caterpillars.R")
} else {
  cat("\n(modelling steps 04-06 skipped - balancer/rstan not installed.",
      "The cohort and table 1 above are the real output; run this on a",
      "machine with the modelling stack for the full run.)\n")
}

cat("\n\nDry run complete. Simulated outputs are in", dir_test,
    "(not the real Data/OG).\n")
cat("For the formal pass/fail logic checks, run run_to_test_logic.R.\n")