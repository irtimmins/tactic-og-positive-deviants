# =============================================================================
# 00  Prepare a simulated run
# -----------------------------------------------------------------------------
# Called by run_master(mode = "simulated"). Builds every input the real stages
# read, then points every pipeline path at Output/sim so a dry run's files never
# land in the real data folders, in the repo's Data/ tree, or in a temporary
# folder that disappears (or is awkward to find) once R exits - these files can
# be large, so they get one clearly-named, keepable home:
#
#   Output/sim/raw            the simulated raw extracts and reference maps
#                              (stands in for dir_raw and dir_ref)
#   Output/sim/patient_level  the patient-level intermediates every stage
#                              writes (stands in for dir_out - the big files:
#                              og_cohort_cwt.rds, pd_cohort.rds, hes_extract.rds)
#   Output/sim/intermediates  the aggregate, non-patient outputs (stands in for
#                              dir_debug - table 1, the flowchart, Stan fits,
#                              the caterpillar plot)
#
# After this returns, run_master runs the ordinary stage code against these
# inputs, file for file the same as a real run.
#
# The one structural difference from a real run: the simulated rapid extract is
# enriched with treatment, stage and covariate columns (and the audit cohort is
# shaped into a clean two-year window over the busier sites) BEFORE the stages
# run, because the raw simulators only make the columns their own stage needs.
# On real data those columns are simply present already.
#
# Nothing here writes to dir_transfer, ever. Output/sim is not the S: transfer
# area and is never meant to be - it is a local, disposable folder; safe to
# delete and rebuild by re-running with mode = "simulated".
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
})

source("R/simulate_test_data/_cohort_extras.R")

prepare_simulated_run <- function(scale = 1, seed = 20260212,
                                  out_root = "Output/sim") {
  
  assign("sim_scale", scale, envir = globalenv())
  assign("sim_seed",  seed,  envir = globalenv())
  
  dir_sim  <- file.path(out_root, "raw")
  dir_pt   <- file.path(out_root, "patient_level")
  dir_agg  <- file.path(out_root, "intermediates")
  assign("dir_sim", dir_sim, envir = globalenv())
  for (d in c(dir_sim, dir_pt, dir_agg))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  
  # point every pipeline path at Output/sim, globally, so the stage scripts
  # (which source R/config/directories.R plainly) leave them alone. dir_debug is
  # included here on purpose - without it, the aggregate outputs (table 1, the
  # flowchart, the Stan fits) fall back to directories.R's ordinary default and
  # land in Data/intermediates instead, splitting a run's output across two
  # unrelated folders.
  assign("allow_repo_relative_out", TRUE, envir = globalenv())
  assign("dir_raw",      dir_sim, envir = globalenv())
  assign("dir_out",      dir_pt,  envir = globalenv())
  assign("dir_ref",      dir_sim, envir = globalenv())
  assign("dir_debug",    dir_agg, envir = globalenv())
  assign("dir_transfer", NULL,    envir = globalenv())
  source("R/config/directories.R")   # re-derive the f_* paths from these
  
  cat("Simulated run output root:", normalizePath(out_root, mustWork = FALSE), "\n")
  
  # -- 1. the raw COSD-stage inputs (rapid + cosd extracts, reference maps) ---
  cat("\n-- simulating the rapid and COSD extracts --\n")
  sys.source("R/simulate_test_data/03b_simulate_cosd_inputs.R", envir = new.env())
  
  # -- 2. enrich the simulated rapid with everything later stages read --------
  # treatment dates and intents, stage, covariates, and the endoscopy fields;
  # then shape the audit cohort into a clean recent two-year window over the
  # busier trusts so the volume floor is clearable and there is real output to
  # inspect. The site build passes unknown columns through, so enriching the
  # rapid here is enough for every stage after it.
  f_rapid <- path_rapid_dta
  rapid <- read_dta(f_rapid)
  rapid <- add_cohort_extras(rapid, seed = seed)
  
  set.seed(seed)
  trust_col <- intersect(c("diag_trust", "true_trust", "diagnosis_trust"),
                         names(rapid))[1]
  busy <- rapid %>% count(.data[[trust_col]], sort = TRUE) %>%
    slice_head(n = 40) %>% pull(1)
  in_window <- rapid[[trust_col]] %in% busy
  rapid$sotn_cohort <- 0
  rapid$sotn_cohort[in_window] <- as.numeric(runif(sum(in_window)) < 0.85)
  
  sotn <- rapid$sotn_cohort == 1
  n_sotn <- sum(sotn)
  new_dx <- as.Date("2023-01-01") + sample(0:729, n_sotn, TRUE)
  shift  <- as.integer(new_dx - as.Date(rapid$diagnosisdate[sotn]))
  shift_date <- function(v) { v[sotn] <- as.Date(v[sotn]) + shift; v }
  for (col in c("diagnosisdate", "endodateHES", "surgery_date",
                "sact_first_date", "rt_first_date", "EMR_ESDdateHES"))
    rapid[[col]] <- shift_date(rapid[[col]])
  rapid$diagnosis_year[sotn] <-
    as.integer(format(as.Date(rapid$diagnosisdate[sotn]), "%Y"))
  
  write_dta(rapid, f_rapid)
  cat(sprintf("enriched the simulated rapid: %d-patient audit cohort over 2023-2024 across %d trusts\n",
              n_sotn, length(busy)))
  
  # -- 3. the CWT and HES extracts, keyed on the same patients ----------------
  cwt <- make_cwt_extract(rapid, seed = seed)
  write_dta(cwt, path_cwt_dta)
  cat(sprintf("built a CWT extract: %d rows, %d patients\n",
              nrow(cwt), length(unique(cwt$patient_pseudo_id))))
  
  hes <- make_hes_extract(rapid, seed = seed)
  saveRDS(hes, f_hes_extract)   # stands in for stage 00's cohort-filtered cut
  cat(sprintf("built a HES-APC extract: %d episodes, %d patients\n",
              nrow(hes), length(unique(hes$patient_pseudo_id))))
  
  invisible(NULL)
}