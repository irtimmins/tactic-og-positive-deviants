# =============================================================================
# run_master.R  -  the whole analysis, on the real data, on the server
# -----------------------------------------------------------------------------
# WHERE: the analysis server, on the real extracts.
# NEEDS: the 20260212 extracts under dir_raw, and the two reference maps plus the
#        valid-trust list under Data/reference (built on the internet machine by
#        run_1_internet.R and copied over). Does NOT need internet.
#
# Runs the three stages end to end:
#   1. derive the five-character site of diagnosis        (site-code build)
#   2. merge the CWT records and derive the waiting times (CWT merge)
#   3. identify the positive-deviant hospitals            (deviance analysis)
#
# Storage, all defined once in R/config/directories.R:
#   dir_out       patient-level data (og_cohort_site.rds, og_cohort_cwt.rds,
#                 pd_cohort.rds, ...) - the restricted, encrypted W: drive
#   dir_debug     aggregate, non-patient intermediates - in the repo
#   dir_transfer  the shareable tables and figures - set below to transfer_root,
#                 the S: results-transfer area, because the W: drive is encrypted
#                 and cannot be transferred off the server. This is set here so a
#                 real run's results actually leave the server; a test run
#                 (run_to_test_logic.R) leaves it NULL and writes to a local
#                 scratch folder instead.
#
# A full console log of the run (everything printed by every stage, in order,
# timestamped) is written alongside the results, as pipeline_log_<time>.txt in
# dir_transfer - so a look back at what happened during a given run does not
# depend on anyone having kept a terminal scrollback. Output still appears on
# screen as normal at the same time.
#
# To override a path on a particular machine, set it before the source() line
# below (e.g. dir_raw <- "..."), and R/config/directories.R will leave it alone.
# =============================================================================

# dir_raw <- "W:/_DATA/IainTimmins/2026 OG SOTN data"   # uncomment to override

source("R/config/directories.R")

# send this run's shareable results to the S: transfer area (see the note above).
dir_transfer <- transfer_root
cat("results for transfer will be written to:\n  ", dir_transfer, "\n\n")

# a small helper so a failure in one stage stops the run with a clear message
# rather than carrying a half-built dataset into the next stage.
run_stage <- function(label, files) {
  cat("\n==================== ", label, " ====================\n", sep = "")
  for (f in files) {
    cat("\n----- ", f, " -----\n", sep = "")
    source(f, local = new.env())
  }
}

# the whole run is wrapped in a function so the log file below is opened and
# closed reliably by on.exit() - which needs a real function call to attach to.
# At the top level of a script, on.exit() only behaves correctly when the script
# is run directly by Rscript; wrapping it here means the log still closes
# properly even if this file is ever sourced from elsewhere, and even if a stage
# below fails partway through (on.exit() fires on an error exit too, so the log
# is still flushed and closed, with whatever ran before the failure captured in
# it - not silently lost or left open).
run_full_pipeline <- function() {
  log_file <- NULL
  if (!is.null(dir_transfer)) {
    dir.create(dir_transfer, recursive = TRUE, showWarnings = FALSE)
    log_file <- file.path(dir_transfer,
                          paste0("pipeline_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
    log_con <- file(log_file, open = "wt")
    sink(log_con, split = TRUE)   # split = TRUE: still visible on screen too
    on.exit({ sink(); close(log_con) }, add = TRUE)
    cat("Pipeline run started:", format(Sys.time()), "\n")
  }
  
  run_stage("1. site of diagnosis", c(
    "R/derive_hospital_code_from_cosd/02_add_site_of_diagnosis.R",
    "R/derive_hospital_code_from_cosd/03_site_diagnostics.R"))
  
  # the endoscopy hospital comes from HES rather than COSD, and is read off the
  # rapid extract directly, so it does not sit in the site -> CWT -> deviance
  # handoff chain. It runs here, before the deviance stage, because that stage
  # joins the lookup it writes. 02a only re-reads the large raw HES file when the
  # cut-down extract is missing.
  run_stage("1b. endoscopy hospital from HES", c(
    "R/derive_hospital_code_from_hes/02a_extract_hes_apc.R",
    "R/derive_hospital_code_from_hes/02_add_endoscopy_site.R",
    "R/derive_hospital_code_from_hes/03_endoscopy_diagnostics.R"))

  run_stage("2. CWT merge and waiting times", c(
    "R/merge_cwt_to_get_dtt/02_derive_pathway.R",
    "R/merge_cwt_to_get_dtt/03_cwt_merge.R"))
  
  run_stage("3. positive-deviance analysis", c(
    "R/identify_positive_deviants/02_build_cohort.R",
    "R/identify_positive_deviants/03_table1_characteristics.R",
    "R/identify_positive_deviants/04_estimation_weights.R",
    "R/identify_positive_deviants/05_shrinkage.R",
    "R/identify_positive_deviants/06_ranks_caterpillars.R"))
  
  cat("\nPipeline complete:", format(Sys.time()), "\n")
  cat("  patient-level data:   ", dir_out, "\n")
  cat("  results for transfer: ", dir_transfer, "\n")
  if (!is.null(log_file)) cat("  run log:              ", basename(log_file), "\n")
}

run_full_pipeline()
