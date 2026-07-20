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
#                 (run_test.R) leaves it NULL and writes to a local scratch folder
#                 instead.
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

run_stage("1. site of diagnosis", c(
  "R/derive_5_digit_site_code/02_add_site_of_diagnosis.R",
  "R/derive_5_digit_site_code/03_site_diagnostics.R"))

run_stage("2. CWT merge and waiting times", c(
  "R/merge_cwt_to_get_dtt/02_derive_pathway.R",
  "R/merge_cwt_to_get_dtt/03_cwt_merge.R"))

run_stage("3. positive-deviance analysis", c(
  "R/identify_positive_deviants/02_build_cohort.R",
  "R/identify_positive_deviants/03_table1_characteristics.R",
  "R/identify_positive_deviants/04_estimation_weights.R",
  "R/identify_positive_deviants/05_shrinkage.R",
  "R/identify_positive_deviants/06_ranks_caterpillars.R"))

cat("\nPipeline complete.\n")
cat("  patient-level data:   ", dir_out, "\n")
cat("  results for transfer: ", dir_transfer, "\n")