# =============================================================================
# 05  Derive the analysis cohort  -  run this step
# -----------------------------------------------------------------------------
# WHERE: the analysis server. Reads the CWT-merged cohort stage 04 wrote and the
# two reference maps from stages 01 and 02; needs no internet and none of the
# modelling packages.
#
# Draws the curative analysis cohort through the patient and hospital funnels,
# writes the attrition flowchart and the per-hospital counts, and describes the
# cohort in Table 1. Everything here is reviewable epidemiology; the modelling
# happens in stage 06.
#
# Needs:    og_cohort_cwt.rds (stage 04), og_endoscopy_hes_lookup.rds (stage
#           03a), site_trust_map.csv and valid_diagnosing_trusts.csv (dir_ref)
# Produces: pd_cohort.rds, pd_flow.csv, pd_hospital_counts.csv, table 1 files
# =============================================================================

source("R/config/directories.R")

run_this_step <- function() {
  dir_build <- "R/05_derive_analysis_cohort"
  step <- function(file) {
    cat("\n========== ", file, " ==========\n", sep = "")
    source(file.path(dir_build, file), local = new.env())
  }
  step("02_build_cohort.R")             # the two funnels and the flowchart
  step("03_table1_characteristics.R")   # describe the cohort
  cat("\nStage 05 complete: pd_cohort.rds and the flowchart written.\n")
  invisible(NULL)
}
run_this_step()
