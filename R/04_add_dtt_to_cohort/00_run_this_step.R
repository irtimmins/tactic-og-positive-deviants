# =============================================================================
# 04  Add the decision-to-treat to the cohort  -  run this step
# -----------------------------------------------------------------------------
# WHERE: the analysis server. Reads the cohort stage 03b wrote and the raw CWT
# extract; needs no internet.
#
# Derives each patient's treatment pathway from the treatment dates and intents,
# then merges the Cancer Waiting Times records to place the decision-to-treat
# and build the waiting-time intervals - including the endoscopy-to-DTT wait the
# analysis runs on.
#
# Needs:    og_cohort_site.rds (stage 03b), the CWT dta
# Produces: og_cohort_cwt.rds  read by stage 05
# =============================================================================

source("R/config/directories.R")

run_this_step <- function() {
  dir_build <- "R/04_add_dtt_to_cohort"
  step <- function(file) {
    cat("\n========== ", file, " ==========\n", sep = "")
    source(file.path(dir_build, file), local = new.env())
  }
  log_file <- NULL
  if (!is.null(dir_transfer)) {
    dir.create(dir_transfer, recursive = TRUE, showWarnings = FALSE)
    log_file <- file.path(dir_transfer,
                          paste0("cwt_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
    log_con <- file(log_file, open = "wt")
    sink(log_con, split = TRUE)
    on.exit({ sink(); close(log_con) }, add = TRUE)
    cat("Build started:", format(Sys.time()), "\n")
  }
  step("02_derive_pathway.R")            # treatment dates -> tx_pathway
  step("03_cwt_merge.R")                 # CWT -> decision-to-treat and waits

  # check_stage_and_dtt.R sits in this folder but is not run here: it is a
  # one-off look at whether the non-positive decision-to-treat intervals are a
  # stage-4 effect rather than a fault in the merge. It writes nothing. Run it
  # by hand when the question comes up.
  #
  # the logic checks for this stage live in R/tests/04_test_cwt_merge.R
  cat("\nStage 04 complete:", f_cohort_cwt, "\n")
  if (!is.null(log_file)) cat("Build log:", basename(log_file), "\n")
  invisible(NULL)
}
run_this_step()
