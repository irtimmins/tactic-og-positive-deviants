# =============================================================================
# 00  List the hospital codes  -  run this step
# -----------------------------------------------------------------------------
# WHERE: the analysis server. Reads the raw extracts; needs no internet.
#
# The one stage that reads the raw HES-APC file, and it reads it once for two
# purposes: the cohort-filtered extract every later HES read reuses, and the
# full list of five-character site codes (from COSD and from HES sitetret
# together) that the internet machine looks up against ODS.
#
# Produces:
#   hes_extract.rds            (dir_out)  reused by stage 03a
#   site_codes_to_lookup.txt   (dir_ref)  copy this ONE file to the internet
#                              machine, then run stages 01 and 02 there.
# =============================================================================

source("R/config/directories.R")

run_this_step <- function() {
  dir_build <- "R/00_list_hospital_codes"
  step <- function(file) {
    cat("\n========== ", file, " ==========\n", sep = "")
    source(file.path(dir_build, file), local = new.env())
  }
  step("01_extract_hes_apc.R")    # raw HES-APC -> cohort-filtered extract
  step("02_list_site_codes.R")    # COSD + HES sitetret -> the code list
  cat("\nStage 00 complete.\n")
  cat("Copy", file.path(dir_ref, "site_codes_to_lookup.txt"),
      "to the internet machine, then run stages 01 and 02 there.\n")
  invisible(NULL)
}
run_this_step()