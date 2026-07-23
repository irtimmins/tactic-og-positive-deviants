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
#   hes_extract.rds            (dir_out)  reused by stage 03a. Patient-level, so
#                              it stays on the restricted drive.
#   site_codes_to_lookup.txt   (dir_ref, and dir_transfer when set) copy this ONE
#                              file to the internet machine, then run stages 01
#                              and 02 there. It holds organisation codes and no
#                              patient rows, which is what makes it safe to
#                              carry off the server.
#
# Set dir_transfer <- transfer_root before running if you want the code list
# placed straight into the transfer area rather than fetched out of dir_ref.
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
  
  # name the copy that is actually reachable: the transfer area if this run has
  # one, since that is the only folder whose contents can leave the server.
  f_codes <- file.path(dir_ref, "site_codes_to_lookup.txt")
  f_out   <- if (!is.null(dir_transfer)) {
    file.path(dir_transfer, "site_codes_to_lookup.txt")
  } else {
    f_codes
  }
  cat("\nStage 00 complete.\n")
  cat("Copy", f_out, "to the internet machine, then run stages 01 and 02 there.\n")
  if (is.null(dir_transfer))
    cat("  (dir_transfer was not set, so this is the copy in the reference",
        "folder.\n   Set dir_transfer <- transfer_root for a copy in the",
        "transfer area.)\n")
  invisible(NULL)
}
run_this_step()