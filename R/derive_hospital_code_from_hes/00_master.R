# =============================================================================
# 00  Master  -  derive the endoscopy hospital site from HES, end to end
# -----------------------------------------------------------------------------
# WHERE: the analysis server, on the real rapid and HES-APC extracts.
#
# Starts from the 20260212 rapid tumour extract and the HES-APC extract. Step 02
# re-finds each patient's diagnostic endoscopy in HES and reads its five-
# character site onto the record; step 03 takes that choice apart so it can be
# reviewed rather than taken on trust.
#
# Storage is defined once in R/config/directories.R (dir_raw, dir_out, dir_ref,
# dir_transfer, transfer_root) and sourced there via 01, so this stage's defaults
# cannot drift from the rest of the pipeline's. Set any of those variables before
# sourcing this file to override them for one run.
#
# dir_transfer is the analysis server's unencrypted output folder - the only
# place a file can be written and then sent off the server. The build log and the
# diagnostics txts from 03 are copied there once the run finishes, if it is set.
# It defaults to NULL (off); set dir_transfer <- transfer_root before sourcing
# this file for a run whose results should reach S:.
#
# The validation script (91) proves the build against worked examples and
# simulated data, so it can be run anywhere - it never reads dir_raw and never
# writes dir_out or dir_transfer. Run it after any change to the rules in 01.
# =============================================================================

source("R/config/directories.R")

run_endoscopy_build <- function() {
  dir_build <- "R/derive_hospital_code_from_hes"
  step <- function(file) {
    cat("\n========== ", file, " ==========\n", sep = "")
    source(file.path(dir_build, file), local = new.env())
  }

  log_file <- NULL
  if (!is.null(dir_transfer)) {
    dir.create(dir_transfer, recursive = TRUE, showWarnings = FALSE)
    log_file <- file.path(dir_transfer,
                          paste0("endoscopy_hes_log_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
    log_con <- file(log_file, open = "wt")
    sink(log_con, split = TRUE)
    on.exit({ sink(); close(log_con) }, add = TRUE)
    cat("Build started:", format(Sys.time()), "\n")
  }

  step("02a_extract_hes_apc.R")       # cut HES-APC down to the cohort, once
  step("02_add_endoscopy_site.R")     # HES -> endoscopy_site on the registry
  step("03_endoscopy_diagnostics.R")  # takes the choice apart for review

  # validation - safe to run anywhere, including off the analysis server
  # step("91_check_endoscopy_logic.R")

  cat("\nBuild complete:", file.path(dir_out, "og_endoscopy_hes_lookup.rds"),
      "\n")

  if (!is.null(dir_transfer)) {
    # .txt rather than .csv: the transfer path appends an encrypted footer to
    # .csv files, which leaves the content readable but the file broken for
    # anything that reads to the end. See the note in 03.
    diag_files <- file.path(dir_out, c("diag_endoscopy_outcome.txt",
                                       "diag_endoscopy_by_window.txt",
                                       "diag_endoscopy_provider_disagree.txt"))
    present <- file.exists(diag_files)
    moved <- if (any(present))
      file.copy(diag_files[present], dir_transfer, overwrite = TRUE) else logical(0)

    cat("\nCopied", sum(moved), "of", sum(present), "diagnostics files to",
        dir_transfer, "\n")
    if (any(!present))
      cat("Not found, so not copied:",
          paste(basename(diag_files[!present]), collapse = ", "), "\n")
    cat("Build log:", basename(log_file), "\n")
  }

  invisible(NULL)
}

run_endoscopy_build()
