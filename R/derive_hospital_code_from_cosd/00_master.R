# =============================================================================
# 00  Master  -  run the OG rapid build end to end
# -----------------------------------------------------------------------------
# Starts from the 20260212 rapid tumour extract and adds the fields the registry
# does not carry. Step 02 brings in the five-character site of diagnosis from the
# COSD extract; the CWT merge follows. Step 03 takes the choice 02 made apart, so
# it can be reviewed rather than taken on trust.
#
# Storage is defined once in R/config/directories.R (dir_raw, dir_out, dir_ref,
# dir_debug, dir_transfer, transfer_root) and sourced here, so this stage's
# defaults can never drift from the rest of the pipeline's. Set any of those
# variables before sourcing this file to override them for one run - e.g.
# dir_out <- "Data/sim" for a trial run against synthetic data, or
# dir_transfer <- NULL to keep a trial's made-up numbers off the S: drive.
#
# dir_transfer is the analysis server's unencrypted output folder - the only
# place a file can be written and then sent off the server (the W: data drive is
# encrypted and cannot be transferred). The build log and the diagnostics csvs
# from 03 are copied there once the run finishes, if it is set. It defaults to
# NULL (off) in R/config/directories.R; set dir_transfer <- transfer_root before
# sourcing this file for a real run that should reach S: (run_master.R does
# this for the whole pipeline in one place).
#
# The validation script proves the build code against fixtures and against
# simulated data, so it can be run anywhere - it never reads dir_raw and never
# writes dir_out or dir_transfer. Run it after any change to the rules in 01.
# =============================================================================

source("R/config/directories.R")

# Everything below runs inside a function rather than at the top level of this
# script. on.exit() only defers to the end of an enclosing function - at the top
# level there is no function for it to attach to, and the sink it sets up below
# would close again almost immediately, catching one line of the log and no
# more. Wrapping the run in run_og_build() gives on.exit something real to hang
# off, so the log stays open for the whole run and still closes cleanly if a
# step errors partway through.
run_og_build <- function() {
  dir_build <- "R/derive_hospital_code_from_cosd"
  step <- function(file) {
    cat("\n========== ", file, " ==========\n", sep = "")
    source(file.path(dir_build, file), local = new.env())
  }
  
  log_file <- NULL
  if (!is.null(dir_transfer)) {
    dir.create(dir_transfer, recursive = TRUE, showWarnings = FALSE)
    log_file <- file.path(dir_transfer,
                          paste0("build_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
    log_con <- file(log_file, open = "wt")
    sink(log_con, split = TRUE)
    on.exit({ sink(); close(log_con) }, add = TRUE)
    cat("Build started:", format(Sys.time()), "\n")
  }
  
  step("02_add_site_of_diagnosis.R")   # COSD -> site_dx_code on the registry
  step("03_site_diagnostics.R")        # takes the choice apart for review
  
  # validation - safe to run anywhere, including off the analysis server
  # step("91_check_site_logic.R")
  
  cat("\nBuild complete:", file.path(dir_out, "og_cohort_site.rds"), "\n")
  
  if (!is.null(dir_transfer)) {
    # .txt rather than .csv: the transfer path appends an encrypted footer to
    # .csv files, which leaves the content readable but the file broken for
    # anything that reads to the end. See the note in 03.
    diag_files <- file.path(dir_out, c("diag_field_effect.txt",
                                       "diag_snomed_topography_conflicts.txt", "diag_trust_vs_tumour_picks.txt"))
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

run_og_build()