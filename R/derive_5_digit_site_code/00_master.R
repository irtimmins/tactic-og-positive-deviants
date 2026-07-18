# =============================================================================
# 00  Master  -  run the OG rapid build end to end
# -----------------------------------------------------------------------------
# Starts from the 20260212 rapid tumour extract and adds the fields the registry
# does not carry. Step 02 brings in the five-character site of diagnosis from the
# COSD extract; the CWT merge follows. Step 03 takes the choice 02 made apart, so
# it can be reviewed rather than taken on trust.
#
# dir_raw, dir_out and dir_transfer are only filled in here if not already set,
# so a script or console session can point them elsewhere first - at Data/sim
# for a trial run, or at a different drive letter - and this file leaves that
# alone rather than overwriting it. dir_ref is not set here since
# 01_define_parameters.R already defaults it, so that default lives in one place.
#
# dir_transfer is the analysis server's unencrypted output folder - the only
# place a file can be written and then sent off the server. The build log and
# the diagnostics csvs from 03 are copied there once the run finishes. For a
# simulated or trial run, set dir_transfer <- NULL (or point it at a scratch
# folder) before sourcing this file, or the trial's made-up numbers land in the
# real folder next to real results.
#
# The validation scripts prove the build code against fixtures and against
# simulated data, so they can be run anywhere - they never read dir_raw and never
# write dir_out or dir_transfer. Run them after any change to the rules in 01.
# =============================================================================

if (!exists("dir_raw")) dir_raw <- "W:/_DATA/IainTimmins/2026 OG SOTN data"
if (!exists("dir_out")) dir_out <- "Data/OG"
if (!exists("dir_transfer"))
  dir_transfer <- paste("S:/NATCAN_Projects/NOGCA/Iain Timmins",
                        "Results transfer out of server/Build dataset", sep = "/")

# Everything below runs inside a function rather than at the top level of this
# script. on.exit() only defers to the end of an enclosing function - at the top
# level there is no function for it to attach to, and the sink it sets up below
# would close again almost immediately, catching one line of the log and no
# more. Wrapping the run in run_og_build() gives on.exit something real to hang
# off, so the log stays open for the whole run and still closes cleanly if a
# step errors partway through.
run_og_build <- function() {
  dir_build <- "R/derive_5_digit_site_code"
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