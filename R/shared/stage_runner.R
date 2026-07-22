# =============================================================================
# stage_runner.R  -  the machinery behind run_master()
# -----------------------------------------------------------------------------
# Sourced by run_master.R at the project root. This file only defines the
# function; it runs nothing. Everything a person needs to decide - which stages,
# on which machine, real or simulated - is in run_master.R, and there is no
# reason to read this file to use the pipeline.
#
# What run_master() does, in order:
#   1  turn the stages asked for into a list of folder names, expanding "3"
#   2  in simulated mode, build the made-up inputs first
#   3  open a log, if this is a real run with a transfer area set
#   4  source each stage's 00_run_this_step.R, in turn
# =============================================================================

source("R/config/directories.R")

# where each stage lives. The names here are the stage labels run_master()
# accepts; keep them in step with the folder names under R/.
stage_dirs <- c(
  "0"  = "R/00_list_hospital_codes",
  "1"  = "R/01_map_hospitals_to_trusts",
  "2"  = "R/02_fetch_snomed_tumour_map",
  "3a" = "R/03a_hospital_from_hes",
  "3b" = "R/03b_hospital_from_cosd",
  "4"  = "R/04_add_dtt_to_cohort",
  "5"  = "R/05_derive_analysis_cohort",
  "6"  = "R/06_identify_positive_deviants")

run_master <- function(stages = 3:6, mode = c("real", "simulated")) {
  mode <- match.arg(mode)
  
  # stage 3 means both halves. 3b (COSD, the hospital of diagnosis) runs first,
  # because stage 4 reads the cohort it writes; then 3a (HES, the endoscopy
  # hospital), whose lookup stage 5 joins on. Asking for "3a" or "3b" alone
  # works too.
  keys <- unlist(lapply(as.character(stages), function(s)
    if (s == "3") c("3b", "3a") else s), use.names = FALSE)
  
  bad <- setdiff(keys, names(stage_dirs))
  if (length(bad))
    stop("unknown stage(s): ", paste(bad, collapse = ", "),
         " - valid: 0, 1, 2, 3, 3a, 3b, 4, 5, 6", call. = FALSE)
  
  if (mode == "simulated") {
    # made-up inputs, test outputs, and no route to the transfer area. The
    # preparation script builds every extract and reference map the chosen
    # stages will read; the identical stage code then runs against them.
    if (any(keys %in% c("1", "2")))
      stop("stages 1 and 2 fetch from the internet and have no simulated ",
           "form - their outputs are simulated directly by the preparation ",
           "script. Ask for stages 0 or 3-6.", call. = FALSE)
    assign("dir_transfer", NULL, envir = globalenv())
    source("R/simulate_test_data/00_prepare_simulated_run.R")
    prepare_simulated_run()
  }
  
  # a full console log, alongside the results, so a look back at what happened
  # does not depend on anyone keeping a terminal scrollback. Only for a real run
  # with a transfer area set - a simulated run has nowhere it should write to.
  log_file <- NULL
  if (!is.null(dir_transfer) && mode == "real") {
    dir.create(dir_transfer, recursive = TRUE, showWarnings = FALSE)
    log_file <- file.path(dir_transfer,
                          paste0("pipeline_log_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
    log_con <- file(log_file, open = "wt")
    sink(log_con, split = TRUE)   # split = TRUE: still visible on screen too
    on.exit({ sink(); close(log_con) }, add = TRUE)
    cat("Pipeline run started:", format(Sys.time()), "| stages:",
        paste(keys, collapse = " "), "\n")
  }
  
  for (k in keys) {
    cat("\n==================== stage ", k, ": ", stage_dirs[[k]],
        " ====================\n", sep = "")
    source(file.path(stage_dirs[[k]], "00_run_this_step.R"), local = new.env())
  }
  
  cat("\nrun_master finished:", format(Sys.time()), "\n")
  if (!is.null(log_file)) cat("  run log:", basename(log_file), "\n")
  invisible(NULL)
}
