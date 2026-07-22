# =============================================================================
# run_to_test_logic.R  -  prove the whole pipeline's logic on synthetic data
# -----------------------------------------------------------------------------
# WHERE: anywhere - the internet machine, a laptop, the server. Needs no real
#        data, no internet, and (for the pure-logic parts) none of the modelling
#        packages. Nothing patient-level is read or written.
#
# Runs each stage's own logic-check script (the 91_check_*.R files, plus the
# reference-fetch checks 11 and 21). Each of those builds its own small
# fixtures - hand-worked examples where the answer is known, plus a simulated
# dataset - and asserts the stage behaves. A green run here means every stage's
# logic passed, not merely that data flowed through without erroring.
#
# Each check is SOURCED directly, in its own fresh environment, rather than run
# as a separate Rscript subprocess. That is deliberate: shelling out (system2 /
# system) to spawn a grandchild Rscript process adds a layer of fragility that
# is not worth it here - PATH resolution for "Rscript" varies across
# RStudio/Windows/Git-Bash setups, and in some environments a parent R session's
# own stdin/output plumbing can deadlock a spawned child's captured output. None
# of that risk exists when the check code just runs in this same process. The
# cost is that each check script's own quit(status=1) (meant for its own
# standalone Rscript use) would try to end THIS session too if it fired here -
# so every 91_check_*.R guards that quit() with interactive(), raising a normal,
# catchable error instead when it is not running as its own top-level script.
# That is what the tryCatch below relies on.
#
# dir_transfer is deliberately left NULL: this run must never write made-up
# numbers to the real S: transfer area. The checks point their own paths at
# temporary folders, so this leaves the real dir_out / dir_ref untouched too.
#
# What this does and does not cover:
#   - the 91_check scripts validate stage LOGIC (covariate coding, the site
#     ranking, the pathway rules, the CWT merge, the standardisation maths, the
#     cohort funnel). That is the part worth trusting.
#   - the modelling checks inside 91_check_positive_deviance.R need balancer and
#     rstan; where those are not installed they skip themselves internally, so
#     this run is still meaningful on a machine without the modelling stack (it
#     just proves less). Run it on the server too, where they are present.
#   - the reference-fetch checks (11, 21) prove the SNOMED and ODS map parsing
#     against fixtures; they need no internet, but do need tidyverse/jsonlite -
#     where those are not installed the check is reported as skipped, not
#     failed, so the suite stays usable on a machine that lacks them.
#
# Run from the repository root, in a batch shell (Rscript run_test.R) or by
# sourcing / pasting into an interactive R session (R console, RStudio) -
# both work the same way.
# =============================================================================

# this run must never write made-up numbers to the real S: transfer area. In an
# interactive session it is easy to have already sourced run_master.R or
# run_4_positive_deviance.R earlier and left dir_transfer set from that, so check
# for it explicitly and say what to do, rather than a bare, unexplained error.
if (!is.null(get0("dir_transfer")))
  stop("dir_transfer is already set (to '", get0("dir_transfer"), "') in this ",
       "session - probably left over from an earlier run_master.R or ",
       "run_4_positive_deviance.R call. Run rm(dir_transfer) (or start a fresh ",
       "R session) before sourcing run_test.R, so a test run can never write to ",
       "the real transfer area.", call. = FALSE)

checks <- c(
  "R/fetch_reference_data/11_check_snomed_map.R",
  "R/fetch_reference_data/21_check_site_trust_map.R",
  "R/derive_hospital_code_from_cosd/91_check_site_logic.R",
  "R/derive_hospital_code_from_hes/91_check_endoscopy_logic.R",
  "R/merge_cwt_to_get_dtt/91_check_cwt_merge.R",
  "R/identify_positive_deviants/91_check_positive_deviance.R")

# these paths are relative to the repository root. Sourcing this file from
# somewhere else (a different working directory in RStudio, or a project opened
# a level up) is a common way for every check to fail with a confusing "cannot
# open file", so catch it here with a clear message instead.
missing_files <- checks[!file.exists(checks)]
if (length(missing_files))
  stop("run_test.R must be run from the repository root (the folder holding ",
       "R/ and this file). Not found from the current working directory (",
       getwd(), "):\n  ", paste(missing_files, collapse = "\n  "),
       "\nsetwd() to the repository root and try again.", call. = FALSE)

results <- data.frame(check = character(), status = character(),
                      note = character(), stringsAsFactors = FALSE)

for (chk in checks) {
  cat("\n########## ", chk, " ##########\n", sep = "")
  # source into a fresh environment each time, so one check's settings/data
  # (e.g. dir_out, in_rds) cannot leak into the next. Also snapshot and restore
  # the working directory around each check: a check that setwd()s into a temp
  # folder and errors before restoring would otherwise leave every later check
  # resolving its relative source() paths from the wrong place - which surfaces
  # as an identical, misleading "cannot open the connection" in every check after
  # the first. Restoring wd here isolates each check from that.
  wd_before <- getwd()
  status <- "pass"; note <- ""
  result <- tryCatch({
    source(chk, local = new.env())
    "pass"
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("there is no package called", msg)) {
      note <<- sub(".*there is no package called ['\"]?([A-Za-z0-9.]+).*",
                   "missing: \\1", msg)
      "skip"
    } else {
      note <<- msg
      "FAIL"
    }
  }, finally = {
    if (getwd() != wd_before) {
      cat("  (note: this check left the working directory changed - restored)\n")
      setwd(wd_before)
    }
  })
  results <- rbind(results, data.frame(check = chk, status = result, note = note))
}

cat("\n\n==================== summary ====================\n")
print(results, row.names = FALSE)
n_fail <- sum(results$status == "FAIL")
n_skip <- sum(results$status == "skip")
if (n_skip)
  cat("\n", n_skip, "check(s) skipped for a missing package - run on a machine",
      "with the full stack (the server) to cover these.\n")
if (n_fail) {
  cat("\n", n_fail, "of", nrow(results), "check scripts FAILED on logic (see above).\n")
  # quit() only makes sense for a batch run (Rscript run_test.R from a shell) -
  # it sets the exit status a CI system or a human would check. Calling it in an
  # interactive session (sourced or pasted into R/RStudio) tries to close that
  # session instead, which can hang or prompt to save the workspace. So only
  # quit when not interactive; otherwise just stop with a message.
  if (!interactive()) quit(status = 1, save = "no")
  stop(n_fail, " check script(s) failed - see the output above.", call. = FALSE)
}
cat("\nAll", nrow(results) - n_skip, "runnable check scripts passed",
    if (n_skip) sprintf("(%d skipped).", n_skip) else ".", "\n")
