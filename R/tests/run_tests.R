# =============================================================================
# run_tests.R  -  prove the whole pipeline's logic on synthetic data
# -----------------------------------------------------------------------------
# WHERE: anywhere - the internet machine, a laptop, the server. Needs no real
#        data, no internet, and (for the pure-logic parts) none of the modelling
#        packages. Nothing patient-level is read or written.
#
# Runs each stage's own logic-check script (one stage-matched test file per
# stage, numbered to match the stage folders). Each of those builds its own small
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
# so every test file guards that quit() with interactive(), raising a normal,
# catchable error instead when it is not running as its own top-level script.
# That is what the tryCatch below relies on.
#
# dir_transfer is deliberately left NULL: this run must never write made-up
# numbers to the real S: transfer area. The checks point their own paths at
# temporary folders, so this leaves the real dir_out / dir_ref untouched too.
#
# What this does and does not cover:
#   - the test scripts validate stage LOGIC (covariate coding, the site
#     ranking, the pathway rules, the CWT merge, the standardisation maths, the
#     cohort funnel). That is the part worth trusting.
#   - the modelling checks inside 06_test_positive_deviance.R need balancer and
#     rstan; where those are not installed they skip themselves internally, so
#     this run is still meaningful on a machine without the modelling stack (it
#     just proves less). Run it on the server too, where they are present.
#   - the reference-map tests (01, 02) prove the SNOMED and ODS map parsing
#     against fixtures; they need no internet, but do need tidyverse/jsonlite -
#     where those are not installed the check is reported as skipped, not
#     failed, so the suite stays usable on a machine that lacks them.
#
# Run from the repository root, in a batch shell (Rscript R/tests/run_tests.R) or by
# sourcing / pasting into an interactive R session (R console, RStudio) -
# both work the same way.
# =============================================================================

# this run must never write made-up numbers to the real S: transfer area. In an
# interactive session it is easy to have already sourced run_master.R or
# run_master.R with dir_transfer set earlier and left dir_transfer set from that, so check
# for it explicitly and say what to do, rather than a bare, unexplained error.
if (!is.null(get0("dir_transfer")))
  stop("dir_transfer is already set (to '", get0("dir_transfer"), "') in this ",
       "session - probably left over from an earlier run_master.R or ",
       "run_master.R with dir_transfer set call. Run rm(dir_transfer) (or start a fresh ",
       "R session) before sourcing R/tests/run_tests.R, so a test run can never write to ",
       "the real transfer area.", call. = FALSE)

checks <- c(
  "R/tests/01_test_site_trust_map.R",
  "R/tests/02_test_snomed_map.R",
  "R/tests/03a_test_hospital_from_hes.R",
  "R/tests/03b_test_hospital_from_cosd.R",
  "R/tests/04_test_cwt_merge.R",
  "R/tests/06_test_positive_deviance.R")

# these paths are relative to the repository root. Sourcing this file from
# somewhere else (a different working directory in RStudio, or a project opened
# a level up) is a common way for every check to fail with a confusing "cannot
# open file", so catch it here with a clear message instead.
missing_files <- checks[!file.exists(checks)]
if (length(missing_files) == length(checks))
  stop("R/tests/run_tests.R must be run from the repository root (the folder ",
       "holding R/ and this file). None of the test files were found from the ",
       "current working directory (", getwd(), ").", call. = FALSE)
if (length(missing_files)) {
  cat("Test file(s) not present - reported as missing, not run:\n  ",
      paste(missing_files, collapse = "\n  "), "\n", sep = "")
  checks <- setdiff(checks, missing_files)
}

results <- data.frame(check = character(), status = character(),
                      note = character(), stringsAsFactors = FALSE)
for (mf in missing_files)
  results <- rbind(results, data.frame(check = mf, status = "MISSING",
                                       note = "test file not present"))

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
  # quit() only makes sense for a batch run (Rscript R/tests/run_tests.R from a shell) -
  # it sets the exit status a CI system or a human would check. Calling it in an
  # interactive session (sourced or pasted into R/RStudio) tries to close that
  # session instead, which can hang or prompt to save the workspace. So only
  # quit when not interactive; otherwise just stop with a message.
  if (!interactive()) quit(status = 1, save = "no")
  stop(n_fail, " check script(s) failed - see the output above.", call. = FALSE)
}
cat("\nAll", nrow(results) - n_skip, "runnable check scripts passed",
    if (n_skip) sprintf("(%d skipped).", n_skip) else ".", "\n")