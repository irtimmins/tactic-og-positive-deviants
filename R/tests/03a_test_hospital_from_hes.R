# =============================================================================
# 91  Check the endoscopy-site logic
# -----------------------------------------------------------------------------
# Checks that the code in 02 does what 01 says it does, without needing the real
# data. Two parts, and neither reads or writes the project's data folder:
#
#   A  worked examples. A handful of patients whose right answer is known because
#      it was built in - a clean match, the code in a later slot, a good site
#      preferred over a blank one, an unusable site, no endoscopy at all, an
#      endoscopy weeks off the date, a patient not in the extract, a dateless
#      code, and the window tightened. Each goes through the real 02 and the
#      answer is checked exactly.
#   B  stand-in data. The full-size made-up extracts from 90, where each
#      patient's true site and intended outcome are known. This says how often
#      the match lands right at a realistic scale, and that the four outcome
#      buckets come out as built.
#
# Run from the project root:
#   Rscript R/tests/03a_test_hospital_from_hes.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(magrittr)   # for %>%
})
dir_build <- "R/03a_hospital_from_hes"

.saved <- mget(c("dir_out", "dir_raw", "dir_ref", "read_rapid", "read_hes",
                 "endo_window_days"),
               ifnotfound = list(NULL), envir = globalenv())
.saved_opt <- options(og_min_input_rows = 1L)   # the examples are tiny on purpose
restore_session <- function() {
  options(.saved_opt)
  for (nm in names(.saved)) {
    if (is.null(.saved[[nm]])) suppressWarnings(rm(list = nm, envir = globalenv()))
    else assign(nm, .saved[[nm]], envir = globalenv())
  }
}

.checks <- new.env(); .checks$rows <- list()
expect <- function(label, cond) {
  ok <- isTRUE(cond)
  .checks$rows[[length(.checks$rows) + 1]] <- list(label = label, ok = ok)
  cat(if (ok) "  pass  " else "  FAIL  ", label, "\n")
}

# run 02 over the given data and hand back the cohort it makes
run_build <- function(rapid, hes, window = 7L) {
  assign("dir_out", tempfile("endo_check_"), envir = globalenv())
  assign("dir_raw", tempdir(), envir = globalenv())
  assign("dir_ref", tempdir(), envir = globalenv())
  assign("read_rapid", function() rapid, envir = globalenv())
  assign("read_hes",   function() hes,   envir = globalenv())
  assign("endo_window_days", window, envir = globalenv())
  env <- new.env()
  invisible(capture.output(
    suppressMessages(sys.source(file.path(dir_build, "02_add_endoscopy_site.R"),
                                envir = env))))
  env$og_cohort
}

# a HES episode row with the 24 operation slots, and at most one code in a slot
hes_row <- function(pid, epistart, epiorder, sitetret, procode3,
                    slot = NA, opcs = NA, opdate = NA) {
  op <- setNames(as.list(rep("-", 24)), sprintf("opertn_%02d", 1:24))
  od <- setNames(as.list(rep("",  24)), sprintf("opdate_%02d", 1:24))
  if (!is.na(slot)) { op[[slot]] <- opcs; od[[slot]] <- opdate }
  bind_cols(
    tibble(patient_pseudo_id = pid, epistart = epistart, epiend = epistart,
           admidate = epistart, epiorder = epiorder, epitype = 1L,
           sitetret = sitetret, procode3 = procode3),
    as_tibble(op), as_tibble(od))
}

# =============================================================================
# A. Worked examples
# =============================================================================
cat("\nA. worked examples\n")

d <- as.Date("2022-03-10")
rapid_eg <- tribble(
  ~patient_pseudo_id, ~tumour_site, ~endoHES, ~endodateHES, ~endotypeHES, ~endotrustHES,
  "w01_clean",        "C15", 1, d,          "G451", "RJ1",
  "w02_slot3",        "C15", 1, d,          "G161", "RJ1",
  "w03_prefer_sited", "C15", 1, d,          "G459", "RJ1",
  "w04_invalid_only", "C15", 1, d,          "G451", "RJ1",
  "w05_no_endo",      "C15", 1, d,          "G451", "RJ1",
  "w06_offwindow",    "C15", 1, d,          "G451", "RJ1",
  "w07_not_in_hes",   "C15", 1, d,          "G451", "RJ1",
  "w08_nodash",       "C15", 1, d,          "G162", "RJ1",
  "w09_no_anchor",    "C15", 0, as.Date(NA), "",    "",
  "w10_edge",         "C15", 1, d,          "G451", "RJ1")

hes_eg <- bind_rows(
  # exact match, slot 1
  hes_row("w01_clean", "2022-03-10", 1, "RJ101", "RJ1", 1, "G451", "2022-03-10"),
  # code in slot 3, a day off
  hes_row("w02_slot3", "2022-03-11", 1, "RJ102", "RJ1", 3, "G161", "2022-03-11"),
  # two endoscopies in the window: a blank site, then a real one
  hes_row("w03_prefer_sited", "2022-03-09", 1, "",     "RJ1", 1, "G459", "2022-03-09"),
  hes_row("w03_prefer_sited", "2022-03-11", 2, "RJ105","RJ1", 1, "G459", "2022-03-11"),
  # the only endoscopy has an unusable site
  hes_row("w04_invalid_only", "2022-03-10", 1, "00000","RJ1", 1, "G451", "2022-03-10"),
  # an APC episode but no endoscopy code
  hes_row("w05_no_endo", "2022-03-10", 1, "RJ108", "RJ1", 1, "G011", "2022-03-10"),
  # an endoscopy, but a month off the reference date
  hes_row("w06_offwindow", "2022-04-09", 1, "RJ109", "RJ1", 1, "G451", "2022-04-09"),
  # w07 has no HES rows at all
  # a dateless-looking date written without dashes
  hes_row("w08_nodash", "2022-03-10", 1, "RJ110", "RJ1", 1, "G162", "20220310"),
  # w09 has an endoscopy but endoHES is 0, so it is never anchored
  hes_row("w09_no_anchor", "2022-03-10", 1, "RJ111", "RJ1", 1, "G451", "2022-03-10"),
  # four days off - inside the default window, outside a window of 1
  hes_row("w10_edge", "2022-03-14", 1, "RJ112", "RJ1", 1, "G451", "2022-03-14"))

# enough rows to clear the row floor is not needed - og_min_input_rows is 1 here
out <- run_build(rapid_eg, hes_eg)
site_of <- function(pid) out$endoscopy_site[out$patient_pseudo_id == pid]
flag_of <- function(pid, col) out[[col]][out$patient_pseudo_id == pid]

expect("a clean match reads the episode's site",
       identical(site_of("w01_clean"), "RJ101"))
expect("a code in a later slot is still found",
       identical(site_of("w02_slot3"), "RJ102"))
expect("a real site is preferred over a blank one",
       identical(site_of("w03_prefer_sited"), "RJ105"))
expect("an unusable site gives no site, but the endoscopy still counts",
       is.na(site_of("w04_invalid_only")) && flag_of("w04_invalid_only", "endo_has_apc"))
# the endoscopy was found at the right date - it is the SITE that failed, so this
# patient must not be reported as an off-window or missing endoscopy
expect("an unusable site is not mistaken for a missed endoscopy",
       flag_of("w04_invalid_only", "endo_matched") &&
         !flag_of("w06_offwindow", "endo_matched"))
expect("a patient with no endoscopy code has none flagged",
       is.na(site_of("w05_no_endo")) && !flag_of("w05_no_endo", "endo_has_apc") &&
         flag_of("w05_no_endo", "endo_in_hes"))
expect("an endoscopy outside the window is not read, but is flagged present",
       is.na(site_of("w06_offwindow")) && flag_of("w06_offwindow", "endo_has_apc"))
expect("a patient not in the extract is flagged out",
       is.na(site_of("w07_not_in_hes")) && !flag_of("w07_not_in_hes", "endo_in_hes"))
expect("a date written without dashes still parses and matches",
       identical(site_of("w08_nodash"), "RJ110"))
expect("a patient with endoHES 0 is never given a site",
       is.na(site_of("w09_no_anchor")) && !flag_of("w09_no_anchor", "endoscopy_site_found"))

out1 <- run_build(rapid_eg, hes_eg, window = 1L)
site_of1 <- function(pid) out1$endoscopy_site[out1$patient_pseudo_id == pid]
expect("four days off is inside the default window",
       identical(site_of("w10_edge"), "RJ112"))
expect("four days off is outside a window of one",
       is.na(site_of1("w10_edge")))
expect("every registry row survives the merge",
       nrow(out) == nrow(rapid_eg))
expect("no site read breaks the five-character rule",
       all(nchar(na.omit(out$endoscopy_site)) == 5))

# =============================================================================
# B. Stand-in data at full size
# =============================================================================
cat("\nB. stand-in data\n")

if (!exists("run_sim_check")) run_sim_check <- TRUE

if (run_sim_check) {
  dir_sim <- tempfile("endo_sim_")
  # a tenth-size run - thousands of patients, enough to exercise the buckets and
  # the volume, but quick. The real run is on the real data, not this.
  sim_scale <- 0.1
  invisible(capture.output(
    sys.source("R/simulate_test_data/03a_simulate_hes_inputs.R", envir = new.env())))

  sim_rapid <- haven::read_dta(file.path(dir_sim,
                  "20260212_Rapidtumour_linked_2026SOTN_clean_OG_postPT.dta"))
  sim_hes <- haven::read_dta(file.path(dir_sim,
                  "20260320_nic709865_hes_apc_202510_OG.dta"))
  truth <- readRDS(file.path(dir_sim, "endoscopy_truth.rds"))

  sim_out <- run_build(sim_rapid, sim_hes)

  expect("every made-up registry row survives", nrow(sim_out) == nrow(sim_rapid))
  expect("one row per patient comes out", !anyDuplicated(sim_out$patient_pseudo_id))
  expect("no site read breaks the five-character rule",
         all(nchar(na.omit(sim_out$endoscopy_site)) == 5))
  expect("no default site escapes",
         !any(na.omit(sim_out$endoscopy_site) %in% c("00000", "89999")))

  scored <- sim_out %>%
    select(patient_pseudo_id, endoscopy_site, endoscopy_site_found,
           endo_in_hes, endo_has_apc) %>%
    left_join(select(truth, patient_pseudo_id, endo_class, true_endo_site),
              by = "patient_pseudo_id")

  cat("\n  outcome against the class each patient was built for:\n")
  scored %>%
    filter(endo_in_hes) %>%
    mutate(outcome = case_when(endoscopy_site_found ~ "site read",
                               endo_has_apc         ~ "apc off-window",
                               TRUE                 ~ "no apc endoscopy")) %>%
    count(endo_class, outcome, name = "patients") %>%
    as.data.frame() %>% print(row.names = FALSE)

  sited <- scored %>% filter(endo_class == "sited")
  expect("sited patients almost all get a site",
         mean(sited$endoscopy_site_found) > 0.95)
  expect("the site read is the true site",
         mean(sited$endoscopy_site == sited$true_endo_site, na.rm = TRUE) > 0.99)

  ow <- scored %>% filter(endo_class == "offwindow")
  expect("off-window patients are not matched",
         mean(!ow$endoscopy_site_found) > 0.95)
  expect("off-window patients still show an APC endoscopy",
         mean(ow$endo_has_apc) > 0.95)

  noapc <- scored %>% filter(endo_class == "no_apc")
  expect("no-apc patients show no APC endoscopy",
         mean(!noapc$endo_has_apc) > 0.95)

  notin <- scored %>% filter(endo_class == "not_in_hes")
  expect("not-in-extract patients are flagged out",
         mean(!notin$endo_in_hes) > 0.99)

  # 03 re-reads what 02 wrote and reports; the check is that it runs to the end
  cat("\n  diagnostics script\n")
  d3 <- tempfile("endo_diag_"); dir.create(d3)
  assign("dir_out", d3, envir = globalenv())
  assign("read_rapid", function() sim_rapid, envir = globalenv())
  assign("read_hes",   function() sim_hes,   envir = globalenv())
  invisible(capture.output(suppressMessages(
    sys.source(file.path(dir_build, "02_add_endoscopy_site.R"), new.env()))))
  invisible(capture.output(suppressMessages(
    sys.source(file.path(dir_build, "03_endoscopy_diagnostics.R"), new.env()))))
  expect("the diagnostics run through and write their outcome file",
         file.exists(file.path(d3, "diag_endoscopy_outcome.txt")))
}

# =============================================================================
# Result
# =============================================================================
res <- bind_rows(lapply(.checks$rows, as_tibble))
n_fail <- sum(!res$ok)
cat("\n", nrow(res), "checks,", n_fail, "failed\n")
if (n_fail) {
  cat("\nfailed:\n"); cat(paste0("  ", res$label[!res$ok], collapse = "\n"), "\n")
  restore_session()
  if (!interactive()) quit(status = 1, save = "no") else
    stop(n_fail, " check(s) failed - see the output above.", call. = FALSE)
}
restore_session()
cat("All checks passed.\n")
