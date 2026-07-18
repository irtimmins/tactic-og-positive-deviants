# =============================================================================
# 91  Check the CWT merge logic
# -----------------------------------------------------------------------------
# Proves the pathway derivation and the CWT merge against small worked examples
# where the right answer is known by hand, and against the simulated data. Needs
# no real data and no internet, so it can be run anywhere after any change to the
# rules in 01.
#
# Run from the project root:
#   Rscript R/merge_cwt_to_get_dtt/91_check_cwt_merge.R
# =============================================================================

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(magrittr)
})

dir_build <- "R/merge_cwt_to_get_dtt"

.saved <- mget(c("dir_out", "dir_raw", "dir_ref", "read_cohort_site",
                 "read_cohort_pathway", "read_cwt", "anchor_prefers_endoscopy"),
               ifnotfound = list(NULL), envir = globalenv())
restore_session <- function() {
  for (nm in names(.saved))
    if (is.null(.saved[[nm]])) suppressWarnings(rm(list = nm, envir = globalenv()))
    else assign(nm, .saved[[nm]], envir = globalenv())
}

.checks <- new.env(); .checks$rows <- list()
expect <- function(label, cond) {
  ok <- isTRUE(cond)
  .checks$rows[[length(.checks$rows) + 1]] <- list(label = label, ok = ok)
  cat(if (ok) "  pass  " else "  FAIL  ", label, "\n")
}

# run 02 then 03 over given registry + cwt data, hand back the merged cohort
run_merge <- function(registry, cwt, endoscopy_anchor = TRUE) {
  d <- tempfile("cwt_check_"); dir.create(d)
  assign("dir_out", d, envir = globalenv())
  assign("dir_raw", tempdir(), envir = globalenv())
  assign("dir_ref", tempdir(), envir = globalenv())
  assign("anchor_prefers_endoscopy", endoscopy_anchor, envir = globalenv())
  assign("read_cohort_site", function() registry, envir = globalenv())
  assign("read_cwt", function() cwt, envir = globalenv())
  invisible(capture.output(suppressMessages(
    sys.source(file.path(dir_build, "02_derive_pathway.R"), new.env()))))
  invisible(capture.output(suppressMessages(
    sys.source(file.path(dir_build, "03_cwt_merge.R"), new.env()))))
  readRDS(file.path(d, "og_cohort_cwt.rds"))
}

# a minimal registry row with everything absent, to be overridden per case
reg_row <- function(id, ...) {
  base <- tibble(
    patient_pseudo_id = id, tumour_site = "C15",
    diagnosisdate = as.Date("2022-01-10"), endodateHES = as.Date(NA),
    EMR_ESDdateHES = as.Date(NA), surgery_date = as.Date(NA),
    sact_first_date = as.Date(NA), rt_first_date = as.Date(NA),
    surgintent = NA_integer_, rt_first_intent = NA_integer_,
    sact_first_intent_pall = NA_integer_,
    surgery_trust = NA_character_, rt_first_trust = NA_character_,
    EMR_ESDtrustHES = NA_character_)
  over <- list(...)
  for (nm in names(over)) base[[nm]] <- over[[nm]]
  base
}
cwt_row <- function(id, dtt, treat, modality, site = "C15") {
  tibble(patient_pseudo_id = id,
         treat_period_start = format(as.Date(dtt), "%Y-%m-%d"),
         treat_start = format(as.Date(treat), "%Y-%m-%d"),
         modality = as.integer(modality), site_icd10 = site,
         org_dec_to_treat = "RAA", org_treat_start = "RAA")
}

cat("\nPathway derivation\n")

reg <- bind_rows(
  reg_row("surg", surgery_date = as.Date("2022-03-01"), surgintent = 1L,
          surgery_trust = "RAA01"),
  reg_row("neoadj", surgery_date = as.Date("2022-05-01"),
          sact_first_date = as.Date("2022-02-01"), surgintent = 1L,
          surgery_trust = "RBB01"),
  reg_row("defcrt", rt_first_date = as.Date("2022-03-01"),
          sact_first_date = as.Date("2022-02-25"), rt_first_intent = 1L,
          rt_first_trust = "RCC01"),
  reg_row("pallrt", rt_first_date = as.Date("2022-03-01"),
          rt_first_intent = 11L, rt_first_trust = "RDD01"),
  reg_row("emr", EMR_ESDdateHES = as.Date("2022-02-10"),
          EMR_ESDtrustHES = "REE"),
  reg_row("none"))
cwt <- bind_rows(
  cwt_row("surg",   "2022-02-20", "2022-03-01", 23),
  cwt_row("neoadj", "2022-01-25", "2022-02-01", 2),
  cwt_row("defcrt", "2022-02-15", "2022-02-25", 4),
  cwt_row("pallrt", "2022-02-20", "2022-03-01", 5),
  cwt_row("emr",    "2022-02-01", "2022-02-10", 23),
  cwt_row("none",   "2022-02-01", "2022-02-05", 7))

out <- run_merge(reg, cwt)
val <- function(id, col) out[[col]][out$patient_pseudo_id == id]

expect("surgery alone gives Surgery only",
       val("surg", "tx_pathway") == "Surgery only")
expect("chemo before surgery gives neoadjuvant chemo",
       val("neoadj", "tx_pathway") == "Surgery + neoadjuvant chemo")
expect("curative RT with chemo, no surgery, gives Definitive chemoRT",
       val("defcrt", "tx_pathway") == "Definitive chemoRT")
expect("palliative RT gives Palliative RT only",
       val("pallrt", "tx_pathway") == "Palliative RT only")
expect("EMR alone gives EMR/ESD only",
       val("emr", "tx_pathway") == "EMR/ESD only")
expect("nothing recorded gives No treatment recorded",
       val("none", "tx_pathway") == "No treatment recorded")

cat("\nThe clock-stop date\n")
expect("surgery pathway stops the clock on the surgery date",
       val("surg", "first_tx_date") == as.Date("2022-03-01"))
# neoadjuvant chemo is itself a first definitive treatment (guidance 3.9.1:
# chemotherapy, including prior to planned surgery), so it stops the clock on the
# chemo date, not the later surgery. This matches the reference pathway script.
expect("neoadjuvant chemo stops the clock on the chemo date",
       val("neoadj", "first_tx_date") == as.Date("2022-02-01"))
expect("definitive chemoRT stops on the earlier of chemo and RT",
       val("defcrt", "first_tx_date") == as.Date("2022-02-25"))
expect("palliative RT leaves the clock-stop unset",
       is.na(val("pallrt", "first_tx_date")))

cat("\nThe decision-to-treat merge\n")
expect("the CWT DTT is attached for a curative patient",
       val("surg", "cwt_dtt_date") == as.Date("2022-02-20"))
expect("DTT-to-treatment is the treatment minus the DTT",
       val("surg", "wt_dtt_to_tx") == 9)
expect("a curative patient with agreeing dates is valid",
       isTRUE(val("surg", "dtt_valid")))

cat("\nModality consistency with the pathway\n")
# a surgery patient whose only CWT row is a palliative modality: not consistent,
# so no DTT should be taken from it (no pathway-consistent record, falls back to
# any in-window - but here the palliative row is the only one, and with no
# consistent record the earliest in-window row still anchors). Check the DTT is
# whatever single row exists, and that group consistency was recorded.
reg2 <- reg_row("s2", surgery_date = as.Date("2022-03-01"), surgintent = 1L,
                surgery_trust = "RAA01")
cwt2 <- bind_rows(
  cwt_row("s2", "2022-02-01", "2022-02-05", 7),      # palliative, inconsistent
  cwt_row("s2", "2022-02-20", "2022-03-01", 23))     # surgery, consistent
out2 <- run_merge(reg2, cwt2)
expect("the pathway-consistent surgery row supplies the DTT, not the palliative one",
       out2$cwt_dtt_date[out2$patient_pseudo_id == "s2"] == as.Date("2022-02-20"))

cat("\nThe surgery transition window\n")
# code 1 ("01") after the window should not count as surgery
reg3 <- reg_row("w1", surgery_date = as.Date("2023-03-01"), surgintent = 1L,
                surgery_trust = "RAA01")
cwt3 <- cwt_row("w1", "2023-02-18", "2023-03-01", 1)   # retired code, after window
out3 <- run_merge(reg3, cwt3)
expect("the retired surgery code after the window is not used",
       is.na(out3$cwt_dtt_date[out3$patient_pseudo_id == "w1"]))

cat("\nThe clock-start anchor\n")
reg4 <- reg_row("a1", diagnosisdate = as.Date("2022-02-01"),
                endodateHES = as.Date("2022-01-15"),
                surgery_date = as.Date("2022-03-01"), surgintent = 1L,
                surgery_trust = "RAA01")
cwt4 <- cwt_row("a1", "2022-02-20", "2022-03-01", 23)
out4a <- run_merge(reg4, cwt4, endoscopy_anchor = TRUE)
out4b <- run_merge(reg4, cwt4, endoscopy_anchor = FALSE)
expect("with the endoscopy anchor, the DTT interval counts from endoscopy",
       out4a$wt_anchor_to_dtt[1] == as.integer(as.Date("2022-02-20") - as.Date("2022-01-15")))
expect("without it, the DTT interval counts from diagnosis",
       out4b$wt_anchor_to_dtt[1] == as.integer(as.Date("2022-02-20") - as.Date("2022-02-01")))
expect("the endoscopy-to-diagnosis lead-in is recorded",
       out4a$wt_endo_to_dx[1] == as.integer(as.Date("2022-02-01") - as.Date("2022-01-15")))

cat("\nValidity reasons\n")
# DTT after the 180-day cap
reg5 <- reg_row("c1", surgery_date = as.Date("2022-12-01"), surgintent = 1L,
                surgery_trust = "RAA01")
cwt5 <- cwt_row("c1", "2022-01-20", "2022-12-01", 23)   # DTT-to-treat > 180
out5 <- run_merge(reg5, cwt5)
expect("an over-long DTT-to-treatment interval is marked over cap",
       out5$dtt_to_tx_reason[1] == "over cap")
expect("and is not counted as valid", !isTRUE(out5$dtt_valid[1]))

cat("\nSimulated data end to end\n")
dir_sim2 <- tempfile("cwt_sim_"); dir.create(dir_sim2)
old_sim <- if (exists("dir_sim")) dir_sim else NULL
assign("dir_sim", dir_sim2, envir = globalenv()); sim_scale <- 1
invisible(capture.output(suppressMessages(
  sys.source(file.path(dir_build, "90_simulate_inputs.R"), new.env()))))
sim_reg <- readRDS(file.path(dir_sim2, "og_cohort_site.rds"))
sim_cwt <- haven::read_dta(file.path(dir_sim2, "20260212_all_cwt_rapid_202601_OG.dta"))
sim_out <- run_merge(sim_reg, sim_cwt)
if (is.null(old_sim)) suppressWarnings(rm("dir_sim", envir = globalenv())) else
  assign("dir_sim", old_sim, envir = globalenv())

expect("every registry patient survives the merge",
       nrow(sim_out) == nrow(sim_reg))
expect("most patients get a DTT",
       mean(!is.na(sim_out$cwt_dtt_date)) > 0.7)
expect("valid DTT-to-treatment intervals are the majority",
       mean(sim_out$dtt_valid, na.rm = TRUE) > 0.6)
expect("no endoscopy date sits after its diagnosis date",
       all(sim_out$wt_endo_to_dx >= 0, na.rm = TRUE))
expect("valid DTT-to-treatment intervals are all positive and within the cap",
       {
         v <- sim_out$wt_dtt_to_tx[sim_out$dtt_valid]
         all(v > 0 & v <= 180, na.rm = TRUE)
       })

# =============================================================================
res <- bind_rows(lapply(.checks$rows, as_tibble))
n_fail <- sum(!res$ok)
cat("\n", nrow(res), "checks,", n_fail, "failed\n")
restore_session()
if (n_fail) {
  cat("\nfailed:\n"); cat(paste0("  ", res$label[!res$ok], collapse = "\n"), "\n")
  quit(status = 1, save = "no")
}
cat("All checks passed.\n")
