# =============================================================================
# 91  Check the site-of-diagnosis logic
# -----------------------------------------------------------------------------
# Checks that the code in 02 does what 01 says it does, without needing the real
# data. Two parts, and neither reads or writes the project's data folder:
#
#   A  worked examples. A handful of patients whose right answer is known
#      because it was built in - a clean code, a trust code pretending to be a
#      site, a GP practice code, the defaults, a row for someone's bowel cancer,
#      a C16 row on a C15 patient, a junction row competing with the patient's
#      own, two codes competing. Each is put through the
#      real build script and the answer checked exactly. This is the part that
#      proves the rules.
#   B  stand-in data. The full-size made-up extracts from 90, where the true site
#      of every patient is known. This cannot prove a rule, but it says how often
#      the ranking lands on the right answer at a realistic scale, and it shows
#      up anything that only appears in volume. C15 is reported on its own,
#      since that is where the analysis ends up.
#
# Run from the project root:
#   Rscript R/derive_5_digit_site_code/91_check_site_logic.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(magrittr)   # for %>%
})
dir_build <- "R/derive_5_digit_site_code"

# This script points dir_out at a temporary folder and hands 02 made-up data.
# Put the session back as it was on the way out, however it exits, so a later
# build in the same session cannot pick up a leftover setting.
.saved <- mget(c("dir_out", "dir_raw", "dir_ref", "read_rapid", "read_cosd",
                 "read_snomed_map", "read_site_trust_map", "site_max_rank",
                 "snomed_from_data", "drop_non_og_rows"),
               ifnotfound = list(NULL), envir = globalenv())
.saved_opt <- options(og_min_input_rows = 1L)   # the examples are tiny on purpose
restore_session <- function() {
  options(.saved_opt)
  for (nm in names(.saved)) {
    if (is.null(.saved[[nm]])) suppressWarnings(rm(list = nm, envir = globalenv()))
    else assign(nm, .saved[[nm]], envir = globalenv())
  }
}
on.exit(restore_session(), add = TRUE)

.checks <- new.env(); .checks$rows <- list()
expect <- function(label, cond) {
  ok <- isTRUE(cond)
  .checks$rows[[length(.checks$rows) + 1]] <- list(label = label, ok = ok)
  cat(if (ok) "  pass  " else "  FAIL  ", label, "\n")
}

# run 02 over the given data and hand back the cohort it makes
run_build <- function(rapid, cosd, max_rank = 3L, use_snomed = TRUE,
                      drop_non_og = TRUE, map = NULL, site_trust = NULL) {
  assign("dir_out", tempfile("og_check_"), envir = globalenv())
  assign("dir_raw", tempdir(), envir = globalenv())
  assign("dir_ref", tempdir(), envir = globalenv())
  assign("read_rapid", function() rapid, envir = globalenv())
  assign("read_cosd",  function() cosd,  envir = globalenv())
  assign("read_snomed_map", function() map, envir = globalenv())
  assign("read_site_trust_map", function() site_trust, envir = globalenv())
  assign("site_max_rank", max_rank, envir = globalenv())
  assign("snomed_from_data", use_snomed, envir = globalenv())
  assign("drop_non_og_rows", drop_non_og, envir = globalenv())
  env <- new.env()
  invisible(capture.output(
    suppressMessages(sys.source(file.path(dir_build, "02_add_site_of_diagnosis.R"),
                                envir = env))))
  env$og_cohort
}

# =============================================================================
# A. Worked examples
# =============================================================================
cat("\nA. worked examples\n")

# One line per patient. The trust column is what the registry believes; the COSD
# rows below are what the site field offers against it.
rapid_eg <- tribble(
  ~patient_pseudo_id, ~tumour_site, ~tumour_morphology_str, ~diagnosis_trust,
  "e01_clean",        "C15",        "8140",                 "RJ1",
  "e02_confirmed",    "C15",        "8140",                 "RJ1",
  "e03_confirm_wins", "C15",        "8140",                 "RJ1",
  "e04_trustcode",    "C15",        "8140",                 "RJ1",
  "e05_gpcode",       "C16",        "8140",                 "RJ1",
  "e06_default",      "C15",        "8140",                 "RJ1",
  "e07_lowercase",    "C15",        "8140",                 "RJ1",
  "e08_bowel_row",    "C15",        "8140",                 "RJ1",
  "e09_c16_on_c15",   "C15",        "8140",                 "RJ1",
  "e10_snomed_other", "C15",        "8140",                 "RJ1",
  "e11_snomed_own",   "C15",        "8140",                 "RZ9",
  "e12_morph_break",  "C15",        "8140",                 "RZ9",
  "e13_count_break",  "C15",        "8140",                 "RZ9",
  "e14_order_break",  "C15",        "8140",                 "RZ9",
  "e15_no_cosd",      "C15",        "8140",                 "RJ1",
  "e16_no_trust",     "C15",        "8140",                 "",
  "e17_topog_text",   "C15",        "8140",                 "RJ1",
  "e18_own_beats_og", "C15",        "8140",                 "RZ9")

# The SNOMED codes have to earn their meaning from the data, so the examples
# include enough plainly labelled rows for the build to learn that 111000001 is
# oesophagus and 222000002 is bowel. These filler patients exist only to teach
# it; they are not checked.
teach <- bind_rows(
  tibble(patient_pseudo_id = sprintf("t_oes_%02d", 1:12),
         site_code_of_diagnosis = "RJ199", topography_icdo3 = "C155",
         morphology_icdo3 = "", diagnosis_snomedct = 111000001,
         snomed_version = 5L),
  tibble(patient_pseudo_id = sprintf("t_bwl_%02d", 1:12),
         site_code_of_diagnosis = "RJ199", topography_icdo3 = "C182",
         morphology_icdo3 = "", diagnosis_snomedct = 222000002,
         snomed_version = 5L))

rapid_eg <- bind_rows(rapid_eg, tibble(
  patient_pseudo_id = teach$patient_pseudo_id,
  tumour_site = "C15", tumour_morphology_str = "8140", diagnosis_trust = "RJ1"))

cosd_eg <- tribble(
  ~patient_pseudo_id, ~site_code_of_diagnosis, ~topography_icdo3, ~morphology_icdo3, ~diagnosis_snomedct, ~snomed_version,
  # one usable code, nothing to support it but the trust
  "e01_clean",        "RJ101",  "",     "",     NA, 5L,
  # the row names the patient's own tumour site
  "e02_confirmed",    "RJ101",  "C155", "",     NA, 5L,
  # a confirmed code beats one that only sits in the right trust
  "e03_confirm_wins", "RJ102",  "",     "",     NA, 5L,
  "e03_confirm_wins", "RA201",  "C155", "",     NA, 5L,
  # a three-character trust code is not a site
  "e04_trustcode",    "RR8",    "",     "",     NA, 5L,
  # nor is a six-character GP practice code
  "e05_gpcode",       "B86012", "",     "",     NA, 5L,
  # nor are the not-known defaults
  "e06_default",      "89997",  "",     "",     NA, 5L,
  "e06_default",      "X99999", "",     "",     NA, 5L,
  # the case the code arrives in is not guaranteed
  "e07_lowercase",    "rj121",  "",     "",     NA, 5L,
  # this row is about the patient's bowel cancer: its site must not be used
  "e08_bowel_row",    "RQQ01",  "C182", "",     NA, 5L,
  # a C16 row on a C15 patient is kept: at the junction the registry and the map
  # can describe one tumour differently, and only non-OG rows are thrown out
  "e09_c16_on_c15",   "RQQ02",  "C161", "",     NA, 5L,
  # a bowel row again, said in SNOMED rather than ICD-O
  "e10_snomed_other", "RQQ03",  "",     "",     222000002, 5L,
  # and SNOMED can confirm the tumour as well as rule it out
  "e11_snomed_own",   "RB201",  "",     "",     111000001, 5L,
  # nothing said about the tumour, trust no help: morphology decides
  "e12_morph_break",  "RB101",  "",     "8140", NA, 5L,
  "e12_morph_break",  "RC101",  "",     "8070", NA, 5L,
  # nothing to separate them but how often each appears
  "e13_count_break",  "RD101",  "",     "",     NA, 5L,
  "e13_count_break",  "RD101",  "",     "",     NA, 5L,
  "e13_count_break",  "RE101",  "",     "",     NA, 5L,
  # nothing at all to separate them
  "e14_order_break",  "RG101",  "",     "",     NA, 5L,
  "e14_order_break",  "RF101",  "",     "",     NA, 5L,
  # the registry has no trust to compare against
  "e16_no_trust",     "RH101",  "",     "",     NA, 5L,
  # topography arrives with its description tacked on
  "e17_topog_text",   "RJ103",  "C155: LOWER THIRD OF OESOPHAGUS", "", NA, 5L,
  # both sides on offer and no trust to help: the row naming the patient's own
  # site has to win, or restricting to C15 later sends them to the wrong hospital
  "e18_own_beats_og", "RQQ04",  "C161", "",     NA, 5L,
  "e18_own_beats_og", "RB301",  "C155", "",     NA, 5L,
  # a COSD patient the registry has never heard of
  "x99_not_in_rapid", "RZZ01",  "",     "",     NA, 5L)

cosd_eg <- bind_rows(cosd_eg, teach)

out <- run_build(rapid_eg, cosd_eg)
val <- function(id, col) out[[col]][out$patient_pseudo_id == id]

expect("every registry row survives the merge", nrow(out) == nrow(rapid_eg))
expect("a COSD-only patient is not added",
       !"x99_not_in_rapid" %in% out$patient_pseudo_id)

expect("a clean code is kept, with the trust supporting it",
       val("e01_clean", "site_dx_code") == "RJ101" &&
         as.character(val("e01_clean", "site_dx_basis")) == "trust matches")
expect("a row naming the patient's own tumour site is recognised",
       as.character(val("e02_confirmed", "site_dx_basis")) == "tumour confirmed")
expect("a confirmed code beats one that only matches the trust",
       val("e03_confirm_wins", "site_dx_code") == "RA201")
expect("a three-character trust code is not taken as a site",
       is.na(val("e04_trustcode", "site_dx_code")))
expect("a six-character GP practice code is not taken as a site",
       is.na(val("e05_gpcode", "site_dx_code")))
expect("the not-known defaults are not taken as a site",
       is.na(val("e06_default", "site_dx_code")))
expect("lower case is tidied up rather than thrown away",
       val("e07_lowercase", "site_dx_code") == "RJ121")

expect("a site from someone's bowel cancer row is not used",
       is.na(val("e08_bowel_row", "site_dx_code")))
expect("a C16 row on a C15 patient is kept, not dropped",
       val("e09_c16_on_c15", "site_dx_code") == "RQQ02")
expect("a SNOMED code meaning a cancer that is not OG rules the row out",
       is.na(val("e10_snomed_other", "site_dx_code")))
expect("a row for the other side of the junction ranks below a confirmed one",
       as.character(val("e09_c16_on_c15", "site_dx_basis")) != "tumour confirmed")
expect("the patient's own site beats the other side when both are offered",
       val("e18_own_beats_og", "site_dx_code") == "RB301" &&
         as.character(val("e18_own_beats_og", "site_dx_basis")) == "tumour confirmed")
expect("a SNOMED code meaning the patient's own cancer confirms the row",
       val("e11_snomed_own", "site_dx_code") == "RB201" &&
         as.character(val("e11_snomed_own", "site_dx_basis")) == "tumour confirmed")

expect("morphology breaks a tie", val("e12_morph_break", "site_dx_code") == "RB101")
expect("the code appearing on more rows breaks a tie",
       val("e13_count_break", "site_dx_code") == "RD101")
expect("an otherwise dead tie is settled the same way every time",
       val("e14_order_break", "site_dx_code") == "RF101")
expect("a patient with no COSD rows gets no site",
       is.na(val("e15_no_cosd", "site_dx_code")) &&
         val("e15_no_cosd", "site_dx_found") == FALSE)
expect("topography with its description tacked on still matches",
       as.character(val("e17_topog_text", "site_dx_basis")) == "tumour confirmed")

expect("a single code is not flagged as a choice",
       val("e01_clean", "site_dx_n_codes") == 1L &&
         val("e01_clean", "site_dx_ambiguous") == FALSE)
expect("competing codes are flagged as a choice",
       val("e03_confirm_wins", "site_dx_n_codes") == 2L &&
         val("e03_confirm_wins", "site_dx_ambiguous") == TRUE)

expect("trust agreement is TRUE when the site sits in the registry's trust",
       isTRUE(val("e01_clean", "site_dx_trust_agrees")))
expect("trust agreement is FALSE when it does not",
       isFALSE(val("e12_morph_break", "site_dx_trust_agrees")))
expect("trust agreement is missing, not FALSE, when there is no trust to check",
       is.na(val("e16_no_trust", "site_dx_trust_agrees")))

# the order the rows sit in the file must not change anything
out_rev <- run_build(rapid_eg, cosd_eg[rev(seq_len(nrow(cosd_eg))), ])
expect("reversing the COSD row order changes nothing",
       identical(out$site_dx_code, out_rev$site_dx_code) &&
         identical(as.character(out$site_dx_basis),
                   as.character(out_rev$site_dx_basis)))

# with the SNOMED reading switched off, the SNOMED-only rows go back to saying
# nothing: the bowel row is no longer ruled out, and the own-cancer row is no
# longer confirmed
out_nos <- run_build(rapid_eg, cosd_eg, use_snomed = FALSE)
val2 <- function(id, col) out_nos[[col]][out_nos$patient_pseudo_id == id]
expect("without SNOMED, the non-OG row is no longer ruled out",
       val2("e10_snomed_other", "site_dx_code") == "RQQ03")
expect("without SNOMED, the own-cancer row is no longer confirmed",
       as.character(val2("e11_snomed_own", "site_dx_basis")) == "no support")

# turning the non-OG rule off should bring those patients back, and change
# nobody else
out_keep <- run_build(rapid_eg, cosd_eg, drop_non_og = FALSE)
val3 <- function(id, col) out_keep[[col]][out_keep$patient_pseudo_id == id]
expect("with the non-OG rule off, the bowel rows come back",
       val3("e08_bowel_row", "site_dx_code") == "RQQ01" &&
         val3("e10_snomed_other", "site_dx_code") == "RQQ03")
expect("with the non-OG rule off, nobody else moves",
       {
         # the teaching patients above carry bowel rows of their own, so they
         # are expected to move too; compare the worked examples only
         moved <- c("e08_bowel_row", "e10_snomed_other")
         same <- function(d) {
           d <- d[startsWith(d$patient_pseudo_id, "e") &
                    !d$patient_pseudo_id %in% moved, ]
           d$site_dx_code[order(d$patient_pseudo_id)]
         }
         identical(same(out), same(out_keep))
       })

# refusing the unsupported codes should drop just those patients
out_r2 <- run_build(rapid_eg, cosd_eg, max_rank = 2L)
expect("site_max_rank = 2 drops the unsupported picks",
       is.na(out_r2$site_dx_code[out_r2$patient_pseudo_id == "e12_morph_break"]) &&
         is.na(out_r2$site_dx_code[out_r2$patient_pseudo_id == "e16_no_trust"]))
expect("site_max_rank = 2 leaves the supported picks alone",
       out_r2$site_dx_code[out_r2$patient_pseudo_id == "e01_clean"] == "RJ101" &&
         out_r2$site_dx_code[out_r2$patient_pseudo_id == "e02_confirmed"] == "RJ101")

# -----------------------------------------------------------------------------
# The published map, and the version rule
# -----------------------------------------------------------------------------
cat("\n  the published map\n")

# what R/reference/10_fetch_snomed_map.R would hand over for these numbers
eg_map <- tibble(snomed = c("111000001", "222000002"),
                 site3 = c("C15", "C18"),
                 icd10_targets = c("C15.9", "C18.9"))

# with the map in place the meanings no longer have to be worked out, so the
# teaching rows are not needed and the switch can be off
out_map <- run_build(rapid_eg, cosd_eg, use_snomed = FALSE, map = eg_map)
valm <- function(id, col) out_map[[col]][out_map$patient_pseudo_id == id]
expect("the map rules out a row for a cancer that is not OG",
       is.na(valm("e10_snomed_other", "site_dx_code")))
expect("the map confirms a row for the patient's own cancer",
       valm("e11_snomed_own", "site_dx_code") == "RB201" &&
         as.character(valm("e11_snomed_own", "site_dx_basis")) == "tumour confirmed")
expect("the map and the worked-out meanings agree",
       identical(out$site_dx_code[out$patient_pseudo_id == "e11_snomed_own"],
                 valm("e11_snomed_own", "site_dx_code")))

cat("\n  the SNOMED version rule\n")

# the same number, but flagged as written in an older SNOMED. It must be
# ignored: in SNOMED II, 111000001 is not oesophageal cancer.
cosd_old <- cosd_eg %>%
  mutate(snomed_version = if_else(patient_pseudo_id == "e11_snomed_own",
                                  1L, snomed_version))
out_old <- run_build(rapid_eg, cosd_old, use_snomed = FALSE, map = eg_map)
expect("a number written in an older SNOMED is not looked up in the map",
       as.character(out_old$site_dx_basis[
         out_old$patient_pseudo_id == "e11_snomed_own"]) == "no support")

cosd_unk <- cosd_eg %>%
  mutate(snomed_version = if_else(patient_pseudo_id == "e11_snomed_own",
                                  99L, snomed_version))
out_unk <- run_build(rapid_eg, cosd_unk, use_snomed = FALSE, map = eg_map)
expect("a number with no version given is still looked up",
       as.character(out_unk$site_dx_basis[
         out_unk$patient_pseudo_id == "e11_snomed_own"]) == "tumour confirmed")

cat("\n  the ODS site-to-trust map\n")

# e01_clean offers site RJ101, whose first three characters RJ1 are the
# registry's diagnosing trust - so without any map it already counts as a trust
# match. Build an ODS map that says RJ101 is actually operated by a DIFFERENT
# trust, RX9. Now the first three characters no longer decide it, and the site
# should no longer count as matching the trust.
map_moved <- tibble(site_code = "RJ101", parent_trust = "RX9",
                    trust_is_prefix = "FALSE", parent_from_ods = "TRUE",
                    status = "Active", is_hospital_site = "TRUE",
                    record_class = "RC2", primary_role = "RO198",
                    predecessor = NA_character_, found = "TRUE")
out_moved <- run_build(rapid_eg, cosd_eg, use_snomed = FALSE,
                       site_trust = map_moved)
vmoved <- function(id, col) out_moved[[col]][out_moved$patient_pseudo_id == id]
expect("the ODS map's trust is used as the site's trust",
       vmoved("e01_clean", "site_dx_trust") == "RX9")
expect("a site ODS moves out of the registry's trust stops matching it",
       as.character(vmoved("e01_clean", "site_dx_basis")) == "no support")

# and the opposite: a site whose first three characters are NOT the registry
# trust, but which ODS says the trust operates, should now count as a match.
# e16_no_trust's registry has no trust, so use a fresh minimal case.
map_rescue <- tibble(site_code = "RQQ02", parent_trust = "RJ1",
                     trust_is_prefix = "FALSE", parent_from_ods = "TRUE",
                     status = "Active", is_hospital_site = "TRUE",
                     record_class = "RC2", primary_role = "RO198",
                     predecessor = NA_character_, found = "TRUE")
# e09_c16_on_c15 is a C15 patient in trust RJ1 whose only code is RQQ02 (prefix
# RQQ). Without the map that code does not match the trust; with the map saying
# RQQ02 is run by RJ1, it should.
out_rescue <- run_build(rapid_eg, cosd_eg, use_snomed = FALSE,
                        site_trust = map_rescue)
vres <- function(id, col) out_rescue[[col]][out_rescue$patient_pseudo_id == id]
expect("a site ODS places inside the registry's trust now matches it",
       as.character(vres("e09_c16_on_c15", "site_dx_basis")) == "trust matches")

expect("with no ODS map, the trust is still the first three characters",
       {
         o <- run_build(rapid_eg, cosd_eg, use_snomed = FALSE, site_trust = NULL)
         o$site_dx_trust[o$patient_pseudo_id == "e01_clean"] == "RJ1"
       })

# the field arriving as text, which is what we have asked NDRS for
cosd_txt <- cosd_eg %>%
  mutate(diagnosis_snomedct = if_else(is.na(diagnosis_snomedct), NA_character_,
                                      sprintf("%.0f", diagnosis_snomedct)))
out_txt <- run_build(rapid_eg, cosd_txt, use_snomed = FALSE, map = eg_map)
expect("a SNOMED field arriving as text works the same as a number",
       identical(out_map$site_dx_code, out_txt$site_dx_code))

# a long number in a numeric field has already lost precision upstream; the
# build cannot fix it but must say so
expect("an over-long SNOMED number in a numeric field is warned about",
       {
         cosd_big <- cosd_eg %>%
           mutate(diagnosis_snomedct = if_else(
             patient_pseudo_id == "e11_snomed_own", 1.653e16, diagnosis_snomedct))
         w <- character()
         withCallingHandlers(
           run_build(rapid_eg, cosd_big, use_snomed = FALSE, map = eg_map),
           warning = function(x) { w <<- c(w, conditionMessage(x))
           invokeRestart("muffleWarning") })
         any(grepl("too large for the numeric field", w))
       })

expect("the build runs with no map file at all",
       {
         o <- run_build(rapid_eg, cosd_eg, map = NULL)
         nrow(o) == nrow(rapid_eg)
       })

# =============================================================================
# B. Stand-in data at full size
# =============================================================================
cat("\nB. stand-in data\n")

if (!exists("run_sim_check")) run_sim_check <- TRUE

if (run_sim_check) {
  dir_sim <- tempfile("og_sim_")
  sim_scale <- 1
  invisible(capture.output(
    sys.source(file.path(dir_build, "90_simulate_inputs.R"), envir = new.env())))
  
  sim_rapid <- haven::read_dta(file.path(dir_sim,
                                         "20260212_Rapidtumour_linked_2026SOTN_clean_OG_postPT.dta"))
  sim_cosd <- haven::read_dta(file.path(dir_sim,
                                        "20260212_all_cosddiagnosis_rapid_202601_OG.dta"))
  truth <- readRDS(file.path(dir_sim, "site_truth.rds"))
  sim_map <- read.csv(file.path(dir_sim, "snomed_og_lookup.csv"),
                      colClasses = "character")
  sim_site_trust <- read.csv(file.path(dir_sim, "site_trust_map.csv"),
                             colClasses = "character")
  
  sim_out <- run_build(sim_rapid, sim_cosd, map = sim_map,
                       site_trust = sim_site_trust)
  
  expect("every made-up registry row survives", nrow(sim_out) == nrow(sim_rapid))
  expect("one row per patient comes out", !anyDuplicated(sim_out$patient_pseudo_id))
  expect("no code escapes the five-character rule",
         all(nchar(na.omit(sim_out$site_dx_code)) == 5))
  expect("no default code escapes",
         !any(na.omit(sim_out$site_dx_code) %in% c("89997", "89999", "X9999")))
  # the trust is the site's first three characters, except where the ODS map
  # deliberately places a site under a different trust - so the rule is that the
  # trust is either the prefix or what the ODS map says
  moved_sites <- sim_site_trust %>%
    filter(trust_is_prefix == "FALSE") %>%
    select(site_dx_code = site_code, ods_trust = parent_trust)
  trust_check <- sim_out %>%
    filter(!is.na(site_dx_code)) %>%
    left_join(moved_sites, by = "site_dx_code") %>%
    mutate(ok = if_else(is.na(ods_trust),
                        substr(site_dx_code, 1, 3) == site_dx_trust,
                        site_dx_trust == ods_trust))
  expect("the trust is the site prefix, or the ODS trust where ODS moved it",
         all(trust_check$ok))
  
  scored <- sim_out %>%
    select(patient_pseudo_id, site_dx_code, site_dx_basis, site_dx_found) %>%
    left_join(truth, by = "patient_pseudo_id") %>%
    mutate(right_site  = site_dx_found & site_dx_code == true_site,
           right_trust = site_dx_found & substr(site_dx_code, 1, 3) == true_trust)
  
  cat("\n  by tumour site:\n")
  scored %>%
    group_by(tumour_site) %>%
    summarise(patients = n(),
              coverage = round(100 * mean(site_dx_found), 1),
              right_site = round(100 * mean(right_site[site_dx_found]), 1),
              right_trust = round(100 * mean(right_trust[site_dx_found]), 1),
              .groups = "drop") %>%
    as.data.frame() %>%
    print(row.names = FALSE)
  
  cat("\n  by how the code was chosen:\n")
  scored %>%
    filter(site_dx_found) %>%
    group_by(site_dx_basis) %>%
    summarise(patients = n(),
              right_site = round(100 * mean(right_site), 1),
              right_trust = round(100 * mean(right_trust), 1),
              .groups = "drop") %>%
    as.data.frame() %>%
    print(row.names = FALSE)
  cat("\n")
  
  c15 <- scored %>% filter(tumour_site == "C15")
  cover <- mean(c15$site_dx_found)
  right_site <- mean(c15$right_site[c15$site_dx_found])
  right_trust <- mean(c15$right_trust[c15$site_dx_found])
  
  # 84.5% of patients are in COSD by construction and about 87% of their rows
  # carry a usable code, so coverage should sit in the low eighties.
  expect("C15 coverage is in the expected range", cover > 0.76 && cover < 0.86)
  expect("the C15 trust is right almost always", right_trust > 0.90)
  expect("the C15 site is right most of the time", right_site > 0.60)
  expect("a confirmed code is more accurate than an unsupported one",
         {
           by_basis <- scored %>% filter(site_dx_found) %>%
             group_by(site_dx_basis) %>% summarise(a = mean(right_site))
           a <- setNames(by_basis$a, as.character(by_basis$site_dx_basis))
           a[["tumour confirmed"]] > a[["no support"]]
         })
  
  # 03 does no work of its own - it re-reads what 02 wrote and reports on it - so
  # the check here is only that it runs to the end and writes its review files,
  # against a build with the maps and one without.
  cat("\n  diagnostics script\n")
  run_diag <- function(map, site_trust) {
    d <- tempfile("og_diag_"); dir.create(d)
    assign("dir_out", d, envir = globalenv())
    assign("read_rapid", function() sim_rapid, envir = globalenv())
    assign("read_cosd",  function() sim_cosd,  envir = globalenv())
    assign("read_snomed_map", function() map, envir = globalenv())
    assign("read_site_trust_map", function() site_trust, envir = globalenv())
    # 03 reads the cohort 02 wrote, so run 02 into the same folder first
    invisible(capture.output(suppressMessages(
      sys.source(file.path(dir_build, "02_add_site_of_diagnosis.R"), new.env()))))
    invisible(capture.output(suppressMessages(
      sys.source(file.path(dir_build, "03_site_diagnostics.R"), new.env()))))
    file.exists(file.path(d, "diag_field_effect.txt"))
  }
  expect("the diagnostics run through with both maps present",
         run_diag(sim_map, sim_site_trust))
  expect("the diagnostics run through with neither map",
         run_diag(NULL, NULL))
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
  quit(status = 1, save = "no")
}
restore_session()
cat("All checks passed.\n")