# =============================================================================
# 01  Define parameters
# -----------------------------------------------------------------------------
# Everything the build needs before any data is read: where the files are, what
# counts as a real site code, how a COSD row is judged to belong to the patient's
# oesophago-gastric tumour rather than to some other cancer, and how one site is
# chosen when a patient offers several. Sourced by every later script; it reads
# and writes nothing itself.
#
# The judgement calls all live here rather than being spread through the build,
# so they can be read, argued with and changed in one place.
# =============================================================================

suppressPackageStartupMessages({
  library(haven)
  library(tidyverse)
})

# -----------------------------------------------------------------------------
# Where the files are
# -----------------------------------------------------------------------------
# dir_raw is the folder holding the 20260212 extracts as supplied. dir_out is
# this project's data folder: everything derived is written there and every later
# script reads from there. Either can be set before this file is sourced (the
# checking script points them at temporary folders), so only unset values are
# filled in with the defaults below.
# dir_ref holds reference files brought in from outside, such as the SNOMED map
# fetched by R/reference/10_fetch_snomed_map.R on a machine with internet.
if (!exists("dir_raw")) dir_raw <- "W:/_DATA/IainTimmins/2026 OG SOTN data"
if (!exists("dir_out")) dir_out <- "Data/OG"
if (!exists("dir_ref")) dir_ref <- "Data/reference"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# dir_transfer is the one folder on the analysis server that is not encrypted at
# rest, so it is the only place output meant to leave the server can be written.
# 00_master.R copies the diagnostics csvs and the build log there. It holds
# nothing patient-identifying beyond the pseudonymised patient_pseudo_id already
# used throughout the project, and only ever small summary tables or a capped
# sample - never a full extract.
#
# Defaults to off (NULL). Set dir_transfer <- "S:/..." before sourcing to turn it
# on for a real run. Leave it NULL, or point it at a scratch folder, for a
# simulated or test run - otherwise made-up numbers land in the real folder next
# to real results, which is exactly the mix-up this default is there to prevent.
if (!exists("dir_transfer")) dir_transfer <- NULL
if (!is.null(dir_transfer))
  dir.create(dir_transfer, recursive = TRUE, showWarnings = FALSE)

path_rapid_dta <- file.path(dir_raw,
                            "20260212_Rapidtumour_linked_2026SOTN_clean_OG_postPT.dta")
path_cosd_dta  <- file.path(dir_raw,
                            "20260212_all_cosddiagnosis_rapid_202601_OG.dta")
path_cwt_dta   <- file.path(dir_raw,
                            "20260212_all_cwt_rapid_202601_OG.dta")   # for the next step, not this one

f_snomed_map   <- file.path(dir_ref, "snomed_og_lookup.csv")

f_site_lookup  <- file.path(dir_out, "og_site_of_diagnosis_lookup.rds")
f_cohort_site  <- file.path(dir_out, "og_cohort_site.rds")

# -----------------------------------------------------------------------------
# What counts as a site code
# -----------------------------------------------------------------------------
# site_code_of_diagnosis is meant to be the five-character code of the hospital
# site where the diagnosis was made, but the field is filled in loosely. Its 608
# distinct values in the January 2026 extract hold four different things:
#
#   RJ121, R0A02, NT216, 8CW14   five-character site codes    - what we want
#   RR8, B86, RRK, rrk           three-character trust codes  - trust, not site
#   B86012, Y03603, XXXXXX       six characters: GP practice codes (the
#                                referrer, not the site) and free text
#   89997, 89999, X99999         the "not known" and "not applicable" defaults
#
# Asking for exactly five characters removes the trust codes and the GP practice
# codes in one step. Codes that start with a digit are kept, because 8xxxx and
# 5xxxx are real historic hospital codes; only an all-digit code is refused, on
# the grounds that no real site code is five digits.
#
# If a list of real site codes is available - Sarah merges one in as
# temp_provider_covariates.dta and uses the non-matches to spot made-up codes -
# put the codes in site_code_allowlist and anything outside it is refused too.
# Left as NULL, only the rules above apply.
site_code_width    <- 5L
site_code_defaults <- c("89997", "89999", "99999", "XXXXX", "ZZZZZ", "X9999")
if (!exists("site_code_allowlist")) site_code_allowlist <- NULL

is_site_code <- function(x) {
  ok <- !is.na(x) &
    str_length(x) == site_code_width &
    str_detect(x, "^[A-Z0-9]{5}$") &
    !str_detect(x, "^[0-9]{5}$") &
    !x %in% site_code_defaults
  if (!is.null(site_code_allowlist)) ok <- ok & x %in% site_code_allowlist
  coalesce(ok, FALSE)
}

# -----------------------------------------------------------------------------
# Which COSD rows belong to the patient's OG tumour
# -----------------------------------------------------------------------------
# COSD is filed by patient, not by tumour. 130,170 rows carry 80,554 patients,
# and the topographies present run well past C15 and C16 - C18, C34, C50, C61 and
# so on - so some rows are plainly a different cancer, diagnosed somewhere else.
# There is no tumour identifier and no usable date in the file to join on, so
# each row has to be read for what it says about itself.
#
# Two fields say what the tumour was, and neither is filled in often:
#   topography_icdo3    on 8.2% of rows. Says it directly, in ICD-O.
#   diagnosis_snomedct  on 38.9% of rows. Says it in SNOMED, which means nothing
#                       to us on its own - a SNOMED number has to be looked up.
#
# So each row's cancer is worked out from up to three things, in this order:
#
#   1 the topography, when it is there. It is already in the language we want.
#   2 the SNOMED number, looked up in NHS England's SNOMED-to-ICD-10 map. The
#     map is published as part of the SNOMED CT UK Clinical Edition;
#     R/reference/10_fetch_snomed_map.R downloads it and writes the C15 and C16
#     codes to Data/reference/snomed_og_lookup.csv. Put that file in place and
#     the build finds it. This is the same job Sarah's hand-written list of
#     SNOMED codes does for bowel, done from the published source instead.
#   3 failing that, the meaning worked out from the data itself: on rows where
#     topography and SNOMED are both filled in, a SNOMED number sits next to a
#     topography, and where it does so consistently and often enough that can be
#     taken as what it means. This was the only route before the map; it is kept
#     because the map has known gaps - NHS England inactivated the ICD-10 map
#     for around 22,000 concepts in 2023 - and because it covers the older
#     SNOMED versions the map does not. The build reports how much each source
#     contributed, so its worth is visible rather than assumed.
#
#   snomed_min_rows    how many labelled rows a SNOMED number needs before its
#                      meaning is taken as settled
#   snomed_min_share   how much of the time it must sit next to the same
#                      topography. Below this it is treated as saying nothing.
if (!exists("snomed_from_data")) snomed_from_data <- TRUE
snomed_min_rows  <- 10L
snomed_min_share <- 0.90

# snomed_version says which SNOMED the number was written in. The NHS Data Model
# and Dictionary gives the codes as: 1 SNOMED II, 2 SNOMED 3, 3 SNOMED 3.5,
# 4 SNOMED RT, 5 SNOMED CT, 99 not known. In this extract 5 covers 82% of the
# rows that have a version, 99 covers 17%, and the older versions about 1%
# between them.
#
# Those older ones are different systems, not older editions of one system, so
# the same number means something else and the map does not apply. Their numbers
# are set aside. Rows that do not say which version was used are kept: an old
# number simply will not be in the map, so there is nothing to gain by dropping
# them in advance.
snomed_old_versions <- c(1L, 2L, 3L, 4L)

# Rows for a cancer that is not OG at all are dropped, not demoted. The site code
# on a row that says C18 is where that patient's bowel cancer was diagnosed; it
# says nothing about where their oesophageal cancer was diagnosed, and using it
# would put the patient at the wrong hospital.
#
# The question asked of each row is only "is this an OG tumour", not "is this the
# same side of the junction as the registry says". A C16 row on a C15 patient is
# kept. It has to be: at the gastro-oesophageal junction the map and the registry
# can describe the same tumour differently - every junction concept in the map
# reads as C16, cardia included - so dropping those rows would throw away the
# patient's own record on a coding technicality. Whether a patient is C15 or C16
# is settled later, from the registry, when the cohort is restricted. That is not
# this step's job.
#
# A row that names the patient's own site still ranks above one that names the
# other side; see the ranking below. Rows that say nothing about the tumour are
# kept: they are unknown, not wrong.
#
# Set drop_non_og_rows to FALSE only to measure what the rule costs.
if (!exists("drop_non_og_rows")) drop_non_og_rows <- TRUE

# what counts as OG. The same pair the map is cut down to in
# R/reference/10_fetch_snomed_map.R.
og_sites <- c("C15", "C16")

# -----------------------------------------------------------------------------
# Choosing between site codes
# -----------------------------------------------------------------------------
# After the wrong-cancer rows have gone, a patient may still offer more than one
# site code. They are ranked on how well supported each one is, best first:
#
#   1 tumour confirmed  a row carrying this code says, in ICD-O or in SNOMED,
#                       that it is about the patient's own tumour site
#   2 trust matches     the row says nothing about the tumour, but the site sits
#                       inside the trust the registry already records as having
#                       made the diagnosis, so a second source agrees with it
#   3 no support        the row says nothing about the tumour, and the trust
#                       disagrees or the registry has no trust to compare with
#
# Where two codes tie: a site inside the registry's trust wins, then a matching
# morphology, then whichever code appears on more rows, and last the code itself.
# That final step is only there so the answer never depends on the order the rows
# happen to sit in the file. Sarah keeps one at random at this point, which is
# the same decision made a different way.
#
# site_max_rank says how far down this list the build will go. 3 accepts
# everything; 2 leaves the unsupported codes out.
site_basis_levels <- c("tumour confirmed", "trust matches", "no support")
if (!exists("site_max_rank")) site_max_rank <- 3L

# -----------------------------------------------------------------------------
# Small shared jobs
# -----------------------------------------------------------------------------

# Tidy an organisation code: upper case, drop anything that is not a letter or a
# digit, and treat an empty string as missing rather than as a value.
tidy_org_code <- function(x) {
  x <- str_to_upper(str_trim(as.character(x)))
  x <- str_replace_all(x, "[^A-Z0-9]", "")
  na_if(x, "")
}

# The three-character cancer site from an ICD-O topography. The field usually
# holds a bare code but sometimes carries the description too ("C155: LOWER THIRD
# OF OESOPHAGUS"), and occasionally holds filler such as XXXXXX. Take the leading
# letter-digit-digit and nothing else; anything else becomes missing.
icd_site3 <- function(x) {
  x <- str_to_upper(str_trim(as.character(x)))
  str_extract(x, "^[A-Z][0-9]{2}")
}

# Morphology as four digits, however it arrives: str4 in COSD, str5 in the
# registry where the fifth character can be the behaviour.
morph4 <- function(x) {
  str_extract(str_trim(as.character(x)), "^[0-9]{4}")
}

# A SNOMED identifier as text. It is only ever compared with itself and with the
# map, never used as a number.
#
# If the field arrives as text, which is what we have asked NDRS for and what
# Sarah's bowel extract already has, it is used as it stands. Older extracts
# have it as a Stata double, which cannot hold the longest SNOMED numbers
# exactly - anything above about 9e15 is already approximate by the time it
# reaches us. Those are converted and a warning is printed, because a number
# that has lost a digit will not match the map and there is nothing this code
# can do to recover it.
snomed_id <- function(x) {
  if (is.character(x)) {
    x <- str_trim(x)
    return(if_else(x == "" | x == "." | x == "0", NA_character_, x))
  }
  x <- suppressWarnings(as.numeric(x))
  too_big <- sum(!is.na(x) & x > 2^53)
  if (too_big)
    warning(too_big, " SNOMED numbers are too large for the numeric field they ",
            "arrived in and may have lost precision. Ask for the field as text.",
            call. = FALSE)
  ifelse(is.na(x) | x <= 0, NA_character_, sprintf("%.0f", x))
}

# Stop early and clearly if a file is not the one we think it is. A missing
# column noticed here saves a puzzling empty result later.
check_input <- function(df, needed, label, path, min_rows = 1000L) {
  absent <- setdiff(needed, names(df))
  if (length(absent))
    stop(label, " has no ", paste(absent, collapse = ", "), " column.",
         "\n  - check that ", basename(path), " is the expected file.",
         call. = FALSE)
  if (nrow(df) < min_rows)
    stop(label, " has only ", nrow(df), " rows, expected at least ", min_rows,
         ".\n  - this looks like a part-read or a test file rather than the real",
         " one: ", basename(path), call. = FALSE)
  invisible(TRUE)
}

cat("01 parameters set: writing to", dir_out,
    "| site codes", site_code_width, "characters",
    "| SNOMED read from the data:", snomed_from_data, "\n")