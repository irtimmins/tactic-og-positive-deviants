# =============================================================================
# common/paths_and_helpers.R  -  shared by every stage of the project
# -----------------------------------------------------------------------------
# Where the data lives and the small jobs that every stage needs. Sourced first
# by each stage's own parameter file, so paths and helpers are defined in one
# place rather than copied between the site-code build, the CWT merge, and the
# deviance analysis.
#
# It defines only things that are genuinely common: the folders, the paths to
# the raw extracts and the derived files, the codes that count as OG, and a
# handful of tidying functions. Anything specific to one stage - what counts as a
# valid site code, how CWT waiting times are derived, and so on - lives in that
# stage's own parameter file, not here.
#
# Reads and writes nothing itself. Every path can be set before this file is
# sourced (the check scripts point them at temporary folders), so only unset
# values take the defaults.
# =============================================================================

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(magrittr)   # for %>%
})

# -----------------------------------------------------------------------------
# Where the files are
# -----------------------------------------------------------------------------
# dir_raw   the folder holding the 20260212 extracts as supplied
# dir_out   this project's data folder: every stage writes its derived files
#           here and the next stage reads them back
# dir_ref   reference files brought in from outside - the SNOMED map and the
#           site-to-trust map, both built on a machine with internet
if (!exists("dir_raw")) dir_raw <- "W:/_DATA/IainTimmins/2026 OG SOTN data"
if (!exists("dir_out")) dir_out <- "Data/OG"
if (!exists("dir_ref")) dir_ref <- "Data/reference"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# dir_transfer is the one folder on the analysis server that is not encrypted at
# rest, so it is the only place output meant to leave the server can be written.
# The master scripts copy their logs and small summary files there. It holds
# nothing patient-identifying beyond the pseudonymised patient_pseudo_id already
# used throughout, and only ever small summaries or a capped sample.
#
# Defaults to off (NULL). Set dir_transfer <- "S:/..." before sourcing to turn it
# on for a real run; leave it NULL for a simulated or test run, so made-up
# numbers never land in the real folder.
if (!exists("dir_transfer")) dir_transfer <- NULL
if (!is.null(dir_transfer))
  dir.create(dir_transfer, recursive = TRUE, showWarnings = FALSE)

# the raw extracts
path_rapid_dta <- file.path(dir_raw,
                            "20260212_Rapidtumour_linked_2026SOTN_clean_OG_postPT.dta")
path_cosd_dta  <- file.path(dir_raw,
                            "20260212_all_cosddiagnosis_rapid_202601_OG.dta")
path_cwt_dta   <- file.path(dir_raw,
                            "20260212_all_cwt_rapid_202601_OG.dta")

# reference files brought in from the internet machine
f_snomed_map     <- file.path(dir_ref, "snomed_og_lookup.csv")
f_site_trust_map <- file.path(dir_ref, "site_trust_map.csv")

# the files each stage derives and hands to the next. The chain runs:
#   site-code build ->  og_cohort_site.rds
#   CWT merge       ->  og_cohort_cwt.rds
#   deviance        ->  its own outputs
# Naming them all here keeps the handoffs in one place, so a later stage never
# has to guess what the earlier one called its output.
f_site_lookup <- file.path(dir_out, "og_site_of_diagnosis_lookup.rds")
f_cohort_site <- file.path(dir_out, "og_cohort_site.rds")
f_cohort_cwt  <- file.path(dir_out, "og_cohort_cwt.rds")

# -----------------------------------------------------------------------------
# What counts as OG
# -----------------------------------------------------------------------------
# The two cancer sites the whole project is about: C15 oesophagus, C16 stomach.
# The same pair the SNOMED map is cut down to in R/reference/10_fetch_snomed_map.
og_sites <- c("C15", "C16")

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