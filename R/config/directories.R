# =============================================================================
# R/config/directories.R  -  where the files live, for the whole project
# -----------------------------------------------------------------------------
# The one place the folders are named. Sourced first by each stage's parameter
# file, so a path is defined once rather than copied between the site-code build,
# the CWT merge and the deviance analysis. It holds directories and file paths
# and nothing else - no functions, no library() calls: those belong to the
# stages (their _load_packages.R and _helpers.R) and to shared/utils.R.
#
# Every path can be set before this file is sourced - the check scripts point
# them at temporary folders - so only unset values take the defaults.
# =============================================================================

# -----------------------------------------------------------------------------
# Roots
# -----------------------------------------------------------------------------
# dir_raw    the folder holding the 20260212 extracts as supplied by NDRS - not
#            this project's to move
# dir_out    every derived, PATIENT-LEVEL file the project produces. Each stage
#            writes its handoff file here and the next reads it back. This is
#            deliberately NOT inside the repository: it lives on the restricted
#            data drive, so patient data can never be committed or synced with
#            the code.
# dir_ref    reference lookups (the SNOMED map, the site-to-trust map, the
#            valid-trust list) - built once, not patient data, so they live in
#            the repo.
# dir_debug  small, AGGREGATE intermediates (no patient_pseudo_id) kept for
#            debugging. Also in the repo.
# dir_sim    the synthetic tree the test run writes to - fixtures for every
#            stage. In the repo but .gitignored, so the fakes are never pushed.
if (!exists("dir_raw")) dir_raw <- "W:/_DATA/IainTimmins/2026 OG SOTN data"
if (!exists("dir_out")) dir_out <- "W:/_DATA/IainTimmins/2026 OG SOTN data/Analysis"
if (!exists("dir_ref")) dir_ref <- "Data/reference"
if (!exists("dir_debug")) dir_debug <- "Data/intermediates"
if (!exists("dir_sim")) dir_sim <- "Data/sim"
for (.dir_to_make in c(dir_out, dir_ref, dir_debug))
  dir.create(.dir_to_make, recursive = TRUE, showWarnings = FALSE)

# a note for anyone adding to this file: every stage's 01 uses plain source(),
# whose default (local = FALSE) evaluates the sourced code in .GlobalEnv - not
# in the caller's environment, however deeply nested the source() calls are.
# That means any variable assigned at the top level here (loop variables
# included) lands in the global environment and can collide with anything else
# using the same short name at the top level of a session. Prefixing throwaway
# names with a dot (as above) keeps them out of the way; a genuine setting
# (dir_out, dir_ref, ...) is meant to be global and is named accordingly.

# a quick, loud check that no stage has been left writing patient data into the
# repo by accident - dir_out should be an absolute path (a drive letter, a root,
# or a UNC path), not something relative that resolves inside the working copy.
# The test scripts set dir_out to a tempfile(), which is absolute, so they pass.
if (!grepl("^([A-Za-z]:[\\\\/]|/|\\\\\\\\)", dir_out))
  stop("dir_out (", dir_out, ") looks like a path relative to the repository, ",
       "not the restricted data drive. Patient-level data must not be written ",
       "inside the repo - check dir_out.", call. = FALSE)

# dir_transfer is the one folder on the analysis server that is not encrypted at
# rest, so it is the only place output meant to leave the server can be written.
# Each stage's results (tables, figures, logs, aggregate summaries) go here when
# it is set. It holds nothing patient-identifying beyond the pseudonymised
# patient_pseudo_id already used throughout, and even that should not normally
# appear - transfer-folder outputs are aggregate results, not patient rows.
#
# Defaults to off (NULL). A run_*.R script sets it to a real S: path for a real
# run; a simulated or test run leaves it NULL, so made-up numbers never land in
# the real folder.
if (!exists("dir_transfer")) dir_transfer <- NULL
if (!is.null(dir_transfer))
  dir.create(dir_transfer, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# The raw extracts
# -----------------------------------------------------------------------------
path_rapid_dta <- file.path(dir_raw,
                            "20260212_Rapidtumour_linked_2026SOTN_clean_OG_postPT.dta")
path_cosd_dta  <- file.path(dir_raw,
                            "20260212_all_cosddiagnosis_rapid_202601_OG.dta")
path_cwt_dta   <- file.path(dir_raw,
                            "20260212_all_cwt_rapid_202601_OG.dta")

# -----------------------------------------------------------------------------
# Reference lookups brought in from the internet machine
# -----------------------------------------------------------------------------
f_snomed_map     <- file.path(dir_ref, "snomed_og_lookup.csv")
f_site_trust_map <- file.path(dir_ref, "site_trust_map.csv")
f_valid_trusts   <- file.path(dir_ref, "valid_diagnosing_trusts.csv")

# -----------------------------------------------------------------------------
# The files each stage derives and hands to the next
# -----------------------------------------------------------------------------
# The chain runs:
#   site-code build ->  og_cohort_site.rds
#   CWT merge       ->  og_cohort_cwt.rds
#   deviance        ->  its own outputs
# Naming them here keeps the handoffs in one place, so a later stage never has to
# guess what the earlier one called its output.
f_site_lookup <- file.path(dir_out, "og_site_of_diagnosis_lookup.rds")
f_cohort_site <- file.path(dir_out, "og_cohort_site.rds")
f_cohort_cwt  <- file.path(dir_out, "og_cohort_cwt.rds")

# -----------------------------------------------------------------------------
# What counts as OG
# -----------------------------------------------------------------------------
# The two cancer sites the whole project is about: C15 oesophagus, C16 stomach.
og_sites <- c("C15", "C16")