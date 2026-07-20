# =============================================================================
# 01  Define parameters  -  the CWT merge and waiting-time derivation
# -----------------------------------------------------------------------------
# The judgement calls specific to this stage: how a CWT treatment period is
# grouped into a broad modality, which modality groups can legitimately stop the
# clock for each treatment pathway, and the tolerances used when a waiting time
# is accepted or set aside. Sourced by 02 and 03 here.
#
# Shared paths and helpers come from common/paths_and_helpers.R. The stage reads
# og_cohort_site.rds (the registry with the site of diagnosis) and the CWT
# extract, and writes og_cohort_cwt.rds.
#
# The coding follows the National Cancer Waiting Times Monitoring Dataset
# Guidance v12.1 (NHS England, July 2025) for the treatment modality codes and
# the meaning of the decision-to-treat date, and follows the bowel-cancer rapid
# CWT script (Creating_CR_CWT_dataset_rapid.do) for the handling that is specific
# to this extract - the YMD date strings, the integer modality codes, and the
# five-day agreement tolerance between the CWT and registry treatment dates.
# =============================================================================

source("R/config/directories.R")
source("R/shared/utils.R")
source("R/merge_cwt_to_get_dtt/_load_packages.R")

# -----------------------------------------------------------------------------
# The CWT extract
# -----------------------------------------------------------------------------
# One row per treatment period, several per patient. The fields the merge uses:
#   treat_period_start  the decision-to-treat date (CANCER TREATMENT PERIOD
#                       START DATE in the guidance), as a DMY... no: YMD string
#   treat_start         the treatment start date, YMD string
#   modality            the treatment modality, as an integer (leading zeros
#                       dropped, so guidance code "02" arrives as 2, "23" as 23)
#   site_icd10          the cancer the treatment period is for, used to set aside
#                       periods that are plainly for a different cancer
if (!exists("path_cwt_dta"))
  path_cwt_dta <- file.path(dir_raw,
                            "20260212_all_cwt_rapid_202601_OG.dta")

# CWT stores its dates as YMD strings. The bowel script parses them the same way;
# note this differs from some older extracts that used DMY.
cwt_date_order <- "%Y-%m-%d"

# -----------------------------------------------------------------------------
# Modality groups
# -----------------------------------------------------------------------------
# The guidance modality codes, grouped as the analysis needs them. Codes arrive
# as integers, so 2 is guidance "02" and so on. Code 1 is the surgical code that
# the guidance retired in 2020 ("01"); it still appears in older rows and counts
# as surgery, which the bowel script also does (modality 1 / 23 / 24).
#
#   surgery       1 (retired), 23 (surgery), 24 (surgery, enabling)
#   chemo         2 (cytotoxic), 14 (other), 15 (immunotherapy)
#   hormone       3
#   chemort       4 (chemoradiotherapy)
#   radiotherapy  5 (teletherapy), 6 (brachytherapy), 13 (proton)
#   palliative    7 (specialist palliative), 8 (active monitoring),
#                 9 (non-specialist palliative)
#   other         97
#   declined      98
# Anything else is left ungrouped and set aside.
modality_groups <- list(
  surgery      = c(1L, 23L, 24L),
  chemo        = c(2L, 14L, 15L),
  hormone      = c(3L),
  chemort      = c(4L),
  radiotherapy = c(5L, 6L, 13L),
  palliative   = c(7L, 8L, 9L),
  other        = c(97L),
  declined     = c(98L))

modality_group_of <- function(code) {
  code <- suppressWarnings(as.integer(code))
  out <- rep(NA_character_, length(code))
  for (g in names(modality_groups))
    out[code %in% modality_groups[[g]]] <- g
  out
}

# The surgical code 1 ("01") was retired from 2020 and 23/24 took its place, with
# a transition window in between where systems were moving over. Within the
# window all three count; outside it, only the code that was live. The bowel
# script does not apply this window (it accepts 1/23/24 throughout), but keeping
# it makes the surgery grouping robust to the changeover, and it costs nothing
# when a row is already unambiguous.
surg_window_start <- as.Date("2020-01-01")
surg_window_end   <- as.Date("2021-06-30")

# -----------------------------------------------------------------------------
# Which modality groups can stop the clock for each pathway
# -----------------------------------------------------------------------------
# A patient's treatment pathway is derived first, from the registry (02). When
# the CWT records are brought in, only a record whose modality group is
# consistent with that pathway should supply the decision-to-treat date - a
# palliative CWT period should not set the DTT for a patient who had curative
# surgery. This table, and the primary group below, follow the OG treatment
# pathways rather than the bowel script's surgery-only approach.
pathway_consistent_groups <- list(
  "EMR/ESD only"                  = c("surgery", "other"),
  "EMR/ESD then surgery"          = c("surgery"),
  "Surgery + neoadjuvant chemoRT" = c("surgery", "chemort", "chemo", "radiotherapy"),
  "Surgery + neoadjuvant chemo"   = c("surgery", "chemo"),
  "Surgery + neoadjuvant RT"      = c("surgery", "radiotherapy", "chemort"),
  "Surgery + adjuvant chemo"      = c("surgery", "chemo"),
  "Surgery only"                  = c("surgery"),
  "Surgery + other"               = c("surgery", "other"),
  "Definitive chemoRT"            = c("chemort", "chemo", "radiotherapy"),
  "Curative RT only"              = c("radiotherapy", "chemort"),
  "Palliative chemo + RT"         = c("chemo", "radiotherapy", "chemort", "palliative"),
  "SACT only"                     = c("chemo", "hormone", "palliative"),
  "Palliative RT only"            = c("radiotherapy", "palliative"),
  "No treatment recorded"         = c("palliative", "other"))

# The single defining modality group for each pathway. When a patient has more
# than one pathway-consistent CWT record, the one matching this group is
# preferred before falling back to the earliest decision-to-treat date.
pathway_primary_group <- c(
  "EMR/ESD only"                  = "surgery",
  "EMR/ESD then surgery"          = "surgery",
  "Surgery + neoadjuvant chemoRT" = "chemort",
  "Surgery + neoadjuvant chemo"   = "chemo",
  "Surgery + neoadjuvant RT"      = "radiotherapy",
  "Surgery + adjuvant chemo"      = "surgery",
  "Surgery only"                  = "surgery",
  "Surgery + other"               = "surgery",
  "Definitive chemoRT"            = "chemort",
  "Curative RT only"              = "radiotherapy",
  "Palliative chemo + RT"         = "chemo",
  "SACT only"                     = "chemo",
  "Palliative RT only"            = "radiotherapy",
  "No treatment recorded"         = "palliative")

# -----------------------------------------------------------------------------
# Interval windows and tolerances
# -----------------------------------------------------------------------------
# tx_window_days       how long after diagnosis a treatment still counts as the
#                      first treatment for these purposes (nine months)
# dtt_min_offset       a DTT is allowed to sit a little before diagnosis, since
#                      the registry diagnosis date and the clinical decision can
#                      be days apart either way; beyond this it is not this
#                      patient's DTT
# treat_agree_days     how far the CWT treatment date and the registry treatment
#                      date may differ and still be taken as the same treatment.
#                      The bowel script uses five days and replaces the interval
#                      with the CWT date where the registry date sits before the
#                      DTT
# dtt_valid_max_days   a diagnosis-to-DTT or DTT-to-treatment interval longer
#                      than this is treated as not a reliable measure of this
#                      pathway's waiting time (the bowel script caps at 180)
tx_window_days     <- 270L
dtt_min_offset     <- -30L
treat_agree_days   <- 5L
dtt_valid_max_days <- 180L

# -----------------------------------------------------------------------------
# The clock-start anchor
# -----------------------------------------------------------------------------
# The OG timed pathway (guidance 2.2.11) can start the clock at the endoscopy
# where an abnormality is seen. So the waiting-time intervals are anchored on the
# endoscopy date where there is one, and fall back to the diagnosis date where
# there is not. Both the endoscopy-anchored and the diagnosis-anchored intervals
# are kept, so the choice is visible rather than hidden.
if (!exists("anchor_prefers_endoscopy")) anchor_prefers_endoscopy <- TRUE

cat("01 parameters set (CWT merge): writing to", dir_out,
    "| DTT window", tx_window_days, "days",
    "| treat agreement", treat_agree_days, "days\n")