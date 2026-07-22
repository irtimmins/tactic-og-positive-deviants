# =============================================================================
# 01  Define parameters  -  the endoscopy-site-from-HES build
# -----------------------------------------------------------------------------
# The judgement calls specific to this stage: which OPCS codes count as a
# diagnostic endoscopy, how close a HES operation date has to sit to the date the
# registry already records before it is taken as the same endoscopy, what counts
# as a real site code, and how one episode is chosen when a patient has several
# endoscopies near the reference date. Sourced by 02 and 03.
#
# The shared paths (where the files are), the generic utilities (tidy_org_code,
# check_input) and this stage's own parsing helpers (norm_opcs, hes_date) are
# brought in here, so they are defined once and 02/03 inherit them by sourcing
# this file.
# =============================================================================

source("R/config/directories.R")
source("R/shared/utils.R")
source("R/derive_hospital_code_from_hes/_load_packages.R")
source("R/derive_hospital_code_from_hes/_helpers.R")

# Every path this stage uses - path_rapid_dta, path_hes_apc_dta, f_hes_extract,
# f_endoscopy_lookup, f_cohort_endoscopy - comes from R/config/directories.R,
# sourced above, so the deviance stage can read f_endoscopy_lookup without
# sourcing this file to find out what it is called.

# -----------------------------------------------------------------------------
# What counts as a diagnostic endoscopy
# -----------------------------------------------------------------------------
# The non-therapeutic (diagnostic) endoscopy OPCS-4 codes from Appendix 6 - the
# same list the State of the Nation methodology and the GOLD pipeline use. A HES
# episode is a diagnostic endoscopy if any of its operation fields holds one of
# these. The rapid record's endoscopy date was itself built by searching
# opertn_01..opertn_24 for exactly these codes, so re-finding the endoscopy this
# way lands back on the same episode.
opcs_diagnostic_endoscopy <- c(
  "G142", "G143", "G145", "G147",
  "G152", "G153", "G154", "G156", "G157", "G158", "G159",
  "G161", "G162", "G168", "G169",
  "G172", "G173", "G188", "G189",
  "G191", "G198", "G199",
  "G201", "G202", "G208", "G209",
  "G214", "G215", "G218", "G219",
  "G422", "G432", "G433", "G435",
  "G441", "G443", "G445", "G446", "G448", "G449",
  "G451", "G452", "G454", "G458", "G459",
  "G462", "G463", "G468", "G469")

# -----------------------------------------------------------------------------
# How close the dates have to be
# -----------------------------------------------------------------------------
# For each endoscopy code found in HES, the endoscopy date is the opdate on the
# same operation slot (Table 4.3), and it is matched against endodateHES, the
# date already on the rapid record. Because endodateHES is itself the earliest
# such opdate, the match is essentially exact - on the January 2026 data all but
# one matched patient sat on day 0. The window is only slack for the odd date
# quirk; it could be tightened to 1 with almost no loss. A patient whose only APC
# endoscopy sits outside the window is left without a site here (it usually means
# their endodateHES was set by an earlier outpatient endoscopy, which this APC
# extract does not carry).
if (!exists("endo_window_days")) endo_window_days <- 7L

# -----------------------------------------------------------------------------
# What counts as a site code
# -----------------------------------------------------------------------------
# sitetret is the five-character site of treatment - the hospital, not the trust.
# It is what we are after, but like any HES org field it carries the odd default
# and test value ("00000", "89999", "12345") and, rarely, a short code. A usable
# site is five characters of letters and digits, is not all digits (no real site
# code is), and is not one of the known not-known defaults. An episode whose site
# fails this is treated as having no site, so a second endoscopy episode with a
# real site is preferred over it.
hes_site_width    <- 5L
hes_site_defaults <- c("00000", "89997", "89998", "89999", "99999",
                       "X99999", "XXXXX", "ZZZZZ")

is_hes_site <- function(x) {
  ok <- !is.na(x) &
    str_length(x) == hes_site_width &
    str_detect(x, "^[A-Z0-9]{5}$") &
    !str_detect(x, "^[0-9]{5}$") &
    !x %in% hes_site_defaults
  coalesce(ok, FALSE)
}

# -----------------------------------------------------------------------------
# Choosing between endoscopy episodes
# -----------------------------------------------------------------------------
# A patient can have more than one endoscopy episode inside the window. They are
# ordered, best first: an episode carrying a real site code before one that does
# not, then the operation date closest to the reference date, then the earliest
# date, then the lowest episode order. The last two steps are only there so the
# answer never depends on the order the rows happen to sit in the file.
if (!exists("prefer_sited_episode")) prefer_sited_episode <- TRUE

cat("01 parameters set: writing to", dir_out,
    "| endoscopy codes", length(opcs_diagnostic_endoscopy),
    "| window", endo_window_days, "days\n")
