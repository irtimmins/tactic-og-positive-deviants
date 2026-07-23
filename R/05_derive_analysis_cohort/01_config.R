# =============================================================================
# 01  configuration for deriving the analysis cohort
# -----------------------------------------------------------------------------
# Sourced at the top of every script in this stage, and by stage 06, which
# inherits all of it and adds the modelling settings on top. Everything here is
# about WHICH PATIENTS AND HOSPITALS are in the analysis; nothing here is about
# how they are then modelled.
#
# What this stage does, in plain terms:
#   Draws the curatively-treated oesophageal analysis cohort from the CWT-merged
#   data through a patient funnel (pathway, tumour site, stage, a valid
#   endoscopy-to-decision-to-treat time) and a hospital funnel (a volume floor
#   on the site of diagnosis), writes the attrition flowchart and the
#   per-hospital counts, and describes the cohort in Table 1.
#
# It reads the CWT-merged cohort (og_cohort_cwt.rds) produced by stage 04, so
# run that stage first.
# ============================================================================

# paths ----------------------------------------------------------------------
# Shared paths - the same file every stage in the pipeline sources. Defines
# dir_raw, dir_out (patient-level data, on the restricted drive), dir_ref
# (reference lookups), dir_debug (aggregate, non-patient intermediates, in the
# repo) and dir_transfer (the results-transfer area, off by default). See
# R/config/directories.R for the full explanation.
source("R/config/directories.R")
source("R/shared/utils.R")   # title_case, for the ODS-derived hospital name
source("R/05_derive_analysis_cohort/_load_packages.R")
# no _helpers.R here: the standardisation, weighting, shrinkage and ranking
# functions are stage 06's, and no script in this stage calls them.

# input: the CWT-merged patient cohort. f_cohort_cwt is defined in the shared
# paths file, so this stage never has to guess what the merge stage called it.
in_rds <- f_cohort_cwt

# reference lookups (site->trust map, valid-trust list) - dir_ref, from the
# shared paths file.
site_trust_map_csv <- file.path(dir_ref, "site_trust_map.csv")
valid_trusts_csv   <- file.path(dir_ref, "valid_diagnosing_trusts.csv")

# --- where this stage's own outputs go, on the same three-way split as the rest
# of the pipeline -------------------------------------------------------------
# patient-level (carries patient_pseudo_id): dir_out, the restricted drive.
cohort_rds <- file.path(dir_out, "pd_cohort.rds")   # 02 -> 03/04, patient rows
# key shareable outputs (aggregate tables and figures): dir_transfer, if it has
# been set for a real run (see block B of run_master.R). If it has not, fall back
# to dir_debug, so nothing here ever requires the S: drive to exist just to run
# the code.
#
# The fallback is dir_debug itself (Output/local for a real run, Output/sim/
# intermediates for a simulated one) rather than a subfolder named for testing.
# An earlier version wrote to "positive_deviance_test_output", which was fine
# while only test runs landed there - but this is also where a REAL run's
# results go when someone forgets to set dir_transfer, and a folder named
# "test_output" holding genuine results is a good way to have them ignored or
# deleted.
if (is.null(dir_transfer)) {
  warning("dir_transfer is not set, so this stage's results are being written ",
          "to ", dir_debug, " rather than the S: transfer area. That is correct ",
          "for a simulated or exploratory run. If this is a real run whose ",
          "results need to leave the server, stop, set dir_transfer <- ",
          "transfer_root, and run again.", call. = FALSE)
  out_dir <- dir_debug
} else {
  out_dir <- dir_transfer
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

flow_csv       <- file.path(out_dir, "pd_flow.csv")            # patient attrition
hosp_excl_csv  <- file.path(out_dir, "pd_hospitals_excluded.csv")

# analysis window ------------------------------------------------------------
# the diagnosis-date window is NOT set here. sotn_cohort (applied in 02) already
# carries the SoTN audit period as part of its own definition, so the analysis
# window is whatever period sotn_cohort represents - no separate window filter
# is applied on top of it. A future, different reporting period (e.g. Jan 2023
# to Mar 2025) would replace the sotn_cohort restriction, not add to it; see the
# commented-out block in 02_build_cohort.R where that would go.

# inclusion: which treatment pathways count as curative --------------------
# The tx_pathway values (from the merge stage) that define the curative cohort.
# Resection, neoadjuvant-then-surgery, adjuvant, and the non-surgical definitive
# pathways are in; EMR/ESD only is included as a curative endoscopic pathway.
# Everything else (palliative, SACT only, no treatment) is excluded as
# non-curative.
pathways_include <- c(
  "Surgery only",
  "Surgery + neoadjuvant chemo",
  "Surgery + neoadjuvant RT",
  "Surgery + neoadjuvant chemoRT",
  "Surgery + adjuvant chemo",
  "Definitive chemoRT",
  "Curative RT only",
  "EMR/ESD only")

# held out, and excluded, pending investigation. Reported in the attrition so
# the count is visible rather than silently dropped.
#
# "EMR/ESD then surgery": these patients had EMR/ESD as their FIRST treatment
# but went on to a major resection. Stage 04 sets this pathway's clock-stop to
# the surgery date, so their decision-to-treat and the treatment-side intervals
# (decision-to-treat-to-treatment, endoscopy-to-treatment) are all measured to
# surgery rather than to the EMR/ESD that actually started their treatment.
# Including them would put a materially different quantity in the same column as
# everyone else's, so they stay out until the anchoring is resolved.
#
# To include them properly, the clock-stop for this pathway has to be reassigned
# to the EMR/ESD date in stage 04 (R/04_add_dtt_to_cohort/03_cwt_merge.R); once
# that is done they can be moved into pathways_include above and grouped with
# "EMR/ESD only" for reporting.
pathways_flagged <- c("EMR/ESD then surgery")

# tumour site restriction ----------------------------------------------------
# the analysis is oesophageal only: keep C15, drop C16 (gastric). tumour_site is
# the 3-character ICD-10 site in the registry. Set include_sites to c("C15","C16")
# to analyse both, or c("C16") for gastric.
include_sites <- c("C15")

# route to diagnosis ---------------------------------------------------------
# tfinal_route is the registry's route-to-diagnosis field. Two of its six values
# are excluded from the analysis:
#
#   Emergency presentation  the patient came in acutely, so there is usually no
#                           orderly endoscopy-to-decision-to-treat pathway to
#                           measure - the wait, where one exists at all, is not
#                           comparable with an elective route and would make a
#                           hospital seeing more emergencies look artificially
#                           different.
#   Unknown                 the route is not recorded, so it cannot be confirmed
#                           as elective.
#
# The four routes kept - TWW, GP referral, Other outpatient, Inpatient elective
# - are all elective pathways where the interval means the same thing. On the
# January 2026 extract this drops roughly a quarter of the registry (Emergency
# presentation 23.9%, Unknown 2.8%), so the step is a large one and is reported
# in the attrition with a breakdown by route.
#
# An exclusion list rather than an inclusion list, so a route added to the
# registry in future is kept and shows up in Table 1 rather than being silently
# dropped; the funnel prints what it removed, so an unexpected value is visible.
routes_exclude <- c("Emergency presentation", "Unknown")

# hospital unit and volume floor ---------------------------------------------
# the unit of analysis is the hospital that did the diagnostic endoscopy - the
# 5-character site code derived from HES in stage 03a (endoscopy_site), joined
# onto the cohort in 02_build_cohort.R. This is where the endoscopy-to-DTT clock
# actually started, which is the more natural hospital for that outcome, and it
# is a real hospital site rather than the trust-level diagnosis_trust the
# registry carries directly.
#
# The alternative is site_dx_code, the COSD-derived site of diagnosis from stage
# 03b - the hospital that made the diagnosis, which need not be the same
# hospital that did the endoscopy. Switch hosp_var to that if the analysis
# question is about the diagnosing hospital rather than the endoscopy hospital.
#
# Either way, its operating trust is resolved from the ODS-derived site->trust
# map (site_trust_map.csv, built in stage 01): parent_trust is what ODS says
# operates the site, which differs from the site's own first three characters
# where ODS relocated it. A site must have at least min_per_year diagnoses in
# EVERY calendar year of the window to be analysed; smaller sites are excluded
# (and listed in pd_hospitals_excluded.csv).
hosp_var       <- "endoscopy_site"
min_per_year   <- 5

# the ODS site->trust map from the reference stage. Keyed on site_code, giving
# parent_trust (the operating trust). If absent, the trust falls back to the
# site's first three characters, with a warning.
site_trust_map_csv <- file.path(dir_ref, "site_trust_map.csv")

# valid diagnosing trusts ----------------------------------------------------
# only sites whose operating trust is a recognised OG-cancer diagnoser are
# analysed. The list is the NOGCA State of the Nation data-quality table (English
# NHS trusts with recorded OG diagnoses), read from valid_trusts_csv - a file
# with a trust_code column, prepared from the SoTN workbook. Kept in dir_ref
# alongside the other reference lookups (it is a general lookup, not specific to
# a single analysis run). If the file is absent the restriction is skipped with
# a warning, so the pipeline still runs.
valid_trusts_csv <- file.path(dir_ref, "valid_diagnosing_trusts.csv")

# outcome --------------------------------------------------------------------
# time from diagnostic endoscopy to the decision to treat, in days. The merge
# stage anchors this strictly on the diagnostic endoscopy; patients with no
# endoscopy date have no endoscopy-anchored wait and are excluded (and counted)
# in 02. Only positive (> 0) waits are analysed, and waits beyond max_wait are
# set aside as implausible.
outcome_var    <- "wt_endo_to_dtt"
max_wait       <- 180
drop_zero_wait <- TRUE

# ----------------------------------------------------------------------------
# The case-mix coding and the modelling settings (balance strictness, shrinkage
# priors, candidate thresholds) are NOT here - no script in this stage uses
# them. They live in R/06_identify_positive_deviants/01_config.R, which sources
# this file and adds them on top.