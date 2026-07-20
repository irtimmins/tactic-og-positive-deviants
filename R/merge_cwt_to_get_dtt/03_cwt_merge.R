# =============================================================================
# 03  Merge the CWT records and derive the waiting times
# -----------------------------------------------------------------------------
# Reads the pathway cohort from 02 and the raw CWT extract, then:
#   1. parse the CWT dates and assign each row a broad modality group, applying
#      the surgery transition-window rule
#   2. keep in-window records whose modality group is consistent with the
#      patient's pathway
#   3. get to one row per patient - the pathway's primary modality first, then
#      the earliest decision-to-treat date
#   4. attach the DTT and derive the waiting-time family, anchored on endoscopy
#      where there is one and diagnosis otherwise
#   5. mark each interval valid or, where not, why not
#
# The decision-to-treat date is the one thing the registry does not carry, so
# that is what the CWT merge exists to supply. The treatment date and the
# pathway are already known from 02; the CWT treatment date is used only to check
# the two sources agree and, where the registry treatment date sits before the
# DTT, to stand in for it (following the bowel rapid script).
#
# Reads : og_cohort_pathway.rds, the CWT .dta
# Writes: og_cohort_cwt.rds
# =============================================================================

source("R/merge_cwt_to_get_dtt/01_define_parameters.R")

f_cohort_pathway <- file.path(dir_out, "og_cohort_pathway.rds")

if (!exists("read_cohort_pathway"))
  read_cohort_pathway <- function() readRDS(f_cohort_pathway)
if (!exists("read_cwt"))
  read_cwt <- function() read_dta(path_cwt_dta, col_select = c(
    patient_pseudo_id, treat_period_start, treat_start, modality, site_icd10,
    org_dec_to_treat, org_treat_start))

pathway <- read_cohort_pathway()
cwt_raw <- read_cwt()
cat("Read", nrow(pathway), "patients and", nrow(cwt_raw), "CWT rows\n")

# -----------------------------------------------------------------------------
# 1. Parse the CWT rows, assign a modality group, apply the surgery window
# -----------------------------------------------------------------------------
# Dates are YMD strings in this extract. The modality is an integer already.
cwt <- cwt_raw %>%
  mutate(
    patient_pseudo_id = as.character(patient_pseudo_id),
    cwt_dtt_date   = as.Date(as.character(treat_period_start), format = cwt_date_order),
    cwt_treat_date = as.Date(as.character(treat_start),        format = cwt_date_order),
    mod_group      = modality_group_of(modality),
    cwt_site3      = str_sub(str_to_upper(as.character(site_icd10)), 1, 3))

# transition-window rule: the retired surgical code 1 ("01") does not count after
# the window, and 23/24 do not count before it, so each surgery code only counts
# where it was the live coding.
cwt <- cwt %>%
  mutate(mod_group = case_when(
    modality == 1L  & cwt_treat_date > surg_window_end   ~ NA_character_,
    modality %in% c(23L, 24L) & cwt_treat_date < surg_window_start ~ NA_character_,
    TRUE ~ mod_group))

# keep usable, dated, non-declined records
cwt <- cwt %>%
  filter(!is.na(mod_group), mod_group != "declined", !is.na(cwt_dtt_date))

cat("Usable CWT rows after grouping:", nrow(cwt), "\n")

# -----------------------------------------------------------------------------
# 2. Candidate DTT rows: in-window, and consistent with the patient's pathway
# -----------------------------------------------------------------------------
anchor_keys <- pathway %>%
  select(patient_pseudo_id, tumour_site, diagnosis_date, endoscopy_date,
         tx_pathway, first_tx_date)

cand <- cwt %>%
  inner_join(anchor_keys, by = "patient_pseudo_id") %>%
  mutate(days_dx_to_dtt = as.integer(cwt_dtt_date - diagnosis_date)) %>%
  filter(days_dx_to_dtt >= dtt_min_offset, days_dx_to_dtt <= tx_window_days)

# is this record's modality group a plausible clock-stop for the pathway?
group_ok <- function(pathway_name, group) {
  ok <- logical(length(group))
  for (i in seq_along(group)) {
    allowed <- pathway_consistent_groups[[pathway_name[i]]]
    ok[i] <- !is.null(allowed) && group[i] %in% allowed
  }
  ok
}
cand <- cand %>%
  mutate(group_ok = group_ok(tx_pathway, mod_group),
         primary_group = pathway_primary_group[tx_pathway],
         is_primary = mod_group == primary_group)

# if a patient has any pathway-consistent record keep only those; otherwise keep
# all their in-window records, so the earliest DTT can still anchor
cand <- cand %>%
  group_by(patient_pseudo_id) %>%
  mutate(any_match = any(group_ok)) %>%
  ungroup() %>%
  filter((any_match & group_ok) | !any_match)

# -----------------------------------------------------------------------------
# 3. One row per patient: primary modality first, then earliest DTT
# -----------------------------------------------------------------------------
anchor <- cand %>%
  group_by(patient_pseudo_id) %>%
  arrange(desc(is_primary), cwt_dtt_date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(patient_pseudo_id, cwt_dtt_date, cwt_treat_date, cwt_modality = modality,
         cwt_mod_group = mod_group)

cat("Patients with a CWT decision-to-treat date:", nrow(anchor), "\n")

# -----------------------------------------------------------------------------
# 4. Attach to the pathway cohort and derive the waiting times
# -----------------------------------------------------------------------------
og <- pathway %>%
  left_join(anchor, by = "patient_pseudo_id")

# the clock-start: the diagnostic endoscopy. The OG timed pathway starts the
# clock at the endoscopy where an abnormality is seen (guidance 2.2.11), so the
# endoscopy date is the anchor. Patients with no endoscopy date have no
# endoscopy-anchored waiting time - wt_endo_to_dtt is left missing for them, and
# the reason column records why, so the loss is explicit rather than silent. The
# diagnosis-anchored intervals (wt_dx_to_dtt, wt_dx_to_tx) are kept alongside for
# any analysis that wants them.
di <- function(a, b) as.integer(a - b)   # date difference in days

og <- og %>%
  mutate(
    # the endoscopy-to-diagnosis lead-in, kept in its own right
    wt_endo_to_dx = di(diagnosis_date, endoscopy_date),
    
    # endoscopy anchor to the decision-to-treat date (missing where no
    # endoscopy). The intervals to the treatment date are derived below, once the
    # treatment date to use has been settled.
    wt_endo_to_dtt = di(cwt_dtt_date, endoscopy_date),
    wt_dtt_to_tx   = di(first_tx_date, cwt_dtt_date),
    wt_dx_to_dtt   = di(cwt_dtt_date, diagnosis_date))

# how well the CWT treatment date and the registry treatment date agree. Where
# they are close the DTT is trustworthy; where the registry treatment date sits
# just before the DTT, the bowel script substitutes the CWT treatment date, so
# the DTT-to-treatment interval is not negative on a few days' disagreement.
#
# The curative pathways have a registry clock-stop date (first_tx_date) to check
# against. The non-curative ones (palliative, SACT-only, no treatment) do not -
# there is no single curative act - so for those the CWT treatment date stands on
# its own as the treatment date, and there is nothing to disagree with.
og <- og %>%
  mutate(
    has_clock_stop = !is.na(first_tx_date),
    treat_gap = di(first_tx_date, cwt_treat_date),
    treat_agrees = case_when(
      !has_clock_stop & !is.na(cwt_treat_date) ~ TRUE,   # CWT date stands alone
      !is.na(treat_gap) & abs(treat_gap) <= treat_agree_days ~ TRUE,
      TRUE ~ FALSE),
    # the treatment date used for the intervals: the registry clock-stop where
    # there is one, else the CWT treatment date
    tx_date_used = coalesce(first_tx_date, cwt_treat_date),
    wt_dtt_to_tx = di(tx_date_used, cwt_dtt_date),
    wt_dtt_to_tx = if_else(!is.na(wt_dtt_to_tx) & wt_dtt_to_tx < 0 &
                             !is.na(cwt_treat_date),
                           di(cwt_treat_date, cwt_dtt_date), wt_dtt_to_tx),
    wt_endo_to_tx = di(tx_date_used, endoscopy_date),
    wt_dx_to_tx   = di(tx_date_used, diagnosis_date))

# -----------------------------------------------------------------------------
# 5. Validity, and where a waiting time is missing, why
# -----------------------------------------------------------------------------
# Rather than a bare missing value, each interval carries a reason it could not
# be measured. The reasons follow the bowel script's scheme, with one addition
# for the endoscopy-anchored interval:
#   ok                  a usable interval
#   no endoscopy date   the patient has no diagnostic endoscopy date, so no
#                       endoscopy-anchored waiting time (endoscopy interval only)
#   no CWT DTT          the patient has no CWT decision-to-treat date
#   treat disagree      the CWT and registry treatment dates are more than a few
#                       days apart, so the DTT is not reliably this treatment's
#   over cap            the interval is longer than the cap, not a reliable measure
#   non-positive        the interval is negative
# start_missing flags, per patient, that the interval's own clock-start date is
# absent; it is passed only for the endoscopy interval, so a missing endoscopy is
# named as such rather than lumped in as a generic "date missing".
# small helper so the case_when condition is a plain logical vector whether or
# not start_missing was supplied (recycled to length 1 when NULL)
isTRUE_vec <- function(x) if (is.null(x)) FALSE else (!is.na(x) & x)

dtt_reason <- function(interval, start_missing = NULL) {
  case_when(
    isTRUE_vec(start_missing)              ~ "no endoscopy date",
    is.na(og$cwt_dtt_date)                 ~ "no CWT DTT",
    !og$treat_agrees                       ~ "treat dates disagree",
    is.na(interval)                        ~ "date missing",
    interval > dtt_valid_max_days          ~ "over cap",
    interval < 0                           ~ "non-positive",
    TRUE                                   ~ "ok")
}

og <- og %>%
  mutate(
    dtt_to_tx_reason     = dtt_reason(wt_dtt_to_tx),
    endo_to_dtt_reason   = dtt_reason(wt_endo_to_dtt,
                                      start_missing = is.na(endoscopy_date)),
    dtt_valid = dtt_to_tx_reason == "ok" & endo_to_dtt_reason == "ok",
    # blank the intervals that are not reliable, keeping the reason column so the
    # loss is explained rather than silent
    wt_endo_to_dtt = if_else(endo_to_dtt_reason == "ok", wt_endo_to_dtt,
                             NA_integer_),
    wt_dtt_to_tx   = if_else(dtt_to_tx_reason == "ok", wt_dtt_to_tx,
                             NA_integer_))

# -----------------------------------------------------------------------------
# Report and save
# -----------------------------------------------------------------------------
n <- nrow(og)
cat("\nMerged cohort:", n, "patients\n")
cat("With a CWT decision-to-treat date:",
    sum(!is.na(og$cwt_dtt_date)),
    sprintf("(%.1f%%)\n", 100 * mean(!is.na(og$cwt_dtt_date))))
cat("With a valid diagnosis/endoscopy-to-DTT-to-treatment chain:",
    sum(og$dtt_valid, na.rm = TRUE),
    sprintf("(%.1f%%)\n", 100 * mean(og$dtt_valid, na.rm = TRUE)))

cat("\nWhy the DTT-to-treatment interval is missing, where it is:\n")
og %>%
  count(dtt_to_tx_reason, name = "patients") %>%
  mutate(pct = round(100 * patients / sum(patients), 1)) %>%
  arrange(desc(patients)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

saveRDS(og, f_cohort_cwt)
cat("\nSaved", f_cohort_cwt, "\n")
cat("03 complete. Next: the positive-deviance analysis.\n")