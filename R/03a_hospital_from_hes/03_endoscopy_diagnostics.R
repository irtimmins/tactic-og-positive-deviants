# =============================================================================
# 03  Endoscopy-site diagnostics
# -----------------------------------------------------------------------------
# 02 reads a site for each patient it can and reports the headline coverage. This
# runs straight after it and takes the same evidence apart, so a single run on
# the server answers the questions we would otherwise need several runs to
# explore:
#
#   A  the outcome for every anchored patient, in one place: site read, an
#      endoscopy matched but its site unusable, an APC endoscopy on the wrong
#      date, no APC endoscopy record, or not in the extract - the breakdown the
#      "do we need HES-OP" question turns on
#   B  the off-window group: how far off they are, which way, and how many more
#      would be caught if the window were widened. A group that is mostly late
#      by weeks is the outpatient-endoscopy-first case, which widening will not
#      fix well
#   C  the provider disagreements: sites whose provider is not the trust the
#      rapid record held. Most are the site-versus-provider difference we expect;
#      the ones where the provider disagrees too are the real ones to look at
#   D  coverage by tumour site, and the busiest sites as a sanity check
#
# It reads the cohort 02 wrote and re-reads the two extracts to rebuild the match
# at wider windows. It writes a few small txt files for review and prints the
# rest. It changes nothing the build produced.
#
# Run as part of 00_master, after 02.
# =============================================================================

source("R/03a_hospital_from_hes/01_define_parameters.R")

# tab-separated .txt, not .csv, for the same reason as the COSD stage: the
# server's transfer path appends an encrypted footer to .csv files, which leaves
# the content readable but breaks anything that reads to the end. .txt goes
# through untouched.
f_diag_outcome    <- file.path(dir_out, "diag_endoscopy_outcome.txt")
f_diag_window     <- file.path(dir_out, "diag_endoscopy_by_window.txt")
f_diag_provider   <- file.path(dir_out, "diag_endoscopy_provider_disagree.txt")

write_diag <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE,
              fileEncoding = "ASCII")
}

if (!exists("read_rapid"))
  read_rapid <- function() read_dta(path_rapid_dta)
if (!exists("read_hes"))
  read_hes <- function() {
    if (file.exists(f_hes_extract)) return(readRDS(f_hes_extract))
    read_dta(path_hes_apc_dta)
  }

section <- function(letter, title) {
  cat("\n", strrep("-", 74), "\n", letter, "  ", title, "\n",
      strrep("-", 74), "\n", sep = "")
}
show <- function(df) { print(as.data.frame(df), row.names = FALSE); invisible(df) }

# -----------------------------------------------------------------------------
# Rebuild the match, the same way 02 does, so the two agree
# -----------------------------------------------------------------------------
rapid     <- read_rapid()
hes       <- read_hes()
og_cohort <- readRDS(f_cohort_endoscopy)

op_cols     <- sprintf("opertn_%02d", 1:24)
opdate_cols <- sprintf("opdate_%02d", 1:24)

anchor <- rapid %>%
  mutate(patient_pseudo_id = as.character(patient_pseudo_id),
         endo_ref_date = hes_date(endodateHES)) %>%
  filter(endoHES == 1, !is.na(endo_ref_date)) %>%
  select(patient_pseudo_id, endo_ref_date, endotrustHES)

hes_cohort <- hes %>%
  mutate(patient_pseudo_id = as.character(patient_pseudo_id)) %>%
  filter(patient_pseudo_id %in% anchor$patient_pseudo_id) %>%
  mutate(epi_id = row_number(), epistart = hes_date(epistart))

endo_ops <- hes_cohort %>%
  select(epi_id, patient_pseudo_id, epiorder, sitetret, procode3,
         all_of(op_cols), all_of(opdate_cols)) %>%
  pivot_longer(cols = c(all_of(op_cols), all_of(opdate_cols)),
               names_to = c(".value", "slot"),
               names_pattern = "(opertn|opdate)_(\\d+)") %>%
  mutate(opcs = norm_opcs(opertn), opdate = hes_date(opdate)) %>%
  filter(opcs %in% opcs_diagnostic_endoscopy, !is.na(opdate))

# every endoscopy near the reference date, with no window applied yet, so the
# window can be varied below
near <- endo_ops %>%
  inner_join(anchor, by = "patient_pseudo_id") %>%
  mutate(days_off = as.integer(opdate - endo_ref_date))

# =============================================================================
section("A", "The outcome for every anchored patient")
# =============================================================================
cat(
  "The task: give every patient the hospital site of their diagnostic\n",
  "endoscopy. The rapid record already dates the endoscopy and names a trust\n",
  "for it, but the trust field mixes three-character provider codes with real\n",
  "five-character site codes, so it is not a site. HES-APC records the site\n",
  "(sitetret) on the episode, so the job is to find that endoscopy episode in\n",
  "HES and read its site.\n\n",
  "Only patients who appear in the HES-APC extract can be judged for an APC\n",
  "endoscopy record, so they are the denominator. Every one of them is placed\n",
  "in exactly one outcome below.\n", sep = "")

assessable <- og_cohort %>% filter(endo_in_hes)
outcome <- assessable %>%
  mutate(outcome = case_when(
    endoscopy_site_found ~ "site read from the endoscopy episode",
    endo_matched         ~ "endoscopy found, but its site was unusable",
    endo_has_apc         ~ "APC endoscopy present, but not within the window",
    TRUE                 ~ "no HES-APC endoscopy record at all")) %>%
  count(outcome, name = "patients") %>%
  mutate(pct = round(100 * patients / sum(patients), 1)) %>%
  arrange(desc(patients))

cat("\n")
show(outcome)
write_diag(outcome, f_diag_outcome)

not_in <- sum(!og_cohort$endo_in_hes & og_cohort$endoHES == 1 &
                !is.na(hes_date(og_cohort$endodateHES)))
cat("\nAnchored patients not in this extract (cannot be judged):", not_in, "\n")

# =============================================================================
section("B", "The off-window group, and whether a wider window helps")
# =============================================================================
# For each patient, the nearest endoscopy to the reference date. If that nearest
# one is outside the current window, how far outside, and which way? A group that
# is mostly late by weeks points at an outpatient endoscopy having set the
# reference date, which the APC extract cannot match well however wide the window.
nearest <- near %>%
  group_by(patient_pseudo_id) %>%
  slice_min(abs(days_off), n = 1, with_ties = FALSE) %>%
  ungroup()

off <- nearest %>% filter(abs(days_off) > endo_window_days)
cat("Patients whose nearest APC endoscopy is outside the",
    endo_window_days, "day window:", nrow(off), "\n")
if (nrow(off)) {
  cat("\nWhich way, and how far:\n")
  off %>%
    mutate(direction = if_else(days_off > 0,
                               "endoscopy after the reference date",
                               "endoscopy before the reference date"),
           band = cut(abs(days_off),
                      breaks = c(endo_window_days, 14, 30, 90, Inf),
                      labels = c("just past window..14", "15..30",
                                 "31..90", "over 90"))) %>%
    count(direction, band, name = "patients") %>%
    show()
}

# how many patients would gain a match as the window widens - the direct read on
# whether widening is worth it
window_gain <- tibble(window_days = c(endo_window_days, 14L, 30L, 60L, 90L)) %>%
  mutate(patients_matched = map_int(window_days, function(w)
    n_distinct(near$patient_pseudo_id[abs(near$days_off) <= w])))
cat("\nPatients matched as the window widens:\n")
show(window_gain)
write_diag(window_gain, f_diag_window)
cat("\nIf the numbers barely move past the current window, the misses are not a",
    "\nwindow problem - they are patients with no usable APC endoscopy date.\n")

# =============================================================================
section("C", "Provider disagreements")
# =============================================================================
# Among the sites read, the ones whose provider does not match the trust the
# rapid record held. sitetret is a site and endotrustHES is (usually) a provider,
# so the first three characters differing is common and expected; the ones where
# the provider code itself disagrees are the few worth an eye.
found <- og_cohort %>% filter(endoscopy_site_found)
disagree <- found %>%
  filter(!is.na(endo_provider_agrees), !endo_provider_agrees)
cat("Sites whose provider does not match endotrustHES:", nrow(disagree),
    "of", sum(!is.na(found$endo_provider_agrees)), "checkable\n")
if (nrow(disagree)) {
  # pattern-level only, no patient id, since this file may leave the server
  disagree %>%
    mutate(hes_provider = str_sub(endoscopy_provider, 1, 3),
           rapid_trust  = str_sub(endotrustHES, 1, 3)) %>%
    count(hes_provider, rapid_trust, name = "patients") %>%
    arrange(desc(patients)) %>%
    write_diag(f_diag_provider)
  cat("Top provider disagreements:\n")
  disagree %>%
    mutate(hes_provider = str_sub(endoscopy_provider, 1, 3),
           rapid_trust  = str_sub(endotrustHES, 1, 3)) %>%
    count(hes_provider, rapid_trust, name = "patients") %>%
    arrange(desc(patients)) %>%
    head(15) %>%
    show()
}

# =============================================================================
section("D", "Coverage by tumour site, and a sanity check")
# =============================================================================
if ("tumour_site" %in% names(og_cohort)) {
  cat("Coverage among anchored patients, by tumour site:\n")
  og_cohort %>%
    filter(endoHES == 1) %>%
    group_by(tumour_site) %>%
    summarise(anchored = n(),
              site_read = sum(endoscopy_site_found),
              pct = round(100 * mean(endoscopy_site_found), 1),
              .groups = "drop") %>%
    show()
}

cat("\nBusiest endoscopy sites (should be the big OG centres):\n")
og_cohort %>%
  filter(endoscopy_site_found) %>%
  count(endoscopy_site, name = "patients", sort = TRUE) %>%
  head(15) %>%
  show()

cat("\n", strrep("=", 74), "\n", sep = "")
cat("Diagnostics complete.\n")
written <- c(f_diag_outcome, f_diag_window,
             if (nrow(disagree)) f_diag_provider)
cat("Written for review:\n")
for (f in written) cat("  ", f, "\n")
if (!nrow(disagree))
  cat("Not written, nothing to put in it:", basename(f_diag_provider), "\n")