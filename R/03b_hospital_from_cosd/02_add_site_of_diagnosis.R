# =============================================================================
# 02  Add the site of diagnosis
# -----------------------------------------------------------------------------
# The registry extract records the trust that made the diagnosis
# (diagnosis_trust, three characters) but not the hospital site within it. COSD
# does carry site_code_of_diagnosis, but it is filed by patient rather than by
# tumour, holds several rows for many patients, and mixes real site codes with
# trust codes, GP practice codes and "not known" defaults.
#
# This step works in four passes:
#   1  read each COSD row and decide which cancer it is about
#   2  throw away the rows that are about a different cancer, and the site codes
#      that are not really site codes
#   3  choose one code per patient from what is left, using the ranking in 01
#   4  put it on the registry extract as site_dx_code
#
# The choice is written out alongside the code - site_dx_basis, site_dx_n_codes,
# site_dx_trust_agrees - so a reader can see how firm each value is rather than
# having to take it on trust.
#
# Reads : the rapid tumour dta, the COSD diagnosis dta
# Writes: og_site_of_diagnosis_lookup.rds  (one row per patient, f_site_lookup)
#         og_cohort_site.rds               (the registry plus the site, f_cohort_site)
# =============================================================================

source("R/03b_hospital_from_cosd/01_define_parameters.R")

# read_rapid() and read_cosd() are here so the checking script can hand this
# script made-up data instead of the real files: define either one before
# sourcing and it will be used. Otherwise they read the dta files. COSD is read
# a few columns at a time - the file is 26 columns wide and we need six.
if (!exists("read_rapid"))
  read_rapid <- function() read_dta(path_rapid_dta)

if (!exists("read_cosd"))
  read_cosd <- function() read_dta(path_cosd_dta, col_select = c(
    patient_pseudo_id, site_code_of_diagnosis, topography_icdo3,
    morphology_icdo3, behaviour_icdo3, diagnosis_snomedct, snomed_version))

# the SNOMED map fetched by stage 02, if it is there
if (!exists("read_snomed_map"))
  read_snomed_map <- function() {
    if (!file.exists(f_snomed_map)) return(NULL)
    read.csv(f_snomed_map, colClasses = "character")
  }

# the ODS site-to-trust map from stage 01, if it is
# there. Without it the build falls back to reading the trust as the first three
# characters of the site code, which is right for most sites but not all.
if (!exists("read_site_trust_map"))
  read_site_trust_map <- function() {
    if (!file.exists(f_site_trust_map)) return(NULL)
    read.csv(f_site_trust_map, colClasses = "character")
  }

# -----------------------------------------------------------------------------
# Read the two files
# -----------------------------------------------------------------------------
rapid <- read_rapid()
check_input(rapid,
            c("patient_pseudo_id", "tumour_site", "diagnosis_trust",
              "tumour_morphology_str"),
            "the rapid tumour extract", path_rapid_dta,
            min_rows = getOption("og_min_input_rows", 1000L))

cosd_raw <- read_cosd()
check_input(cosd_raw,
            c("patient_pseudo_id", "site_code_of_diagnosis",
              "topography_icdo3", "morphology_icdo3", "diagnosis_snomedct",
              "snomed_version"),
            "the COSD diagnosis extract", path_cosd_dta,
            min_rows = getOption("og_min_input_rows", 1000L))

cat("Read", nrow(rapid), "registry rows and", nrow(cosd_raw), "COSD rows\n")

# Everything below joins on the patient alone. That is only safe because the
# registry extract holds one row per patient; if it ever stops doing so the join
# would quietly multiply rows, so stop here instead.
if (anyDuplicated(rapid$patient_pseudo_id))
  stop("the rapid extract has more than one row for some patients, and this ",
       "script assumes one.", call. = FALSE)

# -----------------------------------------------------------------------------
# What the registry already knows about each patient
# -----------------------------------------------------------------------------
# Three things are used to judge a COSD row against: the tumour site (C15 or
# C16), the morphology, and the trust that made the diagnosis. diagnosis_trust is
# three characters in this extract although it is stored as str5, so it is cut to
# three rather than compared whole.
tumour_keys <- rapid %>%
  mutate(patient_pseudo_id = as.character(patient_pseudo_id),
         reg_site3  = icd_site3(tumour_site),
         reg_morph4 = morph4(tumour_morphology_str),
         reg_trust3 = str_sub(tidy_org_code(diagnosis_trust), 1, 3)) %>%
  select(patient_pseudo_id, reg_site3, reg_morph4, reg_trust3)

# -----------------------------------------------------------------------------
# Pass 1: what each COSD row says about itself
# -----------------------------------------------------------------------------
cosd <- cosd_raw %>%
  mutate(patient_pseudo_id = as.character(patient_pseudo_id),
         site_raw    = tidy_org_code(site_code_of_diagnosis),
         topog_site3 = icd_site3(topography_icdo3),
         cosd_morph4 = morph4(morphology_icdo3),
         snomed_ver  = suppressWarnings(as.integer(snomed_version)),
         snomed      = snomed_id(diagnosis_snomedct)) %>%
  # a number written in an older SNOMED means something else, so set it aside
  mutate(snomed = if_else(snomed_ver %in% snomed_old_versions,
                          NA_character_, snomed)) %>%
  select(patient_pseudo_id, site_raw, topog_site3, cosd_morph4, snomed)

# What each SNOMED number means. First the published map, then anything left
# over that the data itself can settle.
from_map <- read_snomed_map()

snomed_meaning <- tibble(snomed = character(), snomed_site3 = character(),
                         meaning_from = character())

if (!is.null(from_map)) {
  snomed_meaning <- from_map %>%
    mutate(snomed = as.character(snomed), snomed_site3 = as.character(site3),
           meaning_from = "map") %>%
    select(snomed, snomed_site3, meaning_from)
  cat("SNOMED map read from", f_snomed_map, ":", nrow(snomed_meaning), "codes\n")
} else {
  cat("No SNOMED map at", f_snomed_map,
      "- run stage 02 on a networked machine and copy",
      "the result over.\n")
}

if (snomed_from_data) {
  learned <- cosd %>%
    filter(!is.na(snomed), !is.na(topog_site3)) %>%
    count(snomed, topog_site3, name = "rows") %>%
    group_by(snomed) %>%
    mutate(labelled = sum(rows), share = rows / labelled) %>%
    slice_max(share, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    filter(labelled >= snomed_min_rows, share >= snomed_min_share) %>%
    mutate(snomed_site3 = topog_site3, meaning_from = "data") %>%
    select(snomed, snomed_site3, meaning_from)
  
  # where both sources name a code, the map wins; note any disagreement, since
  # that is worth knowing about rather than silently resolving
  clash <- inner_join(snomed_meaning, learned, by = "snomed",
                      suffix = c("_map", "_data")) %>%
    filter(snomed_site3_map != snomed_site3_data)
  if (nrow(clash))
    cat("SNOMED codes where the map and the data disagree:", nrow(clash),
        "- the map is used\n")
  
  snomed_meaning <- bind_rows(snomed_meaning,
                              anti_join(learned, snomed_meaning, by = "snomed"))
}

cat("SNOMED numbers with a settled meaning:", nrow(snomed_meaning), "\n")

# How much of the data the map actually reaches. This is the test that matters:
# a map of 145 codes is worth nothing if the numbers in this extract are not
# among them.
if (!is.null(from_map)) {
  in_data <- unique(na.omit(cosd$snomed))
  cat("SNOMED numbers in this extract:", length(in_data), "| in the map:",
      sum(in_data %in% from_map$snomed), "\n")
}

# The row's cancer: what the topography says, or failing that what the SNOMED
# number means. Neither present means the row says nothing, which is not the
# same as saying it is a different cancer.
cosd <- cosd %>%
  left_join(snomed_meaning, by = "snomed") %>%
  mutate(row_site3 = coalesce(topog_site3, snomed_site3),
         said_by = case_when(!is.na(topog_site3) ~ "topography",
                             !is.na(snomed_site3) ~ paste("SNOMED,", meaning_from),
                             TRUE ~ "nothing said"))

cat("\nWhat told us which cancer each COSD row is about:\n")
cosd %>%
  count(said_by, name = "rows") %>%
  mutate(pct = round(100 * rows / sum(rows), 1)) %>%
  arrange(desc(rows)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# -----------------------------------------------------------------------------
# Pass 2: drop the rows we must not use
# -----------------------------------------------------------------------------
cosd <- cosd %>%
  inner_join(tumour_keys, by = "patient_pseudo_id") %>%
  mutate(
    site_ok     = is_site_code(site_raw),
    same_tumour = !is.na(row_site3) & !is.na(reg_site3) & row_site3 == reg_site3,
    other_og    = !is.na(row_site3) & row_site3 %in% og_sites & !same_tumour,
    not_og      = !is.na(row_site3) & !row_site3 %in% og_sites)

# what happened to every COSD row for a patient we have, so the losses are on
# the record rather than left implicit
row_audit <- cosd %>%
  mutate(outcome = case_when(
    not_og & drop_non_og_rows               ~ "dropped: not an OG cancer",
    is.na(site_raw)                         ~ "dropped: site code blank",
    !site_ok & str_length(site_raw) < site_code_width ~
      "dropped: trust code, not a site",
    !site_ok & str_length(site_raw) > site_code_width ~
      "dropped: GP practice code or free text",
    !site_ok                                ~ "dropped: a not-known default",
    same_tumour                             ~ "kept: tumour confirmed",
    other_og                                ~ "kept: OG, other side of the junction",
    TRUE                                    ~ "kept: tumour not stated")) %>%
  count(outcome, name = "rows") %>%
  mutate(pct = round(100 * rows / sum(rows), 1)) %>%
  arrange(desc(rows))

cat("\nEvery COSD row for a patient in the registry extract:\n")
print(as.data.frame(row_audit), row.names = FALSE)

usable <- cosd %>%
  filter(site_ok) %>%
  filter(!(not_og & drop_non_og_rows))

# -----------------------------------------------------------------------------
# The trust each site sits under
# -----------------------------------------------------------------------------
# A site code belongs to a trust. Most of the time that trust is the first three
# characters of the code, but not always - some sites are run by a trust whose
# code is not their own first three characters, usually after a merger. Which
# trust a site really belongs to matters here, because a site is trusted more
# when it agrees with the trust the registry says made the diagnosis.
#
# Stage 01 asks ODS - the NHS body that assigns
# these codes - for the real operating trust of every site code in the data, and
# writes it to a file. If that file is present the real trust is used; if it is
# not, the build falls back to the first three characters and says so, so the
# result is never silently worse than it looks.
site_trust_map <- read_site_trust_map()
if (!is.null(site_trust_map)) {
  cat("\nODS site-to-trust map read from", f_site_trust_map, ":",
      nrow(site_trust_map), "site codes\n")
  trust_of_site <- site_trust_map %>%
    mutate(site_raw = site_code, ods_trust = parent_trust) %>%
    select(site_raw, ods_trust)
  usable <- usable %>%
    left_join(trust_of_site, by = "site_raw") %>%
    mutate(site_trust3 = if_else(is.na(ods_trust),
                                 str_sub(site_raw, 1, 3), ods_trust))
  n_via_ods <- sum(!is.na(usable$ods_trust))
  n_moved   <- sum(!is.na(usable$ods_trust) &
                     usable$ods_trust != str_sub(usable$site_raw, 1, 3))
  cat("  rows whose trust ODS confirmed:", n_via_ods, "| of those, rows where",
      "ODS puts the site under a different trust than its first three",
      "characters:", n_moved, "\n")
} else {
  cat("\nNo ODS site-to-trust map at", f_site_trust_map,
      "- using the first three characters of each site code as its trust.\n")
  usable <- usable %>% mutate(site_trust3 = str_sub(site_raw, 1, 3))
}

# -----------------------------------------------------------------------------
# Pass 3: choose one code per patient
# -----------------------------------------------------------------------------
# One line per patient and site code, carrying the evidence for that code. A code
# counts as confirmed if any of its rows confirms the tumour; the morphology and
# the trust are read the same way.
candidates <- usable %>%
  mutate(morph_hit = !is.na(cosd_morph4) & !is.na(reg_morph4) &
           cosd_morph4 == reg_morph4,
         trust_hit = !is.na(site_trust3) & !is.na(reg_trust3) &
           site_trust3 == reg_trust3) %>%
  group_by(patient_pseudo_id, site_dx_code = site_raw) %>%
  summarise(n_rows    = n(),
            confirmed = any(same_tumour),
            morph_ok  = any(morph_hit),
            trust_ok  = any(trust_hit),
            site_trust3 = first(site_trust3),
            .groups   = "drop") %>%
  mutate(site_rank = case_when(confirmed ~ 1L,
                               trust_ok  ~ 2L,
                               TRUE      ~ 3L))

# site_dx_n_codes is how many codes the patient was offering when the choice was
# made, so a single clear answer can be told apart from a pick between several.
site_lookup <- candidates %>%
  filter(site_rank <= site_max_rank) %>%
  group_by(patient_pseudo_id) %>%
  mutate(site_dx_n_codes = n_distinct(site_dx_code)) %>%
  arrange(site_rank, desc(trust_ok), desc(morph_ok), desc(n_rows), site_dx_code,
          .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(site_dx_trust     = site_trust3,
         site_dx_basis     = factor(site_basis_levels[site_rank],
                                    levels = site_basis_levels),
         site_dx_n_codes   = as.integer(site_dx_n_codes),
         site_dx_ambiguous = site_dx_n_codes > 1L) %>%
  select(patient_pseudo_id, site_dx_code, site_dx_trust, site_dx_basis,
         site_dx_n_codes, site_dx_ambiguous)

saveRDS(site_lookup, f_site_lookup)

# -----------------------------------------------------------------------------
# Pass 4: put it on the registry extract
# -----------------------------------------------------------------------------
# site_dx_trust_agrees compares the chosen site's trust with the trust the
# registry already holds. It is the main read on whether the merge is sound: a
# low rate would mean the two sources are describing different episodes of care.
# It is left missing where either side is missing, so it is never confused with a
# real disagreement.
rows_in <- nrow(rapid)

og_cohort <- rapid %>%
  mutate(patient_pseudo_id = as.character(patient_pseudo_id)) %>%
  left_join(site_lookup, by = "patient_pseudo_id") %>%
  mutate(reg_trust3 = str_sub(tidy_org_code(diagnosis_trust), 1, 3),
         site_dx_found = !is.na(site_dx_code),
         site_dx_trust_agrees = if_else(is.na(site_dx_trust) | is.na(reg_trust3),
                                        NA, site_dx_trust == reg_trust3)) %>%
  select(-reg_trust3)

stopifnot(nrow(og_cohort) == rows_in)

saveRDS(og_cohort, f_cohort_site)

# -----------------------------------------------------------------------------
# What came out
# -----------------------------------------------------------------------------
cat("\nSite of diagnosis found for", sum(og_cohort$site_dx_found), "of",
    nrow(og_cohort), "patients (",
    round(100 * mean(og_cohort$site_dx_found), 1), "% )\n")

cat("\nHow the code was chosen, by tumour site:\n")
og_cohort %>%
  filter(site_dx_found) %>%
  count(tumour_site, site_dx_basis, name = "patients") %>%
  group_by(tumour_site) %>%
  mutate(pct = round(100 * patients / sum(patients), 1)) %>%
  ungroup() %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\nCoverage by tumour site:\n")
og_cohort %>%
  group_by(tumour_site) %>%
  summarise(patients = n(),
            with_a_site = sum(site_dx_found),
            pct = round(100 * mean(site_dx_found), 1),
            more_than_one_offered = sum(site_dx_ambiguous, na.rm = TRUE),
            .groups = "drop") %>%
  as.data.frame() %>%
  print(row.names = FALSE)

agree <- og_cohort$site_dx_trust_agrees
cat("\nChosen site sits in the registry's diagnosing trust:",
    sum(agree, na.rm = TRUE), "of", sum(!is.na(agree)), "that can be checked (",
    round(100 * mean(agree, na.rm = TRUE), 1), "% )\n")
cat("Distinct site codes in the cohort:",
    n_distinct(og_cohort$site_dx_code, na.rm = TRUE), "across",
    n_distinct(og_cohort$site_dx_trust, na.rm = TRUE), "trusts\n")

cat("\nSaved", f_site_lookup, "and", f_cohort_site, "\n")
cat("02 complete. Next: the CWT merge.\n")

# ---- things worth a look, when there is a reason to (uncomment) -------------
# what the SNOMED codes were taken to mean, in full
# print(snomed_meaning %>% count(snomed_site3, sort = TRUE), n = Inf)
# the C15 patients whose site sits outside the trust the registry records
# og_cohort %>% filter(tumour_site == "C15", site_dx_found, !site_dx_trust_agrees) %>%
#   count(diagnosis_trust, site_dx_code, sort = TRUE) %>% print(n = 30)
# how many patients each site code accounts for
# og_cohort %>% count(site_dx_code, sort = TRUE) %>% print(n = 30)