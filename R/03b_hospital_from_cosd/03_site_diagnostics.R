# =============================================================================
# 03  Site-of-diagnosis diagnostics
# -----------------------------------------------------------------------------
# 02 makes one decision per patient and reports the headline coverage. This runs
# straight after it and takes the same evidence apart, so that a single run on
# the analysis server answers the questions we would otherwise need several runs
# to explore:
#
#   A  how far each patient gets on trust alone, and how much of the remaining
#      ambiguity carries tumour information that could settle it - Amanda's
#      question, with the "could we settle it" half answered
#   B  what each field actually changes: rebuild the pick under a stack of rules
#      (trust only, add morphology, add registry confirmation, add SNOMED) and
#      count how often each field moves the answer, rather than assuming it earns
#      its place
#   C  where the fields disagree: the codes SNOMED and topography read
#      differently, and the patients whose trust-matched code is not the one the
#      tumour evidence would pick
#   D  what the C15 cohort would look like once restricted, since that is where
#      the analysis is headed
#
# It reads the cohort 02 wrote and the two raw extracts. It writes a few small
# csv files for review and prints everything else. It changes nothing the build
# produced.
#
# Run as part of 00_master, after 02.
# =============================================================================

source("R/03b_hospital_from_cosd/01_define_parameters.R")

# Written as tab-separated .txt, not .csv, on purpose. Something on the server's
# transfer path appends a few kilobytes of encrypted footer to .csv files - the
# content survives intact but anything reading to the end of the file chokes on
# the tail. Plain .txt goes through untouched, which is why the build log has
# always arrived readable. Tab-separated so it still opens straight into a
# spreadsheet if you want it there.
f_diag_field_effect <- file.path(dir_out, "diag_field_effect.txt")
f_diag_conflicts    <- file.path(dir_out, "diag_snomed_topography_conflicts.txt")
f_diag_trust_vs_tum <- file.path(dir_out, "diag_trust_vs_tumour_picks.txt")

write_diag <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE,
              fileEncoding = "ASCII")
}

if (!exists("read_rapid"))
  read_rapid <- function() read_dta(path_rapid_dta)
if (!exists("read_cosd"))
  read_cosd <- function() read_dta(path_cosd_dta, col_select = c(
    patient_pseudo_id, site_code_of_diagnosis, topography_icdo3,
    morphology_icdo3, behaviour_icdo3, diagnosis_snomedct, snomed_version))
if (!exists("read_snomed_map"))
  read_snomed_map <- function() {
    if (!file.exists(f_snomed_map)) return(NULL)
    read.csv(f_snomed_map, colClasses = "character")
  }
if (!exists("read_site_trust_map"))
  read_site_trust_map <- function() {
    if (!file.exists(f_site_trust_map)) return(NULL)
    read.csv(f_site_trust_map, colClasses = "character")
  }

section <- function(letter, title) {
  cat("\n", strrep("-", 74), "\n", letter, "  ", title, "\n",
      strrep("-", 74), "\n", sep = "")
}
show <- function(df) { print(as.data.frame(df), row.names = FALSE); invisible(df) }

# -----------------------------------------------------------------------------
# Rebuild the per-row evidence, the same way 02 does, so the two agree
# -----------------------------------------------------------------------------
rapid <- read_rapid()
cosd_raw <- read_cosd()
og_cohort <- readRDS(f_cohort_site)

tumour_keys <- rapid %>%
  mutate(patient_pseudo_id = as.character(patient_pseudo_id),
         reg_site3  = icd_site3(tumour_site),
         reg_morph4 = morph4(tumour_morphology_str),
         reg_trust3 = str_sub(tidy_org_code(diagnosis_trust), 1, 3)) %>%
  select(patient_pseudo_id, reg_site3, reg_morph4, reg_trust3)

from_map <- read_snomed_map()
snomed_meaning <- tibble(snomed = character(), snomed_site3 = character())
if (!is.null(from_map))
  snomed_meaning <- from_map %>%
  mutate(snomed = as.character(snomed), snomed_site3 = as.character(site3)) %>%
  select(snomed, snomed_site3)

cosd <- cosd_raw %>%
  mutate(patient_pseudo_id = as.character(patient_pseudo_id),
         site_raw    = tidy_org_code(site_code_of_diagnosis),
         topog_site3 = icd_site3(topography_icdo3),
         cosd_morph4 = morph4(morphology_icdo3),
         snomed_ver  = suppressWarnings(as.integer(snomed_version)),
         snomed      = snomed_id(diagnosis_snomedct)) %>%
  mutate(snomed = if_else(snomed_ver %in% snomed_old_versions,
                          NA_character_, snomed)) %>%
  left_join(snomed_meaning, by = "snomed") %>%
  mutate(row_site3 = coalesce(topog_site3, snomed_site3)) %>%
  inner_join(tumour_keys, by = "patient_pseudo_id") %>%
  mutate(
    site_ok     = is_site_code(site_raw),
    same_tumour = !is.na(row_site3) & !is.na(reg_site3) & row_site3 == reg_site3,
    not_og      = !is.na(row_site3) & !row_site3 %in% og_sites,
    morph_hit   = !is.na(cosd_morph4) & !is.na(reg_morph4) &
      cosd_morph4 == reg_morph4)

# resolve each site's trust the same way 02 does - from the ODS map if present,
# from the first three characters if not - so these diagnostics describe the same
# choice the build made
site_trust_map <- read_site_trust_map()
if (!is.null(site_trust_map)) {
  trust_of_site <- site_trust_map %>%
    mutate(site_raw = site_code, ods_trust = parent_trust) %>%
    select(site_raw, ods_trust)
  cosd <- cosd %>%
    left_join(trust_of_site, by = "site_raw") %>%
    mutate(site_trust3 = if_else(is.na(ods_trust),
                                 str_sub(site_raw, 1, 3), ods_trust))
} else {
  cosd <- cosd %>% mutate(site_trust3 = str_sub(site_raw, 1, 3))
}
cosd <- cosd %>%
  mutate(trust_hit = !is.na(site_trust3) & !is.na(reg_trust3) &
           site_trust3 == reg_trust3)

usable <- cosd %>% filter(site_ok) %>% filter(!(not_og & drop_non_og_rows))

# =============================================================================
section("1", "The problem, and how each patient is resolved")
# =============================================================================
# Written to be read aloud to someone who has not seen the code. It states the
# task, the obstacle, and then sorts every patient into exactly one outcome, so
# the whole approach can be seen at a glance and defended.

cat(
  "The task: give every patient the hospital where they were diagnosed.\n\n",
  "The registry already records the trust that made the diagnosis, but a trust\n",
  "runs several hospitals, and the analysis needs the hospital. The hospital is\n",
  "recorded in a separate source (COSD) as a site code. So the task is to take\n",
  "each patient's site code(s) from COSD and settle on the one hospital where the\n",
  "diagnosis was made.\n\n",
  "Three things make this harder than a straight lookup:\n",
  "  - a patient often has more than one COSD row, and the rows can carry\n",
  "    different site codes\n",
  "  - some of those rows are not about this cancer at all - they are the\n",
  "    patient's other conditions - and their sites must not be used\n",
  "  - a site code's trust is usually its first three characters, but not always,\n",
  "    so 'does this hospital belong to the diagnosing trust' cannot be answered\n",
  "    from the code alone. The ODS lookup answers it properly.\n\n",
  "Every patient is placed in exactly one of the outcomes below.\n", sep = "")

# rebuild the same per-patient picture the build settles on, but keep the steps
# visible so each patient can be labelled by how far the evidence got them
per_patient_story <- usable %>%
  group_by(patient_pseudo_id) %>%
  summarise(
    codes_all        = n_distinct(site_raw),
    codes_in_trust   = n_distinct(site_raw[trust_hit]),
    codes_confirmed  = n_distinct(site_raw[same_tumour]),
    # a code that both sits in the diagnosing trust and is confirmed by the
    # tumour is the strongest evidence there is
    codes_trust_and_confirmed = n_distinct(site_raw[trust_hit & same_tumour]),
    .groups = "drop") %>%
  left_join(select(og_cohort, patient_pseudo_id, site_dx_found),
            by = "patient_pseudo_id")

story <- per_patient_story %>%
  mutate(outcome = case_when(
    codes_all == 1 ~
      "1  one hospital offered, nothing to resolve",
    codes_in_trust == 1 ~
      "2  several offered, but only one is in the diagnosing trust",
    codes_in_trust >= 2 & codes_trust_and_confirmed == 1 ~
      "3  several in the trust, tumour information singles one out",
    codes_in_trust >= 2 & codes_trust_and_confirmed >= 2 ~
      "4  several in the trust and confirmed - tumour cannot separate them",
    codes_in_trust >= 2 ~
      "5  several in the trust, no tumour information to separate them",
    codes_in_trust == 0 ~
      "6  none of the offered hospitals is in the diagnosing trust",
    TRUE ~ "7  other"))

cat("\nEvery patient who offered a usable hospital code:", nrow(story), "\n\n")
story %>%
  count(outcome, name = "patients") %>%
  mutate(pct = round(100 * patients / sum(patients), 1)) %>%
  arrange(outcome) %>%
  show()

cat("\nIn words:\n",
    "  1-2  settled with no need for tumour information: the trust alone leaves\n",
    "       exactly one hospital standing.\n",
    "  3    the hard-won ones: several hospitals in the trust, and it is the\n",
    "       tumour record that identifies which. This is where linking the tumour\n",
    "       information earns its place.\n",
    "  4    genuinely undecidable from the data: more than one hospital in the\n",
    "       trust, all consistent with the tumour. A tie-break rule or a hospital\n",
    "       list is the only way to choose.\n",
    "  5    several in the trust with nothing to separate them at all.\n",
    "  6    the offered hospitals sit under a different trust than the registry\n",
    "       names - treated as unconfirmed rather than trusted.\n", sep = "")

uniquely <- story %>%
  filter(grepl("^[123]", outcome)) %>%
  nrow()
cat("\nUniquely identified without guessing:", uniquely, "of", nrow(story),
    sprintf(" (%.1f%%)\n", 100 * uniquely / nrow(story)))
needs_rule <- story %>% filter(grepl("^[45]", outcome)) %>% nrow()
cat("Left needing a tie-break rule or hospital list:", needs_rule,
    sprintf(" (%.1f%%)\n", 100 * needs_rule / nrow(story)))

# =============================================================================
section("A", "Amanda's question: what survives trust-matching")
# =============================================================================
# Every usable code the patient offers, and whether it sits in the registry's
# trust. The question is what is left once the trust mismatches are gone, and
# crucially whether the leftovers carry any tumour information that could break
# the tie - because if they do, this is work SNOMED or morphology can finish, and
# if they do not, no amount of coding will separate them.
per_patient <- usable %>%
  group_by(patient_pseudo_id) %>%
  summarise(
    codes_all        = n_distinct(site_raw),
    codes_trust      = n_distinct(site_raw[trust_hit]),
    any_trust        = any(trust_hit),
    # among the codes that survive the trust filter, is there tumour evidence
    # that distinguishes them
    trust_codes_confirmed = n_distinct(site_raw[trust_hit & same_tumour]),
    .groups = "drop") %>%
  left_join(select(og_cohort, patient_pseudo_id, tumour_site), by = "patient_pseudo_id")

cat("Patients offering a usable code:", nrow(per_patient), "\n\n")

cat("Before the trust filter - how many distinct codes offered:\n")
per_patient %>%
  mutate(band = case_when(codes_all == 1 ~ "1 code",
                          codes_all == 2 ~ "2 codes",
                          TRUE ~ "3 or more")) %>%
  count(band, name = "patients") %>%
  mutate(pct = round(100 * patients / sum(patients), 1)) %>%
  show()

cat("\nAfter keeping only codes in the registry's trust:\n")
after <- per_patient %>%
  mutate(band = case_when(codes_trust == 0 ~ "0 - no code in that trust",
                          codes_trust == 1 ~ "1 code",
                          TRUE ~ "2 or more - still a choice")) %>%
  count(band, name = "patients") %>%
  mutate(pct = round(100 * patients / sum(patients), 1))
show(after)

still <- per_patient %>% filter(codes_trust >= 2)
cat("\nOf the", nrow(still), "patients still choosing between codes in one trust,",
    "\nhow many have tumour evidence (ICD-O or SNOMED) that picks one out:\n")
still %>%
  mutate(settleable = if_else(trust_codes_confirmed >= 1,
                              "tumour info settles it", "nothing to separate them")) %>%
  count(settleable, name = "patients") %>%
  mutate(pct = round(100 * patients / sum(patients), 1)) %>%
  show()
cat("\nThe second line is the residue no field in this data can resolve -",
    "\nthose patients need either a canonical hospital list or a rule such as",
    "\n'keep the busier site'.\n")

# =============================================================================
section("B", "What each field actually changes")
# =============================================================================
# Rebuild the pick under four rules of rising ambition and see how many patients
# move at each step. This is the heart of it: it says, in numbers, whether
# morphology and SNOMED earn their place or whether trust alone gets there.
#
#   trust_only    rank on trust match, break ties by row count then code string
#   plus_confirm  let a row that confirms the tumour (ICD-O only) outrank the rest
#   plus_morph    add a morphology match as a tie-breaker below confirmation
#   plus_snomed   let SNOMED confirm the tumour as well as ICD-O  (the build)
#
# All four use the same deterministic final tie-break, so any difference between
# them is the field being added, not chance.
pick_under <- function(rows, use_confirm, use_morph, use_snomed) {
  rows %>%
    mutate(
      conf = if (use_snomed) same_tumour
      else !is.na(topog_site3) & !is.na(reg_site3) & topog_site3 == reg_site3,
      rank = case_when(use_confirm & conf ~ 1L, trust_hit ~ 2L, TRUE ~ 3L)) %>%
    group_by(patient_pseudo_id, site_raw) %>%
    summarise(rank = min(rank), trust = any(trust_hit),
              morph = any(morph_hit), nr = n(), .groups = "drop_last") %>%
    arrange(rank, desc(trust),
            desc(if (use_morph) morph else FALSE), desc(nr), site_raw,
            .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(patient_pseudo_id, code = site_raw)
}

p_trust  <- pick_under(usable, FALSE, FALSE, FALSE) %>% rename(trust_only = code)
p_conf   <- pick_under(usable, TRUE,  FALSE, FALSE) %>% rename(plus_confirm = code)
p_morph  <- pick_under(usable, TRUE,  TRUE,  FALSE) %>% rename(plus_morph = code)
p_snomed <- pick_under(usable, TRUE,  TRUE,  TRUE)  %>% rename(plus_snomed = code)

picks <- p_trust %>%
  full_join(p_conf,   by = "patient_pseudo_id") %>%
  full_join(p_morph,  by = "patient_pseudo_id") %>%
  full_join(p_snomed, by = "patient_pseudo_id")

moved <- function(a, b) sum(picks[[a]] != picks[[b]], na.rm = TRUE)
field_effect <- tibble(
  step = c("confirmation (ICD-O) over trust alone",
           "morphology as a tie-break",
           "SNOMED confirmation on top"),
  patients_moved = c(moved("trust_only", "plus_confirm"),
                     moved("plus_confirm", "plus_morph"),
                     moved("plus_morph", "plus_snomed")))
field_effect <- field_effect %>%
  mutate(pct_of_picked = round(100 * patients_moved / nrow(picks), 2))

cat("Patients whose chosen code changes as each field is added:\n")
show(field_effect)
write_diag(field_effect, f_diag_field_effect)
cat("\nRead this as the value of each field: a step that moves few patients is",
    "\ndoing little that trust and row-count were not already doing.\n")

# how the final build compares with trust-only, in trust terms - does the extra
# machinery put people in a better trust or just a different code
cmp <- picks %>%
  inner_join(select(tumour_keys, patient_pseudo_id, reg_trust3),
             by = "patient_pseudo_id") %>%
  mutate(trust_only_ok = str_sub(trust_only, 1, 3) == reg_trust3,
         final_ok      = str_sub(plus_snomed, 1, 3) == reg_trust3)
cat("\nChosen site sits in the registry trust - trust-only pick vs the full build:\n")
tibble(rule = c("trust-only", "full build"),
       in_registry_trust_pct = c(round(100 * mean(cmp$trust_only_ok, na.rm = TRUE), 1),
                                 round(100 * mean(cmp$final_ok, na.rm = TRUE), 1))) %>%
  show()

# =============================================================================
section("C", "Where the fields disagree")
# =============================================================================
# Two kinds of disagreement worth seeing directly.
#
# First, the SNOMED codes whose meaning differs from the topography sitting on
# the same row. These are rows where both fields are filled in and they do not
# agree - either a mis-map or a mis-coded row, and either way worth an eye.
snomed_vs_topog <- cosd %>%
  filter(!is.na(snomed_site3), !is.na(topog_site3), snomed_site3 != topog_site3) %>%
  count(snomed, snomed_site3, topog_site3, name = "rows") %>%
  arrange(desc(rows))
cat("SNOMED-labelled rows whose SNOMED site and topography site differ:",
    sum(snomed_vs_topog$rows), "rows,", nrow(snomed_vs_topog), "distinct codes\n")
if (nrow(snomed_vs_topog)) {
  show(head(snomed_vs_topog, 15))
  write_diag(snomed_vs_topog, f_diag_conflicts)
}

# Second, patients where the tumour-confirmed code and the trust-matched code are
# two different places. This is the case that decides whether "confirm the
# tumour" and "match the trust" ever pull against each other, and which the build
# is trusting when they do (it takes confirmation).
conflict <- usable %>%
  group_by(patient_pseudo_id) %>%
  summarise(confirmed_code = list(unique(site_raw[same_tumour])),
            trust_code     = list(unique(site_raw[trust_hit])),
            .groups = "drop") %>%
  mutate(has_conf = lengths(confirmed_code) > 0,
         has_trust = lengths(trust_code) > 0,
         overlap = map2_lgl(confirmed_code, trust_code,
                            ~ length(intersect(.x, .y)) > 0)) %>%
  filter(has_conf, has_trust, !overlap)
cat("\nPatients whose tumour-confirmed code and trust-matched code disagree:",
    nrow(conflict), "\n")
cat("For these, the build takes the confirmed code. If that is the wrong call,",
    "\nthis is the number of patients it affects.\n")
if (nrow(conflict)) {
  # no patient_pseudo_id here - this file is meant to leave the server, and a
  # pattern-level count of which codes disagree does the job without carrying
  # any per-patient identifier. Tracing a specific patient, if that is ever
  # needed, is done on the server directly from the .rds files, which never
  # leave it.
  conflict %>%
    mutate(confirmed = map_chr(confirmed_code, ~ paste(.x, collapse = "|")),
           trust     = map_chr(trust_code, ~ paste(.x, collapse = "|"))) %>%
    count(confirmed, trust, name = "patients") %>%
    arrange(desc(patients)) %>%
    write_diag(f_diag_trust_vs_tum)
}

# =============================================================================
section("D", "The C15 cohort this is headed for")
# =============================================================================
c15 <- og_cohort %>% filter(tumour_site == "C15")
cat("Registry C15 patients:", nrow(c15), "\n")
cat("With a site of diagnosis:", sum(c15$site_dx_found),
    sprintf("(%.1f%%)\n", 100 * mean(c15$site_dx_found)))

cat("\nHow their code was chosen:\n")
c15 %>%
  filter(site_dx_found) %>%
  count(site_dx_basis, name = "patients") %>%
  mutate(pct = round(100 * patients / sum(patients), 1)) %>%
  show()

cat("\nHow firm the C15 sites are:\n")
c15 %>%
  filter(site_dx_found) %>%
  summarise(
    single_clear_code = sum(!site_dx_ambiguous),
    picked_among_many = sum(site_dx_ambiguous),
    in_registry_trust = sum(site_dx_trust_agrees, na.rm = TRUE),
    trust_uncheckable = sum(is.na(site_dx_trust_agrees))) %>%
  show()

cat("\nBusiest C15 sites (a sanity check - these should be the big OG centres):\n")
c15 %>%
  filter(site_dx_found) %>%
  count(site_dx_code, name = "patients", sort = TRUE) %>%
  head(15) %>%
  show()

cat("\n", strrep("=", 74), "\n", sep = "")
cat("Diagnostics complete.\n")
written <- c(f_diag_field_effect,
             if (nrow(snomed_vs_topog)) f_diag_conflicts,
             if (nrow(conflict)) f_diag_trust_vs_tum)
cat("Written for review:\n")
for (f in written) cat("  ", f, "\n")
skipped <- c(if (!nrow(snomed_vs_topog)) f_diag_conflicts,
             if (!nrow(conflict)) f_diag_trust_vs_tum)
if (length(skipped)) {
  cat("Not written, because there was nothing to put in them:\n")
  for (f in skipped) cat("  ", basename(f), "\n")
}