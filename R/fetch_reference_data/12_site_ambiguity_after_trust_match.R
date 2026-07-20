# =============================================================================
# 12  How much of the site ambiguity survives trust-matching
# -----------------------------------------------------------------------------
# Amanda's question, from her email of 17 July: once COSD rows whose trust does
# not match the registry's diagnosis_trust are thrown out, are the patients left
# with more than one COSD row mostly disagreeing on trust anyway - which no
# amount of SNOMED work can fix - or are they the same trust submitting more
# than one 5-character site code, which is the case SNOMED (or morphology, or
# just picking one) is actually needed for.
#
# This does not change the build. It answers the question on its own, using the
# same registry and COSD files, so it can be run and read independently.
#
# Run from the project root:
#   Rscript R/fetch_reference_data/12_site_ambiguity_after_trust_match.R
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(haven))

source("R/config/directories.R")   # dir_raw, path_rapid_dta, path_cosd_dta
source("R/shared/utils.R")          # tidy_org_code

rapid <- read_dta(path_rapid_dta, col_select = c(
  patient_pseudo_id, diagnosis_trust))
cosd <- read_dta(path_cosd_dta, col_select = c(
  patient_pseudo_id, site_code_of_diagnosis))

reg <- rapid %>%
  transmute(patient_pseudo_id = as.character(patient_pseudo_id),
            reg_trust3 = str_sub(tidy_org_code(diagnosis_trust), 1, 3))

rows <- cosd %>%
  transmute(patient_pseudo_id = as.character(patient_pseudo_id),
            site = tidy_org_code(site_code_of_diagnosis)) %>%
  filter(!is.na(site), str_length(site) == 5) %>%
  inner_join(reg, by = "patient_pseudo_id") %>%
  mutate(site_trust3 = str_sub(site, 1, 3))

cat("COSD rows with a five-character site code:", nrow(rows), "\n")

# Amanda's filter: keep only rows whose trust matches the registry's
matched <- rows %>% filter(site_trust3 == reg_trust3)
cat("Rows left after dropping trust mismatches:", nrow(matched),
    sprintf("(%.1f%%)\n", 100 * nrow(matched) / nrow(rows)))

per_patient <- matched %>%
  distinct(patient_pseudo_id, site) %>%
  group_by(patient_pseudo_id) %>%
  summarise(codes = n(), .groups = "drop")

cat("\nOf the patients with a trust-matched COSD row:\n")
per_patient %>%
  mutate(has_choice = codes > 1) %>%
  count(has_choice, name = "patients") %>%
  mutate(pct = round(100 * patients / sum(patients), 1)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

still_ambiguous <- per_patient %>% filter(codes > 1)
cat("\n", nrow(still_ambiguous),
    "patients still offer more than one site code from the same trust -",
    "different hospital, same trust, or a coding slip within one hospital.\n")

if (nrow(still_ambiguous)) {
  cat("\nHow many different codes these patients are choosing between:\n")
  still_ambiguous %>%
    count(codes, name = "patients") %>%
    as.data.frame() %>%
    print(row.names = FALSE)
  
  cat("\nA sample, to look at directly:\n")
  still_ambiguous %>%
    slice_head(n = 5) %>%
    inner_join(matched, by = "patient_pseudo_id") %>%
    distinct(patient_pseudo_id, site) %>%
    arrange(patient_pseudo_id, site) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
}

cat("\nFor comparison, before Amanda's trust filter:\n")
rows %>%
  distinct(patient_pseudo_id, site) %>%
  count(patient_pseudo_id, name = "codes") %>%
  mutate(has_choice = codes > 1) %>%
  count(has_choice, name = "patients") %>%
  mutate(pct = round(100 * patients / sum(patients), 1)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\nThe gap between the two tables above is what trust-matching alone settles.",
    "\nWhat's left in the first table is the true test of whether SNOMED, or",
    "\nmorphology, or picking the more common code, is pulling its weight.\n")