# =============================================================================
# 03b  Make stand-in COSD and registry data
# -----------------------------------------------------------------------------
# Writes a pair of dta files that stand in for the two real extracts, so the
# build can be run and checked away from the analysis server. The shapes and the
# messy bits are taken from the profiles of the January 2026 data:
#
#   rapid tumour   95,333 rows, one per patient. tumour_site C15 67.9% / C16
#                  32.1%. tumour_morphology_str blank 13.5%. diagnosis_trust
#                  blank 0.5%, 121 distinct three-character codes.
#   COSD diagnosis 130,170 rows carrying 80,554 patients, so 1.62 rows each.
#                  site_code_of_diagnosis blank 8.4%, 608 distinct values mixing
#                  site codes with trust codes, GP practice codes and defaults.
#                  topography_icdo3 blank 91.8%, morphology_icdo3 blank 83.0%,
#                  diagnosis_snomedct blank 61.1%. snomed_version blank 27.2%,
#                  and of the rest 82% say SNOMED CT, 17% say not known, and
#                  about 1% say one of the older SNOMEDs.
#
# Because this script knows which site each patient was really diagnosed at, and
# which cancer each COSD row is really about, it also writes those answers out.
# 91 uses them to ask how often the build gets it right, which no check against
# the real data can do. Nothing here is meant to be analysed - only the code is
# being tested.
#
# Writes: <dir_sim>/*.dta and site_truth.rds
# =============================================================================

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(magrittr)   # for %>%
})

if (!exists("dir_sim")) dir_sim <- "Output/sim/raw"
if (!exists("sim_scale")) sim_scale <- 1        # 0.1 gives a tenth-size run
if (!exists("sim_seed")) sim_seed <- 20260212
dir.create(dir_sim, recursive = TRUE, showWarnings = FALSE)
set.seed(sim_seed)

n_patients <- round(95333 * sim_scale)

# -----------------------------------------------------------------------------
# 121 trusts, each with a handful of hospital sites
# -----------------------------------------------------------------------------
trusts <- unique(paste0("R", sample(LETTERS, 200, TRUE), sample(0:9, 200, TRUE)))
trusts <- sort(trusts[1:121])

sites <- map_dfr(trusts, function(tr) {
  tibble(trust = tr,
         site = paste0(tr, sprintf("%02d", sample(1:20, sample(2:8, 1)))))
})
sites_by_trust <- split(sites$site, sites$trust)

# -----------------------------------------------------------------------------
# How cancers are written down
# -----------------------------------------------------------------------------
# Every cancer site has its own SNOMED codes. Two codes are deliberately sloppy:
# they turn up against more than one cancer, so the build should decide they say
# nothing rather than trust them.
cancer_sites <- c("C15", "C16", "C18", "C20", "C25", "C34", "C50", "C61", "C64")

snomed_by_site <- map(cancer_sites, function(s) {
  base <- 100000000 + 1000 * match(s, cancer_sites)
  sprintf("%.0f", base + c(1, 2, 3))
})
names(snomed_by_site) <- cancer_sites
snomed_sloppy <- c("999000001", "999000002")

# numbers from the older SNOMEDs. They are reused across cancers, so a build
# that reads them without checking the version will place patients wrongly.
snomed_old <- c("81403", "80703", "81409")

morph_pool <- c("8140", "8070", "8010", "8480", "8000", "8144", "8020",
                "8041", "8071", "8046", "8490", "8260")

# -----------------------------------------------------------------------------
# The rapid tumour extract
# -----------------------------------------------------------------------------
truth <- tibble(
  patient_pseudo_id = sprintf("P%09d", seq_len(n_patients)),
  true_trust = sample(trusts, n_patients, TRUE),
  tumour_site = sample(c("C15", "C16"), n_patients, TRUE, c(0.679, 0.321)),
  tumour_morphology_str = sample(morph_pool, n_patients, TRUE)) %>%
  mutate(true_site = map_chr(true_trust, ~ sample(sites_by_trust[[.x]], 1)))

rapid <- truth %>%
  mutate(
    ind_pseudo_id    = seq_len(n()),
    tumour_pseudo_id = seq_len(n()),
    tumourgroup      = "og",
    # 13.5% of morphologies are blank in the real file
    tumour_morphology_str = if_else(runif(n()) < 0.135, "", tumour_morphology_str),
    # 0.5% of diagnosing trusts are blank
    diagnosis_trust  = if_else(runif(n()) < 0.005, "", true_trust),
    diagnosisdate    = as.Date("2018-01-28") + sample(0:2803, n(), TRUE),
    tdiagnosisdate   = format(diagnosisdate, "%Y-%m-%d"),
    gender           = sample(1:2, n(), TRUE, c(0.75, 0.25)),
    age              = as.integer(pmin(95, pmax(20, round(rnorm(n(), 73, 11))))),
    sotn_cohort      = as.numeric(runif(n()) < 0.201),
    audit_eligible   = as.numeric(runif(n()) < 0.800)) %>%
  select(patient_pseudo_id, ind_pseudo_id, tumour_pseudo_id, tumour_site,
         tumourgroup, tumour_morphology_str, diagnosis_trust, diagnosisdate,
         tdiagnosisdate, gender, age, sotn_cohort, audit_eligible)

# -----------------------------------------------------------------------------
# The COSD diagnosis extract
# -----------------------------------------------------------------------------
# 84.5% of registry patients appear, at 1.62 rows each. The first row is the
# patient's own OG tumour; the rest are either that same tumour submitted again
# from a different site in the same trust, or a different cancer altogether -
# which for a C15 patient includes C16.
in_cosd <- truth %>% slice_sample(prop = 0.845)

in_cosd <- in_cosd %>%
  mutate(n_rows = 1L + rbinom(n(), 3, 0.22))   # mean about 1.62 rows each

# expand each row n_rows times, and number the copies 1, 2, 3... within each
# patient - what tidyr::uncount(n_rows, .id = "row_in_patient") would do. Done
# in base R instead, so this script does not need tidyr for the sake of one
# call; sequence() gives exactly that 1..n, 1..n, ... numbering.
idx <- rep(seq_len(nrow(in_cosd)), in_cosd$n_rows)
cosd <- in_cosd[idx, , drop = FALSE]
cosd$row_in_patient <- sequence(in_cosd$n_rows)
rownames(cosd) <- NULL

cosd <- cosd %>%
  mutate(
    row_kind = case_when(
      row_in_patient == 1L ~ "the patient's own tumour",
      runif(n()) < 0.45    ~ "a different cancer",
      TRUE                 ~ "own tumour, another site"),
    # which cancer this row is really about
    row_site3 = if_else(row_kind == "a different cancer",
                        sample(cancer_sites, n(), TRUE),
                        tumour_site),
    # a "different cancer" row that lands back on the patient's own site is not
    # different after all
    row_kind = if_else(row_kind == "a different cancer" & row_site3 == tumour_site,
                       "own tumour, another site", row_kind),
    # where the row says the diagnosis was made
    site_true = case_when(
      row_kind == "the patient's own tumour" ~ true_site,
      row_kind == "own tumour, another site" ~ map_chr(true_trust,
                                                       ~ sample(sites_by_trust[[.x]], 1)),
      TRUE                                   ~ map_chr(sample(trusts, n(), TRUE),
                                                       ~ sample(sites_by_trust[[.x]], 1))))

cosd <- cosd %>%
  mutate(
    # topography: on 8.2% of rows, and a few of those carry the description too
    topography_icdo3 = if_else(runif(n()) < 0.082,
                               paste0(row_site3, sample(0:9, n(), TRUE)), ""),
    topography_icdo3 = if_else(topography_icdo3 != "" & runif(n()) < 0.05,
                               paste0(topography_icdo3, ": LOWER THIRD OF OESOPHAGUS"),
                               topography_icdo3),
    # SNOMED: on 38.9% of rows. Two codes in fifty are one of the sloppy pair,
    # which turn up against any cancer.
    snomed_draw = runif(n()),
    diagnosis_snomedct = case_when(
      snomed_draw >= 0.389 ~ NA_character_,
      snomed_draw < 0.008  ~ sample(snomed_sloppy, n(), TRUE),
      TRUE                 ~ map_chr(row_site3, ~ sample(snomed_by_site[[.x]], 1))),
    # which SNOMED the number was written in. The old versions carry numbers
    # that look like the others but mean something else entirely, so the build
    # must not read them.
    ver_draw = runif(n()),
    snomed_version = case_when(
      ver_draw < 0.272 ~ NA_integer_,
      ver_draw < 0.279 ~ 1L,
      ver_draw < 0.285 ~ 2L,
      ver_draw < 0.288 ~ 4L,
      ver_draw < 0.412 ~ 99L,
      TRUE             ~ 5L),
    diagnosis_snomedct = if_else(snomed_version %in% c(1L, 2L, 4L),
                                 sample(snomed_old, n(), TRUE),
                                 diagnosis_snomedct),
    # morphology: on 17% of rows, and agrees with the registry when the row is
    # the patient's own tumour
    morphology_icdo3 = case_when(
      runif(n()) >= 0.17               ~ "",
      row_kind == "a different cancer" ~ sample(morph_pool, n(), TRUE),
      TRUE                             ~ tumour_morphology_str),
    behaviour_icdo3 = if_else(morphology_icdo3 == "", "", "3"),
    # and now spoil the site code the way the real field is spoiled
    spoil = runif(n()),
    site_code_of_diagnosis = case_when(
      spoil < 0.084 ~ "",                                   # blank
      spoil < 0.104 ~ str_sub(site_true, 1, 3),             # trust, not site
      spoil < 0.116 ~ paste0("B86", sprintf("%03d", sample(1:99, n(), TRUE))),
      spoil < 0.124 ~ sample(c("89997", "89999", "X99999", "XXXXXX", "ZZ201",
                               "5E813"), n(), TRUE),
      spoil < 0.130 ~ str_to_lower(site_true),              # case is not certain
      TRUE          ~ site_true),
    diagnosis_snomedct = as.numeric(diagnosis_snomedct))

cosd_out <- cosd %>%
  select(patient_pseudo_id, morphology_icdo3, behaviour_icdo3,
         topography_icdo3, site_code_of_diagnosis, diagnosis_snomedct,
         snomed_version)

# a few COSD patients the registry does not have: the build must drop them
cosd_out <- bind_rows(cosd_out, tibble(
  patient_pseudo_id = sprintf("Q%09d", 1:200),
  morphology_icdo3 = "", behaviour_icdo3 = "", topography_icdo3 = "",
  site_code_of_diagnosis = "RZZ01", diagnosis_snomedct = NA_real_,
  snomed_version = NA_integer_)) %>%
  slice_sample(prop = 1)                    # the row order must not matter

# -----------------------------------------------------------------------------
# Write
# -----------------------------------------------------------------------------
f_rapid <- file.path(dir_sim,
                     "20260212_Rapidtumour_linked_2026SOTN_clean_OG_postPT.dta")
f_cosd <- file.path(dir_sim,
                    "20260212_all_cosddiagnosis_rapid_202601_OG.dta")

write_dta(rapid, f_rapid)
write_dta(cosd_out, f_cosd)
saveRDS(truth %>% select(patient_pseudo_id, tumour_site, true_trust, true_site),
        file.path(dir_sim, "site_truth.rds"))

# the map the fetch script would have produced from TRUD for these numbers
sim_map <- tibble(snomed = unlist(snomed_by_site[c("C15", "C16")]),
                  site3 = rep(c("C15", "C16"), each = 3)) %>%
  mutate(icd10_targets = paste0(site3, ".9"))
write.csv(sim_map, file.path(dir_sim, "snomed_og_lookup.csv"), row.names = FALSE)

# the ODS site-to-trust map R/reference/20_fetch_site_trust_map.R would produce.
# For most sites the trust is the first three characters; a small number are
# deliberately placed under a different trust, as happens in real ODS data after
# a merger, so the build's handling of that case is exercised. A few codes are
# marked as not hospital sites, and a few as unknown to ODS.
set.seed(sim_seed + 7)
all_sites <- sort(unique(sites$site))
moved <- sample(all_sites, max(1, round(0.02 * length(all_sites))))
not_hosp <- sample(setdiff(all_sites, moved),
                   max(1, round(0.01 * length(all_sites))))

sim_site_trust <- tibble(site_code = all_sites) %>%
  mutate(
    prefix       = str_sub(site_code, 1, 3),
    # a moved site is handed to some other trust from the list
    parent_trust = if_else(site_code %in% moved,
                           sample(trusts, n(), TRUE), prefix),
    parent_trust = if_else(parent_trust == prefix & site_code %in% moved,
                           trusts[(match(prefix, trusts) %% length(trusts)) + 1],
                           parent_trust),
    trust_is_prefix  = parent_trust == prefix,
    parent_from_ods  = !(site_code %in% c(not_hosp)),
    status           = "Active",
    is_hospital_site = !(site_code %in% not_hosp),
    record_class     = if_else(is_hospital_site, "RC2", "RC2"),
    primary_role     = if_else(is_hospital_site, "RO198", "RO177"),
    predecessor      = NA_character_,
    found            = TRUE) %>%
  select(site_code, name = prefix, parent_trust, trust_is_prefix,
         parent_from_ods, status, is_hospital_site, record_class,
         primary_role, predecessor, found)
write.csv(sim_site_trust, file.path(dir_sim, "site_trust_map.csv"),
          row.names = FALSE)

cat("Made", nrow(rapid), "registry rows and", nrow(cosd_out), "COSD rows (",
    n_distinct(cosd_out$patient_pseudo_id), "patients ) in", dir_sim, "\n")
cat("  blank: site code",
    round(100 * mean(cosd_out$site_code_of_diagnosis == ""), 1), "% |",
    "topography", round(100 * mean(cosd_out$topography_icdo3 == ""), 1), "% |",
    "SNOMED", round(100 * mean(is.na(cosd_out$diagnosis_snomedct)), 1), "%\n")
cat("  distinct site code values:",
    n_distinct(cosd_out$site_code_of_diagnosis), "\n")