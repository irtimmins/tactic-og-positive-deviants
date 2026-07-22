# =============================================================================
# 02  Add the endoscopy site
# -----------------------------------------------------------------------------
# The rapid record already carries the diagnostic endoscopy date (endodateHES)
# and a trust for it (endotrustHES), but endotrustHES is a mix of three-character
# provider codes and five-character site codes, so it is not a clean hospital
# site. This step goes back to the HES-APC extract, re-finds the endoscopy
# episode by its operation date, and reads the five-character sitetret off that
# same episode, which gives a consistent hospital site.
#
# It works in four passes:
#   1  unfold the 24 operation slots and keep the diagnostic endoscopy codes
#   2  match each endoscopy's operation date to the date on the rapid record
#   3  choose one episode per patient and read its site
#   4  put it on the registry extract as endoscopy_site
#
# Two flags are carried alongside so a reader can see, for the patients we have
# HES-APC on, what share have no APC endoscopy record at all: endo_in_hes (the
# patient appears in the extract) and endo_has_apc (they have an APC endoscopy
# code, on any date). These are what answers "do we need HES-OP as well".
#
# Reads : the rapid tumour dta, the HES-APC extract (rds or dta)
# Writes: og_endoscopy_hes_lookup.rds  (one row per patient with a site)
#         og_cohort_endoscopy_hes.rds  (the registry plus the site and flags)
# =============================================================================

source("R/derive_hospital_code_from_hes/01_define_parameters.R")

# read_rapid() and read_hes() are here so the checking script can hand this
# script made-up data instead of the real files: define either one before
# sourcing and it will be used. Otherwise read_hes prefers the cohort-filtered
# rds and falls back to reading the raw APC dta.
if (!exists("read_rapid"))
  read_rapid <- function() read_dta(path_rapid_dta)

if (!exists("read_hes"))
  read_hes <- function() {
    if (file.exists(f_hes_extract)) return(readRDS(f_hes_extract))
    cat("No filtered HES extract at", f_hes_extract,
        "- reading the raw APC dta, this is large.\n")
    read_dta(path_hes_apc_dta)
  }

# -----------------------------------------------------------------------------
# Read the two files
# -----------------------------------------------------------------------------
rapid <- read_rapid()
check_input(rapid,
            c("patient_pseudo_id", "endoHES", "endodateHES", "endotrustHES"),
            "the rapid tumour extract", path_rapid_dta,
            min_rows = getOption("og_min_input_rows", 1000L))

op_cols     <- sprintf("opertn_%02d", 1:24)
opdate_cols <- sprintf("opdate_%02d", 1:24)

hes <- read_hes()
check_input(hes,
            c("patient_pseudo_id", "epistart", "epiorder", "sitetret",
              "procode3", op_cols, opdate_cols),
            "the HES-APC extract", f_hes_extract,
            min_rows = getOption("og_min_input_rows", 1000L))

cat("Read", nrow(rapid), "registry rows and", nrow(hes), "HES episodes\n")

# Everything below joins on the patient alone, which is only safe because the
# registry holds one row per patient. Stop rather than let a stray duplicate
# multiply rows quietly.
if (anyDuplicated(rapid$patient_pseudo_id))
  stop("the rapid extract has more than one row for some patients, and this ",
       "script assumes one.", call. = FALSE)

# -----------------------------------------------------------------------------
# The endoscopy date the registry already holds
# -----------------------------------------------------------------------------
# One reference date per patient, kept only where the rapid record says there was
# a HES endoscopy at all. endotypeHES and endotrustHES come along so the result
# can be set against what we had before.
anchor <- rapid %>%
  mutate(patient_pseudo_id = as.character(patient_pseudo_id),
         endo_ref_date = hes_date(endodateHES)) %>%
  filter(endoHES == 1, !is.na(endo_ref_date)) %>%
  select(patient_pseudo_id, endo_ref_date, endotypeHES, endotrustHES)

cat("Patients with a HES endoscopy anchor (endoHES==1):", nrow(anchor), "\n")

# only these patients' HES episodes matter; give each episode a key
hes_cohort <- hes %>%
  mutate(patient_pseudo_id = as.character(patient_pseudo_id)) %>%
  filter(patient_pseudo_id %in% anchor$patient_pseudo_id) %>%
  mutate(epi_id = row_number(),
         epistart = hes_date(epistart))

# -----------------------------------------------------------------------------
# Pass 1: unfold the operation slots and keep the endoscopies
# -----------------------------------------------------------------------------
# Each episode has 24 operation slots and 24 matching date slots. Unfold them
# together - one row per coded operation, its date beside it - and keep the rows
# whose code is a diagnostic endoscopy. The site and provider ride along on the
# episode they belong to.
endo_ops <- hes_cohort %>%
  select(epi_id, patient_pseudo_id, epistart, epiorder, sitetret, procode3,
         all_of(op_cols), all_of(opdate_cols)) %>%
  pivot_longer(cols = c(all_of(op_cols), all_of(opdate_cols)),
               names_to  = c(".value", "slot"),
               names_pattern = "(opertn|opdate)_(\\d+)") %>%
  mutate(opcs   = norm_opcs(opertn),
         opdate = hes_date(opdate)) %>%
  filter(opcs %in% opcs_diagnostic_endoscopy)

# patients with any APC endoscopy code, whatever its date - this is the flag that
# separates "no APC endoscopy record" from "APC endoscopy on the wrong date"
apc_endo_ids <- unique(endo_ops$patient_pseudo_id)

cat("Endoscopy-coded operations found:", nrow(endo_ops), "across",
    length(apc_endo_ids), "patients\n")

# -----------------------------------------------------------------------------
# Pass 2: match each endoscopy to the reference date
# -----------------------------------------------------------------------------
# The endoscopy date is the opdate on the endoscopy's own slot, so an endoscopy
# with no date cannot be placed and drops out here. What remains is matched to
# the rapid reference date and kept if it falls inside the window.
matched <- endo_ops %>%
  filter(!is.na(opdate)) %>%
  inner_join(anchor, by = "patient_pseudo_id") %>%
  mutate(days_off = as.integer(opdate - endo_ref_date)) %>%
  filter(abs(days_off) <= endo_window_days)

# -----------------------------------------------------------------------------
# Pass 3: choose one episode per patient and read its site
# -----------------------------------------------------------------------------
# Order best first (a real site before none, then closest date, then earliest,
# then lowest episode order) and take the top row. Where the chosen episode's
# site fails the site-code rule the site is left missing rather than passed on.
best <- matched %>%
  mutate(site_ok = is_hes_site(sitetret)) %>%
  arrange(patient_pseudo_id, desc(prefer_sited_episode & site_ok),
          abs(days_off), opdate, epiorder) %>%
  distinct(patient_pseudo_id, .keep_all = TRUE) %>%
  mutate(endoscopy_site     = if_else(site_ok, sitetret, NA_character_),
         endoscopy_provider = procode3,
         endoscopy_opcs     = opcs,
         endoscopy_opdate   = opdate,
         endoscopy_days_off = days_off) %>%
  select(patient_pseudo_id, endoscopy_site, endoscopy_provider,
         endoscopy_opcs, endoscopy_opdate, endoscopy_days_off)

# one row per anchored patient, with the three flags and whether a site was read.
# endo_matched is kept apart from endoscopy_site_found on purpose: an endoscopy
# can be found at the right date and still yield no site, because the episode's
# sitetret was blank or a default. Without this the two would be lumped together
# and a missing site would be reported as a missing endoscopy, which it is not.
anchor_out <- anchor %>%
  mutate(endo_in_hes  = patient_pseudo_id %in% hes_cohort$patient_pseudo_id,
         endo_has_apc = patient_pseudo_id %in% apc_endo_ids,
         endo_matched = patient_pseudo_id %in% matched$patient_pseudo_id) %>%
  left_join(best, by = "patient_pseudo_id") %>%
  mutate(endoscopy_site_found = !is.na(endoscopy_site))

# the lookup the main build joins in: patients who got a site, and the evidence
endoscopy_lookup <- anchor_out %>%
  filter(endoscopy_site_found) %>%
  select(patient_pseudo_id, endoscopy_site, endoscopy_provider,
         endoscopy_opcs, endoscopy_opdate, endoscopy_days_off)

saveRDS(endoscopy_lookup, f_endoscopy_lookup)

# -----------------------------------------------------------------------------
# Pass 4: put it on the registry extract
# -----------------------------------------------------------------------------
rows_in <- nrow(rapid)

carry <- anchor_out %>%
  select(patient_pseudo_id, endoscopy_site, endoscopy_provider, endoscopy_opcs,
         endoscopy_opdate, endoscopy_days_off, endoscopy_site_found,
         endo_in_hes, endo_has_apc, endo_matched)

og_cohort <- rapid %>%
  mutate(patient_pseudo_id = as.character(patient_pseudo_id)) %>%
  left_join(carry, by = "patient_pseudo_id") %>%
  mutate(endoscopy_site_found = coalesce(endoscopy_site_found, FALSE),
         endo_in_hes          = coalesce(endo_in_hes, FALSE),
         endo_has_apc         = coalesce(endo_has_apc, FALSE),
         endo_matched         = coalesce(endo_matched, FALSE),
         # does the site's provider agree with the trust the rapid record held?
         # left missing where either side is missing, so it is never confused
         # with a real disagreement
         endo_provider_agrees = if_else(
           is.na(endoscopy_provider) | is.na(endotrustHES) | endotrustHES == "",
           NA,
           str_sub(endoscopy_provider, 1, 3) == str_sub(endotrustHES, 1, 3)))

stopifnot(nrow(og_cohort) == rows_in)

saveRDS(og_cohort, f_cohort_endoscopy)

# -----------------------------------------------------------------------------
# What came out
# -----------------------------------------------------------------------------
# The denominator that matters is the patients we actually have HES-APC on: only
# for those can "no APC endoscopy record" be said at all.
assessable <- anchor_out %>% filter(endo_in_hes)
N <- nrow(assessable)

site_found <- sum(assessable$endoscopy_site_found)
no_site    <- sum(assessable$endo_matched & !assessable$endoscopy_site_found)
apc_offwin <- sum(!assessable$endo_matched & assessable$endo_has_apc)
no_apc     <- sum(!assessable$endo_has_apc)

cat("\nAnchored patients:", nrow(anchor_out),
    "| in this HES-APC extract:", N,
    "| not in extract:", sum(!anchor_out$endo_in_hes), "\n")

cat("\nOf", N, "patients with HES-APC data:\n")
cat(sprintf("  APC endoscopy at the reference date, site read : %6d (%.1f%%)\n",
            site_found, 100 * site_found / N))
cat(sprintf("  endoscopy found, but its site was unusable     : %6d (%.1f%%)\n",
            no_site, 100 * no_site / N))
cat(sprintf("  APC endoscopy present but not within %d days     : %6d (%.1f%%)\n",
            endo_window_days, apc_offwin, 100 * apc_offwin / N))
cat(sprintf("  no HES-APC endoscopy record at all             : %6d (%.1f%%)\n",
            no_apc, 100 * no_apc / N))
cat(sprintf("\n  %.1f%% of the %d patients with HES-APC data have no HES-APC endoscopy record\n",
            100 * no_apc / N, N))

cat("\nDays between the re-found opdate and the rapid endoscopy date:\n")
print(table(anchor_out$endoscopy_days_off, useNA = "ifany"))

# what the HES site gains over endotrustHES, among the sites found
gained <- og_cohort %>%
  filter(endoscopy_site_found) %>%
  mutate(trust_len = str_length(str_trim(endotrustHES)))
cat("\nendotrustHES among patients who now have a HES site:\n")
print(table(gained$trust_len, useNA = "ifany"))
cat("  identical to endotrustHES        :",
    sum(gained$endoscopy_site == gained$endotrustHES, na.rm = TRUE), "\n")
cat("  same provider (first 3 chars)    :",
    sum(gained$endo_provider_agrees, na.rm = TRUE), "\n")
cat("  endotrustHES was <5 char, now 5  :",
    sum(gained$trust_len < 5, na.rm = TRUE), "\n")
cat("  provider disagrees, worth a look :",
    sum(!gained$endo_provider_agrees, na.rm = TRUE), "\n")

cat("\nDistinct endoscopy sites in the cohort:",
    n_distinct(og_cohort$endoscopy_site, na.rm = TRUE), "\n")
cat("\nSaved", f_endoscopy_lookup, "and", f_cohort_endoscopy, "\n")
cat("02 complete. Next: 03_endoscopy_diagnostics.R\n")

# ---- things worth a look, when there is a reason to (uncomment) -------------
# busiest endoscopy sites - should be the big OG centres
# og_cohort %>% filter(endoscopy_site_found) %>%
#   count(endoscopy_site, sort = TRUE) %>% print(n = 20)
# the off-window group, by how far off they are
# matched_all <- endo_ops %>% filter(!is.na(opdate)) %>%
#   inner_join(anchor, by = "patient_pseudo_id") %>%
#   mutate(days_off = as.integer(opdate - endo_ref_date))
# matched_all %>% filter(abs(days_off) > endo_window_days) %>%
#   count(sign(days_off)) %>% print()
