# 02  build the positive-deviance analysis cohort
# -----------------------------------------------------------------------------
# Draws the analysis sample from the CWT-merged cohort through two funnels and
# writes both so the loss is fully accounted for:
#
#   patient funnel   curative pathway -> stage 1-3 -> valid, positive, plausible
#                    endoscopy-to-decision-to-treat time
#   hospital funnel  the site of diagnosis is the unit; a site is kept only if it
#                    has at least min_per_year diagnoses in every year of the
#                    window. Excluded sites are listed with their yearly counts.
#
# The outcome is the time from diagnostic endoscopy to the decision to treat
# (wt_endo_to_dtt from the merge stage), and only positive values are analysed.
#
# Reads : og_cohort_cwt.rds
# Writes: pd_cohort.rds, pd_flow.csv, pd_hospitals_excluded.csv

library(dplyr)
library(lubridate)

source("R/identify_positive_deviants/01_config.R")

og <- readRDS(in_rds)
cat("Read", nrow(og), "patients from", basename(in_rds), "\n")

# a small running flowchart: each step records how many patients remain and how
# many were dropped, so the attrition reads top to bottom.
flow  <- tibble(step = character(), n = integer(), dropped = integer(),
                n_hosp = integer())
.prev <- nrow(og)
fc <- function(step, d, n_hosp = NA_integer_) {
  n <- nrow(d)
  flow <<- bind_rows(flow, tibble(step = step, n = n, dropped = .prev - n,
                                  n_hosp = n_hosp))
  .prev <<- n
}
fc("All patients in the CWT-merged cohort", og)

# main audit cohort -----------------------------------------------------------
# sotn_cohort flags the registry's own audit-eligible cohort (State of the
# Nation exclusions already applied - stage, route, DCO and similar validity
# checks at the registry level). It also carries its own diagnosis-date window
# (the SoTN audit period), so restricting to it here is itself a date
# restriction. Deliberately no separate rolling-window filter is applied on top
# of this: doing so would be a second, different date restriction that can cut
# patients sotn_cohort already includes, for no analytic reason. If a different
# reporting period is wanted later, it should replace sotn_cohort's restriction,
# not stack on top of it - see the note below.
og <- og %>% mutate(endoscopy_date = as.Date(endoscopy_date))
df <- og %>% filter(sotn_cohort == 1)
fc("Main audit cohort (sotn_cohort == 1)", df)

# the date range actually present, used only to build the season / calendar-year
# covariates below (quarter dummies, and an earlier/later half split for the
# yr_late term) - not to exclude anyone. Taken from the diagnostic endoscopy date
# (the clock-start) over the sotn_cohort-restricted data.
start_date <- min(df$endoscopy_date, na.rm = TRUE)
end_date   <- max(df$endoscopy_date, na.rm = TRUE)
cat(sprintf("endoscopy dates present (sotn_cohort): %s to %s\n", start_date, end_date))

# future option: restricting to a specific, different reporting period (e.g.
# Jan 2023 to Mar 2025) instead of relying on sotn_cohort's own window. Not
# applied yet - left here, commented out, as a marker for when that decision is
# made. If uncommented, this should REPLACE the sotn_cohort filter above (or be
# combined deliberately), not sit alongside it as an extra cut.
#   start_date <- as.Date("2023-01-01")
#   end_date   <- as.Date("2025-03-31")
#   df <- df %>% filter(endoscopy_date >= start_date, endoscopy_date <= end_date)
#   fc(sprintf("Endoscopy %s to %s", format(start_date, "%b %Y"),
#              format(end_date, "%b %Y")), df)

# curative pathway -----------------------------------------------------------
# flag the held-out pathway first, so its count is reported before it is dropped,
# then keep only the curative-inclusion pathways.
n_flagged <- sum(df$tx_pathway %in% pathways_flagged)
df <- df %>% filter(!tx_pathway %in% pathways_flagged)
fc(sprintf("Excluding flagged pathway (%s): held pending clinical review",
           paste(pathways_flagged, collapse = ", ")), df)

df <- df %>% filter(tx_pathway %in% pathways_include)
fc("Curative treatment pathway", df)

# tumour site ----------------------------------------------------------------
# oesophageal only by default (C15); see include_sites in the config.
df <- df %>% filter(tumour_site %in% include_sites)
fc(sprintf("Tumour site %s", paste(include_sites, collapse = "/")), df)

# stage 1-3 ------------------------------------------------------------------
# stage_RR is the collapsed 1-4 coding; keep 1-3, drop stage 4 and unstaged.
df <- df %>% filter(stage_RR %in% 1:3)
fc("Stage 1-3", df)

# endoscopy date -------------------------------------------------------------
# the outcome is anchored strictly on the diagnostic endoscopy, so a patient
# with no endoscopy date has no measurable endoscopy-to-DTT wait. Drop them here,
# as their own explicit step, so the number excluded for a missing endoscopy is
# visible in the attrition rather than folded into the outcome-recorded step.
df <- df %>% filter(!is.na(endoscopy_date))
fc("Diagnostic endoscopy date recorded", df)

# valid, positive, plausible endoscopy-to-DTT time ---------------------------
# with an endoscopy date present, a missing outcome here is a missing/invalid
# decision-to-treat date rather than a missing endoscopy.
df <- df %>% filter(!is.na(.data[[outcome_var]]))
df$wait <- df[[outcome_var]]
fc("Endoscopy-to-decision-to-treat time recorded", df)

df <- df %>% filter(wait <= max_wait)
fc(sprintf("Wait <= %d days", max_wait), df)

if (drop_zero_wait) {
  df <- df %>% filter(wait > 0)
  fc("Wait > 0 days (positive endoscopy-to-DTT only)", df)
}

# hospital unit, trust resolution and valid-trust restriction ----------------
# the unit of analysis is the site of diagnosis. Drop patients with no site code
# first (they cannot be assigned to a hospital).
df <- df %>% mutate(diag_hosp_canon = trimws(as.character(.data[[hosp_var]])))
df <- df %>% filter(!is.na(diag_hosp_canon), diag_hosp_canon != "")
fc("Site of diagnosis recorded", df, n_hosp = n_distinct(df$diag_hosp_canon))

# resolve each site's operating trust from the ODS site->trust map. parent_trust
# is authoritative (it handles the sites ODS relocated to a different trust than
# their code prefix); where the map is absent, or a site is not in it, fall back
# to the first three characters of the site code.
if (file.exists(site_trust_map_csv)) {
  site_map <- read.csv(site_trust_map_csv, colClasses = "character") %>%
    distinct(site_code, parent_trust)
  # the hospital's full name, for the per-hospital table further down - kept
  # separate from site_map so df's own columns are unaffected here.
  hosp_name_lookup <- read.csv(site_trust_map_csv, colClasses = "character") %>%
    distinct(site_code, name)
  df <- df %>%
    left_join(site_map, by = c("diag_hosp_canon" = "site_code")) %>%
    mutate(diag_trust = ifelse(is.na(parent_trust) | parent_trust == "",
                               substr(diag_hosp_canon, 1, 3), parent_trust)) %>%
    select(-parent_trust)
} else {
  cat(sprintf("note: no site->trust map at %s - using the site code prefix as the trust\n",
              site_trust_map_csv))
  df <- df %>% mutate(diag_trust = substr(diag_hosp_canon, 1, 3))
  hosp_name_lookup <- tibble(site_code = character(), name = character())
}

# restrict to sites whose operating trust is a recognised OG-cancer diagnoser.
if (file.exists(valid_trusts_csv)) {
  valid_trusts_full <- read.csv(valid_trusts_csv, colClasses = "character")
  valid_trusts <- trimws(valid_trusts_full$trust_code)
  # the trust's full name, for the per-hospital table further down.
  trust_name_lookup <- valid_trusts_full %>%
    mutate(trust_code = trimws(trust_code)) %>%
    distinct(trust_code, trust_name)
  n_before <- nrow(df)
  dropped_trust <- df %>%
    filter(!diag_trust %in% valid_trusts) %>%
    count(diag_trust, name = "patients") %>%
    arrange(desc(patients))
  df <- df %>% filter(diag_trust %in% valid_trusts)
  fc("Site in a recognised OG-diagnosing trust", df,
     n_hosp = n_distinct(df$diag_hosp_canon))
  if (nrow(dropped_trust)) {
    write.csv(dropped_trust, file.path(out_dir, "pd_trusts_excluded.csv"),
              row.names = FALSE)
    cat(sprintf("%d patients at %d trusts not on the valid-diagnoser list (see %s)\n",
                n_before - nrow(df), nrow(dropped_trust), "pd_trusts_excluded.csv"))
  }
} else {
  cat(sprintf("note: no valid-trust list at %s - trust restriction skipped\n",
              valid_trusts_csv))
  trust_name_lookup <- tibble(trust_code = character(), trust_name = character())
}

# per-year diagnosis counts by site, completed over every year in the window so a
# site missing a whole year counts as zero for that year (and so is excluded).
# The year set is taken from the data present; over a clean Jan 2023 - Dec 2024
# cohort that is c(2023, 2024), and both the per-year average and the every-year
# floor divide by / span those two years.
wyears <- sort(unique(df$diagnosis_year))
n_years <- length(wyears)
site_year <- df %>%
  count(diag_hosp_canon, diagnosis_year, name = "n") %>%
  tidyr::complete(diag_hosp_canon, diagnosis_year = wyears, fill = list(n = 0L))

site_min <- site_year %>%
  group_by(diag_hosp_canon) %>%
  summarise(min_year = min(n), total = sum(n), .groups = "drop") %>%
  mutate(per_year = total / n_years)

# a per-hospital table, produced BEFORE the volume floor is applied, so every
# site is visible - including those about to be excluded. One row per hospital
# (site of diagnosis) with its full name and operating trust (code and full
# name), total patients, patients per year (mean over the window), and whether
# it would be kept under the >= min_per_year-in-every-year rule. That rule uses
# the strictest reading: a site must reach min_per_year in each individual year,
# so the yes/no is driven by min_year, not by the average.
site_trust_of <- df %>% distinct(diag_hosp_canon, diag_trust)
hosp_table <- site_min %>%
  left_join(site_trust_of, by = "diag_hosp_canon") %>%
  left_join(hosp_name_lookup, by = c("diag_hosp_canon" = "site_code")) %>%
  left_join(trust_name_lookup, by = c("diag_trust" = "trust_code")) %>%
  group_by(diag_trust) %>%
  mutate(hosps_in_trust = n_distinct(diag_hosp_canon)) %>%
  ungroup() %>%
  mutate(
    # a code with no match in the reference lookup (a decommissioned or
    # miscoded site, or the trust list absent) still shows something, rather
    # than a blank cell.
    name       = ifelse(is.na(name) | name == "", diag_hosp_canon, name),
    trust_name = ifelse(is.na(trust_name) | trust_name == "", diag_trust, trust_name),
    included   = ifelse(min_year >= min_per_year, "Yes", "No")) %>%
  # trusts reporting more than one hospital first, then by trust and hospital
  arrange(desc(hosps_in_trust), diag_trust, desc(total), diag_hosp_canon) %>%
  transmute(
    trust        = diag_trust,
    trust_name   = trust_name,
    hospitals_in_trust = hosps_in_trust,
    hospital     = diag_hosp_canon,
    hospital_name = name,
    patients     = total,
    per_year     = round(per_year, 1),
    min_year_n   = min_year,
    included)

# the included/excluded column gets its full wording only here, right before
# the table is shown and written - the plain "included" name above is what the
# rest of the code works with.
names(hosp_table)[names(hosp_table) == "included"] <-
  sprintf("Included in analysis (>=%d cases per year)", min_per_year)

cat(sprintf("\nPer-hospital patient counts (%d-year window: %s). Multi-hospital",
            n_years, paste(range(wyears), collapse = "-")),
    "trusts first; a site is included only if it reaches",
    sprintf(">= %d diagnoses in every year of the window:\n", min_per_year))
print(as.data.frame(hosp_table), row.names = FALSE)
write.csv(hosp_table, file.path(out_dir, "pd_hospital_counts.csv"),
          row.names = FALSE)
cat(sprintf("(written to %s)\n", "pd_hospital_counts.csv"))

sites_keep <- site_min %>% filter(min_year >= min_per_year) %>% pull(diag_hosp_canon)
sites_drop <- site_min %>% filter(min_year <  min_per_year)

df <- df %>% filter(diag_hosp_canon %in% sites_keep)
fc(sprintf("Site with >= %d diagnoses per year", min_per_year),
   df, n_hosp = n_distinct(df$diag_hosp_canon))

# excluded-hospital output: which sites were dropped for volume, with the yearly
# counts that failed the floor, so the exclusion is auditable.
excl_wide <- site_year %>%
  filter(diag_hosp_canon %in% sites_drop$diag_hosp_canon) %>%
  tidyr::pivot_wider(names_from = diagnosis_year, values_from = n,
                     names_prefix = "y", values_fill = 0) %>%
  left_join(site_min, by = "diag_hosp_canon") %>%
  arrange(desc(total)) %>%
  rename(site = diag_hosp_canon, min_year_count = min_year, total_in_window = total)

write.csv(excl_wide, hosp_excl_csv, row.names = FALSE)
cat(sprintf("\n%d sites excluded on the >= %d/year floor (see %s)\n",
            nrow(sites_drop), min_per_year, basename(hosp_excl_csv)))

# covariates -----------------------------------------------------------------
# recode the rapid-registry fields to the working names the analysis expects.
# gender 1 = Male, 2 = Female; ethnicity_grp 1-5; imd_2019_RR 1 (most deprived)
# to 5; rcs_ch_score is the RCS Charlson score 0/1/2/3.
df <- df %>%
  mutate(
    agediag = as.numeric(age),
    male    = as.integer(gender == 1),
    sexf    = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
    
    cci_n_conditions = as.numeric(rcs_ch_score),
    cci_strata = factor(ifelse(cci_n_conditions >= 2, "2+", "0-1"),
                        levels = c("0-1", "2+")),
    
    stage = factor(stage_RR, levels = c(1, 2, 3), labels = c("1", "2", "3")),
    stage_2 = as.integer(stage == "2"),
    stage_3 = as.integer(stage == "3"),
    
    eth = factor(ethnicity_grp, levels = 1:5,
                 labels = c("White", "Mixed", "Asian", "Black", "Other")),
    
    imd_q = as.integer(imd_2019_RR),
    imd   = factor(imd_q, levels = 1:5),
    
    ydiag = as.integer(diagnosis_year),
    
    canalliance = as.character(diagnosis_alliance),
    
    # calendar-year half (for yr_late) and season dummies, split on the actual
    # midpoint of the endoscopy dates present. Every patient here has an
    # endoscopy date (dropped above otherwise), so the clock-start is complete.
    mid_date = start_date + floor(as.numeric(end_date - start_date) / 2),
    period   = factor(ifelse(endoscopy_date < mid_date, "first", "second"),
                      levels = c("first", "second")),
    yr_late  = as.integer(period == "second"),
    qtr      = quarter(endoscopy_date),
    q2 = as.integer(qtr == 2),
    q3 = as.integer(qtr == 3),
    q4 = as.integer(qtr == 4))

# numeric hospital id for the estimation step (contiguous after the filters)
df <- df %>% mutate(hospid = as.integer(factor(diag_hosp_canon)), hosp = hospid)
df <- df %>% left_join(count(df, diag_hosp_canon, name = "volume"),
                       by = "diag_hosp_canon")

# report ---------------------------------------------------------------------
cat("\npatient attrition:\n")
print(as.data.frame(flow), row.names = FALSE)
cat(sprintf("\nheld-out flagged pathway patients (in window): %d\n", n_flagged))
cat(sprintf("\nanalysis cohort: %d patients, %d diagnosing sites\n",
            nrow(df), n_distinct(df$diag_hosp_canon)))
cat(sprintf("median endoscopy-to-DTT wait %.0f days (IQR %.0f-%.0f)\n",
            median(df$wait), quantile(df$wait, .25), quantile(df$wait, .75)))

saveRDS(df, cohort_rds)
write.csv(as.data.frame(flow), flow_csv, row.names = FALSE)
cat(sprintf("\nsaved %s and %s\n", basename(cohort_rds), basename(flow_csv)))
cat("02 complete. Next: 03 builds Table 1.\n")