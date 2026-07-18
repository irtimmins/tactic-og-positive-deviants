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
# (wt_anchor_to_dtt from the merge stage), and only positive values are analysed.
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

# analysis window ------------------------------------------------------------
# the analysis uses the most recent stretch of diagnoses. window_end defaults to
# the latest endoscopy-anchored diagnosis in the cohort; set it explicitly to fix
# the window (e.g. as.Date("2024-12-31")).
og <- og %>% mutate(anchor_date = as.Date(anchor_date))
end_date   <- if (is.na(window_end)) max(og$anchor_date, na.rm = TRUE) else as.Date(window_end)
start_date <- end_date %m-% months(window_months) %m+% days(1)

# future option: a fixed, precise diagnosis window rather than "the last N
# months". Not applied yet - left here, commented out, as a marker for when the
# analysis moves to a specific reporting period. Uncomment and set window_months
# to NA (or just replace the two lines above) when that decision is made.
#   start_date <- as.Date("2023-01-01")
#   end_date   <- as.Date("2025-03-31")

cat(sprintf("window: %s to %s\n", start_date, end_date))

# main audit cohort -----------------------------------------------------------
# sotn_cohort flags the registry's own audit-eligible cohort (State of the
# Nation exclusions already applied - stage, route, DCO and similar validity
# checks at the registry level). Restricting to it here, before anything else,
# means every later exclusion count is relative to a population already agreed
# to be legitimately diagnosed OG cancer.
df <- og %>% filter(sotn_cohort == 1)
fc("Main audit cohort (sotn_cohort == 1)", df)

# window -----------------------------------------------------------------
df <- df %>% filter(!is.na(anchor_date), anchor_date >= start_date,
                    anchor_date <= end_date)
fc(sprintf("Diagnosed %s to %s", format(start_date, "%b %Y"),
           format(end_date, "%b %Y")), df)

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

# valid, positive, plausible endoscopy-to-DTT time ---------------------------
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
  df <- df %>%
    left_join(site_map, by = c("diag_hosp_canon" = "site_code")) %>%
    mutate(diag_trust = ifelse(is.na(parent_trust) | parent_trust == "",
                               substr(diag_hosp_canon, 1, 3), parent_trust)) %>%
    select(-parent_trust)
} else {
  cat(sprintf("note: no site->trust map at %s - using the site code prefix as the trust\n",
              site_trust_map_csv))
  df <- df %>% mutate(diag_trust = substr(diag_hosp_canon, 1, 3))
}

# restrict to sites whose operating trust is a recognised OG-cancer diagnoser.
if (file.exists(valid_trusts_csv)) {
  valid_trusts <- read.csv(valid_trusts_csv, colClasses = "character")$trust_code
  valid_trusts <- trimws(valid_trusts)
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
}

# per-year diagnosis counts by site, completed over every year in the window so a
# site missing a whole year counts as zero for that year (and so is excluded).
wyears <- sort(unique(df$diagnosis_year))
site_year <- df %>%
  count(diag_hosp_canon, diagnosis_year, name = "n") %>%
  tidyr::complete(diag_hosp_canon, diagnosis_year = wyears, fill = list(n = 0L))

site_min <- site_year %>%
  group_by(diag_hosp_canon) %>%
  summarise(min_year = min(n), total = sum(n), .groups = "drop")

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

    # window halves for the calendar-year term and season dummies
    mid_date = start_date %m+% months(window_months / 2),
    period   = factor(ifelse(anchor_date < mid_date, "first", "second"),
                      levels = c("first", "second")),
    yr_late  = as.integer(period == "second"),
    qtr      = quarter(anchor_date),
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
