# =============================================================================
# 02  Derive the treatment pathway
# -----------------------------------------------------------------------------
# Reads og_cohort_site.rds - the registry cohort with the site of diagnosis - and
# derives the treatment pathway from the treatment dates and intents the registry
# already carries. Unlike the bowel rapid script, which reconstructed treatment
# from a separate HES file, the OG rapid extract holds the surgery, SACT, RT and
# EMR/ESD dates itself, so the pathway can be built here directly.
#
# It writes, per patient:
#   the treatment-presence flags   (had_surgery, had_sact, ...)
#   the sequencing flags           (sact_before_surgery, concurrent_chemo_rt)
#   tx_pathway                     the treatment-pathway classification
#   first_tx_date                  the clock-stop date for that pathway
#   tx_trust                       the provider of the clock-stop treatment
#
# The next script (03) merges the CWT records on to this to get the
# decision-to-treat date and the waiting times.
#
# Reads : og_cohort_site.rds
# Writes: og_cohort_pathway.rds   (one row per patient; input to 03)
# =============================================================================

source("R/merge_cwt_to_get_dtt/01_define_parameters.R")

f_cohort_pathway <- file.path(dir_out, "og_cohort_pathway.rds")

if (!exists("read_cohort_site"))
  read_cohort_site <- function() readRDS(f_cohort_site)

og <- read_cohort_site()
cat("Read", nrow(og), "patients from", basename(f_cohort_site), "\n")

# -----------------------------------------------------------------------------
# Pull the treatment fields into tidy, consistently named columns
# -----------------------------------------------------------------------------
# The registry names are kept as they are elsewhere; here they are copied into
# short working names so the pathway logic reads cleanly. Dates from the registry
# are already proper dates (haven reads the Stata day counts as Date); the intent
# fields are integer lookups documented in the registry.
#
#   surgintent          1 = curative, 2 = palliative, 3 = indeterminate
#   rt_first_intent     1 = CRT definitive, 2 = RT definitive, 3 = CRT neoadj,
#                       5 = CRT adjuvant, 11 = palliative, 19 = indeterminate
#   sact_first_intent_pall  1 = palliative intent
as_date <- function(x) if (inherits(x, "Date")) x else as.Date(x)

pw <- og %>%
  mutate(
    emresd_date   = as_date(EMR_ESDdateHES),
    surgery_date  = as_date(surgery_date),
    sact_date     = as_date(sact_first_date),
    rt_date       = as_date(rt_first_date),
    endoscopy_date = as_date(endodateHES),
    diagnosis_date = as_date(diagnosisdate),

    surg_curative = surgintent == 1,
    rt_curative   = rt_first_intent %in% c(1L, 2L, 3L, 5L),
    rt_definitive = rt_first_intent %in% c(1L, 2L),
    rt_palliative = rt_first_intent == 11L,
    sact_palliative = sact_first_intent_pall == 1L,

    surgery_provider = as.character(surgery_trust),
    rt_provider      = as.character(rt_first_trust),
    emresd_provider  = as.character(EMR_ESDtrustHES))

# -----------------------------------------------------------------------------
# Treatment-presence flags
# -----------------------------------------------------------------------------
pw <- pw %>%
  mutate(
    had_emresd           = !is.na(emresd_date),
    had_surgery          = !is.na(surgery_date),
    had_curative_surgery = had_surgery & coalesce(surg_curative, FALSE),
    had_sact             = !is.na(sact_date),
    had_rt               = !is.na(rt_date),
    had_curative_rt      = had_rt & coalesce(rt_curative, FALSE),
    had_palliative_rt    = had_rt & coalesce(rt_palliative, FALSE),
    # chemo that can define a non-surgical definitive chemoRT pathway. In the OG
    # extract SACT chemo always counts; there is no separate HES-only chemo date
    # to guard against here, unlike the ICON pipeline, so the guard is simply
    # "had SACT".
    had_chemo_for_chemort = had_sact)

# -----------------------------------------------------------------------------
# Sequencing flags
# -----------------------------------------------------------------------------
pw <- pw %>%
  mutate(
    sact_before_surgery = had_sact & had_surgery & sact_date < surgery_date,
    sact_after_surgery  = had_sact & had_surgery & sact_date > surgery_date,
    rt_before_surgery   = had_rt   & had_surgery & rt_date   < surgery_date,
    rt_after_surgery    = had_rt   & had_surgery & rt_date   > surgery_date,
    concurrent_chemo_rt = had_sact & had_curative_rt &
                          abs(as.integer(sact_date - rt_date)) <= 14L)

# -----------------------------------------------------------------------------
# tx_pathway - first matching rule wins
# -----------------------------------------------------------------------------
# The same ladder as the ICON pipeline, adapted to the OG registry's fields. Read
# top to bottom: the first rule a patient satisfies is their pathway.
pw <- pw %>%
  mutate(tx_pathway = case_when(
    had_emresd & !had_surgery & !had_sact & !concurrent_chemo_rt
      ~ "EMR/ESD only",
    had_emresd & had_surgery
      ~ "EMR/ESD then surgery",
    had_surgery & sact_before_surgery & rt_before_surgery
      ~ "Surgery + neoadjuvant chemoRT",
    had_surgery & sact_before_surgery & !rt_before_surgery
      ~ "Surgery + neoadjuvant chemo",
    had_surgery & rt_before_surgery & !sact_before_surgery
      ~ "Surgery + neoadjuvant RT",
    had_surgery & sact_after_surgery & !sact_before_surgery
      ~ "Surgery + adjuvant chemo",
    had_surgery & !had_sact & !concurrent_chemo_rt
      ~ "Surgery only",
    had_surgery
      ~ "Surgery + other",
    !had_surgery & had_curative_rt & had_chemo_for_chemort
      ~ "Definitive chemoRT",
    !had_surgery & had_curative_rt & !had_chemo_for_chemort
      ~ "Curative RT only",
    !had_surgery & had_palliative_rt & had_sact
      ~ "Palliative chemo + RT",
    !had_surgery & had_sact & !had_curative_rt
      ~ "SACT only",
    !had_surgery & had_palliative_rt & !had_sact
      ~ "Palliative RT only",
    TRUE
      ~ "No treatment recorded"))

# -----------------------------------------------------------------------------
# first_tx_date - the clock-stop date for the pathway
# -----------------------------------------------------------------------------
# The clock stops at the first definitive treatment. For the neoadjuvant
# pathways that is the neoadjuvant treatment itself, not the later surgery -
# chemotherapy and radiotherapy given before surgery are both first definitive
# treatments under the guidance (3.9.1, 3.11.1). So neoadjuvant chemo stops the
# clock on the chemo date, neoadjuvant RT on the RT date, and neoadjuvant
# chemoRT on the earlier of the two. The curative act being the surgery does not
# change which treatment started first. Where two dates could apply, the earlier
# is taken.
pmin_date <- function(a, b) as.Date(pmin(as.integer(a), as.integer(b),
                                         na.rm = TRUE), origin = "1970-01-01")

pw <- pw %>%
  mutate(first_tx_date = case_when(
    tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery") ~ emresd_date,
    tx_pathway == "Surgery + neoadjuvant chemoRT" ~ pmin_date(sact_date, rt_date),
    tx_pathway == "Surgery + neoadjuvant RT"      ~ rt_date,
    tx_pathway == "Surgery + neoadjuvant chemo"   ~ sact_date,
    tx_pathway %in% c("Surgery + adjuvant chemo", "Surgery only",
                      "Surgery + other")          ~ surgery_date,
    tx_pathway == "Definitive chemoRT"            ~ pmin_date(sact_date, rt_date),
    tx_pathway == "Curative RT only"              ~ rt_date,
    TRUE ~ as.Date(NA)))   # non-curative pathways leave the clock-stop unset

# -----------------------------------------------------------------------------
# tx_trust - provider of the clock-stop treatment (three characters)
# -----------------------------------------------------------------------------
# Surgical and EMR pathways take the surgery provider; RT-anchored pathways take
# the RT provider. SACT's provider is never the source: neoadjuvant chemo still
# takes surgery's, because the curative act is the surgery.
pw <- pw %>%
  mutate(tx_trust = case_when(
    tx_pathway %in% c("EMR/ESD only", "EMR/ESD then surgery",
                      "Surgery + neoadjuvant chemo", "Surgery + adjuvant chemo",
                      "Surgery only", "Surgery + other")
      ~ str_sub(coalesce(surgery_provider, emresd_provider), 1, 3),
    tx_pathway %in% c("Surgery + neoadjuvant chemoRT", "Surgery + neoadjuvant RT",
                      "Definitive chemoRT", "Curative RT only")
      ~ str_sub(rt_provider, 1, 3),
    TRUE ~ NA_character_))

# -----------------------------------------------------------------------------
# Report the pathway mix and save
# -----------------------------------------------------------------------------
cat("\nDerived pathway mix:\n")
pw %>%
  count(tx_pathway, name = "patients", sort = TRUE) %>%
  mutate(pct = round(100 * patients / sum(patients), 1)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

keep_cols <- c("patient_pseudo_id", "tumour_site",
  "diagnosis_date", "endoscopy_date",
  "emresd_date", "surgery_date", "sact_date", "rt_date",
  "had_emresd", "had_surgery", "had_curative_surgery", "had_sact", "had_rt",
  "had_curative_rt", "had_palliative_rt",
  "sact_before_surgery", "sact_after_surgery", "rt_before_surgery",
  "concurrent_chemo_rt",
  "tx_pathway", "first_tx_date", "tx_trust",
  "surgery_provider", "rt_provider", "emresd_provider")
keep_cols <- intersect(keep_cols, names(pw))

saveRDS(pw[keep_cols], f_cohort_pathway)
cat("\nSaved", f_cohort_pathway, "\n")
cat("02 complete. Next: 03 merges the CWT records for the decision-to-treat date.\n")
