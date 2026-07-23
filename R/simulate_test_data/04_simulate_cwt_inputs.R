# =============================================================================
# 90  Simulate inputs for the CWT merge
# -----------------------------------------------------------------------------
# Writes a stand-in og_cohort_site.rds (the registry cohort with treatment dates
# and intents, as 02 expects) and a stand-in CWT .dta, so the merge can be run
# end to end without the real data. The shapes and the field names match the real
# extracts; the values are made up but plausible - treatment dates a sensible way
# after diagnosis, a mix of pathways, and CWT rows whose decision-to-treat dates
# sit a little before the treatment.
#
# Writes: <dir_sim>/og_cohort_site.rds
#         <dir_sim>/20260212_all_cwt_rapid_202601_OG.dta
# =============================================================================

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(magrittr)
})

if (!exists("dir_sim"))   dir_sim   <- "Output/sim/raw"
if (!exists("sim_scale")) sim_scale <- 1
if (!exists("sim_seed"))  sim_seed  <- 20260212
set.seed(sim_seed)
dir.create(dir_sim, recursive = TRUE, showWarnings = FALSE)

n <- round(9000 * sim_scale)
pid <- sprintf("P%08d", seq_len(n))

trusts <- sort(unique(paste0("R", sample(LETTERS, 400, TRUE),
                             sample(0:9, 400, TRUE))))[1:120]
site_of <- function(tr) paste0(tr, sprintf("%02d", sample(1:30, length(tr), TRUE)))

dx <- as.Date("2018-06-01") + sample(0:2600, n, TRUE)

# assign each patient a latent pathway, then generate dates consistent with it
pathway_type <- sample(
  c("surgery_only", "neoadj_chemo", "neoadj_chemort", "def_chemort",
    "curative_rt", "emr", "sact_only", "pall_rt", "none"),
  n, TRUE,
  prob = c(0.22, 0.14, 0.10, 0.10, 0.05, 0.05, 0.13, 0.06, 0.15))

blankdate <- as.Date(NA)
surgery_date <- rt_date <- sact_date <- emresd_date <- rep(blankdate, n)
surgintent <- rep(NA_integer_, n)
rt_first_intent <- rep(NA_integer_, n)
sact_first_intent_pall <- rep(NA_integer_, n)

is_surg_only  <- pathway_type == "surgery_only"
is_neo_chemo  <- pathway_type == "neoadj_chemo"
is_neo_crt    <- pathway_type == "neoadj_chemort"
is_def_crt    <- pathway_type == "def_chemort"
is_cur_rt     <- pathway_type == "curative_rt"
is_emr        <- pathway_type == "emr"
is_sact_only  <- pathway_type == "sact_only"
is_pall_rt    <- pathway_type == "pall_rt"

surgery_date[is_surg_only] <- dx[is_surg_only] + sample(20:120, sum(is_surg_only), TRUE)
surgintent[is_surg_only] <- 1L

surgery_date[is_neo_chemo] <- dx[is_neo_chemo] + sample(90:200, sum(is_neo_chemo), TRUE)
sact_date[is_neo_chemo]    <- dx[is_neo_chemo] + sample(20:60, sum(is_neo_chemo), TRUE)
surgintent[is_neo_chemo] <- 1L; sact_first_intent_pall[is_neo_chemo] <- 0L

surgery_date[is_neo_crt] <- dx[is_neo_crt] + sample(120:230, sum(is_neo_crt), TRUE)
sact_date[is_neo_crt]    <- dx[is_neo_crt] + sample(20:60, sum(is_neo_crt), TRUE)
rt_date[is_neo_crt]      <- dx[is_neo_crt] + sample(25:65, sum(is_neo_crt), TRUE)
surgintent[is_neo_crt] <- 1L; rt_first_intent[is_neo_crt] <- 3L

sact_date[is_def_crt] <- dx[is_def_crt] + sample(20:70, sum(is_def_crt), TRUE)
rt_date[is_def_crt]   <- dx[is_def_crt] + sample(25:75, sum(is_def_crt), TRUE)
rt_first_intent[is_def_crt] <- 1L; sact_first_intent_pall[is_def_crt] <- 0L

rt_date[is_cur_rt] <- dx[is_cur_rt] + sample(30:90, sum(is_cur_rt), TRUE)
rt_first_intent[is_cur_rt] <- 2L

emresd_date[is_emr] <- dx[is_emr] + sample(15:70, sum(is_emr), TRUE)

sact_date[is_sact_only] <- dx[is_sact_only] + sample(20:80, sum(is_sact_only), TRUE)
sact_first_intent_pall[is_sact_only] <- 1L

rt_date[is_pall_rt] <- dx[is_pall_rt] + sample(20:80, sum(is_pall_rt), TRUE)
rt_first_intent[is_pall_rt] <- 11L

# some patients have an endoscopy a little before diagnosis
endo_date <- dx - sample(5:40, n, TRUE)
endo_date[runif(n) < 0.25] <- blankdate

tr_for <- sample(trusts, n, TRUE)
site_dx <- site_of(tr_for)

registry <- tibble(
  patient_pseudo_id = pid,
  tumour_site   = sample(c("C15", "C16"), n, TRUE, c(0.68, 0.32)),
  diagnosisdate = dx,
  endodateHES   = endo_date,
  EMR_ESDdateHES = emresd_date,
  surgery_date  = surgery_date,
  sact_first_date = sact_date,
  rt_first_date   = rt_date,
  surgintent      = surgintent,
  rt_first_intent = rt_first_intent,
  sact_first_intent_pall = sact_first_intent_pall,
  surgery_trust   = if_else(!is.na(surgery_date), tr_for, NA_character_),
  rt_first_trust  = if_else(!is.na(rt_date), tr_for, NA_character_),
  EMR_ESDtrustHES = if_else(!is.na(emresd_date), tr_for, NA_character_),
  # the fields 02 does not use but that come through from the site-code stage
  site_dx_code  = site_dx,
  site_dx_trust = tr_for,
  diagnosis_trust = tr_for)

saveRDS(registry, file.path(dir_sim, "og_cohort_site.rds"))

# -----------------------------------------------------------------------------
# The CWT extract: one or more treatment periods per patient
# -----------------------------------------------------------------------------
# For most patients, a CWT row whose DTT sits a week or two before the treatment,
# with a modality matching their pathway. Some patients get an extra, earlier
# palliative or other row to exercise the pathway-consistency filter. Some
# patients have no CWT row at all.
modality_for <- function(pt) dplyr::case_when(
  pt %in% c("surgery_only", "neoadj_chemo") ~ 23L,
  pt == "neoadj_chemort" ~ 4L,
  pt == "def_chemort"    ~ 4L,
  pt == "curative_rt"    ~ 5L,
  pt == "emr"            ~ 23L,
  pt == "sact_only"      ~ 2L,
  pt == "pall_rt"        ~ 5L,
  TRUE                   ~ 7L)

# the treatment date each pathway's clock stops on (mirrors 02's first_tx_date)
tx_date_for <- function(i) {
  pt <- pathway_type[i]
  if (pt == "surgery_only")   return(surgery_date[i])
  if (pt == "neoadj_chemo")   return(sact_date[i])
  if (pt == "neoadj_chemort") return(min(sact_date[i], rt_date[i], na.rm = TRUE))
  if (pt == "def_chemort")    return(min(sact_date[i], rt_date[i], na.rm = TRUE))
  if (pt == "curative_rt")    return(rt_date[i])
  if (pt == "emr")            return(emresd_date[i])
  if (pt == "sact_only")      return(sact_date[i])
  if (pt == "pall_rt")        return(rt_date[i])
  as.Date(NA)
}

has_cwt <- runif(n) < 0.88
rows <- list()
for (i in which(has_cwt)) {
  txd <- tx_date_for(i)
  if (is.na(txd)) { txd <- dx[i] + sample(20:120, 1); mg <- 7L } else mg <- modality_for(pathway_type[i])
  # DTT a week or two before treatment; occasionally the treatment date wobbles
  dtt <- txd - sample(3:21, 1)
  treat <- txd + sample(c(-3:3), 1)
  rows[[length(rows) + 1]] <- tibble(
    patient_pseudo_id = pid[i],
    treat_period_start = format(dtt, "%Y-%m-%d"),
    treat_start = format(treat, "%Y-%m-%d"),
    modality = mg,
    site_icd10 = registry$tumour_site[i],
    org_dec_to_treat = tr_for[i],
    org_treat_start = tr_for[i])
  # ~15% get an earlier, non-matching palliative row to test the filter
  if (runif(1) < 0.15) {
    rows[[length(rows) + 1]] <- tibble(
      patient_pseudo_id = pid[i],
      treat_period_start = format(dtt - sample(20:60, 1), "%Y-%m-%d"),
      treat_start = format(dtt - sample(5:15, 1), "%Y-%m-%d"),
      modality = 7L, site_icd10 = registry$tumour_site[i],
      org_dec_to_treat = tr_for[i], org_treat_start = tr_for[i])
  }
}
cwt <- bind_rows(rows)
write_dta(cwt, file.path(dir_sim, "20260212_all_cwt_rapid_202601_OG.dta"))

cat("Simulated", n, "registry patients and", nrow(cwt), "CWT rows in", dir_sim, "\n")
cat("  patients with a CWT row:", sum(has_cwt), "\n")