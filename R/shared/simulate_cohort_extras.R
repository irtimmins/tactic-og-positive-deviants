# =============================================================================
# simulate_cohort_extras.R  -  add treatment, stage and covariate columns to a
#                             simulated registry cohort
# -----------------------------------------------------------------------------
# The site-code stage's simulator writes a registry extract with only the columns
# that stage needs (ids, tumour site, diagnosis date, trust, the site code). The
# CWT merge and positive-deviance stages need more: treatment dates and intents,
# an endoscopy date, stage, and the patient covariates. This helper adds those,
# so a single simulated cohort produced for the site-code stage can be carried
# all the way through the pipeline by the real handoff (run_on_test_data.R does
# exactly this), rather than each stage inventing its own separate fake data.
#
# It takes a data frame that already has, per patient:
#   patient_pseudo_id, tumour_site, a diagnosis date, and the operating trust of
#   the site of diagnosis (any of diag_trust / true_trust / site_dx_trust)
# and returns it with the downstream columns added. Treatment dates and intents
# follow the same scheme as the CWT stage's own simulator, so behaviour is
# consistent wherever the extras come from.
#
# Assumes dplyr is loaded.
# =============================================================================

add_cohort_extras <- function(df, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(df)
  
  # the diagnosis date, whatever it is called in the input
  dx_col <- intersect(c("diagnosisdate", "diagnosis_date"), names(df))[1]
  if (is.na(dx_col)) stop("add_cohort_extras: no diagnosis date column found")
  dx <- as.Date(df[[dx_col]])
  
  # the treating trust, whatever it is called
  trust_col <- intersect(c("diag_trust", "true_trust", "site_dx_trust",
                           "diagnosis_trust"), names(df))[1]
  tr_for <- if (is.na(trust_col)) NA_character_ else as.character(df[[trust_col]])
  
  blankdate <- as.Date(NA)
  pathway_type <- sample(
    c("surgery_only", "neoadj_chemo", "neoadj_chemort", "def_chemort",
      "curative_rt", "emr", "sact_only", "pall_rt", "none"),
    n, TRUE,
    prob = c(0.22, 0.14, 0.10, 0.10, 0.05, 0.05, 0.13, 0.06, 0.15))
  
  surgery_date <- rt_date <- sact_date <- emresd_date <- rep(blankdate, n)
  surgintent <- rt_first_intent <- sact_first_intent_pall <- rep(NA_integer_, n)
  
  is_surg_only <- pathway_type == "surgery_only"
  is_neo_chemo <- pathway_type == "neoadj_chemo"
  is_neo_crt   <- pathway_type == "neoadj_chemort"
  is_def_crt   <- pathway_type == "def_chemort"
  is_cur_rt    <- pathway_type == "curative_rt"
  is_emr       <- pathway_type == "emr"
  is_sact_only <- pathway_type == "sact_only"
  is_pall_rt   <- pathway_type == "pall_rt"
  
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
  
  # an endoscopy a little before diagnosis for most, missing for the rest
  endo_date <- dx - sample(5:40, n, TRUE)
  endo_date[runif(n) < 0.25] <- blankdate
  
  df %>%
    mutate(
      diagnosisdate  = dx,
      endodateHES    = endo_date,
      EMR_ESDdateHES = emresd_date,
      surgery_date   = surgery_date,
      sact_first_date = sact_date,
      rt_first_date   = rt_date,
      surgintent      = surgintent,
      rt_first_intent = rt_first_intent,
      sact_first_intent_pall = sact_first_intent_pall,
      surgery_trust   = ifelse(!is.na(surgery_date), tr_for, NA_character_),
      rt_first_trust  = ifelse(!is.na(rt_date), tr_for, NA_character_),
      EMR_ESDtrustHES = ifelse(!is.na(emresd_date), tr_for, NA_character_),
      
      # stage and the patient covariates the positive-deviance stage needs. Stage
      # 1-3 dominate so the curative cohort is not emptied by the stage filter.
      stage_RR      = sample(1:4, n, TRUE, c(0.30, 0.30, 0.28, 0.12)),
      age           = as.integer(pmin(95, pmax(30, round(rnorm(n, 72, 10))))),
      gender        = sample(1:2, n, TRUE, c(0.72, 0.28)),
      rcs_ch_score  = sample(0:3, n, TRUE, c(0.55, 0.25, 0.14, 0.06)),
      ethnicity_grp = sample(1:5, n, TRUE, c(0.82, 0.04, 0.07, 0.04, 0.03)),
      imd_2019_RR   = sample(1:5, n, TRUE),
      diagnosis_year = as.integer(format(dx, "%Y")),
      diagnosis_alliance = paste0("Alliance ",
                                  sample(LETTERS[1:8], n, TRUE)))
}


# =============================================================================
# make_cwt_extract  -  build a CWT extract keyed on a cohort's own patients
# -----------------------------------------------------------------------------
# The CWT stage merges a CWT extract onto the registry by patient_pseudo_id, so a
# realistic dry run needs the CWT rows to belong to the SAME patients as the
# cohort - not a separately-invented set. This builds that extract from a cohort
# that already carries the treatment dates (as add_cohort_extras leaves it): for
# each patient it works out the pathway's primary treatment date, places a
# decision-to-treat a week or two before it, and writes one row (plus, for some,
# an earlier non-matching palliative row to exercise the consistency filter).
# Most patients get a row; some get none, as in the real data.
#
# Columns match what the CWT stage's 03 reads: patient_pseudo_id,
# treat_period_start, treat_start, modality, site_icd10, org_dec_to_treat,
# org_treat_start. Returns the extract as a data frame (write it out with
# haven::write_dta where the stage expects a .dta).
# =============================================================================

make_cwt_extract <- function(cohort, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(cohort)
  
  trust_col <- intersect(c("diag_trust", "true_trust", "site_dx_trust",
                           "diagnosis_trust"), names(cohort))[1]
  tr <- if (is.na(trust_col)) NA_character_ else as.character(cohort[[trust_col]])
  
  dx  <- as.Date(cohort$diagnosisdate)
  sur <- as.Date(cohort$surgery_date)
  sac <- as.Date(cohort$sact_first_date)
  rt  <- as.Date(cohort$rt_first_date)
  emr <- as.Date(cohort$EMR_ESDdateHES)
  
  # the primary treatment date and its modality group, from the dates present -
  # the same precedence the pathway derivation uses (surgery, then chemo/RT,
  # then EMR), falling back to palliative where only late/other dates exist.
  tx_date  <- as.Date(rep(NA, n), origin = "1970-01-01")
  modality <- rep(7L, n)   # 7 = specialist palliative, the fallback
  has_sur <- !is.na(sur); has_sac <- !is.na(sac)
  has_rt  <- !is.na(rt);  has_emr <- !is.na(emr)
  
  # surgery first
  tx_date[has_sur]  <- sur[has_sur];  modality[has_sur] <- 23L
  # then chemo (neoadjuvant sits before surgery, so surgery already won above;
  # this catches chemo-led pathways)
  pick <- has_sac & is.na(tx_date)
  tx_date[pick] <- sac[pick]; modality[pick] <- 2L
  # then radiotherapy
  pick <- has_rt & is.na(tx_date)
  tx_date[pick] <- rt[pick];  modality[pick] <- 5L
  # then EMR/ESD (counts as surgery)
  pick <- has_emr & is.na(tx_date)
  tx_date[pick] <- emr[pick]; modality[pick] <- 23L
  # anyone still without a date gets a late palliative-style date
  pick <- is.na(tx_date)
  tx_date[pick] <- dx[pick] + sample(20:120, sum(pick), TRUE)
  
  has_cwt <- runif(n) < 0.88
  rows <- vector("list", 0)
  for (i in which(has_cwt)) {
    dtt   <- tx_date[i] - sample(3:21, 1)
    treat <- tx_date[i] + sample(-3:3, 1)
    rows[[length(rows) + 1]] <- data.frame(
      patient_pseudo_id  = as.character(cohort$patient_pseudo_id[i]),
      treat_period_start = format(dtt,   "%Y-%m-%d"),
      treat_start        = format(treat, "%Y-%m-%d"),
      modality           = modality[i],
      site_icd10         = cohort$tumour_site[i],
      org_dec_to_treat   = tr[i],
      org_treat_start    = tr[i],
      stringsAsFactors   = FALSE)
    if (runif(1) < 0.15) {
      rows[[length(rows) + 1]] <- data.frame(
        patient_pseudo_id  = as.character(cohort$patient_pseudo_id[i]),
        treat_period_start = format(dtt - sample(20:60, 1), "%Y-%m-%d"),
        treat_start        = format(dtt - sample(5:15, 1),  "%Y-%m-%d"),
        modality           = 7L,
        site_icd10         = cohort$tumour_site[i],
        org_dec_to_treat   = tr[i],
        org_treat_start    = tr[i],
        stringsAsFactors   = FALSE)
    }
  }
  do.call(rbind, rows)
}