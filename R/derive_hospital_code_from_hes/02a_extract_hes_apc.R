# =============================================================================
# 02a  Cut the HES-APC extract down to the cohort
# -----------------------------------------------------------------------------
# The raw HES-APC file covers far more patients and columns than this stage
# needs, and reading it takes several minutes. This makes the smaller, faster
# copy 02 and 03 read (f_hes_extract): the same episodes, restricted to the
# patients in the rapid extract and to the columns the endoscopy match uses.
#
# It only needs running when the raw HES file changes. If f_hes_extract is
# already there it is left alone; set refresh_hes <- TRUE before sourcing to
# rebuild it anyway.
#
# This replaces doing the same thing by hand at the console, which is how the
# first version was made - and which produced an extract holding only a few
# thousand patients, because it was filtered against a working subset rather
# than the whole cohort. Filtering against the rapid extract here means the
# denominator in 02 is the real one.
#
# Reads : the rapid tumour dta, the raw HES-APC dta
# Writes: hes_extract.rds
# =============================================================================

source("R/derive_hospital_code_from_hes/01_define_parameters.R")

if (!exists("refresh_hes")) refresh_hes <- FALSE

if (file.exists(f_hes_extract) && !refresh_hes) {
  cat("HES extract already at", f_hes_extract,
      "- left alone. Set refresh_hes <- TRUE to rebuild it.\n")
} else {
  op_cols     <- sprintf("opertn_%02d", 1:24)
  opdate_cols <- sprintf("opdate_%02d", 1:24)
  keep_cols <- c("patient_pseudo_id", "epistart", "epiend", "admidate",
                 "epiorder", "epitype", "sitetret", "procode3",
                 op_cols, opdate_cols)

  rapid <- read_dta(path_rapid_dta, col_select = "patient_pseudo_id")
  ids <- unique(as.character(rapid$patient_pseudo_id))
  cat("Cohort patients to keep:", length(ids), "\n")

  cat("Reading", path_hes_apc_dta, "- this is large and takes a few minutes.\n")
  hes <- read_dta(path_hes_apc_dta, col_select = all_of(keep_cols))
  cat("Read", nrow(hes), "episodes for",
      n_distinct(hes$patient_pseudo_id), "patients\n")

  hes <- hes %>%
    mutate(patient_pseudo_id = as.character(patient_pseudo_id)) %>%
    filter(patient_pseudo_id %in% ids)

  saveRDS(hes, f_hes_extract)
  cat("Kept", nrow(hes), "episodes for", n_distinct(hes$patient_pseudo_id),
      "patients ->", f_hes_extract, "\n")
}

cat("02a complete. Next: 02_add_endoscopy_site.R\n")
