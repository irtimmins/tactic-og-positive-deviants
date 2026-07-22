# =============================================================================
# 02  List the site codes to look up
# -----------------------------------------------------------------------------
# Run this where the data is - the analysis server. It reads every source that
# records a hospital site, pulls out the distinct five-character codes, and
# writes them to a plain list. No internet is needed.
#
# Two sources, because the project now uses two different site fields:
#   COSD  site_code_of_diagnosis   the hospital that made the diagnosis
#   HES   sitetret                 the hospital that treated the episode, which
#                                  is where the endoscopy site comes from
#
# They do not hold the same codes. COSD is filled in by cancer services and is
# mostly NHS trust sites; HES sitetret covers everything an episode can be
# delivered at, so it also carries independent-sector sites (the NT4xx, NVCxx and
# NPGxx families), the newer merged-trust codes (R0A02, R1HM0), and a tail of
# placeholders and test values. Looking up only the COSD codes would leave the
# endoscopy sites unmapped, so both go to ODS.
#
# The list is the only thing that then crosses to a machine with internet, where
# stage 01 asks ODS about each code. Splitting it this way keeps the file
# reading on the server and the web lookups off it.
#
# Reads:  the COSD diagnosis extract, the HES-APC extract
# Writes: Data/reference/site_codes_to_lookup.txt   one code per row, with the
#         source(s) it came from
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(haven)
})

source("R/config/directories.R")   # dir_raw, dir_ref, path_cosd_dta, HES paths
source("R/shared/utils.R")         # tidy_org_code

f_site_codes <- file.path(dir_ref, "site_codes_to_lookup.txt")

# how the raw field breaks down, so it is clear what is being kept and what is
# set aside. Only the well-formed five-character codes go to ODS; the shorter and
# longer ones are trust codes, GP practice codes, placeholders and free text,
# which ODS would not recognise as sites.
#
# In HES a three-character sitetret is the provider rather than the site. It is
# still reported here, because knowing how much of the field is provider-only is
# worth having, but it is not sent to ODS as a site.
describe_field <- function(x, label) {
  tibble(code = x) %>%
    mutate(kind = case_when(
      is.na(code)           ~ "blank",
      str_length(code) == 5 ~ "five characters - a site code",
      str_length(code) == 3 ~ "three characters - a trust or provider code",
      str_length(code) == 6 ~ "six characters - a GP practice code",
      TRUE                  ~ "other length")) %>%
    count(kind, name = "rows") %>%
    mutate(pct = round(100 * rows / sum(rows), 1)) %>%
    arrange(desc(rows)) %>%
    mutate(source = label)
}

five_char <- function(x) {
  keep <- !is.na(x) & str_length(x) == 5
  sort(unique(x[keep]))
}

# -----------------------------------------------------------------------------
# COSD: the site of diagnosis
# -----------------------------------------------------------------------------
cosd <- read_dta(path_cosd_dta, col_select = site_code_of_diagnosis)
cosd_tidied <- tidy_org_code(cosd$site_code_of_diagnosis)
cat("COSD rows read:", nrow(cosd), "\n")

# -----------------------------------------------------------------------------
# HES: the site of treatment
# -----------------------------------------------------------------------------
# Prefer the cohort-filtered extract the endoscopy stage builds, since it is far
# quicker to read; fall back to the raw APC file, taking only the one column.
# Reading the whole raw file for a single column is still slow, so a note is
# printed rather than leaving it looking hung.
hes_source <- NULL
if (exists("f_hes_extract") && file.exists(f_hes_extract)) {
  hes <- readRDS(f_hes_extract)
  hes_source <- f_hes_extract
} else if (file.exists(path_hes_apc_dta)) {
  cat("No cohort-filtered HES extract - reading sitetret from the raw APC file,",
      "which takes a few minutes.\n")
  hes <- read_dta(path_hes_apc_dta, col_select = sitetret)
  hes_source <- path_hes_apc_dta
} else {
  hes <- NULL
}

if (is.null(hes)) {
  cat("\nNo HES extract found at either\n  ", f_hes_extract, "\n  ",
      path_hes_apc_dta,
      "\n  - the list will cover the COSD codes only. Run this again once the",
      "HES extract is in place, or the endoscopy sites will be unmapped.\n")
  hes_tidied <- character(0)
} else {
  hes_tidied <- tidy_org_code(hes$sitetret)
  cat("HES episodes read:", nrow(hes), "from", basename(hes_source), "\n")
}

# -----------------------------------------------------------------------------
# What each field holds
# -----------------------------------------------------------------------------
breakdown <- bind_rows(
  describe_field(cosd_tidied, "COSD site_code_of_diagnosis"),
  if (length(hes_tidied)) describe_field(hes_tidied, "HES sitetret"))

cat("\nWhat the site fields hold:\n")
breakdown %>%
  select(source, kind, rows, pct) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# -----------------------------------------------------------------------------
# The union of the five-character codes
# -----------------------------------------------------------------------------
cosd_codes <- five_char(cosd_tidied)
hes_codes  <- five_char(hes_tidied)
codes <- sort(union(cosd_codes, hes_codes))

out <- tibble(site_code = codes) %>%
  mutate(in_cosd = site_code %in% cosd_codes,
         in_hes  = site_code %in% hes_codes,
         source  = case_when(in_cosd & in_hes ~ "both",
                             in_cosd          ~ "COSD only",
                             TRUE             ~ "HES only"))

cat("\nDistinct five-character site codes:\n")
out %>%
  count(source, name = "codes") %>%
  arrange(desc(codes)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
cat("  total to look up:", nrow(out), "\n")

# The HES-only codes are the ones this change is for - without them the
# endoscopy sites have no trust. Show a few so it is clear what they look like.
hes_only <- out %>% filter(source == "HES only")
if (nrow(hes_only)) {
  cat("\nA sample of the codes HES contributes that COSD never had:\n")
  cat("  ", paste(head(hes_only$site_code, 20), collapse = " "), "\n")
  cat("  (", nrow(hes_only), "in total - largely independent-sector sites and",
      "the newer merged-trust codes)\n")
}

# .txt rather than .csv: this file has to leave the server, and something on the
# transfer path appends an encrypted footer to .csv files - the content survives
# but the file breaks for anything reading to the end. Plain tab-separated text
# goes through untouched.
write.table(out, f_site_codes, sep = "\t", quote = FALSE, row.names = FALSE,
            fileEncoding = "ASCII")

cat("\nWrote", f_site_codes, "\n")
cat("Copy it to a machine with internet and run stage 01 there.\n")