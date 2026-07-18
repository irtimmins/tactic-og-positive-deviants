# =============================================================================
# 20a  List the site codes to look up
# -----------------------------------------------------------------------------
# Run this where the COSD file is - the analysis server. It reads the COSD
# extract, pulls out every distinct five-character site code, and writes them to
# a plain list. No internet is needed.
#
# The list is the only thing that then crosses to a machine with internet, where
# 20b_fetch_site_trust_map.R asks ODS about each code. Splitting it this way
# keeps the file-reading on the server and the web lookups off it.
#
# Reads:  the COSD diagnosis extract
# Writes: Data/reference/site_codes_to_lookup.csv   one code per row
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(haven)
})

if (!exists("dir_raw")) dir_raw <- "W:/_DATA/IainTimmins/2026 OG SOTN data"
if (!exists("dir_ref")) dir_ref <- "Data/reference"
dir.create(dir_ref, recursive = TRUE, showWarnings = FALSE)

if (!exists("path_cosd_dta"))
  path_cosd_dta <- file.path(dir_raw,
                             "20260212_all_cosddiagnosis_rapid_202601_OG.dta")

f_site_codes <- file.path(dir_ref, "site_codes_to_lookup.csv")

# the same tidying the build uses, so the codes listed here are the codes the
# build will later try to match: trimmed, upper-cased, and cut at the first run
# of letters and digits.
tidy_org_code <- function(x) {
  x <- str_trim(str_to_upper(as.character(x)))
  x <- str_extract(x, "^[A-Z0-9]+")
  if_else(x == "" | is.na(x), NA_character_, x)
}

cosd <- read_dta(path_cosd_dta, col_select = site_code_of_diagnosis)
cat("COSD rows read:", nrow(cosd), "\n")

tidied <- tidy_org_code(cosd$site_code_of_diagnosis)

# how the raw field breaks down, so it is clear what is being kept and what is
# set aside. Only the well-formed five-character codes go to ODS; the shorter
# and longer ones are trust codes, GP practice codes, placeholders and free
# text, which the build handles itself and ODS would not recognise.
lengths <- tibble(code = tidied) %>%
  mutate(kind = case_when(
    is.na(code)              ~ "blank",
    str_length(code) == 5    ~ "five characters - a site code",
    str_length(code) == 3    ~ "three characters - a trust code",
    str_length(code) == 6    ~ "six characters - a GP practice code",
    TRUE                     ~ "other length")) %>%
  count(kind, name = "rows") %>%
  mutate(pct = round(100 * rows / sum(rows), 1)) %>%
  arrange(desc(rows))

cat("\nWhat the site field holds, across all COSD rows:\n")
print(as.data.frame(lengths), row.names = FALSE)

site_codes <- tidied[!is.na(tidied) & str_length(tidied) == 5]
codes <- sort(unique(site_codes))

out <- tibble(site_code = codes)
write.csv(out, f_site_codes, row.names = FALSE)

cat("\nDistinct five-character site codes:", nrow(out), "\n")
cat("These appear on", length(site_codes), "COSD rows.\n")
cat("\nWrote", f_site_codes, "\n")
cat("Copy it to a machine with internet and run 20b_fetch_site_trust_map.R.\n")