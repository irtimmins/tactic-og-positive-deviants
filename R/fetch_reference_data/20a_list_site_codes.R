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
# Writes: Data/reference/site_codes_to_lookup.txt   one code per row
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(haven)
})

source("R/config/directories.R")   # dir_raw, dir_ref, path_cosd_dta
source("R/shared/utils.R")          # tidy_org_code

f_site_codes <- file.path(dir_ref, "site_codes_to_lookup.txt")

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
# .txt rather than .csv: this file has to leave the server, and something on
# the transfer path appends an encrypted footer to .csv files - the content
# survives but the file breaks for anything reading to the end. Plain
# tab-separated text goes through untouched.
write.table(out, f_site_codes, sep = "\t", quote = FALSE, row.names = FALSE,
            fileEncoding = "ASCII")

cat("\nDistinct five-character site codes:", nrow(out), "\n")
cat("These appear on", length(site_codes), "COSD rows.\n")
cat("\nWrote", f_site_codes, "\n")
cat("Copy it to a machine with internet and run 20b_fetch_site_trust_map.R.\n")