# minimal stand-in for the repo's shared utilities, for offline checking.
suppressPackageStartupMessages({ library(stringr) })

tidy_org_code <- function(x) {
  x <- str_to_upper(str_trim(as.character(x)))
  ifelse(x %in% c("", "."), NA_character_, x)
}

check_input <- function(df, needed, label, extract_path,
                        min_rows = getOption("og_min_input_rows", 1000L)) {
  miss <- setdiff(needed, names(df))
  if (length(miss))
    stop(label, " is missing columns: ", paste(miss, collapse = ", "),
         call. = FALSE)
  if (nrow(df) < min_rows)
    stop(label, " has only ", nrow(df), " rows (expected >= ", min_rows, ").",
         call. = FALSE)
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# title_case()  -  tidy an ALL-CAPS organisation name into natural case
# -----------------------------------------------------------------------------
# ODS returns organisation names in capitals (site_trust_map.csv carries these
# straight through), so any hospital or trust name read from it needs this
# before it goes in a table a person will read. Small joining words stay lower
# case except as the first word; NHS stays upper case.
name_small_words <- c("and", "of", "the", "for", "in", "on", "at", "to", "by", "an", "a", "or")
name_acronyms    <- c("nhs")
title_case <- function(x) {
  vapply(x, function(one) {
    words <- strsplit(tolower(one), " ")[[1]]
    for (i in seq_along(words)) {
      if (words[i] %in% name_acronyms) {
        words[i] <- toupper(words[i])
      } else if (!(words[i] %in% name_small_words) || i == 1) {
        substr(words[i], 1, 1) <- toupper(substr(words[i], 1, 1))
      }
    }
    paste(words, collapse = " ")
  }, character(1), USE.NAMES = FALSE)
}