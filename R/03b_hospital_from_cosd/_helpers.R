# =============================================================================
# _helpers.R  -  functions specific to the five-character site-code build
# -----------------------------------------------------------------------------
# The parsing jobs this stage needs and no other: pulling a three-character site
# from an ICD-O topography, morphology as four digits, and a SNOMED identifier as
# text. The generic org-code tidier and the input guard are not here - they are
# shared, in shared/utils.R.
#
# Assumes _load_packages.R has been sourced first (these use stringr).
# =============================================================================

# The three-character cancer site from an ICD-O topography. The field usually
# holds a bare code but sometimes carries the description too ("C155: LOWER THIRD
# OF OESOPHAGUS"), and occasionally holds filler such as XXXXXX. Take the leading
# letter-digit-digit and nothing else; anything else becomes missing.
icd_site3 <- function(x) {
  x <- str_to_upper(str_trim(as.character(x)))
  str_extract(x, "^[A-Z][0-9]{2}")
}

# Morphology as four digits, however it arrives: str4 in COSD, str5 in the
# registry where the fifth character can be the behaviour.
morph4 <- function(x) {
  str_extract(str_trim(as.character(x)), "^[0-9]{4}")
}

# A SNOMED identifier as text. It is only ever compared with itself and with the
# map, never used as a number.
#
# If the field arrives as text, which is what we have asked NDRS for and what
# Sarah's bowel extract already has, it is used as it stands. Older extracts have
# it as a Stata double, which cannot hold the longest SNOMED numbers exactly -
# anything above about 9e15 is already approximate by the time it reaches us.
# Those are converted and a warning is printed, because a number that has lost a
# digit will not match the map and there is nothing this code can do to recover
# it.
snomed_id <- function(x) {
  if (is.character(x)) {
    x <- str_trim(x)
    return(if_else(x == "" | x == "." | x == "0", NA_character_, x))
  }
  x <- suppressWarnings(as.numeric(x))
  too_big <- sum(!is.na(x) & x > 2^53)
  if (too_big)
    warning(too_big, " SNOMED numbers are too large for the numeric field they ",
            "arrived in and may have lost precision. Ask for the field as text.",
            call. = FALSE)
  ifelse(is.na(x) | x <= 0, NA_character_, sprintf("%.0f", x))
}