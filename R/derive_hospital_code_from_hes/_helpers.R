# =============================================================================
# _helpers.R  -  functions specific to the endoscopy-site-from-HES build
# -----------------------------------------------------------------------------
# The two parsing jobs this stage needs and no other: normalising an OPCS code so
# it matches the code list, and reading a HES date field whatever form it arrives
# in. The generic org-code tidier and the input guard are not here - they are
# shared, in shared/utils.R.
#
# Assumes _load_packages.R has been sourced first (these use stringr/lubridate).
# =============================================================================

# An OPCS-4 code, ready to match against the code list: upper case, no dot, no
# surrounding space. The HES field also uses "-" and "&" as fillers for an empty
# slot; those simply will not be in the code list, so they fall out on their own.
norm_opcs <- function(x) {
  str_replace_all(str_to_upper(str_trim(as.character(x))), "\\.", "")
}

# A HES date, however it reaches us. The APC extract stores its dates as text in
# year-month-day order ("2022-03-10", occasionally "20220310"); a value read from
# the rapid record through haven arrives already as a Date, or as a Stata day
# count. Handle all three so the same code works on the real extract and on the
# stand-in data, and so an unreadable value becomes NA rather than an error.
#
# HES also uses 1800-01-01 and 1801-01-01 as null-date sentinels. This extract
# starts in 2017, so anything on or before 1801 is a sentinel and is returned as
# missing rather than as a real early date.
hes_date <- function(x) {
  d <- if (inherits(x, "Date")) {
    x
  } else if (is.numeric(x)) {
    as.Date(x, origin = "1960-01-01")
  } else {
    # text, usually "2022-03-10" but occasionally "20220310". as.Date needs one
    # format to fit the whole vector, so a mix defeats a single call - parse each
    # form and take whichever hit.
    s <- str_trim(as.character(x))
    coalesce(suppressWarnings(as.Date(s, format = "%Y-%m-%d")),
             suppressWarnings(as.Date(s, format = "%Y%m%d")))
  }
  d[!is.na(d) & d <= as.Date("1801-01-01")] <- NA
  d
}
