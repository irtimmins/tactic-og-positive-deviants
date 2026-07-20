# =============================================================================
# _helpers.R  -  functions specific to the CWT merge and waiting-time derivation
# -----------------------------------------------------------------------------
# The one job this stage needs and no other: turning a CWT modality code into
# the broad group the pathway logic reasons about. The grouping itself - which
# codes belong to which group - is a judgement call and stays in
# 01_define_parameters.R, next to the other judgement calls (the transition
# window, the tolerances); this file holds only the function that applies it.
#
# Assumes _load_packages.R has been sourced first.
# =============================================================================

# The guidance modality codes, grouped as the analysis needs them. Codes arrive
# as integers, so 2 is guidance "02" and so on. Code 1 is the surgical code that
# the guidance retired in 2020 ("01"); it still appears in older rows and counts
# as surgery, which the bowel script also does (modality 1 / 23 / 24).
#
#   surgery       1 (retired), 23 (surgery), 24 (surgery, enabling)
#   chemo         2 (cytotoxic), 14 (other), 15 (immunotherapy)
#   hormone       3
#   chemort       4 (chemoradiotherapy)
#   radiotherapy  5 (teletherapy), 6 (brachytherapy), 13 (proton)
#   palliative    7 (specialist palliative), 8 (active monitoring),
#                 9 (non-specialist palliative)
#   other         97
#   declined      98
# Anything else is left ungrouped and set aside.
modality_groups <- list(
  surgery      = c(1L, 23L, 24L),
  chemo        = c(2L, 14L, 15L),
  hormone      = c(3L),
  chemort      = c(4L),
  radiotherapy = c(5L, 6L, 13L),
  palliative   = c(7L, 8L, 9L),
  other        = c(97L),
  declined     = c(98L))

modality_group_of <- function(code) {
  code <- suppressWarnings(as.integer(code))
  out <- rep(NA_character_, length(code))
  for (g in names(modality_groups))
    out[code %in% modality_groups[[g]]] <- g
  out
}