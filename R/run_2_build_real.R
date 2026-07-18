# =============================================================================
# run_2_build_real.R  -  the real build, on the server
# -----------------------------------------------------------------------------
# WHERE: the analysis server, on the real COSD data.
# NEEDS: the real extracts, plus the two maps built on the internet machine:
#          Data/reference/snomed_og_lookup.csv
#          Data/reference/site_trust_map.csv
#        Does NOT need internet.
#
# This is the real run. It derives the five-character site of diagnosis for
# every patient (02), then reports on how that was done (03). The build log and
# the diagnostics files are copied to the transfer folder at the end.
#
# Nothing is set here on purpose: the defaults in 00_master.R and 01 point at the
# real data path, the real Data/reference, and the S: transfer folder. If any of
# those differ on the server, set the variable before the source() line below.
#
# Before running, it is worth confirming the two maps are in place - 02 prints
# "SNOMED map read from ..." and "ODS site-to-trust map read from ..." early on.
# If either says "No ... map", stop and check Data/reference.
# =============================================================================

# dir_raw <- "W:/_DATA/IainTimmins/2026 OG SOTN data"   # uncomment to override

source("R/derive_5_digit_site_code/00_master.R")
