# =============================================================================
# run_3_cwt.R  -  merge the CWT records and derive the waiting times
# -----------------------------------------------------------------------------
# WHERE: the analysis server, on the real data.
# NEEDS: og_cohort_site.rds (from the site-code build) and the CWT extract.
#        Does NOT need internet.
#
# Derives the treatment pathway from the registry (02), then merges the CWT
# records to get the decision-to-treat date and all the waiting-time variables
# (03). Reads og_cohort_site.rds from Data/OG and writes og_cohort_cwt.rds there.
#
# Nothing is set here: the defaults point at the real data path and Data/OG. Set
# a variable before the source lines only if a path differs on the server.
# =============================================================================

# dir_raw <- "W:/_DATA/IainTimmins/2026 OG SOTN data"   # uncomment to override

source("R/merge_cwt_to_get_dtt/02_derive_pathway.R")
source("R/merge_cwt_to_get_dtt/03_cwt_merge.R")
