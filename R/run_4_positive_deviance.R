# =============================================================================
# run_4_positive_deviance.R  -  identify positive-deviant hospitals
# -----------------------------------------------------------------------------
# WHERE: the analysis server, on the real data.
# NEEDS: og_cohort_cwt.rds (from run_3_cwt.R), the site-to-trust map and the
#        valid-diagnosing-trusts list under Data/reference. Does NOT need
#        internet.
#
# Builds the curative, stage 1-3, C15 positive-deviance cohort (02), Table 1
# (03), the balancing-weights standardisation (04), Bayesian shrinkage (05), and
# the ranking / caterpillar plot / candidate selection (06).
#
# dir_transfer is set here, not in the common paths file, because this is the
# one stage whose whole purpose is producing results to move off the server -
# setting it as a blanket default in paths_and_helpers.R would risk a simulated
# or exploratory run from any stage quietly writing to S: by accident. Every
# other run_N script leaves it NULL.
# =============================================================================

# dir_raw <- "W:/_DATA/IainTimmins/2026 OG SOTN data"   # uncomment to override

# turn on writing to the results-transfer area. transfer_root is the canonical
# S: path, defined once in R/config/directories.R (sourced here first so it is
# available); setting dir_transfer to it is what sends this stage's tables and
# figures off the encrypted W: drive to somewhere they can leave the server.
source("R/config/directories.R")
dir_transfer <- transfer_root

source("R/identify_positive_deviants/00_master.R")