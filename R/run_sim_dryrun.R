# =============================================================================
# run_sim_dryrun.R  -  the whole thing on made-up data
# -----------------------------------------------------------------------------
# WHERE: anywhere - your laptop, the internet machine, the server. It touches no
#        real data and needs no internet.
# NEEDS: nothing.
#
# This proves the pipeline end to end on simulated data that matches the shape of
# the real extracts. Use it to check the code runs on a given machine before the
# real run, and to see what the outputs look like without waiting on the server.
#
# It writes everything to Data/sim and Data/OG_test, well away from the real
# Data/OG, and sets dir_transfer to NULL so made-up numbers never reach the real
# transfer folder.
#
# What it does, in order:
#   1. makes simulated COSD, registry, and both reference maps (90)
#   2. runs the real build scripts (00_master -> 02, 03) against them
# The check scripts (91, 11, 21) are the formal tests; run those separately if
# you want the pass/fail counts.
# =============================================================================

dir_sim <- "Data/sim"
sim_scale <- 1                      # 1 = full size; 0.1 for a quick small run

# 90 builds the simulated .dta files and the two reference maps into dir_sim
source("R/derive_5_digit_site_code/90_simulate_inputs.R")

# hand the build the simulated data instead of the real extracts
rm(list = setdiff(ls(), c("dir_sim", "sim_scale")))

dir_raw      <- "Data/sim"          # read the simulated extracts
dir_out      <- "Data/OG_test"      # write well away from the real Data/OG
dir_ref      <- "Data/sim"          # the simulated maps live here
dir_transfer <- NULL                # never touch the real transfer folder

source("R/derive_5_digit_site_code/00_master.R")

message("\nSimulated dry run complete. Outputs in Data/OG_test (not the real ",
        "Data/OG). For the formal pass/fail tests, run ",
        "R/derive_5_digit_site_code/91_check_site_logic.R.")
