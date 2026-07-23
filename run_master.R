# =============================================================================
# run_master.R  -  run the pipeline, a block at a time
# -----------------------------------------------------------------------------
# Sourcing this file runs nothing. Find the block below that matches where you
# are sitting, and run the line you want.
#
# It has to be used a block at a time because the pipeline crosses an airgap:
# stages 1 and 2 fetch from the internet and cannot run on the analysis server,
# and everything else reads patient data and cannot run anywhere else. Stages 0
# and 2 each end by printing which file to carry to the other machine.
#
# The stages, in order:
#   0   list the hospital codes             server    reads COSD + HES
#   1   map the hospitals to their trusts   internet  ODS API
#   2   fetch the SNOMED tumour map         internet  TRUD download
#   3b  the hospital of diagnosis, from COSD  server
#   3a  the endoscopy hospital, from HES      server
#   4   add the decision-to-treat           server    CWT merge
#   5   derive the analysis cohort          server    funnel, flowchart, table 1
#   6   identify the positive deviants      server    weights, shrinkage, ranks
#
# 3b runs before 3a: stage 4 reads the cohort 3b writes, and stage 5 joins the
# lookup 3a writes. Asking for stage 3 does both, in that order.
#
# The machinery is in R/shared/stage_runner.R. There is no need to read it.
# =============================================================================

source("R/shared/stage_runner.R")


# -----------------------------------------------------------------------------
# A.  Dry run  -  anywhere, on made-up data
# -----------------------------------------------------------------------------
# Touches no real data and needs no internet. It builds a full set of simulated
# extracts and writes everything - raw fixtures, patient-level intermediates,
# and the aggregate outputs (table 1, the flowchart, Stan fits) - under
# Output/sim, then runs the same stage code a real run does, file for file.
# Nothing can reach the transfer area. Output/sim can be large; it is safe to
# delete and rebuild any time by running this block again.
# Do this first on a new machine, to check the whole thing runs.
#
# For pass/fail logic checks rather than output to look at, run instead:
#   Rscript R/tests/run_tests.R

# Generate synthetic list of hospital 5-digit codes.
#run_master(stages = 0,                        mode = "simulated")
# Build analysis cohort using hes/cosd info.
# run_master(stages = c("3b", "3a", 4, 5),      mode = "simulated")
# Test run of positive deviance
# run_master(stages = 6,                        mode = "simulated")   # needs balancer + rstan


# -----------------------------------------------------------------------------
# B.  The real run  -  on the analysis server
# -----------------------------------------------------------------------------
# Stage 0 only needs running when the raw extracts change: it reads the large
# HES file once and writes both the cohort-filtered copy and the site-code list
# to carry to the internet machine.
#
# Stages 3-6 are the analysis proper, and need site_trust_map.csv and
# snomed_og_lookup.csv already copied back into Data/reference from block C.
#
# Set dir_transfer first if the results should leave the server - without it,
# tables and figures are written to a local folder instead of the S: area.

# dir_transfer <- transfer_root

# run_master(stages = 0)        # first time, or after a new data cut
# run_master(stages = 3:6)      # the usual run: 3b, 3a, 4, 5, 6

# or a stage at a time, to inspect the output as it goes:
# run_master(stages = "3b")
# run_master(stages = "3a")
# run_master(stages = 4)
# run_master(stages = 5)
# run_master(stages = 6)


# -----------------------------------------------------------------------------
# C.  The reference lookups  -  on a machine with internet
# -----------------------------------------------------------------------------
# Neither stage reads patient data: stage 1 takes only the code list from stage
# 0, and stage 2 downloads a public NHS release. That is what makes them safe to
# run off the server.
#
# Stage 2 needs a TRUD account, a subscription to item 101, and TRUD_API_KEY in
# .Renviron. Both stages print which file to copy back to Data/reference when
# they finish.

# run_master(stages = 1:2)