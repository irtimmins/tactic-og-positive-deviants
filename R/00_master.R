# =============================================================================
# 00  Master  -  run the OG rapid build end to end
# -----------------------------------------------------------------------------
# Starts from the 20260212 rapid tumour extract and adds the fields the registry
# does not carry. Step 02 brings in the five-character site of diagnosis from the
# COSD extract; the CWT merge follows.
#
# Set dir_raw below to wherever the extracts sit. dir_out is where everything
# derived is written.
#
# The validation scripts prove the build code against fixtures and against
# simulated data, so they can be run anywhere - they never read dir_raw and never
# write dir_out. Run them after any change to the rules in 01.
# =============================================================================

dir_raw <- "W:/_DATA/IainTimmins/2026 OG SOTN data"
dir_out <- "Data/OG"

dir_build <- "R/build_og_site"
step <- function(file) {
  message("\n========== ", file, " ==========")
  source(file.path(dir_build, file), local = new.env())
}

step("02_add_site_of_diagnosis.R")   # COSD -> site_dx_code on the registry

# validation - safe to run anywhere, including off the analysis server
# step("91_check_site_logic.R")

message("\nBuild complete: ", file.path(dir_out, "og_cohort_site.rds"))
