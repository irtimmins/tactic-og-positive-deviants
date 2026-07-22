# =============================================================================
# 02  Fetch the SNOMED tumour map  -  run this step
# -----------------------------------------------------------------------------
# WHERE: a machine with internet access. NOT the analysis server.
#
# Nothing here reads patient data. The input is a public NHS release and the
# output is a code lookup, so this stage is safe to run entirely off the server -
# the same rule as stage 01.
#
# NHS England publishes a SNOMED-to-ICD-10 map as part of the SNOMED CT UK
# Clinical Edition. This downloads that release and cuts it to the codes that
# mean an oesophageal (C15) or gastric (C16) tumour, so the COSD site build
# (stage 03b) can read a tumour's meaning off a lookup rather than inferring it.
#
# Before running, three things are needed, and an API key alone is not enough:
#   1. a TRUD account            https://isd.digital.nhs.uk/trud
#   2. a subscription to item 101, "SNOMED CT UK Clinical Edition, RF2: Full,
#      Snapshot & Delta", with its licence accepted
#   3. the API key in a .Renviron file as TRUD_API_KEY=your-key, then restart R
# The full notes are at the top of 01_fetch_snomed_map.R.
#
# The download is about 1 GB and only one file inside it is used, so this takes
# a while. TRUD asks that it be run on a weekday between 8am and 6pm, or between
# midnight and 6am, to avoid their maintenance windows.
#
# Produces, for the build:
#   snomed_og_lookup.csv              the codes meaning C15 or C16
#   snomed_map_source.txt             which release they came from
# Produces, for a clinical colleague to check:
#   snomed_og_for_review.csv          the same codes, named, arguable ones first
#   snomed_og_dropped_for_review.csv  codes that could not be placed
#   snomed_map_survey.csv             what the release contained
#
# Copy snomed_og_lookup.csv back to the server under Data/reference. The three
# review files stay here - they are for reading, not for the build.
# =============================================================================

source("R/config/directories.R")

run_this_step <- function() {
  dir_build <- "R/02_fetch_snomed_tumour_map"
  step <- function(file) {
    cat("\n========== ", file, " ==========\n", sep = "")
    source(file.path(dir_build, file), local = new.env())
  }
  
  if (!nzchar(Sys.getenv("TRUD_API_KEY")))
    stop("TRUD_API_KEY is not set, so the release cannot be downloaded.\n",
         "  - put TRUD_API_KEY=your-key in .Renviron and restart R.\n",
         "  - the account must also be subscribed to TRUD item 101; a key on ",
         "its own is not enough.", call. = FALSE)
  
  step("01_fetch_snomed_map.R")
  
  cat("\nStage 02 complete.\n")
  cat("Copy", file.path(dir_ref, "snomed_og_lookup.csv"),
      "back to the server under Data/reference.\n")
  cat("The review files stay here for a clinical colleague to read:\n")
  cat("  snomed_og_for_review.csv, snomed_og_dropped_for_review.csv,",
      "snomed_map_survey.csv\n")
  cat("The logic checks for this stage live in R/tests/02_test_snomed_map.R.\n")
  invisible(NULL)
}
run_this_step()