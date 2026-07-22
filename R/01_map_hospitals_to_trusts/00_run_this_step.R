# =============================================================================
# 01  Map the hospitals to their trusts  -  run this step
# -----------------------------------------------------------------------------
# WHERE: a machine with internet access. NOT the analysis server.
#
# Nothing here reads patient data. The only input is the code list stage 00
# wrote, which holds organisation codes and nothing else - that is what makes it
# safe to carry off the server. Keep it that way: if a script in this folder
# ever needs a patient-level file, it belongs in a server stage instead.
#
# Asks the NHS ODS API about every five-character site code in the list - both
# the COSD sites of diagnosis and the HES sites of treatment - and records, for
# each, the trust that operates it, whether that trust matches the code's own
# first three characters, the site's status, and what kind of organisation ODS
# says it is. The API is open: no key, no account. It asks callers to stay under
# five requests a second, which the fetch respects, so roughly 1,700 codes takes
# about six minutes.
#
# Needs:    site_codes_to_lookup.txt in dir_ref (copied from the server)
# Produces: site_trust_map.csv, site_trust_lookup_source.txt
#           Copy site_trust_map.csv back to the server under Data/reference.
# =============================================================================

source("R/config/directories.R")

run_this_step <- function() {
  dir_build <- "R/01_map_hospitals_to_trusts"
  step <- function(file) {
    cat("\n========== ", file, " ==========\n", sep = "")
    source(file.path(dir_build, file), local = new.env())
  }
  
  step("01_fetch_site_trust_map.R")   # ODS lookup, one row per site code
  
  cat("\nStage 01 complete.\n")
  cat("Copy", file.path(dir_ref, "site_trust_map.csv"),
      "back to the server under Data/reference.\n")
  cat("The logic checks for this stage live in R/tests/01_test_site_trust_map.R.\n")
  invisible(NULL)
}
run_this_step()