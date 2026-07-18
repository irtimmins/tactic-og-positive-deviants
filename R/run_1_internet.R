# =============================================================================
# run_1_internet.R  -  build the two reference maps
# -----------------------------------------------------------------------------
# WHERE: a machine with internet access. Not the analysis server.
# NEEDS: internet. Does NOT need the COSD data or any patient data.
#
# This builds the two small reference files the airgap build later reads:
#   - the SNOMED-to-C15/C16 map, from NHS TRUD
#   - the site-to-trust map, from the NHS ODS API
#
# The site map needs the list of site codes that run_0_list_codes.R produced on
# the server. Copy that list (site_codes_to_lookup.txt) into Data/reference here
# before running the second half.
#
# Run the whole file, or a block at a time. The check scripts need no internet
# and can be run again anywhere.
#
# Afterwards, copy these to the server under Data/reference:
#   snomed_og_lookup.csv
#   site_trust_map.csv
# =============================================================================

dir_ref <- "Data/reference"

# -- the SNOMED map, from TRUD ------------------------------------------------
# Needs a TRUD account and the API key in TRUD_API_KEY. See the notes at the top
# of 10_fetch_snomed_map.R.
source("R/reference/10_fetch_snomed_map.R")
source("R/reference/11_check_snomed_map.R")   # proves 10, no internet needed

# -- the site-to-trust map, from ODS -----------------------------------------
# Reads Data/reference/site_codes_to_lookup.txt (from run_0 on the server) and
# asks ODS about each code. Open API, no key.
source("R/reference/20b_fetch_site_trust_map.R")
source("R/reference/21_check_site_trust_map.R")   # proves 20b, no internet needed

message("\nDone. Copy snomed_og_lookup.csv and site_trust_map.csv to the ",
        "server under Data/reference.")