# =============================================================================
# run_0_list_codes.R  -  list the site codes to look up
# -----------------------------------------------------------------------------
# WHERE: the analysis server, on the real COSD data.
# NEEDS: the real COSD file. Does NOT need internet.
#
# This is the first thing to run. It reads the COSD extract and writes the list
# of distinct five-character site codes to
# Data/reference/site_codes_to_lookup.txt. That list carries no patient data -
# it is just codes - so it is safe to take off the server to the internet
# machine, where run_1_internet.R turns it into the site-to-trust map.
#
# dir_raw points at the real extracts. Edit it if the path on the server differs
# from the default in the script.
# =============================================================================

# dir_raw <- "W:/_DATA/IainTimmins/2026 OG SOTN data"   # uncomment to override
dir_ref <- "Data/reference"

source("R/reference/20a_list_site_codes.R")

message("\nDone. Copy Data/reference/site_codes_to_lookup.txt to the internet ",
        "machine and run run_1_internet.R.")