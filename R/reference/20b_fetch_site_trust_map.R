# =============================================================================
# 20b  Look up the parent trust of every site code, from ODS
# -----------------------------------------------------------------------------
# Run this on a machine with internet access. It is the second half of the site
# lookup: 20a_list_site_codes.R runs on the server and writes the list of codes;
# this reads that list and asks ODS about each one. Neither the COSD file nor the
# analysis server is touched here.
#
# The COSD extract records where a patient was diagnosed as a five-character site
# code - a specific hospital, not the trust that runs it. Most of the time the
# trust is simply the first three characters of the site code, but not always:
# some sites have been through mergers, some sit under a trust whose code is not
# their own first three characters, and some codes are not hospital sites at all
# (a GP practice, or a placeholder). Guessing the trust from the first three
# characters is right most of the time and quietly wrong the rest, and there is
# no way to tell which from the code alone.
#
# The Organisation Data Service (ODS) is the NHS body that assigns these codes
# and records, for each one, what kind of organisation it is, whether it is still
# open, which trust operates it, and what it was before any merger. This script
# asks ODS about every code in the list and writes what it says to a small file
# the build reads.
#
# What it records for each site code:
#   parent_trust     the trust ODS says operates the site (not the first three
#                    characters - the actual operator)
#   trust_is_prefix  whether that parent happens to match the first three
#                    characters, so the scale of the "not always" is visible
#   status           whether the site is open or closed
#   is_hospital_site whether it is a hospital site at all, as opposed to a GP
#                    practice, a trust code, or a placeholder
#   predecessor      the code this one succeeded, where a merger is recorded
#   name             the site's name, so the file can be read by a person
#
# Before running: nothing. The ODS ORD API is open - no key, no sign-up. It asks
# that callers stay below five requests a second, which the pause below respects.
#
# Reads:  Data/reference/site_codes_to_lookup.csv   the list 20a wrote
# Writes: Data/reference/site_trust_map.csv         one row per site code
#         Data/reference/site_trust_lookup_source.txt   when it was built
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(magrittr)   # for %>%
})

if (!exists("dir_ref")) dir_ref <- "Data/reference"
dir.create(dir_ref, recursive = TRUE, showWarnings = FALSE)

f_site_codes        <- file.path(dir_ref, "site_codes_to_lookup.csv")
f_site_trust_map    <- file.path(dir_ref, "site_trust_map.csv")
f_site_trust_source <- file.path(dir_ref, "site_trust_lookup_source.txt")

ods_base <- "https://directory.spineservices.nhs.uk/ORD/2-0-0/organisations"

# the codes ODS uses that matter here. RO197 is an NHS trust, RO198 an NHS trust
# site; RC1 is an organisation record (a trust), RC2 a site record. RE6 is the
# "is operated by" relationship that points a site at its trust.
role_trust      <- "RO197"
role_trust_site <- "RO198"
rel_operated_by <- "RE6"

# -----------------------------------------------------------------------------
# the codes to look up
# -----------------------------------------------------------------------------
# These come straight from the list 20a wrote on the server - already tidied and
# already restricted to distinct five-character codes - so there is nothing to
# clean here, only to read.
read_site_codes <- function() {
  if (!file.exists(f_site_codes))
    stop("no code list at ", f_site_codes,
         "\n  - run 20a_list_site_codes.R on the server first and copy the ",
         "result here.", call. = FALSE)
  codes <- read.csv(f_site_codes, colClasses = "character")$site_code
  sort(unique(codes[!is.na(codes) & nchar(codes) == 5]))
}

# -----------------------------------------------------------------------------
# asking ODS about one code
# -----------------------------------------------------------------------------
# Returns one row. A code ODS has never heard of comes back marked not found
# rather than stopping the run - an unknown code is a finding, not an error.
ask_ods <- function(code) {
  url <- file.path(ods_base, code)
  res <- tryCatch(fromJSON(url, simplifyVector = FALSE),
                  error = function(e) NULL)
  
  if (is.null(res) || is.null(res$Organisation))
    return(derive_fields(tibble(site_code = code, found = FALSE,
                                name = NA_character_, record_class = NA_character_,
                                status = NA_character_, primary_role = NA_character_,
                                parent_trust = NA_character_, predecessor = NA_character_)))
  
  org <- res$Organisation
  
  roles <- org$Roles$Role
  primary <- NA_character_
  if (length(roles)) {
    is_primary <- map_lgl(roles, ~ isTRUE(.x$primaryRole))
    primary <- if (any(is_primary)) roles[[which(is_primary)[1]]]$id
    else roles[[1]]$id
  }
  
  # the parent trust is the active "operated by" relationship pointing at a
  # trust. Take the active one if there is a choice.
  parent <- NA_character_
  rels <- org$Rels$Rel
  if (length(rels)) {
    op <- keep(rels, ~ identical(.x$id, rel_operated_by))
    active <- keep(op, ~ identical(.x$Status, "Active"))
    pick <- if (length(active)) active[[1]] else if (length(op)) op[[1]] else NULL
    if (!is.null(pick)) parent <- pick$Target$OrgId$extension
  }
  
  # the code this site succeeded, where a merger left a predecessor
  predecessor <- NA_character_
  succs <- org$Succs$Succ
  if (length(succs)) {
    pre <- keep(succs, ~ identical(.x$Type, "Predecessor"))
    if (length(pre)) predecessor <- pre[[1]]$Target$OrgId$extension
  }
  
  derive_fields(tibble(site_code = code, found = TRUE,
                       name = org$Name %||% NA_character_,
                       record_class = org$orgRecordClass %||% NA_character_,
                       status = org$Status %||% NA_character_,
                       primary_role = primary,
                       parent_trust = parent %||% NA_character_,
                       predecessor = predecessor))
}

# turn the raw answer into the fields the build wants. Shared so that a single
# lookup and the whole batch describe a code the same way.
derive_fields <- function(row) {
  row %>%
    mutate(
      is_hospital_site = found & record_class == "RC2" &
        primary_role %in% c(role_trust_site, role_trust),
      # where ODS gives no parent - a site with no recorded operator, or a code
      # that is itself a trust - fall back to the first three characters so the
      # column is never empty, and mark that this was assumed
      parent_from_ods  = !is.na(parent_trust),
      parent_trust     = if_else(is.na(parent_trust),
                                 str_sub(site_code, 1, 3), parent_trust),
      trust_is_prefix  = parent_trust == str_sub(site_code, 1, 3),
      status           = if_else(found, status, "not in ODS"))
}

# -----------------------------------------------------------------------------
# do it
# -----------------------------------------------------------------------------
if (!exists("skip_ods_run") || !isTRUE(skip_ods_run)) {
  codes <- read_site_codes()
  cat("Distinct five-character site codes to look up:", length(codes), "\n")
  cat("At five a second this is about",
      ceiling(length(codes) / 5 / 60), "minutes.\n\n")
  
  raw <- vector("list", length(codes))
  for (i in seq_along(codes)) {
    raw[[i]] <- ask_ods(codes[i])
    if (i %% 100 == 0) cat("  looked up", i, "of", length(codes), "\n")
    Sys.sleep(0.21)   # stay under five requests a second
  }
  # ask_ods already derived the fields, so the batch only needs the wanted
  # column order
  site_trust <- bind_rows(raw) %>%
    select(site_code, name, parent_trust, trust_is_prefix, parent_from_ods,
           status, is_hospital_site, record_class, primary_role, predecessor,
           found)
  
  write.csv(site_trust, f_site_trust_map, row.names = FALSE)
  
  # -- what the lookup found, in plain terms ----------------------------------
  cat("\nOf", nrow(site_trust), "site codes looked up:\n")
  summ <- tibble(
    finding = c(
      "known to ODS as a hospital site",
      "known to ODS, but not a hospital site (GP practice, trust code, other)",
      "not found in ODS at all"),
    codes = c(
      sum(site_trust$is_hospital_site),
      sum(site_trust$found & !site_trust$is_hospital_site),
      sum(!site_trust$found)))
  print(as.data.frame(summ), row.names = FALSE)
  
  hosp <- site_trust %>% filter(is_hospital_site)
  cat("\nOf the", nrow(hosp), "hospital sites, does the parent trust match the",
      "first three characters of the code:\n")
  hosp %>%
    mutate(answer = if_else(trust_is_prefix,
                            "yes - first three characters are the trust",
                            "no - ODS puts it under a different trust")) %>%
    count(answer, name = "codes") %>%
    as.data.frame() %>%
    print(row.names = FALSE)
  
  differ <- hosp %>% filter(!trust_is_prefix)
  if (nrow(differ)) {
    cat("\nThe sites whose trust is not their first three characters -",
        "the ones a simple rule would place wrongly:\n")
    differ %>%
      mutate(prefix = str_sub(site_code, 1, 3)) %>%
      select(site_code, name, prefix, parent_trust) %>%
      head(30) %>%
      as.data.frame() %>%
      print(row.names = FALSE)
  }
  
  closed <- hosp %>% filter(status != "Active")
  cat("\nHospital sites recorded as closed:", nrow(closed),
      "- kept in the file, since a patient can have been diagnosed at a site",
      "that has since shut.\n")
  
  writeLines(c(
    "Site to parent-trust lookup, from the ODS ORD API",
    paste("Endpoint:", ods_base),
    paste("Site codes looked up:", nrow(site_trust)),
    paste("  hospital sites:", sum(site_trust$is_hospital_site)),
    paste("  not a hospital site:",
          sum(site_trust$found & !site_trust$is_hospital_site)),
    paste("  not found in ODS:", sum(!site_trust$found)),
    paste("  parent trust differs from first three characters:",
          sum(site_trust$is_hospital_site & !site_trust$trust_is_prefix)),
    paste("Built:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))),
    f_site_trust_source)
  
  cat("\nWrote", f_site_trust_map, "\n")
  cat("Copy it to the analysis server under Data/reference.\n")
}