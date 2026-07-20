# =============================================================================
# 21  Check the ODS parsing
# -----------------------------------------------------------------------------
# The lookup in 20 turns each ODS record into one row. The turning is the part
# that can go wrong - the records nest the parent trust and the merger history a
# few levels down - so this checks the parsing against real ODS payloads (a live
# hospital site and a trust, fetched from the API and pasted in verbatim) and
# against made-up ones for the awkward cases: a site under an unexpected trust, a
# closed site, a GP practice, a code ODS has never heard of.
#
# It does not touch the network. ask_ods() reads a URL, so the checks replace the
# reading step and feed it the stored json instead.
#
# Run from the project root:
#   Rscript R/reference/21_check_site_trust_map.R
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(magrittr)
})

skip_ods_run <- TRUE
dir_ref <- tempfile("ref_")
source("R/reference/20_fetch_site_trust_map.R")

.checks <- new.env(); .checks$rows <- list()
expect <- function(label, cond) {
  ok <- isTRUE(cond)
  .checks$rows[[length(.checks$rows) + 1]] <- list(label = label, ok = ok)
  cat(if (ok) "  pass  " else "  FAIL  ", label, "\n")
}

# feed ask_ods() stored json instead of a live fetch, by pointing fromJSON at a
# file. ask_ods builds the url itself, so the seam is fromJSON: override it in
# the function's environment for the length of one call.
ask_from_json <- function(code, json) {
  f <- tempfile(fileext = ".json"); writeLines(json, f)
  environment(ask_ods)$fromJSON <- function(url, ...) jsonlite::fromJSON(f, ...)
  on.exit(rm("fromJSON", envir = environment(ask_ods)), add = TRUE)
  ask_ods(code)
}

# -----------------------------------------------------------------------------
# real payloads, fetched from the ODS ORD API and pasted verbatim
# -----------------------------------------------------------------------------
# a live hospital site: Kent & Canterbury, operated by trust RVV, with a
# predecessor code from a past merger
json_site <- '{"Organisation":{"Name":"KENT & CANTERBURY HOSPITAL","Date":[{"Type":"Operational","Start":"1999-04-01"}],"OrgId":{"root":"2.16.840.1.113883.2.1.3.2.4.18.48","assigningAuthorityName":"HSCIC","extension":"RVVKC"},"Status":"Active","LastChangeDate":"2021-07-21","orgRecordClass":"RC2","GeoLoc":{"Location":{"AddrLn1":"ETHELBERT ROAD","Town":"CANTERBURY","County":"KENT","PostCode":"CT1 3NG","Country":"ENGLAND","UPRN":100062280103}},"Roles":{"Role":[{"id":"RO198","uniqueRoleId":13437,"primaryRole":true,"Date":[{"Type":"Operational","Start":"1999-04-01"}],"Status":"Active"}]},"Rels":{"Rel":[{"Date":[{"Type":"Operational","Start":"1999-04-01"}],"Status":"Active","Target":{"OrgId":{"root":"2.16.840.1.113883.2.1.3.2.4.18.48","assigningAuthorityName":"HSCIC","extension":"RVV"},"PrimaryRoleId":{"id":"RO197","uniqueRoleId":41756}},"id":"RE6","uniqueRelId":15806}]},"Succs":{"Succ":[{"uniqueSuccId":28060,"Date":[{"Type":"Legal","Start":"1999-04-01"}],"Type":"Predecessor","Target":{"OrgId":{"root":"2.16.840.1.113883.2.1.3.2.4.18.48","assigningAuthorityName":"HSCIC","extension":"RGWKC"},"PrimaryRoleId":{"id":"RO198","uniqueRoleId":3501}}}]}}}'

# a trust record: an organisation, not a site, and inactive
json_trust <- '{"Organisation":{"Name":"WIGAN AND LEIGH HEALTH SERVICES NHS TRUST","Date":[{"Type":"Operational","Start":"1993-04-01","End":"2001-03-31"}],"OrgId":{"root":"2.16.840.1.113883.2.1.3.2.4.18.48","assigningAuthorityName":"HSCIC","extension":"RJY"},"Status":"Inactive","LastChangeDate":"2020-03-19","orgRecordClass":"RC1","refOnly":true,"GeoLoc":{"Location":{"AddrLn1":"ROYAL ALBERT EDWARD INFIRMARY","AddrLn2":"WIGAN LANE","Town":"WIGAN","County":"LANCASHIRE","PostCode":"WN1 2NN","Country":"ENGLAND"}},"Roles":{"Role":[{"id":"RO197","uniqueRoleId":37515,"primaryRole":true,"Date":[{"Type":"Operational","Start":"1993-04-01","End":"2001-03-31"}],"Status":"Inactive"}]}}}'

cat("\nReal ODS payloads\n")

site <- ask_from_json("RVVKC", json_site)
expect("a hospital site is recognised as one", site$is_hospital_site)
expect("its parent trust is read from the relationship, not the code",
       site$parent_trust == "RVV")
expect("here the trust does match the first three characters",
       site$trust_is_prefix)
expect("the parent came from ODS, not a fallback", site$parent_from_ods)
expect("the site is active", site$status == "Active")
expect("the predecessor code is read from the merger history",
       site$predecessor == "RGWKC")
expect("the name is kept", site$name == "KENT & CANTERBURY HOSPITAL")

trust <- ask_from_json("RJY", json_trust)
expect("a trust record is not counted as a hospital site",
       !trust$is_hospital_site)
expect("a trust record is still found", trust$found)

# -----------------------------------------------------------------------------
# made-up payloads for the awkward cases
# -----------------------------------------------------------------------------
cat("\nThe awkward cases\n")

# a site whose operating trust is NOT its first three characters - the case the
# whole lookup exists for
json_moved <- '{"Organisation":{"Name":"SOME TREATMENT CENTRE","OrgId":{"extension":"NT301"},"Status":"Active","orgRecordClass":"RC2","Roles":{"Role":[{"id":"RO198","primaryRole":true,"Status":"Active"}]},"Rels":{"Rel":[{"Status":"Active","Target":{"OrgId":{"extension":"RXK"},"PrimaryRoleId":{"id":"RO197"}},"id":"RE6"}]}}}'
moved <- ask_from_json("NT301", json_moved)
expect("a site under an unexpected trust is picked up correctly",
       moved$parent_trust == "RXK")
expect("and is flagged as not matching its first three characters",
       !moved$trust_is_prefix)

# a closed site
json_closed <- '{"Organisation":{"Name":"OLD HOSPITAL","OrgId":{"extension":"RAA01"},"Status":"Inactive","orgRecordClass":"RC2","Roles":{"Role":[{"id":"RO198","primaryRole":true,"Status":"Inactive"}]},"Rels":{"Rel":[{"Status":"Inactive","Target":{"OrgId":{"extension":"RAA"},"PrimaryRoleId":{"id":"RO197"}},"id":"RE6"}]}}}'
closed <- ask_from_json("RAA01", json_closed)
expect("a closed site is still a hospital site", closed$is_hospital_site)
expect("its closed status is recorded", closed$status == "Inactive")
expect("an inactive relationship is still read for the parent",
       closed$parent_trust == "RAA")

# a GP practice - a real code, but not a hospital site
json_gp <- '{"Organisation":{"Name":"A GP SURGERY","OrgId":{"extension":"A81001"},"Status":"Active","orgRecordClass":"RC2","Roles":{"Role":[{"id":"RO177","primaryRole":true,"Status":"Active"}]}}}'
gp <- ask_from_json("A81001", json_gp)
expect("a GP practice is found but not a hospital site",
       gp$found && !gp$is_hospital_site)
expect("a GP practice with no trust relationship falls back to its prefix",
       gp$parent_trust == "A81" && !gp$parent_from_ods)

# a code ODS returns nothing for
json_empty <- '{}'
missing <- ask_from_json("ZZ999", json_empty)
expect("a code ODS does not know is marked not found", !missing$found)
expect("an unknown code is not a hospital site", !missing$is_hospital_site)
expect("an unknown code still gets a prefix fallback for its trust",
       missing$parent_trust == "ZZ9")

# two operated-by relationships, one active and one not: take the active one
json_two <- '{"Organisation":{"Name":"MOVED SITE","OrgId":{"extension":"RBB02"},"Status":"Active","orgRecordClass":"RC2","Roles":{"Role":[{"id":"RO198","primaryRole":true,"Status":"Active"}]},"Rels":{"Rel":[{"Status":"Inactive","Target":{"OrgId":{"extension":"RCC"},"PrimaryRoleId":{"id":"RO197"}},"id":"RE6"},{"Status":"Active","Target":{"OrgId":{"extension":"RBB"},"PrimaryRoleId":{"id":"RO197"}},"id":"RE6"}]}}}'
two <- ask_from_json("RBB02", json_two)
expect("with an old and a current trust, the current one is taken",
       two$parent_trust == "RBB")

# =============================================================================
res <- bind_rows(lapply(.checks$rows, as_tibble))
n_fail <- sum(!res$ok)
cat("\n", nrow(res), "checks,", n_fail, "failed\n")
if (n_fail) {
  cat("\nfailed:\n"); cat(paste0("  ", res$label[!res$ok], collapse = "\n"), "\n")
  quit(status = 1, save = "no")
}
cat("All checks passed. The live download is the only part not covered.\n")