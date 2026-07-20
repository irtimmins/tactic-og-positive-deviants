# =============================================================================
# 11  Check the TRUD map handling
# -----------------------------------------------------------------------------
# Builds a stand-in TRUD release on disk and puts it through the reading and
# filtering in 10, so the fiddly parts - finding the file inside the zip, coping
# with a zip inside a zip, Full versus Snapshot, codes with two ICD-10 targets -
# are tested without needing the network or a TRUD account.
#
# The stand-in release carries the same twenty-odd map files the real one does,
# and the chosen file carries several maps stacked in it, because both of those
# have already caught this code out once:
#   - the first version read the International map instead of the UK one, since
#     the fake release only had one map file and there was no choice to get wrong
#   - the second read every refset in the right file, including the OPCS-4
#     procedure map, whose C15 and C16 are eye operations, and counted them as
#     oesophageal and gastric cancer
#
# The download itself is the only part not covered here. Everything after it is.
#
# Run from the project root:
#   Rscript R/fetch_reference_data/11_check_snomed_map.R
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))

# load the functions from 10 without letting it try to reach TRUD. local = TRUE
# sources 10 into THIS environment, so the skip_trud_run set just above is the
# one its download guard sees. Without it, a plain source() would read 10 in the
# global environment instead - and when this check is itself sourced into a child
# environment (as run_test.R does), skip_trud_run would not be visible there and
# the real 1 GB TRUD download would run. A test must never reach the network.
skip_trud_run <- TRUE
dir_ref <- tempfile("ref_")
source("R/fetch_reference_data/10_fetch_snomed_map.R", local = TRUE)

.checks <- new.env(); .checks$rows <- list()
expect <- function(label, cond) {
  ok <- isTRUE(cond)
  .checks$rows[[length(.checks$rows) + 1]] <- list(label = label, ok = ok)
  cat(if (ok) "  pass  " else "  FAIL  ", label, "\n")
}

# -----------------------------------------------------------------------------
# A stand-in map file
# -----------------------------------------------------------------------------
# The columns are the ones the real extended map refset uses. The rows cover:
#   a plain oesophageal code, a plain gastric code
#   a bowel code, which must not come through
#   a code retired in a later row, which must not come through
#   a code that was retired and then brought back, which must
#   a code with two targets that agree, which must come through
#   a code with two targets that disagree, which must be left out
#   a code with a blank target, which says nothing
map_cols <- c("id", "effectiveTime", "active", "moduleId", "refsetId",
              "referencedComponentId", "mapGroup", "mapPriority", "mapRule",
              "mapAdvice", "mapTarget", "correlationId", "mapCategoryId")

icd10_uk <- "999002271000000101"    # the map we want
opcs4    <- "1126441000000105"      # OPCS-4.9, sitting in the same file

# names, as the release supplies them. The International file and the UK file
# both have to be read: the UK one does not repeat the International content.
desc_row <- function(concept, term, type = "900000000000003001") {
  tibble(id = paste0("d", concept), effectiveTime = "20240101", active = "1",
         moduleId = "900000000000207008", conceptId = concept,
         languageCode = "en", typeId = type, term = term,
         caseSignificanceId = "900000000000448009")
}

fake_descriptions <- bind_rows(
  desc_row("111000001", "Adenocarcinoma of oesophagus (disorder)"),
  desc_row("111000002", "Malignant tumour of stomach (disorder)"),
  # one that should catch a clinician's eye
  desc_row("555000001", "Adenocarcinoma of gastro-oesophageal junction (disorder)"),
  # named without the word junction, which is how the UK concepts read and how
  # the flag came to miss them
  desc_row("444000001", "Siewert type I adenocarcinoma (disorder)"),
  # named to sort alphabetically between the two junction terms above, so an
  # alphabetical ordering would put a benign overlapping-site row between two
  # rows that need a real decision
  desc_row("777000002", "Overlapping malignant neoplasm of esophagus (disorder)"),
  desc_row("666000001", "Malignant neoplasm of overlapping site (disorder)"),
  # a synonym, which must not be used in place of the fully specified name
  desc_row("111000001", "Oesophageal adenocarcinoma", "900000000000013009"),
  desc_row(icd10_uk, "ICD-10 5th edition complex map reference set (foundation metadata concept)"),
  desc_row(opcs4, "OPCS-4.9 complex map reference set (foundation metadata concept)"))

# the UK file carries only UK content, as in the real release
fake_descriptions_uk <- desc_row("901000001",
                                 "Operation on eyelid (procedure)")

row_of <- function(id, eff, active, concept, target, group = 1, priority = 1,
                   refset = icd10_uk) {
  tibble(id = id, effectiveTime = eff, active = active,
         moduleId = "999000031000000106", refsetId = refset,
         referencedComponentId = concept, mapGroup = as.character(group),
         mapPriority = as.character(priority), mapRule = "TRUE",
         mapAdvice = "ALWAYS C15.9", mapTarget = target,
         correlationId = "447561005", mapCategoryId = "447637006")
}

fake_map <- bind_rows(
  row_of("m01", "20240101", "1", "111000001", "C15.9"),   # oesophagus
  row_of("m02", "20240101", "1", "111000002", "C16.9"),   # stomach
  row_of("m03", "20240101", "1", "222000001", "C18.7"),   # bowel
  row_of("m04", "20200101", "1", "333000001", "C15.5"),   # retired later
  row_of("m04", "20240101", "0", "333000001", "C15.5"),
  row_of("m05", "20200101", "0", "444000001", "C15.3"),   # brought back later
  row_of("m05", "20240101", "1", "444000001", "C15.3"),
  row_of("m06", "20240101", "1", "555000001", "C15.4"),   # two targets, agree
  row_of("m07", "20240101", "1", "555000001", "C15.9", group = 2),
  # two alternatives within group 1 that disagree about the site: genuinely
  # ambiguous, so set aside
  row_of("m08", "20240101", "1", "666000001", "C15.9"),
  row_of("m09", "20240101", "1", "666000001", "C18.9", priority = 2),
  # group 1 is the diagnosis, group 2 adds the infectious agent. This is a
  # gastric cancer with a virus code attached, not an ambiguous one, and must be
  # kept as C16.
  row_of("m15", "20240101", "1", "904000001", "C16.9"),
  row_of("m16", "20240101", "1", "904000001", "B97.8", group = 2),
  # an overlapping-site code, to sit alongside the junction one in the ordering
  # check below
  row_of("m17", "20240101", "1", "777000002", "C15.8"),
  row_of("m10", "20240101", "1", "777000001", ""),        # nothing to say
  # a marker, not a code
  row_of("m11", "20240101", "1", "888000001", "#NC"),
  # the OPCS-4 map, in the same file. C15 and C16 here are eye operations, and
  # must not be read as oesophageal or gastric cancer.
  row_of("m12", "20240101", "1", "901000001", "C151", refset = opcs4),
  row_of("m13", "20240101", "1", "902000001", "C162", refset = opcs4),
  row_of("m14", "20240101", "1", "903000001", "G451", refset = opcs4))

stopifnot(identical(names(fake_map), map_cols))

# -----------------------------------------------------------------------------
# A stand-in release: the map inside a zip inside a zip, as some releases have
# -----------------------------------------------------------------------------
# The real release carries every combination of source, release type and refset
# shape. Only one of them is the UK Clinical ICD-10 map; the rest are decoys,
# and the code has to pick out the right one by name rather than take whatever
# it stumbles on first.
map_file_names <- function() {
  out <- character()
  for (shape in c("iissscc", "iisssci"))
    for (who in c("UKCL", "UKCR", "UKED"))
      for (type in c("Snapshot", "Full", "Delta"))
        out <- c(out, sprintf("der2_%sRefset_ExtendedMap%s%s_GB1000000_20240101.txt",
                              shape, who, type))
  c(out, "der2_iisssccRefset_ExtendedMapSnapshot_INT_20240101.txt",
    "der2_iisssccRefset_ExtendedMapFull_INT_20240101.txt")
}

# the file we want, and everything else holding something different, so a wrong
# pick shows up as a wrong answer rather than passing quietly
build_release <- function(map, nested = TRUE, decoy_map = NULL,
                          only_file = NULL) {
  work <- tempfile("trud_fake_"); dir.create(work)
  old <- setwd(work); on.exit(setwd(old), add = TRUE)
  
  root <- "SnomedCT_UKClinicalRF2_PRODUCTION_20240101"
  wanted <- "der2_iisssciRefset_ExtendedMapUKCLSnapshot_GB1000000_20240101.txt"
  if (is.null(decoy_map)) decoy_map <- map %>% mutate(mapTarget = "C34.9")
  
  files <- if (is.null(only_file)) map_file_names() else only_file
  for (nm in files) {
    type <- if (grepl("Snapshot", nm)) "Snapshot" else
      if (grepl("Delta", nm)) "Delta" else "Full"
    d <- file.path(root, type, "Refset", "Map")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    content <- if (identical(nm, wanted) || !is.null(only_file)) map else decoy_map
    write.table(content, file.path(d, nm), sep = "\t", quote = FALSE,
                row.names = FALSE, na = "")
  }
  
  # the names live in the release too, in their own files
  td <- file.path(root, "Snapshot", "Terminology")
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  write.table(fake_descriptions, file.path(td,
                                           "sct2_Description_Snapshot-en_INT_20240101.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  write.table(fake_descriptions_uk, file.path(td,
                                              "sct2_Description_UKCLSnapshot-en_GB1000000_20240101.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  # a text definition file. It is shaped like a description file and would be
  # read as one if the pattern were loose, so it carries a row that must not
  # turn up in the names.
  write.table(desc_row("decoy", "Long free-text definition (disorder)"),
              file.path(td, "sct2_TextDefinition_Snapshot-en_INT_20240101.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  
  dir.create("Documentation", showWarnings = FALSE)
  writeLines("release notes", "Documentation/readme.txt")
  
  if (nested) {
    zip("SnomedCT_UKClinical.zip", root, flags = "-qr")
    unlink(root, recursive = TRUE)
    zip("uk_sct2cl_42.1.0_20240101000001Z.zip",
        c("SnomedCT_UKClinical.zip", "Documentation"), flags = "-qr")
  } else {
    zip("uk_sct2cl_42.1.0_20240101000001Z.zip", c(root, "Documentation"),
        flags = "-qr")
  }
  file.path(work, "uk_sct2cl_42.1.0_20240101000001Z.zip")
}

# =============================================================================
cat("\nReading the map out of a release\n")

n_icd10 <- sum(fake_map$refsetId == icd10_uk)

flat_zip <- build_release(fake_map, nested = FALSE)
map_flat <- read_extended_map(flat_zip)
expect("the map is found in a plain zip, and only the ICD-10 refset comes back",
       nrow(map_flat) == n_icd10 && all(map_flat$refsetId == icd10_uk))

nested_zip <- build_release(fake_map, nested = TRUE)
map_nested <- read_extended_map(nested_zip)
expect("the map is found in a zip inside a zip",
       nrow(map_nested) == n_icd10)
expect("both routes give the same map",
       identical(map_flat$id, map_nested$id))

# the bug this file exists to stop coming back
expect("the UK Clinical ICD-10 map is chosen, not one of the twenty others",
       all(map_flat$mapTarget[map_flat$mapTarget != ""] != "C34.9"))
expect("the International map is not what gets read",
       !any(grepl("C34", map_nested$mapTarget)))
expect("a renamed map file stops the run rather than reading the wrong one",
       {
         msg <- tryCatch(read_extended_map(flat_zip, want = "NotAFile"),
                         error = conditionMessage)
         grepl("expected exactly one map file", msg) &&
           grepl("ExtendedMapUKCLSnapshot", msg)
       })

expect("a release with no map file gives a clear error",
       {
         empty <- tempfile(fileext = ".zip")
         d <- tempfile(); dir.create(d)
         writeLines("nothing here", file.path(d, "readme.txt"))
         # local() gives on.exit a real function frame, so the working directory
         # is put back even if zip() errors (e.g. no zip utility on PATH) rather
         # than being left changed for every later check to trip over.
         local({
           old <- setwd(d); on.exit(setwd(old))
           zip(empty, "readme.txt", flags = "-q")
         })
         msg <- tryCatch(read_extended_map(empty), error = conditionMessage)
         grepl("no ExtendedMap file", msg)
       })

expect("a map missing a column gives a clear error",
       {
         bad <- build_release(fake_map %>% select(-mapTarget))
         msg <- tryCatch(read_extended_map(bad), error = conditionMessage)
         grepl("no mapTarget column", msg)
       })

# =============================================================================
cat("\nSurveying what is in a release\n")

sv <- survey_maps(flat_zip)
expect("the survey reports one line per file and refset, not per file",
       nrow(sv) > 7 && all(c("file", "refset") %in% names(sv)))
expect("the survey shows the OPCS-4 map sitting in the same file as the ICD-10 one",
       {
         same <- sv %>% filter(grepl("iisssciRefset_ExtendedMapUKCLSnapshot", file))
         all(c(icd10_uk, opcs4) %in% same$refset)
       })
expect("the survey shows OPCS-4 contributing bogus C15/C16 codes",
       sv$og_codes_rough[sv$refset == opcs4 &
                           grepl("iisssciRefset_ExtendedMapUKCLSnapshot", sv$file)][1] == 2)
# the survey is a rough count, meant only for telling the files apart: it does
# not do the newest-row-per-entry reduction, so a code retired in a later row
# still shows up here. og_codes_from_map() is what gives the real answer, and
# gets 4 where the survey says 6.
expect("the survey counts the OG codes per map, roughly",
       sv$og_codes_rough[sv$refset == icd10_uk & grepl(
         "iisssciRefset_ExtendedMapUKCLSnapshot", sv$file)][1] == 8)
expect("the survey shows the decoy files have no OG codes",
       all(sv$og_codes_rough[!grepl("iisssciRefset_ExtendedMapUKCLSnapshot",
                                    sv$file)] == 0))
# the report at the bottom of 10 is skipped by these checks (skip_trud_run), so
# run it here on purpose: it once named a column that had been renamed, and
# nothing found out until it hit the real release
expect("the survey report prints without error",
       {
         out <- capture.output(report_survey(sv))
         length(out) > 1
       })
expect("without the names to hand, the report says so rather than inventing",
       {
         out <- paste(capture.output(report_survey(sv)), collapse = " ")
         grepl("not named", out)
       })
expect("a mis-read file is caught rather than producing odd counts",
       {
         shifted <- fake_map %>% mutate(active = "999002271000000101")
         msg <- tryCatch(check_map_parse(shifted, "x.txt"), error = conditionMessage)
         grepl("looks mis-read", msg) && grepl("active column", msg)
       })
expect("a shifted identifier column is caught too",
       {
         shifted <- fake_map %>% mutate(referencedComponentId = "C15.9")
         msg <- tryCatch(check_map_parse(shifted, "x.txt"), error = conditionMessage)
         grepl("looks mis-read", msg) && grepl("referencedComponentId", msg)
       })
expect("a well-formed map passes the parse check",
       isTRUE(check_map_parse(fake_map, "x.txt")))
expect("the report says when a file is simply empty",
       {
         out <- paste(capture.output(report_survey(
           sv %>% mutate(refset = NA_character_))), collapse = " ")
         grepl("file is empty", out)
       })
expect("the report flags a refset the release does not name",
       {
         out <- paste(capture.output(report_survey(
           sv %>% mutate(refset = "999999"), read_descriptions(flat_zip))),
           collapse = " ")
         grepl("not named in this release", out)
       })
expect("every column the report asks for is in the survey",
       all(c("file", "refset", "live_codes", "og_codes_rough",
             "example_targets") %in% names(sv)))

expect("naming a refset that is not there stops the run",
       {
         msg <- tryCatch(read_extended_map(flat_zip, refset = "123"),
                         error = conditionMessage)
         grepl("is not in", msg) && grepl(opcs4, msg)
       })
expect("the survey shows example targets so the file can be recognised",
       any(grepl("C15", sv$example_targets)))

# =============================================================================
cat("\nCutting it down to the OG codes\n")

res <- og_codes_from_map(map_flat)
og <- res$kept

expect("an oesophageal code is kept",
       "111000001" %in% og$snomed &&
         og$site3[og$snomed == "111000001"] == "C15")
expect("a gastric code is kept",
       "111000002" %in% og$snomed &&
         og$site3[og$snomed == "111000002"] == "C16")
expect("a bowel code is left out", !"222000001" %in% og$snomed)
expect("a code retired in a later row is left out", !"333000001" %in% og$snomed)
expect("a code retired and then brought back is kept", "444000001" %in% og$snomed)
expect("a code with two targets that agree is kept", "555000001" %in% og$snomed)
expect("a second map group carrying an infectious agent does not lose the code",
       "904000001" %in% og$snomed &&
         og$site3[og$snomed == "904000001"] == "C16")
expect("a code with two targets that disagree is left out",
       !"666000001" %in% og$snomed)
expect("a code with no target is left out", !"777000001" %in% og$snomed)
expect("a #NC marker is not read as a code", !"888000001" %in% og$snomed)

# the bug that produced 186 codes instead of 147 on the real release
expect("an OPCS-4 eye operation at C15 is not read as oesophageal cancer",
       !"901000001" %in% og$snomed)
expect("an OPCS-4 eye operation at C16 is not read as gastric cancer",
       !"902000001" %in% og$snomed)
expect("nothing else sneaks in", nrow(og) == 6)
expect("the sites are only ever C15 or C16", all(og$site3 %in% c("C15", "C16")))

# a Full file, where the same entry appears many times, must give the same
# answer as a Snapshot
map_full <- read_extended_map(
  build_release(fake_map, only_file =
                  "der2_iisssciRefset_ExtendedMapUKCLSnapshot_GB1000000_20240101.txt"))
og_full <- og_codes_from_map(map_full)$kept
expect("a Full-style file, with a row per version, gives the same answer",
       identical(og$snomed, og_full$snomed))

# and the shape of what we hand on to the build
expect("the written file has the columns 02 expects",
       identical(names(og), c("snomed", "site3", "icd10_targets")))
expect("one row per SNOMED code", !anyDuplicated(og$snomed))

# =============================================================================
cat("\nNaming the codes from the release\n")

terms <- read_descriptions(flat_zip)
expect("names are read from the release",
       nrow(terms) > 0 && "term" %in% names(terms))
expect("the International and the UK description files are both read",
       all(c("111000001", "901000001") %in% terms$snomed))
expect("a text definition file is not mistaken for a description file",
       nrow(terms) == 9 && !"decoy" %in% terms$snomed)
expect("the fully specified name is used, not a synonym",
       terms$term[terms$snomed == "111000001"] == "Adenocarcinoma of oesophagus")
expect("the semantic tag is pulled out of the name",
       terms$semantic_tag[terms$snomed == "111000001"] == "disorder")
expect("a code with no name says so rather than going blank",
       name_of("999999999", terms) == "not named in this release")

expect("the refsets are named from the release, not from a hard-coded list",
       {
         out <- paste(capture.output(report_survey(sv, terms)), collapse = " ")
         grepl("ICD-10 5th edition complex map", out) &&
           grepl("OPCS-4.9 complex map", out)
       })

cat("\nFlagging what a clinician should look at\n")

review <- og %>%
  mutate(term = name_of(snomed, terms)) %>%
  flag_for_review()
expect("a junction code is flagged for a decision",
       review$look_at[review$snomed == "555000001"] == "junction - C15 or C16?")
expect("a plain oesophageal adenocarcinoma is not flagged",
       review$look_at[review$snomed == "111000001"] == "")
expect("a Siewert code with no junction in its name is still flagged",
       review$look_at[review$snomed == "444000001"] == "junction - C15 or C16?")
expect("the flag patterns are lowercase, since the terms are lowercased first",
       identical(review_flags$pattern, str_to_lower(review_flags$pattern)))
expect("the flags do not drop any codes",
       nrow(review) == nrow(og))

expect("codes we could not place are written out, not silently lost",
       nrow(res$dropped) == 1 && res$dropped$snomed == "666000001")
expect("a set-aside code carries a reason",
       grepl("disagree", res$dropped$reason))


# =============================================================================
cat("\nRefusing to run without a key\n")
expect("no API key gives a clear error",
       {
         msg <- withr_env <- local({
           old <- Sys.getenv("TRUD_API_KEY")
           Sys.setenv(TRUD_API_KEY = "")
           on.exit(Sys.setenv(TRUD_API_KEY = old))
           tryCatch(trud_latest_release(), error = conditionMessage)
         })
         grepl("no TRUD API key", msg)
       })

# =============================================================================
res <- bind_rows(lapply(.checks$rows, as_tibble))
n_fail <- sum(!res$ok)
cat("\n", nrow(res), "checks,", n_fail, "failed\n")
if (n_fail) {
  cat("\nfailed:\n"); cat(paste0("  ", res$label[!res$ok], collapse = "\n"), "\n")
  quit(status = 1, save = "no")
}
cat("All checks passed. The download itself is the only part not covered.\n")