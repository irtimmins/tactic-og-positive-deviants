# =============================================================================
# 01  Fetch the SNOMED to ICD-10 map from TRUD
# -----------------------------------------------------------------------------
# Run this on a machine with internet access. The analysis server never runs it,
# and nothing here reads patient data - the release comes from TRUD and the
# output is a code lookup, so it is safe to run entirely off the server.
#
# NHS England publishes a map from SNOMED codes to ICD-10 codes as part of the
# SNOMED CT UK Clinical Edition. This script downloads that release, pulls out
# the map, keeps the codes that mean oesophageal (C15) or gastric (C16) cancer,
# and writes a small file - a few hundred rows - to carry across to the server.
# The COSD site-of-diagnosis build (stage 03b) reads that file. It replaces
# having to work the meanings out from the data.
#
# Before running:
#   1. Create a TRUD account at https://isd.digital.nhs.uk/trud
#   2. Subscribe to item 101, "SNOMED CT UK Clinical Edition, RF2: Full,
#      Snapshot & Delta", and accept its licence. An API key alone is not
#      enough - the account has to be subscribed to the item.
#      https://isd.digital.nhs.uk/trud/users/guest/filters/0/categories/26/items/101/releases
#   3. Put the API key from the account page in a .Renviron file as
#      TRUD_API_KEY=your-key, and restart R.
#
# The download is about 1 GB and only one file inside it is used, so this takes
# a while. Run it on a weekday between 8am and 6pm, or between midnight and 6am,
# which is what TRUD asks for to avoid their maintenance windows.
#
# Writes, for the build:
#   Data/reference/snomed_og_lookup.csv           the codes meaning C15 or C16
#   Data/reference/snomed_map_source.txt          which release they came from
#
# Writes, for a clinical colleague to check:
#   Data/reference/snomed_og_for_review.csv       the same codes, named, with the
#                                                 arguable ones flagged first
#   Data/reference/snomed_og_dropped_for_review.csv  codes we could not place
#   Data/reference/snomed_map_survey.csv          what was in the release
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(tidyverse)
})

if (!exists("dir_ref")) dir_ref <- "Data/reference"
dir.create(dir_ref, recursive = TRUE, showWarnings = FALSE)

trud_item <- 101L
og_sites  <- c("C15", "C16")

# -----------------------------------------------------------------------------
# Which map to use
# -----------------------------------------------------------------------------
# The release holds about twenty files with ExtendedMap in the name. They vary
# in three ways, and every combination is present:
#
#   INT / UKCL / UKCR / UKED   whose map it is - the International edition, the
#                              UK Clinical extension, the UK clinical refsets,
#                              or the merged UK edition
#   Snapshot / Full / Delta    Snapshot is the map as it stands now. Full is
#                              every version of every row ever published. Delta
#                              is only what changed in this release.
#   iissscc / iisssci          two different refset shapes. NHS England's own
#                              release notes name der2_iisssciRefset_
#                              ExtendedMapUKCL... as the SNOMED CT to ICD-10
#                              mapping file, so that is the one we want.
#
# We want the UK Clinical ICD-10 map as a Snapshot.
#
# Picking the file is not enough. One file holds several maps stacked on top of
# each other, told apart only by refsetId, and they do not all map to ICD-10:
#
#   999002271000000101  ICD-10 5th ed, 5-character, UK    <- the one we want
#   1891651000000103    OPCS-4.11, the surgical procedure classification. Live:
#                       72,880 codes, 59 of which carry a C15 or C16 target.
#   1126441000000105    OPCS-4.9, the superseded version. No live rows.
#   1382401000000109    no live rows either.
#   447562003           ICD-10, International edition (in the INT file)
#
# The OPCS-4 map is the one that bites. OPCS-4 chapter C is operations on the
# eye, so OPCS-4 has its own C15 and C16 codes meaning eyelid and orbit
# procedures. Read the whole file and keep anything starting C15 or C16 and
# those come through as oesophageal and gastric cancer: 186 codes, against 153
# from the ICD-10 map alone.
#
# The names above are not hard-coded here - report_survey() reads them out of the
# release, because a refset id is itself a concept. That is how OPCS-4.11 was
# identified after two wrong guesses at what it might be.
snomed_map_file   <- "der2_iisssciRefset_ExtendedMapUKCLSnapshot_GB"
snomed_map_refset <- "999002271000000101"

f_snomed_lookup <- file.path(dir_ref, "snomed_og_lookup.csv")
f_snomed_survey <- file.path(dir_ref, "snomed_map_survey.csv")
f_snomed_source <- file.path(dir_ref, "snomed_map_source.txt")

# for clinical colleagues, not for the build
f_review_kept    <- file.path(dir_ref, "snomed_og_for_review.csv")
f_review_dropped <- file.path(dir_ref, "snomed_og_dropped_for_review.csv")

# -----------------------------------------------------------------------------
# Things worth a clinician's eye
# -----------------------------------------------------------------------------
# The map decides which codes are C15 and which are C16. That is a coding
# decision, not a clinical one, and it will not always be the decision this audit
# wants - most obviously at the gastro-oesophageal junction, where the cohort
# ends up C15-only and the map's choice moves numbers.
#
# These patterns pick out the rows worth arguing about, so a colleague reads
# twenty lines rather than a hundred and fifty. Everything is kept either way;
# the flag only says "look here first".
# The terms are lowercased before matching, so the patterns must be lowercase
# too. They were not at first, and "Barrett" silently matched nothing.
#
# "siewert" earns its place on its own rather than riding on "junction". The UK
# concepts are named plainly - "Siewert type I adenocarcinoma", no junction, no
# cardia - so a junction-only pattern leaves them unflagged at the bottom of the
# file, which is precisely where the disagreement hides: the International
# concept for Siewert I maps to C16 and the UK one to C15.
review_flags <- tribble(
  ~pattern,                                        ~note,
  "junction|cardia|siewert|gastro.?oesophageal|gastro.?esophageal",
  "junction - C15 or C16?",
  "secondary|metasta",                             "secondary, not a primary?",
  "in situ|dysplas",                               "not invasive cancer?",
  "overlapping",                                   "overlaps subsites within the organ",
  "uncertain|unknown|unspecified site",            "site not certain")

stopifnot(identical(review_flags$pattern, str_to_lower(review_flags$pattern)))

flag_for_review <- function(df) {
  notes <- map_chr(str_to_lower(df$term), function(tm) {
    hit <- review_flags$note[str_detect(tm, review_flags$pattern)]
    if (!length(hit)) "" else paste(hit, collapse = "; ")
  })
  df %>% mutate(look_at = notes)
}

# -----------------------------------------------------------------------------
# Talking to TRUD
# -----------------------------------------------------------------------------
# Two calls: ask which releases exist, then download the newest. The API key
# goes in the URL, which is how TRUD built it, so keep these URLs out of logs.

trud_latest_release <- function(item = trud_item,
                                key = Sys.getenv("TRUD_API_KEY")) {
  if (!nzchar(key))
    stop("no TRUD API key. Put TRUD_API_KEY=your-key in a .Renviron file and ",
         "restart R.", call. = FALSE)
  url <- sprintf("https://isd.digital.nhs.uk/trud/api/v1/keys/%s/items/%d/releases?latest",
                 key, item)
  res <- fromJSON(url, simplifyDataFrame = FALSE)
  if (!identical(res$httpStatus, 200L) || !length(res$releases))
    stop("TRUD would not list item ", item, ": ", res$message,
         "\n  - a 400 means the key is wrong; a 404 usually means the account ",
         "is not subscribed to the item.", call. = FALSE)
  res$releases[[1]]
}

trud_download <- function(release, dest_dir = tempdir()) {
  dest <- file.path(dest_dir, release$archiveFileName)
  if (file.exists(dest) && file.size(dest) == release$archiveFileSizeBytes) {
    cat("Already downloaded:", basename(dest), "\n")
    return(dest)
  }
  cat("Downloading", release$archiveFileName,
      sprintf("(%.1f GB)\n", release$archiveFileSizeBytes / 1e9))
  download.file(release$archiveFileUrl, dest, mode = "wb", quiet = FALSE)
  dest
}

# TRUD gives a checksum with every release. Worth using: a download this size
# can fail quietly, and a half-written file would otherwise surface as a
# confusing parse error.
check_download <- function(path, release) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    cat("Note: the digest package is not installed, so the checksum was not",
        "verified.\n")
    return(invisible(NA))
  }
  got <- toupper(digest::digest(path, algo = "sha256", file = TRUE))
  if (!identical(got, toupper(release$archiveFileSha256)))
    stop("the downloaded file does not match TRUD's checksum. Delete it and ",
         "try again.", call. = FALSE)
  cat("Checksum matches.\n")
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Finding the map files inside the release
# -----------------------------------------------------------------------------
# Nothing here assumes where they are: it lists what is in the zip and looks.
# Some releases nest a zip inside the zip, so it looks one level down as well.

list_release_files <- function(zip_path, pattern) {
  entries <- unzip(zip_path, list = TRUE)$Name
  hits <- entries[grepl(pattern, basename(entries))]
  if (length(hits)) return(list(zip = zip_path, entries = hits))
  
  for (z in entries[grepl("\\.zip$", entries)]) {
    tmp <- file.path(tempdir(), "trud_inner")
    dir.create(tmp, showWarnings = FALSE)
    unzip(zip_path, files = z, exdir = tmp, junkpaths = TRUE)
    got <- list_release_files(file.path(tmp, basename(z)), pattern)
    if (!is.null(got)) return(got)
  }
  NULL
}

list_map_files <- function(zip_path)
  list_release_files(zip_path, "ExtendedMap.*\\.txt$")

read_map_entry <- function(zip_path, entry) {
  map <- read.delim(unz(zip_path, entry), sep = "\t", quote = "",
                    colClasses = "character", stringsAsFactors = FALSE)
  needed <- c("id", "effectiveTime", "active", "refsetId",
              "referencedComponentId", "mapTarget")
  absent <- setdiff(needed, names(map))
  if (length(absent))
    stop(basename(entry), " has no ", paste(absent, collapse = ", "), " column.",
         "\n  - columns found: ", paste(names(map), collapse = ", "),
         call. = FALSE)
  check_map_parse(map, entry)
  map
}

# Is the file actually being read the way we think it is? These columns have
# known shapes: active is 0 or 1, and the two identifiers are plain numbers. If a
# stray tab in a free-text column ever shifted the fields, this is where it would
# show, rather than as a strange count much later on.
check_map_parse <- function(map, entry = "the map") {
  if (!nrow(map)) return(invisible(TRUE))
  bad_active <- setdiff(unique(map$active), c("0", "1"))
  if (length(bad_active))
    stop(basename(entry), " looks mis-read: the active column should only hold ",
         "0 or 1 but also holds ", paste(head(bad_active, 3), collapse = ", "),
         call. = FALSE)
  ids <- head(unique(map$referencedComponentId), 10000)
  if (any(!str_detect(ids, "^[0-9]{6,18}$")))
    stop(basename(entry), " looks mis-read: referencedComponentId should be a ",
         "plain SNOMED number but holds values such as ",
         paste(head(ids[!str_detect(ids, "^[0-9]{6,18}$")], 3), collapse = ", "),
         call. = FALSE)
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# What is actually in the release
# -----------------------------------------------------------------------------
# Opens every Snapshot map file and reports what it holds. This is the check
# that the right file was chosen: the one we want should have ICD-10-looking
# targets and a sensible number of codes landing on C15 and C16. It is slow -
# a few files of a couple of hundred thousand rows - but it only runs once, and
# it is the difference between knowing and hoping.

survey_maps <- function(zip_path, sites = og_sites) {
  found <- list_map_files(zip_path)
  if (is.null(found))
    stop("no ExtendedMap file in ", basename(zip_path),
         ". Check this is the SNOMED CT UK Clinical Edition release (item 101).",
         call. = FALSE)
  
  snaps <- found$entries[grepl("Snapshot", found$entries, ignore.case = TRUE)]
  cat("Surveying", length(snaps), "Snapshot map files (of",
      length(found$entries), "map files in the release)\n")
  
  map_dfr(snaps, function(entry) {
    map <- read_map_entry(found$zip, entry)
    if (!nrow(map))
      return(tibble(file = basename(entry), refset = NA_character_, rows = 0L,
                    codes = 0L, live_codes = 0L, og_codes_rough = 0L,
                    example_targets = ""))
    map %>%
      filter(active == "1", !is.na(mapTarget), mapTarget != "") %>%
      group_by(refset = refsetId) %>%
      summarise(
        rows = n(),
        codes = n_distinct(referencedComponentId),
        live_codes = n_distinct(referencedComponentId),
        # a rough count, for telling the maps apart. It skips the newest-row
        # reduction that og_codes_from_map() does, so it runs a little high.
        og_codes_rough = n_distinct(referencedComponentId[
          str_sub(mapTarget, 1, 3) %in% sites]),
        example_targets = paste(head(unique(mapTarget), 4), collapse = " "),
        .groups = "drop") %>%
      mutate(file = basename(entry), .before = 1)
  })
}

# -----------------------------------------------------------------------------
# What every code is called
# -----------------------------------------------------------------------------
# The release carries the names as well as the maps, so nothing has to be looked
# up by hand. Reading them means the output can go straight to a clinician as a
# list of names rather than a list of numbers, and it names the refsets too - a
# refset id is itself a concept, so the same lookup answers "what is
# 1891651000000103" without guessing.
#
# The UK release does not repeat the International descriptions, so both sets of
# description files are read and stacked.
#
# The Fully Specified Name is used rather than a synonym: it is the unambiguous
# one, and it carries the semantic tag in brackets - (disorder), (morphologic
# abnormality) - which is worth seeing.
fsn_type_id <- "900000000000003001"

read_descriptions <- function(zip_path) {
  found <- list_release_files(zip_path, "^sct2_Description.*Snapshot.*\\.txt$")
  if (is.null(found)) {
    cat("No description files in the release, so codes cannot be named.\n")
    return(tibble(snomed = character(), term = character(),
                  semantic_tag = character()))
  }
  cat("Reading names from", length(found$entries), "description files\n")
  
  out <- map_dfr(found$entries, function(entry) {
    d <- read.delim(unz(found$zip, entry), sep = "\t", quote = "",
                    colClasses = "character", stringsAsFactors = FALSE)
    need <- c("active", "conceptId", "typeId", "term")
    if (length(setdiff(need, names(d))))
      stop(basename(entry), " is not shaped like a description file.",
           call. = FALSE)
    d %>%
      filter(active == "1", typeId == fsn_type_id) %>%
      select(snomed = conceptId, term)
  })
  
  out %>%
    distinct(snomed, .keep_all = TRUE) %>%
    mutate(semantic_tag = str_extract(term, "\\(([^()]+)\\)$"),
           semantic_tag = str_remove_all(semantic_tag, "[()]"),
           term = str_trim(str_remove(term, "\\s*\\([^()]+\\)$")))
}

name_of <- function(ids, terms) {
  hit <- terms$term[match(ids, terms$snomed)]
  if_else(is.na(hit), "not named in this release", hit)
}

# Print the survey. This is a function rather than a few lines at the bottom so
# the checks can run it: the first version of this report referred to a column
# that had been renamed, and nothing noticed until it ran against the real
# release, because the checks skip the download block entirely.
report_survey <- function(survey, terms = NULL) {
  out <- survey %>%
    mutate(what = if (is.null(terms)) "not named"
           else name_of(refset, terms),
           what = if_else(is.na(refset), "file is empty", what),
           file = str_remove(file, "^der2_"),
           file = str_remove(file, "_GB.*|_INT.*"))
  out %>%
    select(file, refset, what, live_codes, og_codes_rough) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
  invisible(out)
}

# -----------------------------------------------------------------------------
# Reading the one we want
# -----------------------------------------------------------------------------
# The file is named explicitly rather than guessed at. If the naming ever
# changes this stops with a list of what was there, which is much better than
# quietly reading the wrong map - which is exactly what happened the first time
# this ran, when it picked the International map instead of the UK one.

read_extended_map <- function(zip_path, want = snomed_map_file,
                              refset = snomed_map_refset) {
  found <- list_map_files(zip_path)
  if (is.null(found))
    stop("no ExtendedMap file in ", basename(zip_path), ".", call. = FALSE)
  
  pick <- found$entries[grepl(want, basename(found$entries), fixed = TRUE)]
  if (length(pick) != 1)
    stop("expected exactly one map file matching '", want, "', found ",
         length(pick), ".\n  - map files in the release:\n    ",
         paste(basename(found$entries), collapse = "\n    "),
         call. = FALSE)
  
  cat("Using:", basename(pick), "\n")
  all_rows <- read_map_entry(found$zip, pick)
  
  present <- unique(all_rows$refsetId)
  if (!refset %in% present)
    stop("refset ", refset, " is not in ", basename(pick),
         ".\n  - refsets in that file: ", paste(present, collapse = ", "),
         "\n  - check the survey and set snomed_map_refset accordingly.",
         call. = FALSE)
  
  map <- all_rows %>% filter(refsetId == refset)
  dropped <- setdiff(present, refset)
  # A Snapshot keeps one row per map entry ever created, including entries that
  # have since been switched off, so the total is much larger than the live
  # count. Only the live ones are used. Both are printed because the total on its
  # own looks alarming - it is bigger than the number of SNOMED concepts that
  # exist - and the live count is the one to sanity-check.
  live <- map %>% filter(active == "1")
  cat("  refset", refset, "-", nrow(map), "rows |",
      n_distinct(map$referencedComponentId), "codes, of which",
      n_distinct(live$referencedComponentId), "have a live map\n")
  if (length(dropped))
    cat("  other maps in the same file, dropped:",
        paste(dropped, collapse = " "), "\n")
  map
}

# -----------------------------------------------------------------------------
# Cutting it down to the OG codes
# -----------------------------------------------------------------------------
# Keep the newest row per entry, keep the ones still in use, and take the codes
# whose ICD-10 target is C15 or C16.
#
# A SNOMED code can have more than one ICD-10 target - a cancer with a named
# secondary site gets both, and some targets only apply under a condition set
# out in the map rule. So a code is only accepted when every one of its live
# targets agrees on the same three-character site. Anything pointing at both C15
# and C18 is thrown out rather than guessed at, and reported.

og_codes_from_map <- function(map, sites = og_sites) {
  live <- map %>%
    group_by(id) %>%
    slice_max(effectiveTime, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    filter(active == "1", !is.na(mapTarget), mapTarget != "") %>%
    mutate(target = str_replace_all(mapTarget, "[^A-Z0-9]", "")) %>%
    # the map uses markers such as #NC for "cannot be classified". They are not
    # codes and must not be read as though they were.
    filter(str_detect(target, "^[A-Z][0-9]{2}")) %>%
    # A concept can need more than one ICD-10 code, and mapGroup says so: group 1
    # is the diagnosis, later groups add things like the infectious agent. So
    # "Primary carcinoma of stomach due to Epstein-Barr virus disease" is group 1
    # C169, group 2 B978 - a gastric cancer with a virus code attached, not an
    # ambiguous one. Reading every group and asking them all to agree threw that
    # code away. Only group 1 says where the tumour is.
    filter(mapGroup == "1") %>%
    mutate(site3 = str_sub(target, 1, 3))
  
  by_code <- live %>%
    group_by(snomed = referencedComponentId) %>%
    summarise(sites_seen = n_distinct(site3),
              site3 = first(site3),
              targets = paste(sort(unique(mapTarget)), collapse = " "),
              .groups = "drop")
  
  # a code whose targets disagree about the site is not guessed at. It is set
  # aside and written out, because "we threw these away" is a claim a reader
  # should be able to check rather than take on faith.
  dropped <- by_code %>%
    filter(sites_seen > 1, str_detect(targets, paste(sites, collapse = "|"))) %>%
    select(snomed, icd10_targets = targets) %>%
    mutate(reason = "targets disagree about the site")
  
  kept <- by_code %>%
    filter(sites_seen == 1, site3 %in% sites) %>%
    select(snomed, site3, icd10_targets = targets) %>%
    arrange(site3, snomed)
  
  cat("SNOMED codes meaning", paste(sites, collapse = " or "), ":",
      nrow(kept), "| set aside as unclear:", nrow(dropped), "\n")
  print(as.data.frame(count(kept, site3, name = "codes")), row.names = FALSE)
  list(kept = kept, dropped = dropped)
}

# -----------------------------------------------------------------------------
# Do it
# -----------------------------------------------------------------------------
if (!exists("skip_trud_run") || !isTRUE(skip_trud_run)) {
  release <- trud_latest_release()
  cat("Latest release:", release$name, "dated", release$releaseDate, "\n")
  
  zip_path <- trud_download(release)
  check_download(zip_path, release)
  
  terms <- read_descriptions(zip_path)
  cat("Names available for", nrow(terms), "concepts\n")
  
  survey <- survey_maps(zip_path)
  cat("\nWhat is in the release, one line per map rather than per file:\n")
  survey_named <- report_survey(survey, terms)
  write.csv(survey_named, f_snomed_survey, row.names = FALSE)
  
  cat("\n")
  map <- read_extended_map(zip_path)
  og  <- og_codes_from_map(map)
  
  # the build only needs the code and the site; the names are for people
  lookup <- og$kept
  write.csv(lookup, f_snomed_lookup, row.names = FALSE)
  
  # -- the two files to send a clinical colleague ------------------------------
  review <- og$kept %>%
    mutate(term = name_of(snomed, terms),
           semantic_tag = terms$semantic_tag[match(snomed, terms$snomed)],
           site = if_else(site3 == "C15", "C15 oesophagus", "C16 stomach")) %>%
    flag_for_review() %>%
    select(snomed, term, semantic_tag, site, icd10_targets, look_at) %>%
    arrange(desc(look_at != ""), site, term)
  write.csv(review, f_review_kept, row.names = FALSE)
  
  dropped <- og$dropped %>%
    mutate(term = name_of(snomed, terms)) %>%
    select(snomed, term, icd10_targets, reason) %>%
    arrange(term)
  write.csv(dropped, f_review_dropped, row.names = FALSE)
  
  writeLines(c(
    paste("TRUD item:", trud_item, "(SNOMED CT UK Clinical Edition, RF2)"),
    paste("Release:", release$name),
    paste("Release date:", release$releaseDate),
    paste("Archive:", release$archiveFileName),
    paste("SHA-256:", release$archiveFileSha256),
    paste("Map file:", snomed_map_file),
    paste("Map refset:", snomed_map_refset, "-",
          name_of(snomed_map_refset, terms)),
    paste("Fetched:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("Codes written:", nrow(lookup)),
    paste("Codes set aside as unclear:", nrow(dropped))), f_snomed_source)
  
  # -- what to do next --------------------------------------------------------
  cat("\nFor the build:\n")
  cat("  ", f_snomed_lookup, "-", nrow(lookup), "codes\n")
  cat("  ", f_snomed_source, "- where they came from\n")
  cat("   Copy both to the analysis server under Data/reference.\n")
  
  cat("\nFor a clinical colleague to check:\n")
  cat("  ", f_review_kept, "-", nrow(review), "codes with their names\n")
  cat("  ", f_review_dropped, "-", nrow(dropped),
      "codes we could not place\n")
  cat("   ", sum(review$look_at != ""),
      "of the kept codes are flagged as worth a second look:\n")
  if (any(review$look_at != "")) {
    review %>%
      filter(look_at != "") %>%
      count(look_at, name = "codes") %>%
      as.data.frame() %>%
      print(row.names = FALSE)
  }
}