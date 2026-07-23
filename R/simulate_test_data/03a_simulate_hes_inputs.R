# =============================================================================
# 90  Make stand-in data
# -----------------------------------------------------------------------------
# Writes a pair of dta files that stand in for the rapid tumour extract and the
# HES-APC extract, so the build can be run and checked away from the analysis
# server. The messy bits are taken from the profiles of the 2026 data: sitetret
# 2% blank with the odd default and test code, procode3 three characters,
# operation dates as year-month-day text that is often blank, and an endotrustHES
# that is almost always the three-character provider rather than a site.
#
# Because this script knows each patient's true endoscopy site and which of four
# outcomes they were built to fall into, it writes those answers out too. 91 uses
# them to ask how often the build gets it right, which no check against the real
# data can do.
#
# The number of HES episodes per patient is kept low here on purpose - a handful,
# not the real dozen-plus - because the logic is exercised by the endoscopy
# episode and a little filler, and a faithful episode count would only make the
# reshape slow with nothing gained. Nothing here is meant to be analysed.
#
# Writes: <dir_sim>/*.dta and endoscopy_truth.rds
# =============================================================================

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(magrittr)   # for %>%
})

if (!exists("dir_sim")) dir_sim <- "Output/sim/raw"
if (!exists("sim_scale")) sim_scale <- 1
if (!exists("sim_seed")) sim_seed <- 20260212
dir.create(dir_sim, recursive = TRUE, showWarnings = FALSE)
set.seed(sim_seed)

n_patients <- round(95333 * sim_scale)

# the diagnostic endoscopy codes the build looks for, and a few non-endoscopy
# codes for the filler episodes (a resection, an attendance) so not every episode
# is an endoscopy
endo_codes  <- c("G451", "G459", "G161", "G169", "G162", "G433", "G443", "G191")
other_codes <- c("G011", "G274", "H221", "X403", "U051")

# -----------------------------------------------------------------------------
# Trusts and their hospital sites
# -----------------------------------------------------------------------------
trusts <- unique(paste0("R", sample(LETTERS, 200, TRUE), sample(0:9, 200, TRUE)))
trusts <- sort(trusts[1:121])
sites <- map_dfr(trusts, function(tr)
  tibble(trust = tr,
         site = paste0(tr, sprintf("%02d", sample(1:20, sample(2:8, 1))))))
sites_by_trust <- split(sites$site, sites$trust)

# -----------------------------------------------------------------------------
# The rapid tumour extract, one row per patient
# -----------------------------------------------------------------------------
# 81.1% of patients have a HES endoscopy anchor (endoHES==1); the rest have no
# endoscopy date and never enter the match. Every anchored patient is assigned
# one of four outcomes, in the rough proportions the real within-extract data
# showed, plus a slice who are simply not in the HES extract at all.
truth <- tibble(
  patient_pseudo_id = sprintf("P%09d", seq_len(n_patients)),
  true_trust  = sample(trusts, n_patients, TRUE),
  tumour_site = sample(c("C15", "C16"), n_patients, TRUE, c(0.679, 0.321))) %>%
  mutate(true_site = map_chr(true_trust, ~ sample(sites_by_trust[[.x]], 1)),
         has_anchor = runif(n()) < 0.811,
         u = runif(n()),
         endo_class = case_when(
           !has_anchor        ~ "no_anchor",
           u < 0.930          ~ "sited",           # APC endoscopy at the date
           u < 0.955          ~ "offwindow",       # APC endoscopy, weeks late
           u < 0.965          ~ "no_apc",          # APC, but no endoscopy code
           TRUE               ~ "not_in_hes"))      # not in the APC extract

rapid <- truth %>%
  mutate(
    ind_pseudo_id    = seq_len(n()),
    tumour_pseudo_id = seq_len(n()),
    ref_date         = as.Date("2018-01-23") + sample(0:2800, n(), TRUE),
    endoHES          = if_else(has_anchor, 1, 0),
    endodateHES      = as.Date(if_else(has_anchor, ref_date, as.Date(NA))),
    endotypeHES      = if_else(has_anchor, sample(endo_codes, n(), TRUE), ""),
    # endotrustHES is almost always the three-character provider; a few carry the
    # full five-character site, as in the real field
    endotrustHES     = case_when(
      !has_anchor       ~ "",
      runif(n()) < 0.02 ~ true_site,
      TRUE              ~ str_sub(true_site, 1, 3)),
    tumour_morphology_str = sample(c("8140", "8070", "8010"), n(), TRUE),
    diagnosis_trust  = true_trust,
    gender           = sample(1:2, n(), TRUE, c(0.75, 0.25)),
    age              = as.integer(pmin(95, pmax(20, round(rnorm(n(), 73, 11)))))) %>%
  select(patient_pseudo_id, ind_pseudo_id, tumour_pseudo_id, tumour_site,
         endoHES, endodateHES, endotypeHES, endotrustHES,
         tumour_morphology_str, diagnosis_trust, gender, age)

# -----------------------------------------------------------------------------
# The HES-APC extract, several episodes per patient
# -----------------------------------------------------------------------------
# Build the episodes as a plain list of rows first (patient, dates, site, and at
# most one coded operation with its slot), then lay the 24 operation and 24 date
# columns out as matrices. That keeps it fast at full size, where binding a row
# at a time would not be.
ref_of <- setNames(rapid$endodateHES, rapid$patient_pseudo_id)

make_rows <- function(pid, class, true_site) {
  ref <- ref_of[[pid]]
  prov <- str_sub(true_site, 1, 3)
  rows <- list()
  add <- function(epistart, epiorder, site, code, opdate) {
    rows[[length(rows) + 1]] <<- tibble(
      patient_pseudo_id = pid, epistart = epistart, epiorder = epiorder,
      sitetret = site, procode3 = prov, endo_code = code, endo_date = opdate)
  }
  # a site value spoiled the way the real field is: usually the true site, but
  # sometimes blank, a default, lower case, or the three-character provider
  spoil_site <- function(s) {
    r <- runif(1)
    if (r < 0.02) "" else if (r < 0.03) "00000"
    else if (r < 0.04) str_to_lower(s) else if (r < 0.05) str_sub(s, 1, 3) else s
  }
  if (class == "sited") {
    # the endoscopy episode: a code in one of the first slots, on the reference
    # date give or take a day. A minority also get a second endoscopy the same
    # week whose site is unusable, to check the good site is preferred.
    d <- ref + sample(0:1, 1)
    add(as.character(d), 1L, true_site, sample(endo_codes, 1), fmt_date(d))
    if (runif(1) < 0.15)
      add(as.character(ref - 1), 2L, sample(c("", "00000"), 1),
          sample(endo_codes, 1), fmt_date(ref - 1))
    if (runif(1) < 0.5)
      add(as.character(ref + 60), 1L, true_site, sample(other_codes, 1),
          fmt_date(ref + 60))
  } else if (class == "offwindow") {
    d <- ref + sample(20:40, 1)
    add(as.character(d), 1L, spoil_site(true_site), sample(endo_codes, 1),
        fmt_date(d))
  } else if (class == "no_apc") {
    add(as.character(ref + sample(0:30, 1)), 1L, spoil_site(true_site),
        sample(other_codes, 1), fmt_date(ref))
    if (runif(1) < 0.5)
      add(as.character(ref + 90), 1L, true_site, sample(other_codes, 1),
          fmt_date(ref + 90))
  }
  # "not_in_hes" and "no_anchor" patients get no rows here
  if (length(rows)) bind_rows(rows) else NULL
}

# operation dates are text; a fraction are written with no dashes, and the real
# field is often blank, so leave the filler episodes' dates empty
fmt_date <- function(d) {
  if (is.na(d)) return("")
  if (runif(1) < 0.1) format(d, "%Y%m%d") else format(d, "%Y-%m-%d")
}

sited <- truth %>% filter(endo_class %in% c("sited", "offwindow", "no_apc"))
episodes <- pmap(list(sited$patient_pseudo_id, sited$endo_class,
                      sited$true_site), make_rows)
episodes <- bind_rows(episodes)

# lay out the 24 operation and 24 date columns
n <- nrow(episodes)
OP <- matrix("-", n, 24, dimnames = list(NULL, sprintf("opertn_%02d", 1:24)))
OD <- matrix("",  n, 24, dimnames = list(NULL, sprintf("opdate_%02d", 1:24)))
slot <- sample(1:4, n, TRUE)                     # which slot the code sits in
has_code <- !is.na(episodes$endo_code)
for (i in which(has_code)) {
  OP[i, slot[i]] <- episodes$endo_code[i]
  OD[i, slot[i]] <- episodes$endo_date[i]
}

hes <- bind_cols(
  tibble(patient_pseudo_id = episodes$patient_pseudo_id,
         epistart = episodes$epistart,
         epiend   = episodes$epistart,
         admidate = episodes$epistart,
         epiorder = episodes$epiorder,
         epitype  = 1L,
         sitetret = episodes$sitetret,
         procode3 = episodes$procode3),
  as_tibble(OP), as_tibble(OD))

# a few HES patients the registry does not have, in a shuffled order, so the
# build's patient filter and the row order are both exercised
extra <- tibble(patient_pseudo_id = sprintf("Q%09d", 1:200),
                epistart = "2022-01-01", epiend = "2022-01-01",
                admidate = "2022-01-01", epiorder = 1L, epitype = 1L,
                sitetret = "RZZ01", procode3 = "RZZ")
extra <- bind_cols(extra,
                   as_tibble(matrix("-", 200, 24,
                                    dimnames = list(NULL, sprintf("opertn_%02d", 1:24)))),
                   as_tibble(matrix("", 200, 24,
                                    dimnames = list(NULL, sprintf("opdate_%02d", 1:24)))))
hes <- bind_rows(hes, extra) %>% slice_sample(prop = 1)

# -----------------------------------------------------------------------------
# Write
# -----------------------------------------------------------------------------
f_rapid <- file.path(dir_sim,
                     "20260212_Rapidtumour_linked_2026SOTN_clean_OG_postPT.dta")
f_hes   <- file.path(dir_sim, "20260320_nic709865_hes_apc_202510_OG.dta")

write_dta(rapid, f_rapid)
write_dta(hes, f_hes)
saveRDS(truth %>%
          mutate(true_endo_site = if_else(endo_class == "sited", true_site,
                                          NA_character_)) %>%
          select(patient_pseudo_id, tumour_site, endo_class, true_endo_site),
        file.path(dir_sim, "endoscopy_truth.rds"))

cat("Made", nrow(rapid), "registry rows and", nrow(hes), "HES episodes (",
    n_distinct(hes$patient_pseudo_id), "patients ) in", dir_sim, "\n")
cat("  anchored patients:", sum(rapid$endoHES == 1), "\n")
cat("  outcome classes among anchored:\n")
truth %>% filter(has_anchor) %>% count(endo_class) %>% as.data.frame() %>%
  print(row.names = FALSE)
cat("  sitetret blank/default in HES:",
    round(100 * mean(hes$sitetret == "" | hes$sitetret == "00000"), 1), "%\n")