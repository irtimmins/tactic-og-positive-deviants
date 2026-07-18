library(dplyr)
library(tidyr)
library(stringr)

# reads the cohort saved at the end of 03_cwt_merge.R - has tx_pathway,
# dtt_to_tx_reason and the stage fields carried through from the rapid
# extract untouched, since nothing upstream has dropped or recoded them
if (!exists("f_cohort_cwt")) f_cohort_cwt <- file.path(dir_out, "og_cohort_cwt.rds")
og <- readRDS(f_cohort_cwt)

stage_var <- if ("stage_RR" %in% names(og)) "stage_RR" else "stage"
cat("using", stage_var, "for stage\n\n")

# how much of the cohort is stage 4, and how tx_pathway splits by stage
cat("tx_pathway by stage (row percentages within stage)\n")
og %>%
  mutate(stage_grp = case_when(
    .data[[stage_var]] %in% c(1, 10:13) ~ "1",
    .data[[stage_var]] %in% c(2, 20:23) ~ "2",
    .data[[stage_var]] %in% c(3, 30:33) ~ "3",
    .data[[stage_var]] %in% c(4, 40:43) ~ "4",
    is.na(.data[[stage_var]])           ~ "missing",
    TRUE                                ~ "other")) %>%
  count(stage_grp, tx_pathway) %>%
  group_by(stage_grp) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  arrange(stage_grp, desc(n)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# same split but for the DTT-to-treatment reason, since that's the one that
# looked out of line with the mature cohort
cat("\ndtt_to_tx_reason by stage\n")
og %>%
  mutate(stage_grp = case_when(
    .data[[stage_var]] %in% c(1, 10:13) ~ "1",
    .data[[stage_var]] %in% c(2, 20:23) ~ "2",
    .data[[stage_var]] %in% c(3, 30:33) ~ "3",
    .data[[stage_var]] %in% c(4, 40:43) ~ "4",
    is.na(.data[[stage_var]])           ~ "missing",
    TRUE                                ~ "other")) %>%
  count(stage_grp, dtt_to_tx_reason) %>%
  group_by(stage_grp) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  arrange(stage_grp, desc(n)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# the bit that actually answers the question: is "non-positive" still common
# once stage 4 is set aside, or does it go away
cat("\nnon-positive rate, stage 1-3 vs stage 4/missing\n")
og %>%
  filter(dtt_to_tx_reason %in% c("ok", "non-positive")) %>%
  mutate(is_stage4 = .data[[stage_var]] %in% c(4, 40:43) | is.na(.data[[stage_var]])) %>%
  count(is_stage4, dtt_to_tx_reason) %>%
  group_by(is_stage4) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# a sample of the non-positive rows themselves, split the same way, so the
# raw dates can be eyeballed rather than trusted on the strength of a rate
cols_wanted <- intersect(
  c("patient_pseudo_id", stage_var, "tx_pathway", "first_tx_date",
    "cwt_dtt_date", "cwt_treat_date", "cwt_modality", "wt_dtt_to_tx"),
  names(og))

cat("\nsample non-positive rows, stage 1-3\n")
og %>%
  filter(dtt_to_tx_reason == "non-positive",
         .data[[stage_var]] %in% c(1:3, 10:13, 20:23, 30:33)) %>%
  select(all_of(cols_wanted)) %>%
  slice_sample(n = min(15, nrow(.))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)


# an external check, independent of anything in 02/03: the registry carries
# its own treatment classification (Maintxtype, with curtx flagging whether it
# was curative). worth seeing how closely our tx_pathway lines up with it
# before trusting either
if ("Maintxtype" %in% names(og)) {
  
  maintx_lab <- c(`1` = "Surgery", `2` = "Neoadj CT/RT + surgery",
                  `3` = "EMR/ESD", `4` = "Neoadj FLOT, no surgery",
                  `5` = "Neoadj CRT/RT, no surgery", `6` = "Curative C/RT",
                  `7` = "Other Chemo/RT")
  
  cat("\ntx_pathway (ours) against Maintxtype (registry's own classification)\n")
  og %>%
    mutate(maintx = if_else(is.na(Maintxtype), "missing",
                            maintx_lab[as.character(Maintxtype)])) %>%
    count(maintx, tx_pathway) %>%
    group_by(maintx) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup() %>%
    arrange(maintx, desc(n)) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
  
  cat("\nno treatment (ours) vs Maintxtype missing, by stage\n")
  og %>%
    mutate(stage_grp = case_when(
      .data[[stage_var]] %in% c(4, 40:43) ~ "4",
      is.na(.data[[stage_var]])           ~ "missing",
      TRUE                                ~ "1-3"),
      ours_none = tx_pathway == "No treatment recorded",
      registry_none = is.na(Maintxtype)) %>%
    count(stage_grp, ours_none, registry_none) %>%
    arrange(stage_grp) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
}

if ("curtx" %in% names(og)) {
  curative_pathways <- c("EMR/ESD only", "EMR/ESD then surgery",
                         "Surgery + neoadjuvant chemoRT", "Surgery + neoadjuvant chemo",
                         "Surgery + neoadjuvant RT", "Surgery + adjuvant chemo", "Surgery only",
                         "Surgery + other", "Definitive chemoRT", "Curative RT only")
  
  cat("\nour curative-pathway flag against curtx\n")
  og %>%
    mutate(curative_ours = tx_pathway %in% curative_pathways) %>%
    count(curative_ours, curtx) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
}

# Endototx is the registry's own endoscopy-to-first-treatment interval - a
# second, independent check, this time on the waiting-time side rather than
# the pathway classification
if ("Endototx" %in% names(og) && "endoscopy_date" %in% names(og)) {
  cat("\nendoscopy-to-treatment: ours vs the registry's Endototx\n")
  tx_col <- if ("tx_date_used" %in% names(og)) "tx_date_used" else "first_tx_date"
  cmp <- og %>%
    filter(!is.na(endoscopy_date), !is.na(.data[[tx_col]]), !is.na(Endototx)) %>%
    mutate(ours = as.integer(.data[[tx_col]] - endoscopy_date),
           diff = ours - Endototx)
  cat("rows compared:", nrow(cmp), "\n")
  cat("median difference (ours minus registry):", median(cmp$diff), "days\n")
  cat("within 7 days:", round(100 * mean(abs(cmp$diff) <= 7), 1), "%\n")
}


# the reason function blanks any interval <= 0 as "non-positive", which
# lumps together two very different things: a genuinely impossible negative
# gap, and a same-day CWT record (treat_start == treat_period_start), which
# is ordinary for palliative/monitoring periods and costs nothing to tell
# apart, since tx_date_used and cwt_dtt_date are both still on the object
if (all(c("tx_date_used", "cwt_dtt_date") %in% names(og))) {
  cat("\nsplitting non-positive into exactly-zero vs genuinely negative\n")
  og %>%
    filter(dtt_to_tx_reason == "non-positive") %>%
    mutate(raw_gap = as.integer(tx_date_used - cwt_dtt_date),
           kind = if_else(raw_gap == 0, "same day (0)", "negative")) %>%
    count(kind) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
  
  cat("\nsame-day share of non-positive, by stage\n")
  og %>%
    filter(dtt_to_tx_reason == "non-positive") %>%
    mutate(raw_gap = as.integer(tx_date_used - cwt_dtt_date),
           stage_grp = if_else(.data[[stage_var]] %in% c(4, 40:43) |
                                 is.na(.data[[stage_var]]), "4/missing", "1-3"),
           kind = if_else(raw_gap == 0, "same day (0)", "negative")) %>%
    count(stage_grp, kind) %>%
    group_by(stage_grp) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup() %>%
    as.data.frame() %>%
    print(row.names = FALSE)
}

# the Maintxtype gap: patients the registry's own classification found a
# treatment for, that our tx_pathway calls "No treatment recorded". Pull the
# raw registry date fields for a sample, to see which one is actually blank
if (all(c("Maintxtype", "surgery_date", "sact_date", "rt_date",
          "emresd_date") %in% names(og))) {
  cat("\nsample: Maintxtype found something, ours says no treatment\n")
  og %>%
    filter(tx_pathway == "No treatment recorded", !is.na(Maintxtype)) %>%
    select(patient_pseudo_id, Maintxtype, surgery_date, sact_date, rt_date,
           emresd_date, cwt_treat_date, cwt_modality) %>%
    slice_sample(n = min(15, nrow(.))) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
}

cat("\nsample non-positive rows, stage 4 or missing\n")
og %>%
  filter(dtt_to_tx_reason == "non-positive",
         .data[[stage_var]] %in% c(4, 40:43) | is.na(.data[[stage_var]])) %>%
  select(all_of(cols_wanted)) %>%
  slice_sample(n = min(15, nrow(.))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)