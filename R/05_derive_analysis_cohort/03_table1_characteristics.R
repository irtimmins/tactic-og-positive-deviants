# 03  Table 1: cohort characteristics with the pathway intervals
# -----------------------------------------------------------------------------
# One row per characteristic level: number and percent of patients, then the
# mean (SD) of each of the three pathway intervals in that group -
#   endoscopy -> decision-to-treat   (the interval the cohort is defined on)
#   decision-to-treat -> treatment
#   endoscopy -> treatment           (the whole pathway; the sum of the two)
# Written to Word with flextable where available, and always to CSV. Also a
# short representativeness check by cancer alliance.

library(dplyr)

source("R/05_derive_analysis_cohort/01_config.R")
df <- readRDS(cohort_rds)

# route to diagnosis: guard for a cohort built before this field was carried
# through, so an older pd_cohort.rds still produces a table rather than erroring
# on a column that is not there. The block is left out of the table below when
# nothing is recorded.
if (!"tfinal_route" %in% names(df)) df$tfinal_route <- NA_character_
has_route <- any(!is.na(df$tfinal_route))

# pathway grouping for reporting ----------------------------------------------
# The cohort-membership decision (which raw tx_pathway values are curative) is
# made in 02_build_cohort.R / 01_config.R; this is a purely presentational
# regrouping of pathways that are already in the cohort, so it lives here rather
# than upstream:
#   - the three neoadjuvant routes (chemo / RT / chemoRT) are combined into one
#     "Surgery + neoadjuvant chemo/RT" row, since the neoadjuvant regimen detail
#     is not the axis this table reports on
#   - "EMR/ESD only" and "EMR/ESD then surgery" are combined into "EMR/ESD".
#     NEEDS CLINICAL REVIEW: the "then surgery" patients' wait intervals are
#     still anchored on the surgery date rather than the EMR/ESD date - see the
#     printed flag in 02_build_cohort.R and its comment in 01_config.R. Counted
#     again here so the review need is visible in this table's own output too.
pathway_group_map <- c(
  "Surgery only"                   = "Surgery only",
  "Surgery + neoadjuvant chemo"    = "Surgery + neoadjuvant chemo/RT",
  "Surgery + neoadjuvant RT"       = "Surgery + neoadjuvant chemo/RT",
  "Surgery + neoadjuvant chemoRT"  = "Surgery + neoadjuvant chemo/RT",
  "Surgery + adjuvant chemo"       = "Surgery + adjuvant chemo",
  "Definitive chemoRT"             = "Definitive chemoRT",
  "Curative RT only"               = "Curative RT only",
  "EMR/ESD only"                   = "EMR/ESD",
  "EMR/ESD then surgery"           = "EMR/ESD")
pathway_levels <- c("Surgery only", "Surgery + neoadjuvant chemo/RT",
                    "Surgery + adjuvant chemo", "Definitive chemoRT",
                    "Curative RT only", "EMR/ESD")

n_emr_then_surg <- sum(df$tx_pathway == "EMR/ESD then surgery", na.rm = TRUE)
if (n_emr_then_surg > 0)
  cat(sprintf(paste0(
    "NEEDS CLINICAL REVIEW: the 'EMR/ESD' row below includes %d patient(s) ",
    "whose pathway was 'EMR/ESD then surgery' - their wait intervals are ",
    "anchored on the surgery date, not the EMR/ESD date. See 02_build_cohort.R.\n"),
    n_emr_then_surg))

d <- df %>% mutate(
  age_grp = cut(agediag, c(-Inf, 50, 60, 70, 80, Inf),
                labels = c("<50", "50-59", "60-69", "70-79", "80+")),
  sex     = factor(sexf, levels = c("Male", "Female")),
  stage_f = factor(as.character(stage), levels = c("1", "2", "3")),
  path_f  = factor(unname(pathway_group_map[tx_pathway]), levels = pathway_levels),
  eth_f   = factor(as.character(eth),
                   levels = c("White", "Asian", "Black", "Mixed", "Other")),
  # IMD quintile 1 is the MOST deprived and 5 the LEAST deprived in the
  # registry's own coding (see imd_2019_RR in 02_build_cohort.R) - the same
  # convention the registry's own quintile_2019 field uses ("1 - most
  # deprived" ... "5 - least deprived"). The table is ordered from least to
  # most deprived (5 down to 1), with those two labels made explicit so the
  # direction cannot be misread from the number alone.
  imd_f   = factor(dplyr::recode(as.character(imd),
                                 "1" = "1 - most deprived",
                                 "5" = "5 - least deprived"),
                   levels = c("5 - least deprived", "4", "3", "2",
                              "1 - most deprived")),
  cci_f   = factor(ifelse(cci_n_conditions >= 2, "2+", as.character(cci_n_conditions)),
                   levels = c("0", "1", "2+")),
  season_f = factor(case_when(q2 == 1 ~ "Apr-Jun",
                              q3 == 1 ~ "Jul-Sep",
                              q4 == 1 ~ "Oct-Dec",
                              TRUE    ~ "Jan-Mar"),
                    levels = c("Jan-Mar", "Apr-Jun", "Jul-Sep", "Oct-Dec")),
  year_f   = factor(ydiag),
  # route to diagnosis, in a fixed clinical order rather than alphabetical. Only
  # the elective routes remain by this point (02 drops emergency and unknown),
  # but the factor is built from whatever is present so a route added to the
  # registry later still appears rather than becoming a silent NA.
  route_f  = factor(tfinal_route,
                    levels = intersect(
                      c("TWW", "GP referral", "Other outpatient",
                        "Inpatient elective", "Emergency presentation",
                        "Unknown"),
                      unique(tfinal_route)))
)

# The three pathway intervals reported alongside each characteristic. The first
# defines the cohort, so it is always present; the other two can be missing for a
# patient whose treatment date was unusable, so each is summarised over the
# patients who have it and the mean is shown as "-" where a group has none.
# n_col records how many contributed, so a thin cell is visible as thin rather
# than looking like a confident estimate.
wait_vars <- c(wait          = "Endoscopy to decision-to-treat",
               wt_dtt_to_tx  = "Decision-to-treat to treatment",
               wt_endo_to_tx = "Endoscopy to treatment")
for (v in names(wait_vars))
  if (!v %in% names(d)) d[[v]] <- NA_integer_

msd <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("-")
  if (length(x) == 1) return(sprintf("%.1f (-)", x))
  sprintf("%.1f (%.1f)", mean(x), sd(x))
}

# one block per characteristic: a heading row (the label, no data alongside),
# then one indented row per level with n (%) and the mean (SD) of each interval.
cat_block <- function(d, var, label) {
  levels_tab <- d %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(Level = .data[[var]]) %>%
    summarise(np = n(),
              w1 = msd(wait),
              w2 = msd(wt_dtt_to_tx),
              w3 = msd(wt_endo_to_tx),
              .groups = "drop") %>%
    transmute(item = as.character(Level),
              patients = sprintf("%s (%.1f%%)", formatC(np, format = "d", big.mark = ","),
                                 100 * np / nrow(d)),
              endo_to_dtt = w1, dtt_to_tx = w2, endo_to_tx = w3,
              is_head = FALSE, indent = TRUE)
  bind_rows(tibble(item = label, patients = "", endo_to_dtt = "",
                   dtt_to_tx = "", endo_to_tx = "",
                   is_head = TRUE, indent = FALSE),
            levels_tab)
}

cat(sprintf("Age at diagnosis: mean %.1f years (SD %.1f)\n", mean(d$agediag), sd(d$agediag)))
for (v in names(wait_vars))
  cat(sprintf("%-31s mean %s days, recorded for %d of %d (%.1f%%)\n",
              wait_vars[[v]], msd(d[[v]]), sum(!is.na(d[[v]])), nrow(d),
              100 * mean(!is.na(d[[v]]))))

tab1 <- bind_rows(
  tibble(item = "Patients, total",
         patients = formatC(nrow(d), format = "d", big.mark = ","),
         endo_to_dtt = msd(d$wait),
         dtt_to_tx   = msd(d$wt_dtt_to_tx),
         endo_to_tx  = msd(d$wt_endo_to_tx),
         is_head = FALSE, indent = FALSE),
  cat_block(d, "age_grp",  "Age group"),
  cat_block(d, "sex",      "Sex"),
  cat_block(d, "stage_f",  "Stage at diagnosis"),
  cat_block(d, "path_f",   "Treatment pathway"),
  cat_block(d, "eth_f",    "Ethnicity"),
  cat_block(d, "imd_f",    "Deprivation quintile"),
  cat_block(d, "cci_f",    "RCS Charlson comorbidity score"),
  if (has_route) cat_block(d, "route_f", "Route to diagnosis"),
  cat_block(d, "year_f",   "Calendar year"),
  cat_block(d, "season_f", "Season of diagnosis")
)

tab_cols <- c("item", "patients", "endo_to_dtt", "dtt_to_tx", "endo_to_tx")
write.csv(tab1[, tab_cols],
          file.path(out_dir, "table1_characteristics.csv"), row.names = FALSE)
cat("Table 1 written to", file.path(out_dir, "table1_characteristics.csv"), "\n")

# Word version, if flextable / officer are available. Wrapped so a missing
# package cannot lose the CSV already written.
if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer",   quietly = TRUE)) {
  library(flextable); library(officer)
  
  head_rows   <- which(tab1$is_head)
  indent_rows <- which(tab1$indent)
  total_row   <- which(tab1$item == "Patients, total")
  
  disp <- tab1[, tab_cols]
  ft <- flextable(disp)
  ft <- set_header_labels(ft, item = "Characteristic",
                          patients = "Patients, n (%)",
                          endo_to_dtt = "Endoscopy to\ndecision-to-treat",
                          dtt_to_tx   = "Decision-to-treat\nto treatment",
                          endo_to_tx  = "Endoscopy to\ntreatment")
  ft <- add_header_row(ft, top = TRUE,
                       values = c("", "", "Days, mean (SD)"),
                       colwidths = c(1, 1, 3))
  ft <- bold(ft, part = "header")
  ft <- bold(ft, j = "item", i = c(head_rows, total_row), part = "body")
  ft <- align(ft, j = "item", align = "left", part = "all")
  ft <- align(ft, j = setdiff(tab_cols, "item"), align = "center", part = "all")
  ft <- align(ft, i = 1, align = "center", part = "header")
  ft <- valign(ft, valign = "top", part = "body")
  ft <- add_footer_lines(ft, paste(
    "Patients as n (%); age at diagnosis as mean (SD). The three intervals are",
    "mean (SD) days. The cohort is defined on a positive endoscopy-to-decision-to-treat",
    "time, so that interval is present for every patient; the two treatment-side",
    "intervals are shown for the patients with a usable treatment date, and a group",
    "with none is marked \"-\"."))
  ft <- fontsize(ft, size = 8, part = "all")
  ft <- padding(ft, padding.top = 1, padding.bottom = 1, part = "all")
  ft <- padding(ft, j = "item", padding.left = 3, part = "body")
  ft <- padding(ft, i = indent_rows, j = "item", padding.left = 18, part = "body")
  ft <- line_spacing(ft, space = 1.3, part = "all")
  
  ft <- border_remove(ft)
  box   <- fp_border(color = "black",  width = 1)
  soft  <- fp_border(color = "grey65", width = 0.5)
  faint <- fp_border(color = "grey85", width = 0.5)
  ft <- border_outer(ft, border = box, part = "all")
  ft <- hline_bottom(ft, border = soft, part = "header")
  ft <- hline(ft, i = seq_len(nrow(disp) - 1), border = faint, part = "body")
  
  ft <- autofit(ft)
  ft <- width(ft, j = "item",     width = 1.75)
  ft <- width(ft, j = "patients", width = 1.0)
  ft <- width(ft, j = c("endo_to_dtt", "dtt_to_tx", "endo_to_tx"), width = 1.15)
  ft <- set_table_properties(ft, layout = "fixed", align = "left")
  
  save_as_docx(ft, path = file.path(out_dir, "table1_characteristics.docx"))
  cat("Table 1 also written to", file.path(out_dir, "table1_characteristics.docx"), "\n")
} else {
  cat("flextable / officer not installed - Word Table 1 skipped (CSV written).\n")
}

# representativeness: patient share by cancer alliance -----------------------
alliance <- df %>%
  count(canalliance, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(desc(n))
print(as.data.frame(alliance), row.names = FALSE)
write.csv(alliance, file.path(out_dir, "table1_alliance.csv"), row.names = FALSE)
cat("03 complete. Next: stage 06 estimates the standardised hospital means.\n")