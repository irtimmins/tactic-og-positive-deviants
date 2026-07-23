# 04  posterior ranks, caterpillar plot, and positive-deviant candidates
# -----------------------------------------------------------------------------
# Turn the shrinkage fit into shrunk site means with credible intervals, expected
# ranks, and the probability of sitting in the fastest X%. Draws the caterpillar
# plots, writes the ranking tables, and selects the positive-deviant candidates.
#
# Both estimands are processed the same way, because in both a LOWER number is
# better and so rank 1 is the best performer:
#   SUSTAINED  a lower standardised mean = consistently faster.
#   IMPROVED   a lower (more negative) change = a bigger reduction in waiting
#              time from period 1 to period 2.
# The two differ only in what they are ranked against: the sustained plot sits
# against the average site, the improvement plot against zero (no change).
#
# The improvement outputs are produced only when 02 found the estimand
# estimable; when it did not, this script writes the sustained results alone and
# says so.
#
# Reads : stan_sustained.rds, stan_improve.rds (if present)
# Writes: ranks_sustained.csv, caterpillar_sustained.pdf, top_ranking.csv,
#         candidates.csv, and the _improve counterparts

library(rstan)
library(dplyr)

source("R/06_identify_positive_deviants/01_config.R")

# optional: a site code to highlight in the caterpillar plot
if (!exists("highlight_hosp")) highlight_hosp <- NA

process_fit <- function(obj) {
  draws <- rstan::extract(obj$fit, pars = "y_site_true")$y_site_true  # draws x J
  site  <- obj$site
  rm <- rank_metrics(draws)
  site %>%
    mutate(
      post_mean = colMeans(draws),
      post_sd   = apply(draws, 2, sd),
      ci_lo     = apply(draws, 2, quantile, probs = 0.025),
      ci_hi     = apply(draws, 2, quantile, probs = 0.975)
    ) %>%
    bind_cols(rm)
}

sus <- process_fit(readRDS(file.path(out_dir, "stan_sustained.rds")))
write.csv(sus, file.path(out_dir, "ranks_sustained.csv"), row.names = FALSE)

imp <- if (file.exists(stan_improve_rds)) {
  process_fit(readRDS(stan_improve_rds))
} else {
  NULL
}
if (!is.null(imp)) {
  write.csv(imp, file.path(out_dir, "ranks_improve.csv"), row.names = FALSE)
  cat(sprintf("improvement ranks written for %d sites\n", nrow(imp)))
} else {
  cat("no stan_improve.rds - sustained results only\n")
}

# caterpillar plot, if ggplot2 is available ----------------------------------
# the plot follows the house style used in the summary figure: theme_classic,
# points and credible-interval segments in a single house colour, a heavier
# semi-transparent reference line, and a y-axis that can be floored at zero. The
# style is tunable in the small block below.
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("scales", quietly = TRUE)) {
  library(ggplot2)
  library(scales)
  
  # tunable style (adjust here) ----------------------------------------------
  axis_title_size <- 11      # axis titles
  axis_text_size  <- 10       # axis tick labels
  cat_ci_alpha    <- 0.35    # credible-interval opacity - lighter than the points
  cat_ci_lwd      <- 0.7    # credible-interval line width
  cat_pt_size     <- 1.2     # point size
  cat_xtitle_gap  <- 2.5     # pt gap from the x-axis line to the "Hospitals" title
  col_base        <- "darkblue"    # points and intervals
  col_high        <- "darkorange3" # a highlighted hospital, if one is named
  
  # hospitals ordered by expected posterior rank; segment = credible interval; a
  # grey reference line at the mean. from_zero floors the y-axis at 0.
  caterpillar <- function(d, ylab, title, ref_line = mean(d$post_mean),
                          from_zero = FALSE, highlight = NA) {
    d <- d %>% arrange(exp_rank) %>% mutate(rank = row_number())
    hl <- if (!is.na(highlight)) d$diag_hosp == highlight else rep(FALSE, nrow(d))
    p <- ggplot(d, aes(rank, post_mean)) +
      geom_hline(yintercept = ref_line, linewidth = 1, alpha = 0.6,
                 colour = "gray30") +
      geom_segment(aes(xend = rank, y = ci_lo, yend = ci_hi),
                   colour = col_base, alpha = cat_ci_alpha, linewidth = cat_ci_lwd) +
      geom_point(shape = 16, colour = col_base, size = cat_pt_size) +
      labs(title = title) +
      theme_classic(base_size = 11) +
      theme(axis.title  = element_text(size = axis_title_size),
            axis.text   = element_text(size = axis_text_size),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.length.x = unit(0, "pt"),
            axis.title.x = element_text(size = axis_title_size,
                                        margin = margin(t = cat_xtitle_gap)),
            legend.position = "none") +
      scale_x_continuous("Hospitals") +
      scale_y_continuous(ylab, breaks = breaks_width(10))
    if (from_zero) p <- p + coord_cartesian(ylim = c(0, max(d$ci_hi) * 1.02))
    if (any(hl))
      p <- p + geom_point(data = d[hl, ], colour = col_high, size = cat_pt_size * 2)
    p
  }
  ggsave(file.path(out_dir, "caterpillar_sustained.pdf"),
         caterpillar(sus, "Mean waiting time, endoscopy to decision-to-treat (days)",
                     sprintf("Sustained performance, %s", window_label),
                     from_zero = TRUE, highlight = highlight_hosp),
         width = 7, height = 6)
  cat("caterpillar plot written to caterpillar_sustained.pdf\n")
  
  # the improvement plot is referenced against zero, not the average site: the
  # question is whether a site got faster, not whether it improved more than its
  # peers did. The axis is not floored at zero, since a change is signed.
  if (!is.null(imp)) {
    ggsave(file.path(out_dir, "caterpillar_improve.pdf"),
           caterpillar(imp,
                       sprintf("Change in mean waiting time, %s vs %s (days)",
                               period_2_label, period_1_label),
                       sprintf("Improvement, %s vs %s",
                               period_2_label, period_1_label),
                       ref_line = 0, from_zero = FALSE,
                       highlight = highlight_hosp),
           width = 7, height = 6)
    cat("caterpillar plot written to caterpillar_improve.pdf\n")
  }
} else {
  cat("ggplot2 / scales not installed - caterpillar plot skipped (CSV written).\n")
}

# ranking tables -------------------------------------------------------------
# one row per site, best first: shrunk estimate, effective sample size, and the
# posterior probability of sitting in the top 10 / 20 / 25 percent. The two
# estimands share the layout; only the sample-size columns differ, because the
# improvement estimand has a count in each period rather than one overall.
ranking_table <- function(d, improve = FALSE) {
  d <- d %>% arrange(exp_rank)
  out <- data.frame(
    rank     = seq_len(nrow(d)),
    exp_rank = round(d$exp_rank, 1),
    site     = d$diag_hosp,
    stringsAsFactors = FALSE)
  if (improve) {
    out$n_period_1 <- d$n1
    out$n_period_2 <- d$n2
    out$ess        <- round(d$ess1 + d$ess2, 0)
  } else {
    out$n   <- d$n
    out$ess <- round(d$n_eff, 0)
  }
  out$estimate <- sprintf("%.1f (%.1f to %.1f)", d$post_mean, d$ci_lo, d$ci_hi)
  out$p_top10  <- round(d$p_top10, 2)
  out$p_top20  <- round(d$p_top20, 2)
  out$p_top25  <- round(d$p_top25, 2)
  out
}

tab <- ranking_table(sus)
write.csv(tab, file.path(out_dir, "top_ranking.csv"), row.names = FALSE)
cat("\nfastest 10 sites (sustained):\n")
print(as.data.frame(head(tab, 10)), row.names = FALSE)

if (!is.null(imp)) {
  tab_imp <- ranking_table(imp, improve = TRUE)
  write.csv(tab_imp, file.path(out_dir, "top_ranking_improve.csv"), row.names = FALSE)
  cat(sprintf("\nmost improved 10 sites (%s vs %s; negative = faster):\n",
              period_2_label, period_1_label))
  print(as.data.frame(head(tab_imp, 10)), row.names = FALSE)
}

# formatted Word ranking table -----------------------------------------------
# The same rankings as the CSVs above, laid out for reading: one block per
# estimand in a single document, with a banner naming each, the column titles
# repeated before each block, and a heading spanning the three probability
# columns. Caption and footnote sit outside the table as ordinary paragraphs so
# Word keeps them with it.
#
# Built as a function taking the two processed data frames, rather than inline,
# so the layout can be exercised without a Stan fit to hand.
#
# top_n caps each block; a block with fewer sites than that simply shows all of
# them. When there is no improvement estimand the document is written with the
# sustained block alone rather than being skipped.
ranking_docx <- function(sus, imp, path, top_n = 20) {
  
  hosp_name <- hospital_names()
  col_names <- c("Rank", "ExpRank", "Hospital", "SiteCode", "N", "ESS",
                 "Estimate", "p10", "p20", "p25")
  
  # one block's worth of rows. n_sample / ess_eff differ by estimand: the
  # sustained figures are the window's count and Kish effective n, the
  # improvement ones are the two periods added together.
  top_block <- function(d, improve = FALSE) {
    d <- d %>% arrange(exp_rank)
    d <- d[seq_len(min(top_n, nrow(d))), , drop = FALSE]
    n_sample <- if (improve) d$n1 + d$n2 else d$n
    ess_eff  <- if (improve) d$ess1 + d$ess2 else d$n_eff
    nm <- if (is.null(hosp_name)) rep(NA_character_, nrow(d))
    else unname(hosp_name[as.character(d$diag_hosp)])
    data.frame(
      Rank     = as.character(seq_len(nrow(d))),
      ExpRank  = sprintf("%.1f", d$exp_rank),
      Hospital = ifelse(is.na(nm) | nm == "", as.character(d$diag_hosp), nm),
      SiteCode = as.character(d$diag_hosp),
      N        = as.character(n_sample),
      ESS      = sprintf("%.0f", ess_eff),
      Estimate = sprintf("%.1f (%.1f to %.1f)", d$post_mean, d$ci_lo, d$ci_hi),
      p10      = sprintf("%.2f", d$p_top10),
      p20      = sprintf("%.2f", d$p_top20),
      p25      = sprintf("%.2f", d$p_top25),
      stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  blank_row <- function(first_cell = "") {
    r <- as.data.frame(as.list(rep("", length(col_names))), stringsAsFactors = FALSE)
    names(r) <- col_names
    r[[1]] <- first_cell
    r
  }
  titles_row <- function(est_label) {
    r <- as.data.frame(as.list(c("Rank", "Expected rank", "Hospital", "Site code",
                                 "Sample size (N)", "Effective sample size (ESS)",
                                 est_label, "10%", "20%", "25%")),
                       stringsAsFactors = FALSE)
    names(r) <- col_names
    r
  }
  
  group_row <- blank_row("")
  group_row[["p10"]] <- "Probability in the top X% of performers"
  
  sus_tab <- top_block(sus)
  # the banners carry the period labels from the config, so the caption can
  # never disagree with the window the numbers came from
  parts <- list(group_row,
                blank_row(sprintf("Sustained performance (average waiting time, %s)",
                                  window_label)),
                titles_row("Standardised mean waiting time (days)"),
                sus_tab)
  if (!is.null(imp)) {
    imp_tab <- top_block(imp, improve = TRUE)
    parts <- c(parts,
               list(blank_row(sprintf("Improvement over the period (change, %s vs %s)",
                                      period_2_label, period_1_label)),
                    titles_row("Change in standardised mean waiting time (days)"),
                    imp_tab))
  }
  combined <- do.call(rbind, parts)
  
  # row positions, computed from the block sizes rather than assumed, so a short
  # block (fewer sites than top_n) cannot push the rules onto the wrong rows
  row_group   <- 1L
  row_sus_div <- 2L
  row_titles1 <- 3L
  sus_rows    <- row_titles1 + seq_len(nrow(sus_tab))
  has_imp     <- !is.null(imp)
  if (has_imp) {
    row_imp_div <- max(sus_rows) + 1L
    row_titles2 <- row_imp_div + 1L
    imp_rows    <- row_titles2 + seq_len(nrow(imp_tab))
  }
  
  ft <- flextable(combined)
  ft <- delete_part(ft, part = "header")   # every heading row is already in the body
  
  left_cols <- match(c("Rank", "ExpRank", "Hospital", "SiteCode", "N", "ESS",
                       "Estimate"), col_names)
  prob_cols <- match(c("p10", "p20", "p25"), col_names)
  ft <- merge_at(ft, i = row_group, j = left_cols, part = "body")
  ft <- merge_at(ft, i = row_group, j = prob_cols, part = "body")
  ft <- merge_at(ft, i = row_sus_div, j = seq_along(col_names), part = "body")
  if (has_imp)
    ft <- merge_at(ft, i = row_imp_div, j = seq_along(col_names), part = "body")
  
  bold_rows <- c(row_group, row_sus_div, row_titles1,
                 if (has_imp) c(row_imp_div, row_titles2))
  ft <- bold(ft, i = bold_rows, part = "body")
  ft <- align(ft, part = "body", align = "center")
  ft <- align(ft, j = match(c("Hospital", "SiteCode"), col_names),
              align = "left", part = "body")
  banner_rows <- c(row_sus_div, if (has_imp) row_imp_div)
  ft <- align(ft, i = banner_rows, align = "left", part = "body")
  
  rule  <- fp_border(color = "black",  width = 1)
  faint <- fp_border(color = "grey85", width = 0.5)
  ft <- border_outer(ft, border = rule, part = "body")
  for (i in c(row_group, row_sus_div, row_titles1, max(sus_rows),
              if (has_imp) c(row_imp_div, row_titles2)))
    ft <- hline(ft, i = i, border = rule, part = "body")
  
  # very faint rules between individual hospitals; interior rows only, since each
  # block's last row already carries a black rule or the outer frame
  if (length(sus_rows) > 1)
    ft <- hline(ft, i = sus_rows[-length(sus_rows)], border = faint, part = "body")
  if (has_imp && length(imp_rows) > 1)
    ft <- hline(ft, i = imp_rows[-length(imp_rows)], border = faint, part = "body")
  
  # faint dotted verticals bracketing the sample-size pair and the estimate, on
  # the structured rows only so they do not cut through the merged banners
  vrule <- fp_border(color = "grey60", width = 0.75, style = "dotted")
  vcols <- match(c("ExpRank", "SiteCode", "ESS", "Estimate"), col_names)
  structured <- setdiff(seq_len(nrow(combined)), c(row_group, banner_rows))
  ft <- vline(ft, i = structured, j = vcols, border = vrule, part = "body")
  
  ft <- fontsize(ft, size = 7, part = "all")
  ft <- padding(ft, padding.top = 1, padding.bottom = 1, part = "all")
  ft <- autofit(ft)
  ft <- width(ft, j = "Rank",     width = 0.42)
  ft <- width(ft, j = "ExpRank",  width = 0.62)
  ft <- width(ft, j = "Hospital", width = 1.7)
  ft <- width(ft, j = "SiteCode", width = 0.6)
  ft <- width(ft, j = "N",        width = 0.5)
  ft <- width(ft, j = "ESS",      width = 0.6)
  ft <- width(ft, j = "Estimate", width = 1.35)   # holds "mean (lo to hi)" on one line
  ft <- width(ft, j = c("p10", "p20", "p25"), width = 0.4)
  ft <- set_table_properties(ft, layout = "fixed", align = "left")
  
  caption_text <- paste(
    sprintf("Top %d hospitals by the shrinkage model, ordered by expected posterior rank",
            top_n),
    if (has_imp) ": sustained (upper block) and improvement (lower block)." else ".",
    "Rank is the position in that ordering; Expected rank is the posterior mean",
    "rank (1 = best). The estimate is the case-mix-standardised endoscopy-to-",
    "decision-to-treat time from the Bayesian shrinkage model, with its 95%",
    "posterior credible interval in brackets; the probabilities are the posterior",
    "probability the hospital sits in the top X% of performers.")
  footnote_text <- paste(
    "Sustained: standardised to the whole window's patient mix.",
    if (has_imp) paste(
      "Improvement: the change in standardised waiting time between",
      sprintf("%s and %s, both standardised to %s's patient mix -", period_1_label,
              period_2_label, period_1_label),
      "a negative value means faster, and 'top X%' means the most improved X% of",
      "hospitals.") else "")
  
  doc <- read_docx()
  doc <- body_set_default_section(doc, prop_section(
    page_size    = page_size(orient = "portrait"),
    page_margins = page_mar(top = 0.7, bottom = 0.7, left = 0.6, right = 0.6)))
  doc <- body_add_par(doc, caption_text, style = "Normal")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_fpar(doc, fpar(ftext(footnote_text, prop = fp_text(font.size = 8))))
  print(doc, target = path)
  invisible(path)
}

if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer",   quietly = TRUE)) {
  library(flextable); library(officer)
  docx_path <- file.path(out_dir, "top20_ranking.docx")
  ranking_docx(sus, imp, docx_path)
  cat("\nformatted ranking table written to", basename(docx_path),
      if (is.null(imp)) "(sustained block only)" else "(sustained and improvement)", "\n")
} else {
  cat("\nflextable / officer not installed - Word ranking table skipped",
      "(the CSVs above are unaffected).\n")
}

# candidate selection --------------------------------------------------------
# a candidate has at least prob_cut posterior probability of sitting in the top
# 20% on that estimand: the fastest fifth for sustained, the most improved fifth
# for improvement. The two lists are written separately and then combined, so a
# site can be flagged on one, the other, or both - being consistently fast and
# having got faster are different claims and worth being able to tell apart.
candidate_list <- function(d, estimand) {
  d %>%
    filter(p_top20 >= prob_cut) %>%
    arrange(exp_rank) %>%
    transmute(estimand  = estimand,
              site      = diag_hosp,
              exp_rank  = round(exp_rank, 1),
              post_mean = round(post_mean, 1),
              p_top20   = round(p_top20, 2))
}

candidates <- candidate_list(sus, "sustained")
write.csv(candidates, file.path(out_dir, "candidates.csv"), row.names = FALSE)
cat(sprintf("\nsustained candidates (P(top 20%%) >= %.2f): %d\n",
            prob_cut, nrow(candidates)))
cat("candidates written to candidates.csv\n")

if (!is.null(imp)) {
  candidates_improve <- candidate_list(imp, "improved")
  write.csv(candidates_improve, file.path(out_dir, "candidates_improve.csv"),
            row.names = FALSE)
  cat(sprintf("improved candidates (P(top 20%%) >= %.2f): %d\n",
              prob_cut, nrow(candidates_improve)))
  cat("candidates written to candidates_improve.csv\n")
  
  both <- bind_rows(candidates, candidates_improve)
  write.csv(both, file.path(out_dir, "candidates_all.csv"), row.names = FALSE)
  
  on_both <- intersect(candidates$site, candidates_improve$site)
  cat(sprintf("\nflagged on BOTH estimands: %d site(s)%s\n",
              length(on_both),
              if (length(on_both)) paste0(" - ", paste(on_both, collapse = ", ")) else ""))
  cat("combined list written to candidates_all.csv\n")
}
cat("04 complete. This is the last step of stage 06.\n")