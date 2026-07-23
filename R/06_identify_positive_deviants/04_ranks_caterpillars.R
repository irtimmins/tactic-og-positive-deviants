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
  axis_title_size <- 10      # axis titles
  axis_text_size  <- 9       # axis tick labels
  cat_ci_alpha    <- 0.30    # credible-interval opacity - lighter than the points
  cat_ci_lwd      <- 0.25    # credible-interval line width
  cat_pt_size     <- 0.8     # point size
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