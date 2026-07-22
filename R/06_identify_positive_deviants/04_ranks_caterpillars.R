# 04  posterior ranks, caterpillar plot, and positive-deviant candidates
# -----------------------------------------------------------------------------
# Turn the shrinkage fit into shrunk site means with credible intervals, expected
# ranks, and the probability of sitting in the fastest X%. Draws the caterpillar
# plot, writes the ranking table, and selects the sustained positive-deviant
# candidates. Shorter waits are better, so rank 1 is the best performer.
#
# Reads : stan_sustained.rds
# Writes: ranks_sustained.csv, caterpillar_sustained.pdf, top_ranking.csv,
#         candidates.csv

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
                     "Sustained performance", from_zero = TRUE,
                     highlight = highlight_hosp),
         width = 7, height = 6)
  cat("caterpillar plot written to caterpillar_sustained.pdf\n")
} else {
  cat("ggplot2 / scales not installed - caterpillar plot skipped (CSV written).\n")
}

# ranking table --------------------------------------------------------------
# one row per site, fastest first: shrunk estimate, effective sample size, and
# the posterior probability of sitting in the top 10 / 20 / 25 percent.
tab <- sus %>%
  arrange(exp_rank) %>%
  transmute(
    rank       = row_number(),
    exp_rank   = round(exp_rank, 1),
    site       = diag_hosp,
    n          = n,
    ess        = round(n_eff, 0),
    estimate   = sprintf("%.1f (%.1f to %.1f)", post_mean, ci_lo, ci_hi),
    p_top10    = round(p_top10, 2),
    p_top20    = round(p_top20, 2),
    p_top25    = round(p_top25, 2))
write.csv(tab, file.path(out_dir, "top_ranking.csv"), row.names = FALSE)
cat("\nfastest 10 sites:\n")
print(as.data.frame(head(tab, 10)), row.names = FALSE)

# candidate selection --------------------------------------------------------
# a sustained positive-deviant candidate has at least prob_cut posterior
# probability of sitting in the fastest 20% of sites.
sus_cand <- sus %>% filter(p_top20 >= prob_cut) %>% pull(diag_hosp)
cat(sprintf("\nsustained candidates (P(top 20%%) >= %.2f): %d\n",
            prob_cut, length(sus_cand)))

candidates <- sus %>%
  filter(diag_hosp %in% sus_cand) %>%
  arrange(exp_rank) %>%
  transmute(site = diag_hosp, exp_rank = round(exp_rank, 1),
            post_mean = round(post_mean, 1), p_top20 = round(p_top20, 2))
write.csv(candidates, file.path(out_dir, "candidates.csv"), row.names = FALSE)
cat("candidates written to candidates.csv\n")
cat("04 complete. This is the last step of stage 06.\n")