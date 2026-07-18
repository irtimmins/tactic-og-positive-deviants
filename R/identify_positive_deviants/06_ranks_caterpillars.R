# 06  posterior ranks, caterpillar plot, and positive-deviant candidates
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

source("R/identify_positive_deviants/01_config.R")

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
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  caterpillar <- function(d, ylab, title, highlight = NA) {
    d <- d %>% arrange(exp_rank) %>% mutate(rank_order = row_number())
    hl <- if (!is.na(highlight)) d$diag_hosp == highlight else rep(FALSE, nrow(d))
    p <- ggplot(d, aes(rank_order, post_mean)) +
      geom_hline(yintercept = mean(d$post_mean), colour = "grey50") +
      geom_linerange(aes(ymin = ci_lo, ymax = ci_hi), colour = "grey70") +
      geom_point(size = 0.9) +
      labs(x = "Hospital (ordered by rank)", y = ylab, title = title) +
      theme_bw() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    if (any(hl)) p <- p + geom_point(data = d[hl, ], colour = "firebrick", size = 2)
    p
  }
  ggsave(file.path(out_dir, "caterpillar_sustained.pdf"),
         caterpillar(sus, "Standardised days, endoscopy to decision-to-treat",
                     "Sustained performance", highlight_hosp),
         width = 7, height = 6)
  cat("caterpillar plot written to caterpillar_sustained.pdf\n")
} else {
  cat("ggplot2 not installed - caterpillar plot skipped (CSV written).\n")
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
cat("06 complete.\n")
