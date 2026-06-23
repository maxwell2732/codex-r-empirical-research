# ------------------------------------------------------------------------------
# File:     explorations/staggered_did_demo/R/01_demo.R
# Project:  Staggered DiD demo — TWFE bias and heterogeneity-robust fixes
# Author:   [Instructor]
# Purpose:  Show, on a simulated staggered-treatment panel where the true ATT
#           is heterogeneous in cohort AND event time:
#             (1) why the naive TWFE estimator is biased
#             (2) how Sun-Abraham (fixest::sunab) recovers the truth
#             (3) how Callaway-Sant'Anna (did::att_gt + aggte) does the same
#                 via a different aggregation
#             (4) how Borusyak-Jaravel-Spiess (did2s) does it as a third check
#           The figure overlays all four estimators on one event-study plot
#           against the true effect path.
# Inputs:   (none — data simulated inside the script)
# Outputs:  output/figures/event_study_comparison.{pdf,png}
#           output/figures/data_visualisation.{pdf,png}
#           output/tables/pooled_att_comparison.{tex,csv}
# Log:      explorations/staggered_did_demo/logs/01_demo.log
# ------------------------------------------------------------------------------

if (getRversion() < "4.3.0") stop("Requires R >= 4.3.0; you have ", R.version.string)
options(warn = 1, scipen = 999, stringsAsFactors = FALSE)

source("R/_utils/paths.R")
source("R/_utils/logging.R")
source("R/_utils/theme_journal.R")

demo_dir <- proj_path("explorations", "staggered_did_demo")
log_dir <- file.path(demo_dir, "logs")
fig_dir <- file.path(demo_dir, "output", "figures")
tab_dir <- file.path(demo_dir, "output", "tables")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

start_log("01_demo", dir = log_dir)
on.exit(stop_log(), add = TRUE)

set.seed(20260512)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(fixest)
  library(did)            # Callaway-Sant'Anna
  library(did2s)          # Borusyak-Jaravel-Spiess two-stage DiD
  library(modelsummary)
})


# --- 1. Simulate a staggered-treatment panel --------------------------------
# Setup:
#   - N = 600 units, T = 10 periods
#   - 4 groups defined by treatment timing (cohort g):
#       g = 3 : treated starting in period 3 (early)
#       g = 5 : treated starting in period 5 (mid)
#       g = 7 : treated starting in period 7 (late)
#       g = 0 : never-treated (clean control)
#   - True ATT is HETEROGENEOUS in cohort and event time:
#       early-cohort effect grows fast (intercept 2.0, slope 0.5 per period)
#       mid-cohort effect grows slower (1.0, 0.3)
#       late-cohort effect grows slowest (0.5, 0.2)
#     So the estimator that simply averages "treated minus untreated" will mix
#     these heterogeneous effects with already-treated units used as controls
#     for not-yet-treated ones — that's the TWFE-staggered bias.
#   - Add unit and period fixed effects + Gaussian noise.

N <- 600
T <- 10
n_per_cohort <- N / 4
cohorts <- c(0, 3, 5, 7)              # 0 = never-treated

dat <- expand_grid(
  id     = seq_len(N),
  period = seq_len(T)
) %>%
  mutate(
    g        = rep(cohorts, each = n_per_cohort * T),  # cohort assignment
    treat    = as.integer(g != 0 & period >= g),        # current treatment
    e        = ifelse(g == 0, NA_integer_, period - g)  # event time (NA for never)
  )

# Per-cohort heterogeneous true effects.
true_eff <- function(g, e) {
  if (g == 0 || e < 0) return(0)
  switch(as.character(g),
         "3" = 2.0 + 0.5 * e,
         "5" = 1.0 + 0.3 * e,
         "7" = 0.5 + 0.2 * e,
         0)
}

# Outcome = unit FE + period FE + true heterogeneous treatment effect + noise.
unit_fe   <- rnorm(N, sd = 0.5)
period_fe <- rnorm(T, sd = 0.3)

dat <- dat %>%
  rowwise() %>%
  mutate(true_te = true_eff(g, ifelse(is.na(e), -1L, as.integer(e)))) %>%
  ungroup() %>%
  mutate(
    y = unit_fe[id] + period_fe[period] + true_te + rnorm(n(), sd = 0.8)
  )

cat("\n*** Panel structure ***\n")
print(
  dat %>%
    group_by(g) %>%
    summarise(units      = n_distinct(id),
              periods    = n_distinct(period),
              n_treated  = sum(treat),
              .groups = "drop")
)

# Truth we want every estimator to recover (event-time average ATT).
# Average only over cohorts that actually HAVE observations at this event time:
# with T = 10 and treatment at g in {3, 5, 7}, cohort 7 only reaches e = 3,
# cohort 5 reaches e = 5, cohort 3 reaches e = 7. The estimators average over
# the cohort set that survives at each e — the truth line must do the same to
# be a fair benchmark.
truth_event <- dat %>%
  filter(!is.na(e), e >= -2, e <= 5) %>%
  distinct(g, e) %>%
  rowwise() %>%
  mutate(true_te = true_eff(g, e)) %>%
  ungroup() %>%
  group_by(e) %>%
  summarise(true_atte = mean(true_te), .groups = "drop")

cat("\n*** True average ATT(e) by event time (target) ***\n")
print(truth_event)


# --- 2. Visualise the data ---------------------------------------------------
# A spaghetti plot of y over time, coloured by cohort, with vertical dashes at
# each cohort's treatment date. This makes the heterogeneity visible to the
# eye before any model is fit.

cohort_label <- function(g) ifelse(g == 0, "Never-treated", paste0("Treated at t = ", g))

dat_means <- dat %>%
  group_by(g, period) %>%
  summarise(y = mean(y), .groups = "drop") %>%
  mutate(cohort = factor(cohort_label(g),
                         levels = c("Never-treated",
                                    "Treated at t = 3",
                                    "Treated at t = 5",
                                    "Treated at t = 7")))

cohort_palette <- c("Never-treated"     = pal_journal[["slate"]],
                    "Treated at t = 3"  = pal_journal[["blue"]],
                    "Treated at t = 5"  = pal_journal[["teal"]],
                    "Treated at t = 7"  = pal_journal[["lilac"]])

p_data <- ggplot(dat_means, aes(x = period, y = y, colour = cohort)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4, shape = 21, fill = "white", stroke = 0.8) +
  geom_vline(xintercept = c(3, 5, 7), linetype = "dotted", colour = "grey60") +
  scale_x_continuous(breaks = 1:T, expand = expansion(mult = 0.02)) +
  scale_colour_manual(values = cohort_palette, name = NULL) +
  labs(title    = "Simulated panel: cohort-mean outcomes over time",
       subtitle = "Treatment dates marked at t = 3, 5, 7; never-treated as the slate-grey baseline",
       x        = "Period",
       y        = "Mean outcome",
       caption  = "Source: simulated; true ATT heterogeneous in cohort and event time.") +
  theme_journal(base_size = 11) +
  theme(legend.position = "right")

ggsave(proj_path("explorations", "staggered_did_demo", "output", "figures",
                 "data_visualisation.pdf"),
       plot = p_data, width = 7, height = 4)
ggsave(proj_path("explorations", "staggered_did_demo", "output", "figures",
                 "data_visualisation.png"),
       plot = p_data, width = 7, height = 4, dpi = 300)


# --- 3. Estimator A: naive TWFE event-study --------------------------------
# Build leads and lags of treatment status by event time. The naive TWFE event
# study regresses y on event-time dummies plus unit and period FE. Under
# homogeneous effects this would be unbiased; under heterogeneous-in-cohort
# effects, the implicit "controls" include already-treated units, biasing
# the ATT estimates toward zero (or even flipping sign in extreme cases).

# fixest::i() expands a factor variable into dummies and lets you set a
# reference category. We use e (event time) with reference = -1 (the period
# right before treatment — the standard normalisation).
dat_es <- dat %>%
  mutate(e_safe = ifelse(is.na(e), -1000L, e))    # never-treated to deep negative

m_twfe <- feols(
  y ~ i(e_safe, ref = c(-1, -1000)) | id + period,
  cluster = ~id,
  data    = dat_es
)
cat("\n*** Estimator A: naive TWFE event-study coefficients ***\n")
print(summary(m_twfe))


# --- 4. Estimator B: Sun-Abraham (interaction-weighted, native in fixest) --
# Sun-Abraham (2021) re-weights cohort-period interaction terms so that the
# pooled event-study coefficient reflects a properly-weighted average ATT.
# fixest implements it via the sunab() helper:
#   sunab(cohort_var, period_var)
# Cohort 0 (never-treated) is automatically used as the control.

m_sa <- feols(
  y ~ sunab(g, period) | id + period,
  cluster = ~id,
  data    = dat
)
cat("\n*** Estimator B: Sun-Abraham (sunab) ***\n")
print(summary(m_sa))


# --- 5. Estimator C: Callaway-Sant'Anna (did::att_gt + aggte) ---------------
# CS21 estimates ATT(g, t) for every (cohort, calendar-time) pair using the
# never-treated as the control group, then aggregates them.
# - att_gt() : per (g, t) ATTs
# - aggte(type = "dynamic") : event-study aggregation
# - aggte(type = "simple")  : single pooled ATT

cs <- att_gt(
  yname  = "y",
  tname  = "period",
  idname = "id",
  gname  = "g",
  data   = dat,
  control_group = "nevertreated",
  panel  = TRUE
)
cs_dyn <- aggte(cs, type = "dynamic", min_e = -2, max_e = 5)
cs_simple <- aggte(cs, type = "simple")
cat("\n*** Estimator C: Callaway-Sant'Anna dynamic (event-study) ***\n")
print(summary(cs_dyn))
cat("\n*** Estimator C: Callaway-Sant'Anna pooled simple ATT ***\n")
print(summary(cs_simple))


# --- 6. Estimator D: Borusyak-Jaravel-Spiess (did2s) ------------------------
# BJS / Gardner: a two-stage estimator. Stage 1 absorbs FE on the untreated
# observations only; stage 2 regresses residuals on event-time dummies.
# event_study() is a thin wrapper; we point it at our event-time variable.

bjs <- did2s(
  data       = dat_es,
  yname      = "y",
  first_stage = ~ 0 | id + period,
  second_stage = ~ i(e_safe, ref = c(-1, -1000)),
  treatment  = "treat",
  cluster_var = "id"
)
cat("\n*** Estimator D: Borusyak-Jaravel-Spiess (did2s) ***\n")
print(summary(bjs))


# --- 7. Stitch the four estimators into one tidy event-study frame ---------
# Each estimator gives a coefficient at every observed event time. We tidy
# them into a long frame with (estimator, e, estimate, se), then plot one
# line per estimator together with the truth.

tidy_fixest_es <- function(model, est_label) {
  co <- as.data.frame(coeftable(model))
  co$term <- rownames(co)
  co %>%
    filter(grepl("e_safe::", term)) %>%
    mutate(
      e        = as.integer(sub("e_safe::", "", term)),
      estimate = Estimate,
      se       = `Std. Error`,
      estimator = est_label
    ) %>%
    select(estimator, e, estimate, se)
}

twfe_es <- tidy_fixest_es(m_twfe, "TWFE (naive)")
bjs_es  <- tidy_fixest_es(bjs,    "Borusyak-Jaravel-Spiess")

# Sun-Abraham via fixest already produces per-event-time aggregated effects
# in the default summary(). Subtle gotcha: `summary(m_sa, agg = "ATT")`
# would collapse to a single pooled-ATT row, NOT per-event-time. Plain
# `summary(m_sa)` keeps the event-time rows labelled `period::-2`, `period::0`,
# `period::1`, etc.
sa_co <- as.data.frame(coeftable(summary(m_sa)))
sa_co$term <- rownames(sa_co)
sa_es <- sa_co %>%
  filter(grepl("^period::", term)) %>%
  mutate(
    e         = as.integer(sub("period::", "", term)),
    estimate  = Estimate,
    se        = `Std. Error`,
    estimator = "Sun-Abraham"
  ) %>%
  select(estimator, e, estimate, se)

cs_es <- tibble(
  estimator = "Callaway-Sant'Anna",
  e         = cs_dyn$egt,
  estimate  = cs_dyn$att.egt,
  se        = cs_dyn$se.egt
)

es_long <- bind_rows(twfe_es, sa_es, cs_es, bjs_es) %>%
  filter(e >= -2, e <= 5) %>%
  mutate(
    estimator = factor(estimator,
                       levels = c("TWFE (naive)",
                                  "Sun-Abraham",
                                  "Callaway-Sant'Anna",
                                  "Borusyak-Jaravel-Spiess")),
    ci_lo = estimate - 1.96 * se,
    ci_hi = estimate + 1.96 * se
  )

cat("\n*** Tidy event-study frame (long) ***\n")
print(es_long)


# --- 8. Event-study plot: four estimators vs truth -------------------------

est_palette <- c("TWFE (naive)"            = pal_journal[["slate"]],
                 "Sun-Abraham"             = pal_journal[["blue"]],
                 "Callaway-Sant'Anna"      = pal_journal[["teal"]],
                 "Borusyak-Jaravel-Spiess" = pal_journal[["lilac"]])

p_es <- ggplot(es_long, aes(x = e, y = estimate, colour = estimator)) +
  # Truth as a thin black step-line.
  geom_line(data = truth_event,
            aes(x = e, y = true_atte),
            inherit.aes = FALSE,
            colour = "grey15", linewidth = 0.5, linetype = "dashed") +
  geom_point(data = truth_event,
             aes(x = e, y = true_atte),
             inherit.aes = FALSE,
             colour = "grey15", size = 1.6) +
  # Per-estimator point + error bar, dodged so they don't overlap.
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey60") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                position = position_dodge(width = 0.4),
                width = 0.18, linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.4),
             size = 2.2, shape = 21, fill = "white", stroke = 0.8) +
  scale_colour_manual(values = est_palette, name = "Estimator") +
  scale_x_continuous(breaks = -2:5) +
  labs(title    = "Event-study estimates vs truth",
       subtitle = paste("Dashed black line + black dots = true ATT(e);",
                        "coloured points = estimator with 95% CI"),
       x        = "Event time (periods relative to treatment)",
       y        = "Estimated ATT(e)",
       caption  = "Source: simulated panel (N = 600, T = 10).") +
  theme_journal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(proj_path("explorations", "staggered_did_demo", "output", "figures",
                 "event_study_comparison.pdf"),
       plot = p_es, width = 7.5, height = 5.2)
ggsave(proj_path("explorations", "staggered_did_demo", "output", "figures",
                 "event_study_comparison.png"),
       plot = p_es, width = 7.5, height = 5.2, dpi = 300)


# --- 9. Pooled-ATT comparison table -----------------------------------------
# The single-number summary each estimator reports for "the" treatment effect.
# True pooled ATT (averaged over all post-treatment cohort-period cells):
truth_pooled <- dat %>%
  filter(g != 0, period >= g) %>%
  rowwise() %>%
  mutate(true_te = true_eff(g, period - g)) %>%
  ungroup() %>%
  summarise(truth = mean(true_te)) %>%
  pull(truth)

cat(sprintf("\n*** True pooled ATT (post-treatment average): %.3f ***\n",
            truth_pooled))

# Naive TWFE pooled = coefficient on the binary `treat` indicator in a
# separate fit (the event-study m_twfe has dummies, not a single coefficient).
m_twfe_pooled <- feols(y ~ treat | id + period, cluster = ~id, data = dat)

# Sun-Abraham pooled = the "ATT" aggregation of sunab.
sa_pooled <- summary(m_sa, agg = "att")

# BJS pooled: a SEPARATE did2s fit with the binary `treat` as the
# second-stage regressor. The event-study `bjs` model above has one
# coefficient per event time; coef(bjs)[1] would be the e=-2 placebo
# lead, NOT a pooled ATT — easy mistake to make.
bjs_pooled <- did2s(
  data        = dat,
  yname       = "y",
  first_stage = ~ 0 | id + period,
  second_stage = ~ treat,
  treatment   = "treat",
  cluster_var = "id"
)

pooled_tbl <- tibble(
  Estimator = c("Truth",
                "TWFE (naive)",
                "Sun-Abraham",
                "Callaway-Sant'Anna",
                "Borusyak-Jaravel-Spiess"),
  Pooled_ATT = c(truth_pooled,
                 coef(m_twfe_pooled)[["treat"]],
                 coef(sa_pooled)[1],
                 cs_simple$overall.att,
                 coef(bjs_pooled)[["treat"]]),
  SE         = c(NA,
                 se(m_twfe_pooled)[["treat"]],
                 se(sa_pooled)[1],
                 cs_simple$overall.se,
                 se(bjs_pooled)[["treat"]])
) %>%
  mutate(across(c(Pooled_ATT, SE), ~ round(.x, 3)),
         Bias = round(Pooled_ATT - truth_pooled, 3))

cat("\n*** Pooled-ATT comparison ***\n")
print(pooled_tbl)

write.csv(pooled_tbl,
          file.path(tab_dir, "pooled_att_comparison.csv"),
          row.names = FALSE)

# A modelsummary-rendered LaTeX table grouping the heterogeneity-robust models
# side-by-side (truth and naive TWFE rendered separately as a CSV column above).
ms_models <- list("TWFE (naive)" = m_twfe_pooled,
                  "Sun-Abraham"  = m_sa,
                  "BJS"          = bjs_pooled)

modelsummary(
  ms_models,
  output = file.path(tab_dir, "pooled_att_comparison.tex"),
  stars  = c("*" = 0.10, "**" = 0.05, "***" = 0.01),
  fmt    = 3,
  gof_omit = "AIC|BIC|Log.Lik|RMSE|Std.Errors|FE",
  title  = sprintf("Pooled ATT estimators (true value = %.3f)", truth_pooled),
  notes  = c("Cluster-robust SE at unit level (cluster = ~id).",
             "Heterogeneity-robust estimators (Sun-Abraham, BJS) recover the truth;",
             "the naive TWFE estimator is biased toward zero by the staggered design.")
)


# --- 10. Done ---------------------------------------------------------------

cat("\nPipeline finished. Inspect:\n")
cat("  log:     explorations/staggered_did_demo/logs/01_demo.log\n")
cat("  figures: explorations/staggered_did_demo/output/figures/\n")
cat("  tables:  explorations/staggered_did_demo/output/tables/\n")

stop_log()
