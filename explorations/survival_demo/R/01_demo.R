# ------------------------------------------------------------------------------
# File:     explorations/survival_demo/R/01_demo.R
# Project:  Survival analysis demo - Kaplan-Meier + Cox PH on survival::lung
# Author:   [Instructor]
# Purpose:  A complete, heavily-commented R walk-through for someone new to
#           survival analysis:
#             (1) inspect the lung-cancer survival data + report censoring
#             (2) Kaplan-Meier curves by sex with a log-rank test
#             (3) Cox proportional-hazards model with age + ECOG + sex
#             (4) PH-assumption diagnostic (cox.zph)
#             (5) HR forest plot
#             (6) publication-style HR table
# Inputs:   survival::lung   (built-in: NCCTG lung-cancer trial; 228 patients)
# Outputs:  output/figures/km_by_sex.{pdf,png}
#           output/figures/hr_forest.{pdf,png}
#           output/figures/ph_diagnostic.{pdf,png}
#           output/tables/cox_hr_table.{tex,csv}
# Log:      explorations/survival_demo/logs/01_demo.log
# ------------------------------------------------------------------------------

if (getRversion() < "4.3.0") stop("Requires R >= 4.3.0; you have ", R.version.string)
options(warn = 1, scipen = 999, stringsAsFactors = FALSE)

source("R/_utils/paths.R")
source("R/_utils/logging.R")
source("R/_utils/theme_journal.R")

demo_dir <- proj_path("explorations", "survival_demo")
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
  library(ggplot2)
  library(patchwork)
  library(survival)         # Surv(), survfit(), coxph(), cox.zph(), survdiff()
  library(broom)            # tidy(coxph, exponentiate = TRUE)
})


# --- 1. Load + inspect the data ----------------------------------------------
# survival::lung is the NCCTG advanced-lung-cancer trial: 228 patients,
# followed up to death or censoring. Variables we care about:
#   time   - survival time in days
#   status - 1 = censored, 2 = dead    (note: this convention is non-standard;
#                                       we recode below to the canonical
#                                       0 = censored, 1 = event)
#   sex    - 1 = male, 2 = female
#   age    - years
#   ph.ecog - ECOG performance score: 0 (asymptomatic) to 5 (dead); higher = worse

data("lung", package = "survival")

dat <- lung %>%
  filter(!is.na(ph.ecog)) %>%                  # drop the 1 row missing ECOG
  mutate(
    event = as.integer(status == 2),           # 1 = died, 0 = censored
    sex   = factor(sex, levels = c(1, 2),
                   labels = c("Male", "Female")),
    ecog  = factor(ph.ecog, levels = c(0, 1, 2, 3),
                   labels = c("0 Asymptomatic",
                              "1 Symptomatic, ambulatory",
                              "2 In bed < 50% of day",
                              "3 In bed > 50% of day"))
  )

cat("\n*** Sample characteristics ***\n")
print(dat %>%
        summarise(
          N             = dplyr::n(),
          events        = sum(event),
          censored      = sum(1 - event),
          censoring_pct = round(100 * mean(1 - event), 1),
          median_age    = median(age, na.rm = TRUE),
          median_time   = median(time, na.rm = TRUE)
        ))

cat("\n*** Sex distribution ***\n")
print(dat %>% count(sex) %>%
        mutate(pct = round(100 * n / sum(n), 1)))

cat("\n*** ECOG performance distribution ***\n")
print(dat %>% count(ecog) %>%
        mutate(pct = round(100 * n / sum(n), 1)))


# --- 2. Kaplan-Meier curves by sex + log-rank test ---------------------------
# The Kaplan-Meier (KM) estimator gives a non-parametric estimate of the
# survival function S(t) = Pr(T > t). It treats each event time as a step
# down; censored observations are removed from the risk set but do not
# cause a step.
#
# `Surv(time, event)` constructs the outcome with right-censoring.
# `survfit(formula, data)` fits the KM; with `formula = Surv(...) ~ group`
# you get one curve per group.
#
# `survdiff()` runs the log-rank test of the null "all groups have the
# same survival function". chi-square statistic with df = number of
# groups - 1.

km_fit <- survfit(Surv(time, event) ~ sex, data = dat)
cat("\n*** Kaplan-Meier summary by sex ***\n")
print(summary(km_fit, times = c(180, 365, 540, 720)))     # 6, 12, 18, 24 mo

logrank <- survdiff(Surv(time, event) ~ sex, data = dat)
cat("\n*** Log-rank test (sex) ***\n")
print(logrank)
logrank_pval <- 1 - pchisq(logrank$chisq, df = length(logrank$n) - 1)
cat(sprintf("log-rank p-value = %.4g\n", logrank_pval))


# --- 3. KM figure with risk table (journal style) ----------------------------
# Build the Kaplan-Meier display directly from survfit output so every visual
# element is an explicit ggplot2 layer. The lower panel is a compact risk table
# drawn with geom_text and aligned to the same x-axis.

km_palette <- c(pal_journal[["blue"]], pal_journal[["teal"]])
sex_levels <- c("Male", "Female")

km_summary <- summary(km_fit)
km_df <- tibble(
  time = km_summary$time,
  surv = km_summary$surv,
  lower = km_summary$lower,
  upper = km_summary$upper,
  sex = sub("^sex=", "", km_summary$strata)
)

start_rows <- dat %>%
  count(sex, name = "n.risk") %>%
  transmute(time = 0, surv = 1, lower = 1, upper = 1, sex = as.character(sex))

km_plot_df <- bind_rows(start_rows, km_df) %>%
  mutate(sex = factor(sex, levels = sex_levels)) %>%
  arrange(sex, time)

median_df <- km_plot_df %>%
  filter(surv <= 0.5) %>%
  group_by(sex) %>%
  slice_min(time, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(sex, median = time)

risk_times <- seq(0, 720, by = 180)
risk_summary <- summary(km_fit, times = risk_times, extend = TRUE)
risk_df <- tibble(
  time = risk_summary$time,
  sex = factor(sub("^sex=", "", risk_summary$strata), levels = rev(sex_levels)),
  n_risk = risk_summary$n.risk
)

p_km_curve <- ggplot(km_plot_df, aes(x = time, y = surv, colour = sex)) +
  geom_step(linewidth = 0.8) +
  geom_step(aes(y = lower), linewidth = 0.35, alpha = 0.45, linetype = "dashed") +
  geom_step(aes(y = upper), linewidth = 0.35, alpha = 0.45, linetype = "dashed") +
  geom_segment(data = median_df,
               aes(x = 0, xend = median, y = 0.5, yend = 0.5, colour = sex),
               linewidth = 0.35, linetype = "dotted", inherit.aes = FALSE) +
  geom_segment(data = median_df,
               aes(x = median, xend = median, y = 0, yend = 0.5, colour = sex),
               linewidth = 0.35, linetype = "dotted", inherit.aes = FALSE) +
  annotate("text", x = 30, y = 0.12,
           label = sprintf("Log-rank p = %.4g", logrank_pval),
           hjust = 0, family = "serif", size = 3.3, colour = "grey20") +
  scale_colour_manual(values = setNames(km_palette, sex_levels), name = "Sex") +
  scale_x_continuous(breaks = risk_times, limits = c(0, 720),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) +
  labs(title    = "Kaplan-Meier survival curves by sex",
       subtitle = sprintf("NCCTG lung-cancer trial (N = %d, %d events / %d censored)",
                          nrow(dat), sum(dat$event), sum(1 - dat$event)),
       x        = NULL,
       y        = "Survival probability S(t)") +
  theme_journal(base_size = 11) +
  theme(legend.position = c(0.82, 0.82),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.title.x = element_blank())

p_risk <- ggplot(risk_df, aes(x = time, y = sex, label = n_risk, colour = sex)) +
  geom_text(family = "serif", size = 3.4, show.legend = FALSE) +
  scale_colour_manual(values = setNames(km_palette, sex_levels)) +
  scale_x_continuous(breaks = risk_times, limits = c(0, 720),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(title = "Number at risk",
       x = "Days from enrollment",
       y = NULL) +
  theme_journal(base_size = 10) +
  theme(plot.title = element_text(size = 10),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(face = "bold"))

p_km_combined <- p_km_curve / p_risk +
  plot_layout(heights = c(4, 1))

ggsave(file.path(fig_dir, "km_by_sex.pdf"),
       plot = p_km_combined, width = 7.5, height = 6.5)
ggsave(file.path(fig_dir, "km_by_sex.png"),
       plot = p_km_combined, width = 7.5, height = 6.5, dpi = 300)

# --- 4. Cox proportional-hazards model ---------------------------------------
# The Cox PH model assumes
#     h(t | X) = h0(t) * exp(beta' * X)
# where h0(t) is a non-parametric baseline hazard (estimated as a nuisance)
# and exp(beta_j) is the HAZARD RATIO for a one-unit increase in X_j,
# holding the rest fixed.
#
# Reporting: HR + 95% CI (NOT the raw log-HR). broom::tidy() does this in
# one line with `exponentiate = TRUE, conf.int = TRUE`.

cox_fit <- coxph(Surv(time, event) ~ sex + age + ecog, data = dat)
cat("\n*** Cox PH model ***\n")
print(summary(cox_fit))

hr_tbl <- broom::tidy(cox_fit,
                      exponentiate = TRUE,
                      conf.int = TRUE) %>%
  transmute(
    term,
    HR        = round(estimate, 3),
    CI_lo     = round(conf.low, 3),
    CI_hi     = round(conf.high, 3),
    z         = round(statistic, 2),
    p_value   = signif(p.value, 3)
  )

cat("\n*** HR table ***\n")
print(hr_tbl)

# Concordance index (C-index): probability that, in a random pair of
# patients, the one predicted to die first actually died first. 0.5 = chance;
# 1 = perfect ordering. Report it as a model-fit metric.
c_index <- summary(cox_fit)$concordance
cat(sprintf("\nC-index: %.3f (SE %.3f)\n",
            c_index[["C"]], c_index[["se(C)"]]))

# Pretty labels shared by the PH diagnostic and forest plot.
term_label <- c(
  "sexFemale"                     = "Female (vs Male)",
  "age"                           = "Age (per year)",
  "ecog1 Symptomatic, ambulatory" = "ECOG 1 (vs 0)",
  "ecog2 In bed < 50% of day"     = "ECOG 2 (vs 0)",
  "ecog3 In bed > 50% of day"     = "ECOG 3 (vs 0)"
)


# --- 5. Proportional-hazards diagnostic --------------------------------------
# Cox's PH assumption - that hazard ratios are CONSTANT over time - must
# be checked. cox.zph() correlates the scaled Schoenfeld residuals with
# (a function of) time; a small p means the HR for that covariate is
# actually time-varying.
#
# If you reject for any covariate, the fix is either:
#   (a) stratify by it: coxph(... + strata(covariate))
#   (b) allow a time-varying coefficient: coxph(... + tt(covariate))
# Do NOT just throw in a time*X interaction inside coxph() - that does
# not give you a time-varying coefficient; you need tt() or to expand
# the data into start-stop format via survSplit/tmerge.

zph <- cox.zph(cox_fit)
cat("\n*** Proportional-hazards test (cox.zph) ***\n")
print(zph)

# Plot the scaled Schoenfeld residuals directly with ggplot2. A flat smooth
# around zero supports the proportional-hazards assumption.
zph_df <- as.data.frame(zph$y) %>%
  mutate(transformed_time = zph$x) %>%
  tidyr::pivot_longer(cols = -transformed_time,
                      names_to = "term",
                      values_to = "scaled_residual") %>%
  mutate(term = factor(term, levels = names(term_label), labels = term_label))

p_zph <- ggplot(zph_df, aes(x = transformed_time, y = scaled_residual)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(colour = pal_journal[["navy"]], alpha = 0.45, size = 1.1) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = 0.75,
              colour = pal_journal[["blue"]], fill = pal_journal[["blue"]],
              alpha = 0.18, linewidth = 0.7) +
  facet_wrap(~ term, scales = "free_y", ncol = 2) +
  labs(title    = "Proportional-hazards diagnostic",
       subtitle = sprintf("Scaled Schoenfeld residual smooths from cox.zph; global p = %.4g",
                          zph$table["GLOBAL", "p"]),
       x        = "Transformed event time",
       y        = "Scaled Schoenfeld residual") +
  theme_journal(base_size = 10) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "ph_diagnostic.pdf"),
       plot = p_zph, width = 8, height = 6)
ggsave(file.path(fig_dir, "ph_diagnostic.png"),
       plot = p_zph, width = 8, height = 6, dpi = 300)


# --- 6. Hazard-ratio forest plot ---------------------------------------------
# A picture of the HR table: one row per covariate, HR + 95% CI on a log
# x-axis, dashed vertical line at HR = 1 (the "no effect" mark). Forest
# plots are the standard way clinical trials report Cox results.

forest_df <- hr_tbl %>%
  mutate(
    label = term_label[term],
    label = factor(label, levels = rev(term_label))
  ) %>%
  filter(!is.na(label)) %>%
  mutate(label_x = max(CI_hi) * 1.05)

p_hr <- ggplot(forest_df, aes(x = HR, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey25", linewidth = 0.5) +
  geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi),
                orientation = "y", width = 0.18,
                colour = pal_journal[["navy"]], linewidth = 0.6) +
  geom_point(size = 3, shape = 21, fill = pal_journal[["blue"]],
             colour = pal_journal[["navy"]], stroke = 0.9) +
  geom_text(aes(x = label_x,
                label = sprintf("%.2f (%.2f, %.2f)  p = %.3g",
                                HR, CI_lo, CI_hi, p_value)),
            hjust = 0, family = "serif", size = 3.5, colour = "grey15") +
  scale_x_log10(breaks = c(0.5, 1, 2, 5),
                limits = c(min(forest_df$CI_lo) * 0.9,
                           max(forest_df$CI_hi) * 5)) +
  labs(title    = "Cox proportional-hazards estimates",
       subtitle = sprintf("NCCTG lung-cancer trial (N = %d, %d events); C-index = %.3f",
                          nrow(dat), sum(dat$event), c_index[["C"]]),
       x        = "Hazard ratio (log scale, 95% CI)",
       y        = NULL,
       caption  = "Source: survival::lung. Dashed line = HR of 1 (no effect).") +
  theme_journal(base_size = 11)

ggsave(file.path(fig_dir, "hr_forest.pdf"),
       plot = p_hr, width = 8, height = 4.5)
ggsave(file.path(fig_dir, "hr_forest.png"),
       plot = p_hr, width = 8, height = 4.5, dpi = 300)


# --- 7. Publication-style HR table -------------------------------------------

write.csv(hr_tbl,
          file.path(tab_dir, "cox_hr_table.csv"),
          row.names = FALSE)

stars <- function(p) ifelse(p < 0.01, "***",
                     ifelse(p < 0.05, "**",
                     ifelse(p < 0.10, "*", "")))

tex_rows <- hr_tbl %>%
  mutate(label = term_label[term],
         hr_s  = sprintf("%.3f%s", HR, stars(p_value)),
         ci_s  = sprintf("(%.3f, %.3f)", CI_lo, CI_hi),
         p_s   = sprintf("%.3g", p_value))

tex_lines <- c(
  "\\begin{tabular}{lccc}",
  "  \\toprule",
  "  Covariate & HR & 95\\% CI & $p$ \\\\",
  "  \\midrule",
  paste0("  ", tex_rows$label, " & ", tex_rows$hr_s, " & ",
         tex_rows$ci_s, " & ", tex_rows$p_s, " \\\\"),
  "  \\midrule",
  sprintf(paste0("  \\multicolumn{4}{l}{\\small $N = %d$ \\quad ",
                 "Events = %d \\quad C-index = %.3f \\quad ",
                 "cox.zph global $p$ = %.3g} \\\\"),
          nrow(dat), sum(dat$event), c_index[["C"]], zph$table["GLOBAL", "p"]),
  "  \\bottomrule",
  "\\end{tabular}"
)
writeLines(tex_lines,
           file.path(tab_dir, "cox_hr_table.tex"))


# --- 8. Done -----------------------------------------------------------------

cat("\nPipeline finished. Inspect:\n")
cat("  log:     explorations/survival_demo/logs/01_demo.log\n")
cat("  figures: explorations/survival_demo/output/figures/\n")
cat("  tables:  explorations/survival_demo/output/tables/\n")

stop_log()
