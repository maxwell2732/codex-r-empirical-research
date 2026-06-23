# Exploration: Survival analysis demo (KM + Cox PH on survival::lung)

A self-contained R demo of standard survival analysis on a real
clinical-trial dataset: the NCCTG advanced lung-cancer trial (228 patients
followed to death or censoring).

## What it shows

1. **Sample characteristics** - N, event count, censoring rate, basic demographics.
2. **Kaplan-Meier curves by sex** with 95 % CI ribbons, median-survival
   reference lines, log-rank test, and a risk table below.
3. **Cox proportional-hazards model** for sex + age + ECOG performance score,
   reporting hazard ratios (HR) with 95 % CI - NOT raw log-HRs.
4. **Proportional-hazards diagnostic** via `cox.zph` (Schoenfeld residuals
   correlated with time) for every covariate plus the global test.
5. **HR forest plot** on a log x-axis with HR / CI / p annotated to the right.
6. **Publication-style HR table** (CSV + LaTeX), with N / events / C-index /
   global PH-test p in the footer.

## How to replicate

```bash
Rscript scripts/setup_r.R                # one-time; survival + ggplot2 stack
bash scripts/run_r.sh explorations/survival_demo/R/01_demo.R
```

## Headline results

- N = 227 (1 row dropped for missing ECOG); 164 events / 63 censored (27.8 % censoring)
- Log-rank test by sex: chi-square = 9.69, p = 0.0016 - female survival significantly longer
- Cox PH model HRs (95 % CI):

  | Covariate | HR | 95 % CI | p |
  |---|---|---|---|
  | Female (vs Male) | **0.58** | 0.42-0.81 | 0.001 |
  | Age (per year) | 1.01 | 0.99-1.03 | 0.25 |
  | ECOG 1 (vs 0) | 1.51 | 1.02-2.23 | 0.04 |
  | ECOG 2 (vs 0) | **2.47** | 1.58-3.86 | < 0.001 |
  | ECOG 3 (vs 0) | 7.06 | 0.94-53.1 | 0.06 |

- C-index = 0.637 (modest discrimination, typical for an advanced-cancer cohort
  with limited covariates)
- cox.zph global test p = 0.17 - proportional-hazards assumption not rejected
The textbook clinical reading: **female sex is protective** (hazard ratio ~0.58),
**ECOG performance score is the strongest single mortality predictor**, age
adds no independent signal in this small sample.

## Outputs

| Path | Contents |
|---|---|
| `output/figures/km_by_sex.{pdf,png}` | KM curves + 95 % CI + log-rank + risk table |
| `output/figures/hr_forest.{pdf,png}` | HR + 95 % CI forest plot on log x-axis |
| `output/figures/ph_diagnostic.{pdf,png}` | Schoenfeld residuals per covariate from `cox.zph` |
| `output/tables/cox_hr_table.csv` | HR / CI / z / p per covariate |
| `output/tables/cox_hr_table.tex` | LaTeX version with N / events / C-index footer |

## Status

Exploration (60 / 100 quality threshold per `.claude/rules/exploration-fast-track.md`).
This is a clean, didactic walkthrough - not a re-analysis of the trial.

## Scope deliberately omitted

- Time-varying covariates (`tt(covariate)` or `tmerge` start-stop reshape).
  Would extend if a covariate failed the cox.zph test.
- Stratified Cox model when the PH assumption fails (alternative fix to time-varying).
- Functional form of continuous covariates (martingale residuals via `ggcoxfunctional`)
  and influential observations (`ggcoxdiagnostics`). Both are next-week extensions.
- Competing risks (`cmprsk` / `survival::finegray`); we treat death as the
  only event, which is correct for overall survival.
- Multiple imputation for missing data - we drop the 1 NA-ECOG row instead.

## Notes for the curious

Two API gotchas worth remembering:

1. **`survival::lung$status` is coded `1 = censored, 2 = dead`**, NOT the
   canonical `0 = censored, 1 = event`. The script recodes via
   `event = as.integer(status == 2)`. Forgetting this flips every HR.
2. The KM curve and risk table are built directly from `survfit` summaries
   with `ggplot2` layers, then composed via `patchwork::/`. This keeps the
   exported object explicit and avoids wrapper-specific save behavior.
