# REDCap Retention Analytics
### Survival Analysis of Participant Attrition in Multi-Site Clinical Trials

**Author:** Romario Joseph, MPH · BU SPH (Epidemiology & Biostatistics)
**Stack:** R 4.x · tidyverse · survival · survminer · gtsummary · ggplot2 · Quarto
**Data:** Synthetic REDCap export modeled on a multi-site recruitment registry (mirrors my work at the Broad Institute).

---

## Epidemiological Objective

In multi-site clinical research, participant attrition is not a nuisance — it is a *time-dependent failure process* that silently erodes statistical power, biases treatment-effect estimates, and corrupts equity of representation when drop-out is differential across demographic strata. Conventional REDCap reporting tells investigators **who** has left the study; it cannot tell them **when** drop-out concentrates along the visit timeline, **which covariates** elevate the instantaneous hazard of attrition, or **whether** the hazard structure is proportional across sites and sub-populations. This repository reframes participant retention as a **time-to-event problem** and asks three coupled questions:

1. **Time-to-attrition.** What is the probability that a participant remains active at time *t* after enrollment, and at which visit does the marginal hazard of drop-out peak?
2. **Covariate-adjusted hazard.** Holding site, age, sex, race/ethnicity, preferred language, and visit burden fixed, which covariates are independently associated with elevated hazard of attrition?
3. **Equity-stratified retention.** Are retention curves statistically separable across demographic strata, and does the proportional-hazards assumption hold within each stratum?

The applied motivation is direct: in production, this exact workflow drove a **25% reduction in participant drop-off** at the Broad Institute by surfacing visit-2 as the dominant inflection point in the hazard function and enabling a targeted reminder intervention.

---

## Methodological Framework

This repository is explicitly designed as a **survival analysis showcase** anchored in semi-parametric regression theory. The biostatistical models implemented are:

- **Kaplan–Meier (KM) product-limit estimator.** Non-parametric estimation of the marginal survival function S(t) = Pr(T > t), where T is time from enrollment to drop-out and the time origin is the enrollment date. Censoring is administrative (active participants at the end of follow-up) and treated as non-informative under the standard independent-censoring assumption, which is documented and stress-tested via a sensitivity analysis that re-codes the latest censoring tertile as events.
- **Log-rank test (Mantel–Haenszel).** Used to test the global null of equal survival functions across demographic strata (site, race/ethnicity, language, age band). The log-rank statistic is reported alongside stratified KM curves so the visual and inferential evidence are co-located.
- **Cox proportional hazards (PH) regression.** A semi-parametric multivariable model of the form h(t │ X) = h₀(t) · exp(βᵀX), where the baseline hazard h₀(t) is left unspecified and only the log-hazard ratios β are estimated by partial likelihood. Predictors include site (factor), age (continuous, centered), sex, race/ethnicity (OMB categories), preferred language, and a *visit burden* index (cumulative scheduled visits in the prior 30 days).
- **Diagnostics for the proportional-hazards assumption.** Because the entire interpretive validity of a Cox model collapses if proportionality fails, the pipeline runs three explicit checks:
  1. **Schoenfeld residual test.** Global and per-covariate χ² tests of the null that scaled Schoenfeld residuals are uncorrelated with a transformation of time (`survival::cox.zph`). A non-significant global p-value is required for the headline model to stand.
  2. **Log–log survival plots.** For each categorical covariate, log(−log S(t)) is plotted against log(t); parallel curves are taken as visual support for proportionality.
  3. **Time-by-covariate interaction terms.** Where Schoenfeld diagnostics flag a violation, the offending covariate is re-fit with an explicit `tt()` time-varying coefficient and the time-stratified hazard ratios are reported in the appendix rather than buried.
- **Bootstrapped equity rate ratios.** Stratum-specific retention rate ratios are accompanied by 1,000-replicate non-parametric bootstrap 95% CIs to avoid asymptotic normal-approximation failure in small demographic cells.

This combination — KM for the marginal picture, log-rank for stratified inference, and Cox PH with formal proportionality diagnostics for adjusted hazard ratios — mirrors the standard reporting template of *JAMA*-tier survival analyses and explicitly signals competence in **semi-parametric regression theory** rather than mere REDCap operationalization.

---

## Data Architecture

The pipeline is a five-stage DAG, with each R script owning a single, testable transformation. REDCap exports are messy by design — they are wide for instruments, long for repeating events, and mix coded factors with raw text — so the ingest layer carries most of the data-engineering burden.

**Stage 1 — REDCap ingest (`R/01_ingest_redcap.R`).** Reads both the long-format event log and the wide-format instrument exports natively (`REDCapR` API path supported, with CSV fallback for the synthetic case). Coded factor levels are restored from the REDCap data dictionary so no downstream script depends on numeric codes. A typed manifest of row counts, primary keys, and event-name uniqueness is written to `logs/ingest_manifest.json` and any deviation halts the pipeline.

**Stage 2 — Per-participant timeline construction (`R/02_build_retention.R`).** Each participant’s record is collapsed into a *survival object* with three fields: enrollment date (time origin), event date (first missed scheduled visit beyond grace window, or last completed visit + administrative end), and a 0/1 event indicator. Grace windows are configurable (default: 14 days). The construction explicitly handles four data-engineering edge cases that distinguish a defensible survival dataset from a brittle one:

- **Missing visit dates.** A scheduled visit without a completion record is treated as a candidate event only if (a) the scheduled date is more than *grace_window* days in the past and (b) no later visit in the protocol has been completed; otherwise the visit is flagged as *late but not lost*.
- **Out-of-order events.** Participants occasionally complete visit 3 before visit 2 (real-world site behavior). The timeline collapses by *protocol visit number*, not by chronological order, and the discrepancy is logged.
- **Re-engagement.** A participant who returns after the grace window is coded as left-censored at the first return and treated as a competing process; the headline analysis treats only the *first* attrition event, with a sensitivity analysis using a recurrent-event Andersen–Gill formulation.
- **Administrative censoring.** Active participants at the data freeze are right-censored at the freeze date.

**Stage 3 — KPI tables (`R/03_kpis.R`).** Emits enrollment rate, visit-1 completion, 30-day retention, KM-estimated median time-to-attrition, and the equity retention ratio as long-format Tableau-ready CSVs.

**Stage 4 — Survival modeling (`R/04_survival.R`).** Fits the KM estimator, runs the log-rank test, fits the Cox PH model, and emits `cox.zph` diagnostics, log–log plots, and time-varying-coefficient appendices into `outputs/diagnostics/`. The script fails loudly if the global Schoenfeld p < 0.05 unless an explicit `allow_ph_violation = TRUE` flag is set in the config, so a reviewer can trust that any committed headline hazard ratio has passed proportionality screening.

**Stage 5 — Equity stratification (`R/05_equity.R`).** Computes retention rate ratios across demographic strata with 1,000-replicate bootstrap CIs, with the lowest-vs-highest stratum ratio surfaced as the headline equity statistic.

```
redcap-retention-analytics/
├── README.md
├── R/
│   ├── 01_ingest_redcap.R       ← read & normalize REDCap export (long + wide)
│   ├── 02_build_retention.R     ← per-participant survival object
│   ├── 03_kpis.R                ← retention KPI tables
│   ├── 04_survival.R            ← KM + log-rank + Cox PH + PH diagnostics
│   └── 05_equity.R              ← stratified retention ratios + bootstrap CIs
├── data/
│   ├── synthetic/               ← fake REDCap export (committed)
│   └── README.md
├── reports/
│   └── retention_dashboard.qmd  ← Quarto dashboard
├── outputs/
│   └── diagnostics/             ← cox.zph tables, log–log plots
└── LICENSE
```

---

## KPI definitions

| KPI | Definition |
|-----|------------|
| Enrollment rate | n_enrolled / n_screened |
| Visit-1 completion | n_visit1_complete / n_enrolled |
| 30-day retention | participants still active 30 days after enrollment |
| Median time-to-dropoff | Kaplan–Meier point estimate of the median of T, censored at end of follow-up |
| Equity retention ratio | retention rate in lowest stratum ÷ highest stratum (bootstrap 95% CI) |

---

## Quick start

```r
install.packages(c("tidyverse","survival","survminer","gtsummary","REDCapR","ggplot2"))
source("R/01_ingest_redcap.R")    # data/synthetic/ → tidy tibbles
source("R/02_build_retention.R")  # per-participant survival object
source("R/03_kpis.R")             # KPI summary
source("R/04_survival.R")         # KM + log-rank + Cox PH + cox.zph diagnostics
source("R/05_equity.R")           # retention rate ratios + bootstrap CIs

quarto::quarto_render("reports/retention_dashboard.qmd")
```

---

## Why this matters

In production, the same toolkit:

- Identified **visit 2** as the highest-hazard timepoint via the KM survival curve and confirmed it as an independent risk factor in the Cox model → added a 7-day reminder → retention up 25%.
- Surfaced site-level variance hidden by aggregate reporting → targeted operational support to underperforming sites.
- Documented equity gaps in retention across language preference, enabling translated outreach.

---

## Disclaimer

All data committed here is synthetic. Not derived from any real REDCap project.

## Contact

Romario Joseph · rjoseph3@bu.edu · [LinkedIn](https://www.linkedin.com/in/romariojosephpublichealth/)
