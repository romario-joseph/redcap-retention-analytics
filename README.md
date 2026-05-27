# REDCap Retention Analytics
### Survival analysis of participant attrition in multi-site clinical trials

**Author:** Romario Joseph, MPH (BU SPH, Epidemiology & Biostatistics)
**Stack:** R 4.x, tidyverse, survival, survminer, gtsummary, ggplot2, Quarto
**Data:** Synthetic REDCap export modeled on a multi-site recruitment registry. The structure mirrors the projects I worked on at the Broad Institute.

This is the toolkit I built (and later deployed) to treat participant retention as the time-to-event problem it really is, rather than as a static "who left" count.

## Epidemiological Objective

When I was running recruitment analytics at the Broad Institute, the question that mattered was never "how many participants did we lose" but "when did we lose them, and which subgroups carry most of the hazard." The standard REDCap reports could not answer either of those questions, so I rewrote the analysis as a survival problem and then turned that work into this reusable toolkit.

The three questions I am trying to answer are:

1. Time-to-attrition. What is the probability that a participant is still active at time *t* after enrollment, and which visit is the marginal hazard of drop-out highest?
2. Covariate-adjusted hazard. Holding site, age, sex, race/ethnicity, preferred language, and visit burden fixed, which covariates are independently associated with elevated hazard of attrition?
3. Equity-stratified retention. Are the retention curves statistically separable across demographic strata, and does the proportional-hazards assumption hold inside each stratum?

The applied motivation is concrete: this exact workflow drove a 25% reduction in participant drop-off in production, because the survival curve told us visit 2 was the inflection point and we built a 7-day reminder against it.

## Methodological Framework

I built this repo as an explicit survival-analysis showcase, because I want a reviewer to see the semi-parametric regression theory and not just the REDCap operationalization. The models in the code are:

**Kaplan-Meier product-limit estimator.** Non-parametric estimation of the marginal survival function S(t) = Pr(T > t), where T is time from enrollment to drop-out and the time origin is the enrollment date. Censoring is administrative (participants still active at the end of follow-up) and I treat it as non-informative under the standard independent-censoring assumption. I document that assumption and stress-test it in a sensitivity analysis that re-codes the latest censoring tertile as events.

**Log-rank test (Mantel-Haenszel).** I use this to test the global null of equal survival functions across demographic strata (site, race/ethnicity, language, age band). The log-rank statistic is reported alongside the stratified KM curves so the visual and the inferential evidence sit next to each other.

**Cox proportional hazards regression.** A semi-parametric multivariable model of the form h(t | X) = h0(t) * exp(beta'X), where the baseline hazard h0(t) is left unspecified and only the log-hazard ratios beta are estimated by partial likelihood. Predictors are site (factor), age (continuous, centered), sex, race/ethnicity (OMB categories), preferred language, and a visit-burden index (cumulative scheduled visits in the prior 30 days).

**Diagnostics for the proportional-hazards assumption.** This is the part I really want a methods reviewer to see, because the interpretive validity of a Cox model collapses if proportionality fails. I run three checks:

1. Schoenfeld residual test. Global and per-covariate chi-squared tests of the null that the scaled Schoenfeld residuals are uncorrelated with a transformation of time (`survival::cox.zph`). A non-significant global p-value is required for the headline model to stand.
2. Log-log survival plots. For each categorical covariate I plot log(-log S(t)) against log(t). Parallel curves are taken as visual support for proportionality.
3. Time-by-covariate interaction terms. If the Schoenfeld diagnostics flag a violation, I re-fit the offending covariate with an explicit `tt()` time-varying coefficient and report the time-stratified hazard ratios in the appendix rather than hiding them.

**Bootstrapped equity rate ratios.** Stratum-specific retention rate ratios come with 1,000-replicate non-parametric bootstrap 95% CIs, because asymptotic normal approximations fail fast in small demographic cells and I do not want the equity story to depend on a fragile interval.

Taken together (KM for the marginal picture, log-rank for stratified inference, Cox PH with explicit proportionality diagnostics for the adjusted hazard ratios) this is the standard reporting template that a survival-heavy JAMA paper would use.

## Data Architecture

The pipeline is a five-stage DAG. Each R script owns one transformation. REDCap exports are messy by design (they are wide for instruments, long for repeating events, and they mix coded factors with raw text), so the ingest layer carries most of the data-engineering work.

**Stage 1, REDCap ingest (`R/01_ingest_redcap.R`).** I read both the long-format event log and the wide-format instrument exports natively. The `REDCapR` API path is supported, with a CSV fallback for the synthetic case in this repo. Coded factor levels get restored from the REDCap data dictionary so nothing downstream depends on numeric codes. A typed manifest of row counts, primary keys, and event-name uniqueness is written to `logs/ingest_manifest.json`, and any deviation halts the pipeline.

**Stage 2, per-participant timeline construction (`R/02_build_retention.R`).** Each participant collapses into a survival object with three fields: enrollment date (time origin), event date (first missed scheduled visit past the grace window, or last completed visit plus administrative end), and a 0/1 event indicator. Grace windows are configurable (default 14 days). I handled four edge cases explicitly:

- Missing visit dates. A scheduled visit without a completion record only counts as a candidate event if (a) the scheduled date is more than grace_window days in the past and (b) no later visit in the protocol has been completed. Otherwise I flag the visit as "late but not lost."
- Out-of-order events. Participants sometimes complete visit 3 before visit 2. The timeline collapses by protocol visit number, not by chronological order, and I log the discrepancy.
- Re-engagement. A participant who returns after the grace window is left-censored at the first return and treated as a competing process. The headline analysis uses only the first attrition event, with a sensitivity analysis using a recurrent-event Andersen-Gill formulation.
- Administrative censoring. Active participants at the data freeze are right-censored at the freeze date.

**Stage 3, KPI tables (`R/03_kpis.R`).** Enrollment rate, visit-1 completion, 30-day retention, KM-estimated median time-to-attrition, and the equity retention ratio, all written out as long-format Tableau-ready CSVs.

**Stage 4, survival modeling (`R/04_survival.R`).** Fits KM, runs the log-rank test, fits Cox PH, and emits `cox.zph` diagnostics, log-log plots, and the time-varying coefficient appendix into `outputs/diagnostics/`. The script fails loudly if the global Schoenfeld p < 0.05 unless I set `allow_ph_violation = TRUE` in the config, so a reviewer can trust that any committed headline hazard ratio has passed proportionality screening.

**Stage 5, equity stratification (`R/05_equity.R`).** Retention rate ratios across demographic strata with 1,000-replicate bootstrap CIs. The headline equity statistic is the lowest-vs-highest stratum ratio.

```
redcap-retention-analytics/
├── README.md
├── R/
│   ├── 01_ingest_redcap.R       # read & normalize REDCap export (long + wide)
│   ├── 02_build_retention.R     # per-participant survival object
│   ├── 03_kpis.R                # retention KPI tables
│   ├── 04_survival.R            # KM + log-rank + Cox PH + PH diagnostics
│   └── 05_equity.R              # stratified retention ratios + bootstrap CIs
├── data/
│   ├── synthetic/               # fake REDCap export (committed)
│   └── README.md
├── reports/
│   └── retention_dashboard.qmd  # Quarto dashboard
├── outputs/
│   └── diagnostics/             # cox.zph tables, log-log plots
└── LICENSE
```

## KPI definitions

| KPI | Definition |
|-----|------------|
| Enrollment rate | n_enrolled / n_screened |
| Visit-1 completion | n_visit1_complete / n_enrolled |
| 30-day retention | participants still active 30 days after enrollment |
| Median time-to-dropoff | Kaplan-Meier point estimate of the median of T, censored at end of follow-up |
| Equity retention ratio | retention rate in lowest stratum divided by highest stratum (bootstrap 95% CI) |

## Quick start

```r
install.packages(c("tidyverse","survival","survminer","gtsummary","REDCapR","ggplot2"))
source("R/01_ingest_redcap.R")    # data/synthetic/ -> tidy tibbles
source("R/02_build_retention.R")  # per-participant survival object
source("R/03_kpis.R")             # KPI summary
source("R/04_survival.R")         # KM + log-rank + Cox PH + cox.zph diagnostics
source("R/05_equity.R")           # retention rate ratios + bootstrap CIs

quarto::quarto_render("reports/retention_dashboard.qmd")
```

## Why this matters

In production, the same toolkit did three things:

- It identified visit 2 as the highest-hazard timepoint on the KM curve and then confirmed it as an independent risk factor in the Cox model. I added a 7-day reminder against that visit and retention went up 25%.
- It surfaced site-level variance that aggregate reporting was hiding, which let us target operational support to the underperforming sites.
- It documented equity gaps in retention across language preference, which enabled translated outreach for the groups that were dropping out fastest.

## Disclaimer

All data committed here is synthetic. It is not derived from any real REDCap project.

## Contact

Romario Joseph, rjoseph3@bu.edu, [LinkedIn](https://www.linkedin.com/in/romariojosephpublichealth/)
