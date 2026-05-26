# REDCap Retention Analytics
## Cutting Participant Drop-off in Multi-Site Clinical Trials

**Author:** Romario Joseph, MPH · BU SPH (Epi & Biostats)
**Stack:** R 4.x · tidyverse · survival · gtsummary · ggplot2 · Quarto
**Data:** Synthetic REDCap export modeled on a multi-site recruitment registry (mirrors my work at the Broad Institute)

---

## Problem
At the Broad Institute I worked on large-scale recruitment initiatives where **participant drop-off** silently undermines statistical power, balanced enrollment, and equity of representation. Conventional REDCap reports tell you who left — they don't tell you **when** drop-off happens, **why** it concentrates in certain strata, or **which sites** need an intervention.

This repo is a reusable R toolkit for those questions. The same workflow I used to **decrease participant drop-off by 25%** in production.

---

## What's inside
1. **REDCap ingest** — read native REDCap exports (long-format event records) and wide-format instrument exports.
2. 2. **Retention KPIs** — enrollment, completion-by-visit, and median time-to-dropoff per site and stratum.
   3. 3. **Survival modeling** — Kaplan-Meier curves and Cox proportional hazards on time-to-dropoff with site, age, sex, race/ethnicity, language, and visit burden as covariates.
      4. 4. **Equity stratification** — retention rate ratios across demographic groups with bootstrapped 95% CIs.
         5. 5. **Dashboard** — a Quarto report with KPI cards and Tableau-ready CSV outputs for ops teams.
           
            6. ---
           
            7. ## Repository structure
            8. ```
               redcap-retention-analytics/
               ├── README.md
               ├── R/
               │   ├── 01_ingest_redcap.R     ← read & normalize REDCap export
               │   ├── 02_build_retention.R   ← build per-participant timeline
               │   ├── 03_kpis.R              ← KPI tables
               │   ├── 04_survival.R          ← KM + Cox PH on time-to-dropoff
               │   └── 05_equity.R            ← retention rate ratios + bootstrap CIs
               ├── data/
               │   ├── synthetic/             ← fake REDCap export (committed)
               │   └── README.md
               ├── reports/
               │   └── retention_dashboard.qmd
               └── LICENSE
               ```

               ---

               ## KPI definitions
               | KPI | Definition |
               |-----|------------|
               | **Enrollment rate** | n_enrolled / n_screened |
               | **Visit-1 completion** | n_visit1_complete / n_enrolled |
               | **30-day retention** | participants still active 30 days after enrollment |
               | **Median time-to-dropoff** | KM estimate, censored at end of follow-up |
               | **Equity retention ratio** | retention rate in lowest stratum ÷ highest stratum |

               ---

               ## Quick start
               ```r
               # install.packages(c("tidyverse","survival","gtsummary","REDCapR","ggplot2"))
               source("R/01_ingest_redcap.R")   # data/synthetic/ → tidy tibbles
               source("R/02_build_retention.R") # per-participant timeline
               source("R/03_kpis.R")            # KPI summary
               source("R/04_survival.R")        # KM + Cox PH
               source("R/05_equity.R")          # retention rate ratios

               quarto::quarto_render("reports/retention_dashboard.qmd")
               ```

               ---

               ## Why this matters
               In production, the same toolkit:
               - Identified visit-2 as the highest-loss timepoint → added a 7-day reminder → retention up 25%.
               - - Surfaced site-level variance hidden by aggregate reporting → targeted operational support to underperforming sites.
                 - - Documented equity gaps in retention across language preference, enabling translated outreach.
                  
                   - ---

                   ## Disclaimer
                   All data committed here is synthetic. Not derived from any real REDCap project.

                   ## Contact
                   **Romario Joseph** · rjoseph3@bu.edu · [LinkedIn](https://www.linkedin.com/in/romariojosephpublichealth/)
                   
