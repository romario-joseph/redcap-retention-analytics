# 03_kpis.R
# Build retention KPI tables overall, by site, and by stratum.
# Author: Romario Joseph | BU SPH

library(tidyverse)
library(lubridate)

ret <- readRDS("data/processed/retention_timeline.rds")

# Overall KPIs
overall <- ret |>
  summarise(
        n_enrolled = n(),
        n_visit1_complete = sum(visit1_complete, na.rm = TRUE),
        n_active_at_30d   = sum(days_followed >= 30 & !dropped, na.rm = TRUE),
        visit1_rate       = n_visit1_complete / n_enrolled,
        retention_30d     = n_active_at_30d / n_enrolled,
        median_days_to_drop = median(days_followed[dropped == 1], na.rm = TRUE)
      )

write_csv(overall, "outputs/tables/kpi_overall.csv")

# By site
by_site <- ret |>
  group_by(site) |>
  summarise(
        n_enrolled = n(),
        visit1_rate    = mean(visit1_complete, na.rm = TRUE),
        retention_30d  = mean(days_followed >= 30 & !dropped, na.rm = TRUE),
        median_drop    = median(days_followed[dropped == 1], na.rm = TRUE),
        .groups = "drop"
      )

write_csv(by_site, "outputs/tables/kpi_by_site.csv")

# By preferred language (equity-relevant)
by_lang <- ret |>
  group_by(preferred_language) |>
  summarise(
        n = n(),
        retention_30d = mean(days_followed >= 30 & !dropped, na.rm = TRUE),
        .groups = "drop"
      ) |>
  arrange(retention_30d)

write_csv(by_lang, "outputs/tables/kpi_by_language.csv")

print(overall)
message("Wrote KPI tables to outputs/tables/")
