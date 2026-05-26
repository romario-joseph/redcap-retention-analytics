# 01_ingest_redcap.R
# Read a REDCap export (long event format) and normalize into tidy frames.
# Handles both API pulls (via REDCapR) and exported CSVs.
# Author: Romario Joseph | BU SPH

library(tidyverse)
# library(REDCapR)  # uncomment for live API pulls

raw_path <- "data/synthetic/redcap_export.csv"

raw <- read_csv(raw_path, show_col_types = FALSE)

# Expected REDCap columns: record_id, redcap_event_name, redcap_repeat_instance,
# enrollment_date, dropout_date, site, age_at_enroll, sex, race, preferred_language,
# visit_status (1=complete, 2=incomplete, 0=missed)

participants <- raw |>
  filter(redcap_event_name == "baseline_arm_1") |>
  select(record_id, enrollment_date, dropout_date, site,
                  age_at_enroll, sex, race, preferred_language) |>
  distinct()

visits <- raw |>
  filter(str_detect(redcap_event_name, "visit")) |>
  transmute(
        record_id,
        visit = redcap_event_name,
        visit_status = visit_status,
        visit_date = as.Date(visit_date)
      ) |>
  arrange(record_id, visit_date)

saveRDS(participants, "data/processed/participants.rds")
saveRDS(visits,       "data/processed/visits.rds")

message("Ingested ", nrow(participants), " participants and ", nrow(visits), " visit records")
