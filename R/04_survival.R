# 04_survival.R
# Kaplan-Meier + Cox PH model on time-to-dropoff
# Author: Romario Joseph | BU SPH

library(tidyverse)
library(survival)
library(survminer)

ret <- readRDS("data/processed/retention_timeline.rds")

# Build the survival object: time = days from enrollment to dropoff/end-of-followup
# event = 1 if participant dropped, 0 if still active (censored)
surv_obj <- Surv(time = ret$days_followed, event = ret$dropped)

# --- Kaplan-Meier overall + by site ---
km_overall <- survfit(surv_obj ~ 1, data = ret)
km_site    <- survfit(surv_obj ~ site, data = ret)

p_km <- ggsurvplot(
    km_site, data = ret,
    pval = TRUE, conf.int = TRUE,
    risk.table = TRUE,
    xlab = "Days since enrollment",
    ylab = "Retention probability",
    title = "Time-to-dropoff by site"
  )
ggsave("outputs/figures/km_by_site.png", p_km$plot, width = 8, height = 5, dpi = 200)

# --- Cox proportional hazards ---
cox_fit <- coxph(
    surv_obj ~ age_at_enroll + sex + race + preferred_language + visit_burden + site,
    data = ret
  )
summary(cox_fit)

# Test PH assumption
ph_test <- cox.zph(cox_fit)
print(ph_test)

# Tidy hazard ratios
hr_table <- broom::tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

write_csv(hr_table, "outputs/tables/cox_hr.csv")
message("Saved KM plot and Cox HR table to outputs/")
