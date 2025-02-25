# load packages
pacman::p_load(dplyr, arsenal, survival, stats)

## HCV analysis ##

# baseline characteristics sex work
analysis_data_hcv_clean <- analysis_data_hcv_clean %>%
  group_by(id) %>%
  mutate(id_seq = row_number())

analysis_data_hcv_bl <- subset(analysis_data_hcv_clean, id_seq == 1)

# convert numeric to factor variables
analysis_data_hcv_bl$SEXBRTH <- factor(analysis_data_hcv_bl$SEXBRTH)
analysis_data_hcv_bl$Homeless <- factor(analysis_data_hcv_bl$Homeless)
analysis_data_hcv_bl$TPrisJail6Mo <- factor(analysis_data_hcv_bl$TPrisJail6Mo)
analysis_data_hcv_bl$MetBupPrg6M <- factor(analysis_data_hcv_bl$MetBupPrg6M)
analysis_data_hcv_bl$SexWork <- factor(analysis_data_hcv_bl$SexWork)
analysis_data_hcv_bl$SexWMen <- factor(analysis_data_hcv_bl$SexWMen)
analysis_data_hcv_bl$DGINJFQB <- factor(analysis_data_hcv_bl$DGINJFQB)
analysis_data_hcv_clean$days_risk <- as.numeric(analysis_data_hcv_clean$days_risk)


# male and female
tab_bl_sw_all_hcv <- tableby(~ SEXBRTH + AGE + SexWork + AGE1INJ + Homeless + TPrisJail6Mo + MetBupPrg6M + DGINJFQB, data=analysis_data_hcv_bl)
summary(tab_bl_sw_all_hcv, text=TRUE)
 # stratified
tab_bl_sw_sex_hcv <- tableby(SEXBRTH ~ AGE + SexWork + SexWMen + AGE1INJ + Homeless + TPrisJail6Mo + MetBupPrg6M + DGINJFQB, data=analysis_data_hcv_bl)
summary(tab_bl_sw_sex_hcv, text=TRUE)

### overall incidence rate calculations ###

total_days_hcv <- sum(analysis_data_hcv_clean$days_risk)
total_cases <- sum(analysis_data_hcv_clean$hcv_rslt)
incidence_rate <- (total_cases / total_days_hcv) * 365.25 * 100

# Calculate 95% confidence intervals
incidence_rate_se <- sqrt(total_cases) / total_days_hcv * 365.25 * 100
ci_lower <- incidence_rate - 1.96 * incidence_rate_se
ci_upper <- incidence_rate + 1.96 * incidence_rate_se

cat("Incidence rate of HCV per 100 person years:", incidence_rate, "\n")
cat("95% CI:", ci_lower, "-", ci_upper, "\n")

### sex work HCV incidence rate calculations ###

# selling sex work incidence rate
analysis_data_hcv_clean$sw_time_bin <- analysis_data_hcv_clean$SexWork6Mo
analysis_data_hcv_clean$sw_time_bin[is.na(analysis_data_hcv_clean$sw_time_bin)] <- 0
total_days_hcv_sw <- sum(analysis_data_hcv_clean$days_risk[analysis_data_hcv_clean$sw_time_bin == 1])
total_cases_sw <- sum(analysis_data_hcv_clean$hcv_rslt[analysis_data_hcv_clean$sw_time_bin == 1])
incidence_rate_sw <- (total_cases_sw / total_days_hcv_sw) * 365.25 * 100

# Calculate 95% confidence intervals for sex workers
incidence_rate_sw_se <- sqrt(total_cases_sw) / total_days_hcv_sw * 365.25 * 100
ci_lower_sw <- incidence_rate_sw - 1.96 * incidence_rate_sw_se
ci_upper_sw <- incidence_rate_sw + 1.96 * incidence_rate_sw_se

cat("Incidence rate of HCV per 100 person years among sex workers:", incidence_rate_sw, "\n")
cat("95% CI:", ci_lower_sw, "-", ci_upper_sw, "\n")

# no sex work incidence rate
analysis_data_hcv_clean$sw_time_bin <- analysis_data_hcv_clean$SexWork6Mo
analysis_data_hcv_clean$sw_time_bin <- ifelse(is.na(analysis_data_hcv_clean$SexWork6Mo), 1, analysis_data_hcv_clean$SexWork6Mo)
total_days_hcv_nosw <- sum(analysis_data_hcv_clean$days_risk[analysis_data_hcv_clean$sw_time_bin == 0])
total_cases_nosw <- sum(analysis_data_hcv_clean$hcv_rslt[analysis_data_hcv_clean$sw_time_bin == 0])
incidence_rate_nosw <- (total_cases_nosw / total_days_hcv_nosw) * 365.25 * 100

# Calculate 95% confidence intervals for non-sex workers
incidence_rate_nosw_se <- sqrt(total_cases_nosw) / total_days_hcv_nosw * 365.25 * 100
ci_lower_nosw <- incidence_rate_nosw - 1.96 * incidence_rate_nosw_se
ci_upper_nosw <- incidence_rate_nosw + 1.96 * incidence_rate_nosw_se

cat("Incidence rate of HCV per 100 person years among non-sex workers:", incidence_rate_nosw, "\n")
cat("95% CI:", ci_lower_nosw, "-", ci_upper_nosw, "\n")

# Calculate rate ratio and its 95% confidence interval
rate_ratio <- incidence_rate_sw / incidence_rate_nosw
rate_ratio_se <- sqrt((1 / total_cases_sw) + (1 / total_cases_nosw))
ci_lower_rr <- exp(log(rate_ratio) - 1.96 * rate_ratio_se)
ci_upper_rr <- exp(log(rate_ratio) + 1.96 * rate_ratio_se)

cat("Rate ratio of HCV (sex workers vs non-sex workers):", rate_ratio, "\n")
cat("95% CI:", ci_lower_rr, "-", ci_upper_rr, "\n")

# Calculate person-time for each group
analysis_data_hcv_clean <- analysis_data_hcv_clean %>%
  mutate(sw_time_bin = ifelse(is.na(SexWork6Mo), 0, SexWork6Mo))

# Create a summary dataset for Poisson regression
summary_data <- analysis_data_hcv_clean %>%
  group_by(sw_time_bin, Homeless, TPrisJail6Mo, DGINJFQB, AGE1INJ) %>%
  summarise(
    total_cases = sum(hcv_rslt),
    total_days = sum(days_risk)
  ) %>%
  mutate(
    rate = total_cases / total_days * 365.25 * 100
  )

# Fit Poisson regression model
poisson_model <- glm(total_cases ~ sw_time_bin + offset(log(total_days)), 
                     family = poisson(link = "log"), 
                     data = summary_data)

# Extract rate ratio and confidence intervals
rate_ratio <- exp(coef(poisson_model)[2])
ci <- exp(confint(poisson_model)[2, ])

cat("Rate ratio of HCV (sex workers vs non-sex workers):", rate_ratio, "\n")
cat("95% CI:", ci[1], "-", ci[2], "\n")

# Fit Poisson regression model controlling for Homeless and TPrisJail6Mo
poisson_model <- glm(total_cases ~ sw_time_bin + Homeless + TPrisJail6Mo + offset(log(total_days)), 
                     family = poisson(link = "log"), 
                     data = summary_data)

# Extract rate ratio and confidence intervals
rate_ratio <- exp(coef(poisson_model)[2])
ci <- exp(confint(poisson_model)[2, ])

cat("Rate ratio of HCV (sex workers vs non-sex workers) controlling for Homeless and TPrisJail6Mo:", rate_ratio, "\n")
cat("95% CI:", ci[1], "-", ci[2], "\n")

# Fit Poisson regression model controlling for Homeless, TPrisJail6Mo, DGINJFQB, and AGE1INJ
poisson_model <- glm(total_cases ~ sw_time_bin + Homeless + TPrisJail6Mo + DGINJFQB + AGE1INJ + offset(log(total_days)), 
                     family = poisson(link = "log"), 
                     data = summary_data)

# Extract rate ratio and confidence intervals
rate_ratio <- exp(coef(poisson_model)[2])
ci <- exp(confint(poisson_model)[2, ])

cat("Rate ratio of HCV (sex workers vs non-sex workers) controlling for Homeless, TPrisJail6Mo, DGINJFQB, and AGE1INJ:", rate_ratio, "\n")
cat("95% CI:", ci[1], "-", ci[2], "\n")

# Calculate person-time for each group
analysis_data_hcv_clean <- analysis_data_hcv_clean %>%
  mutate(sw_time_bin = ifelse(is.na(SexWork6Mo), 0, SexWork6Mo))

# Function to calculate incidence rates and rate ratios
calculate_incidence_and_rate_ratio <- function(data, group_label) {
  # selling sex work incidence rate
  total_days_hcv_sw <- sum(data$days_risk[data$sw_time_bin == 1])
  total_cases_sw <- sum(data$hcv_rslt[data$sw_time_bin == 1])
  incidence_rate_sw <- (total_cases_sw / total_days_hcv_sw) * 365.25 * 100

  # Calculate 95% confidence intervals for sex workers
  incidence_rate_sw_se <- sqrt(total_cases_sw) / total_days_hcv_sw * 365.25 * 100
  ci_lower_sw <- incidence_rate_sw - 1.96 * incidence_rate_sw_se
  ci_upper_sw <- incidence_rate_sw + 1.96 * incidence_rate_sw_se

  cat("Incidence rate of HCV per 100 person years among sex workers (", group_label, "):", incidence_rate_sw, "\n")
  cat("95% CI:", ci_lower_sw, "-", ci_upper_sw, "\n")

  # no sex work incidence rate
  total_days_hcv_nosw <- sum(data$days_risk[data$sw_time_bin == 0])
  total_cases_nosw <- sum(data$hcv_rslt[data$sw_time_bin == 0])
  incidence_rate_nosw <- (total_cases_nosw / total_days_hcv_nosw) * 365.25 * 100

  # Calculate 95% confidence intervals for non-sex workers
  incidence_rate_nosw_se <- sqrt(total_cases_nosw) / total_days_hcv_nosw * 365.25 * 100
  ci_lower_nosw <- incidence_rate_nosw - 1.96 * incidence_rate_nosw_se
  ci_upper_nosw <- incidence_rate_nosw + 1.96 * incidence_rate_nosw_se

  cat("Incidence rate of HCV per 100 person years among non-sex workers (", group_label, "):", incidence_rate_nosw, "\n")
  cat("95% CI:", ci_lower_nosw, "-", ci_upper_nosw, "\n")

  # Calculate rate ratio and its 95% confidence interval
  rate_ratio <- incidence_rate_sw / incidence_rate_nosw
  rate_ratio_se <- sqrt((1 / total_cases_sw) + (1 / total_cases_nosw))
  ci_lower_rr <- exp(log(rate_ratio) - 1.96 * rate_ratio_se)
  ci_upper_rr <- exp(log(rate_ratio) + 1.96 * rate_ratio_se)

  cat("Rate ratio of HCV (sex workers vs non-sex workers) (", group_label, "):", rate_ratio, "\n")
  cat("95% CI:", ci_lower_rr, "-", ci_upper_rr, "\n")

  # Create a summary dataset for Poisson regression
  summary_data <- data %>%
    group_by(sw_time_bin, Homeless, TPrisJail6Mo, DGINJFQB, AGE1INJ) %>%
    summarise(
      total_cases = sum(hcv_rslt),
      total_days = sum(days_risk)
    ) %>%
    mutate(
      rate = total_cases / total_days * 365.25 * 100
    )

  # Fit Poisson regression model controlling for Homeless and TPrisJail6Mo
  poisson_model1 <- glm(total_cases ~ sw_time_bin + Homeless + TPrisJail6Mo + offset(log(total_days)), 
                        family = poisson(link = "log"), 
                        data = summary_data)

  # Extract rate ratio and confidence intervals
  rate_ratio1 <- exp(coef(poisson_model1)[2])
  ci1 <- exp(confint(poisson_model1)[2, ])

  cat("Rate ratio of HCV (sex workers vs non-sex workers) controlling for Homeless and TPrisJail6Mo (", group_label, "):", rate_ratio1, "\n")
  cat("95% CI:", ci1[1], "-", ci1[2], "\n")

  # Fit Poisson regression model controlling for Homeless, TPrisJail6Mo, DGINJFQB, and AGE1INJ
  poisson_model2 <- glm(total_cases ~ sw_time_bin + Homeless + TPrisJail6Mo + DGINJFQB + AGE1INJ + offset(log(total_days)), 
                        family = poisson(link = "log"), 
                        data = summary_data)

  # Extract rate ratio and confidence intervals
  rate_ratio2 <- exp(coef(poisson_model2)[2])
  ci2 <- exp(confint(poisson_model2)[2, ])

  cat("Rate ratio of HCV (sex workers vs non-sex workers) controlling for Homeless, TPrisJail6Mo, DGINJFQB, and AGE1INJ (", group_label, "):", rate_ratio2, "\n")
  cat("95% CI:", ci2[1], "-", ci2[2], "\n")
}

# Calculate for men
analysis_data_hcv_men <- subset(analysis_data_hcv_clean, SEXBRTH == 1)
calculate_incidence_and_rate_ratio(analysis_data_hcv_men, "men")

# Calculate for women
analysis_data_hcv_women <- subset(analysis_data_hcv_clean, SEXBRTH == 2)
calculate_incidence_and_rate_ratio(analysis_data_hcv_women, "women")













### MSM HCV incidence rate calculations ###

# Calculate person-time for each group
analysis_data_hcv_clean <- analysis_data_hcv_clean %>%
  mutate(msm_time_bin = ifelse(is.na(SexWMen), 0, SexWMen))

# Function to calculate incidence rates and rate ratios
calculate_incidence_and_rate_ratio <- function(data, group_label) {
  # MSM incidence rate
  total_days_hcv_msm <- sum(data$days_risk[data$msm_time_bin == 1])
  total_cases_msm <- sum(data$hcv_rslt[data$msm_time_bin == 1])
  incidence_rate_msm <- (total_cases_msm / total_days_hcv_msm) * 365.25 * 100

  # Calculate 95% confidence intervals for MSM
  incidence_rate_msm_se <- sqrt(total_cases_msm) / total_days_hcv_msm * 365.25 * 100
  ci_lower_msm <- incidence_rate_msm - 1.96 * incidence_rate_msm_se
  ci_upper_msm <- incidence_rate_msm + 1.96 * incidence_rate_msm_se

  cat("Incidence rate of HCV per 100 person years among MSM (", group_label, "):", incidence_rate_msm, "\n")
  cat("95% CI:", ci_lower_msm, "-", ci_upper_msm, "\n")

  # non-MSM incidence rate
  total_days_hcv_nonmsm <- sum(data$days_risk[data$msm_time_bin == 0])
  total_cases_nonmsm <- sum(data$hcv_rslt[data$msm_time_bin == 0])
  incidence_rate_nonmsm <- (total_cases_nonmsm / total_days_hcv_nonmsm) * 365.25 * 100

  # Calculate 95% confidence intervals for non-MSM
  incidence_rate_nonmsm_se <- sqrt(total_cases_nonmsm) / total_days_hcv_nonmsm * 365.25 * 100
  ci_lower_nonmsm <- incidence_rate_nonmsm - 1.96 * incidence_rate_nonmsm_se
  ci_upper_nonmsm <- incidence_rate_nonmsm + 1.96 * incidence_rate_nonmsm_se

  cat("Incidence rate of HCV per 100 person years among non-MSM (", group_label, "):", incidence_rate_nonmsm, "\n")
  cat("95% CI:", ci_lower_nonmsm, "-", ci_upper_nonmsm, "\n")

  # Calculate rate ratio and its 95% confidence interval
  rate_ratio <- incidence_rate_msm / incidence_rate_nonmsm
  rate_ratio_se <- sqrt((1 / total_cases_msm) + (1 / total_cases_nonmsm))
  ci_lower_rr <- exp(log(rate_ratio) - 1.96 * rate_ratio_se)
  ci_upper_rr <- exp(log(rate_ratio) + 1.96 * rate_ratio_se)

  cat("Rate ratio of HCV (MSM vs non-MSM) (", group_label, "):", rate_ratio, "\n")
  cat("95% CI:", ci_lower_rr, "-", ci_upper_rr, "\n")

  # Create a summary dataset for Poisson regression
  summary_data <- data %>%
    group_by(msm_time_bin, Homeless, TPrisJail6Mo, DGINJFQB, AGE1INJ) %>%
    summarise(
      total_cases = sum(hcv_rslt),
      total_days = sum(days_risk)
    ) %>%
    mutate(
      rate = total_cases / total_days * 365.25 * 100
    )

  # Fit Poisson regression model controlling for Homeless and TPrisJail6Mo
  poisson_model1 <- glm(total_cases ~ msm_time_bin + Homeless + TPrisJail6Mo + offset(log(total_days)), 
                        family = poisson(link = "log"), 
                        data = summary_data)

  # Extract rate ratio and confidence intervals
  rate_ratio1 <- exp(coef(poisson_model1)[2])
  ci1 <- exp(confint(poisson_model1)[2, ])

  cat("Rate ratio of HCV (MSM vs non-MSM) controlling for Homeless and TPrisJail6Mo (", group_label, "):", rate_ratio1, "\n")
  cat("95% CI:", ci1[1], "-", ci1[2], "\n")

  # Fit Poisson regression model controlling for Homeless, TPrisJail6Mo, DGINJFQB, and AGE1INJ
  poisson_model2 <- glm(total_cases ~ msm_time_bin + Homeless + TPrisJail6Mo + DGINJFQB + AGE1INJ + offset(log(total_days)), 
                        family = poisson(link = "log"), 
                        data = summary_data)

  # Extract rate ratio and confidence intervals
  rate_ratio2 <- exp(coef(poisson_model2)[2])
  ci2 <- exp(confint(poisson_model2)[2, ])

  cat("Rate ratio of HCV (MSM vs non-MSM) controlling for Homeless, TPrisJail6Mo, DGINJFQB, and AGE1INJ (", group_label, "):", rate_ratio2, "\n")
  cat("95% CI:", ci2[1], "-", ci2[2], "\n")
}

# Calculate for men
analysis_data_hcv_men <- subset(analysis_data_hcv_clean, SEXBRTH == 1)
calculate_incidence_and_rate_ratio(analysis_data_hcv_men, "men")













######## Ever HCV estimates #########

# unadjusted hazard ratio sw HCV recent - men and women
sw_model_hcv_crude_ever = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork, 
  data = analysis_data_hcv_clean
)

summary(sw_model_hcv_crude_ever)

# adjusted hazard ratio sw HCV recent - men and women
sw_model_hcv_adj1_ever = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork + Homeless + TPrisJail6Mo, 
  data = analysis_data_hcv_clean
)

summary(sw_model_hcv_adj1_ever)

sw_model_hcv_adj2_ever = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork + Homeless + TPrisJail6Mo + DGINJFQB + AGE1INJ, 
  data = analysis_data_hcv_clean
)

summary(sw_model_hcv_adj2_ever)

# unadjusted hazard ratio sw HCV recent - men
analysis_data_hcv_men <- subset(analysis_data_hcv_clean, SEXBRTH == 1)

sw_model_hcv_crude_men_ever = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork, 
  data = analysis_data_hcv_men
)

summary(sw_model_hcv_crude_men_ever)

# adjusted hazard ratio sw HCV recent - men 
sw_model_hcv_adj1_men_ever = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork + Homeless + TPrisJail6Mo, 
  data = analysis_data_hcv_men
)

summary(sw_model_hcv_adj1_men_ever)

sw_model_hcv_adj2_men_ever = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork + Homeless + TPrisJail6Mo + DGINJFQB + AGE1INJ, 
  data = analysis_data_hcv_men
)

summary(sw_model_hcv_adj2_men_ever)

# unadjusted hazard ratio sw HCV recent - women
analysis_data_hcv_women <- subset(analysis_data_hcv_clean, SEXBRTH == 2)

sw_model_hcv_crude_women_ever = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork, 
  data = analysis_data_hcv_women
)

summary(sw_model_hcv_crude_women_ever)

# adjusted hazard ratio sw HCV recent - women 
sw_model_hcv_adj1_women_ever = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork + Homeless + TPrisJail6Mo, 
  data = analysis_data_hcv_women
)

summary(sw_model_hcv_adj1_women_ever)

sw_model_hcv_adj2_women_ever = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork + Homeless + TPrisJail6Mo + DGINJFQB + AGE1INJ, 
  data = analysis_data_hcv_women
)

summary(sw_model_hcv_adj2_women_ever)

# unadjusted hazard ratio MSM HCV recent
msm_model_hcv_crude_men_ever = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWMen, 
  data = analysis_data_hcv_men
)

summary(msm_model_hcv_crude_men_ever)

# adjusted hazard ratio MSM HCV recent
msm_model_hcv_adj1_men_ever = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWMen + Homeless + TPrisJail6Mo, 
  data = analysis_data_hcv_men
)

summary(msm_model_hcv_adj1_men_ever)

msm_model_hcv_adj2_men_ever = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWMen + Homeless + TPrisJail6Mo + DGINJFQB + AGE1INJ, 
  data = analysis_data_hcv_men
)

summary(msm_model_hcv_adj2_men_ever)

