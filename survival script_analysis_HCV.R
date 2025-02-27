# load packages
pacman::p_load(dplyr, arsenal, survival, stats, readxl)

## set wd
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Publications/Sex work and risk of HIV and HCV/Emails to authors/Tijuana data/Data")

# load data
analysis_data_hcv_clean <- read_excel("HCV_data_clean.xlsx")

# baseline characteristics sex work
analysis_data_hcv_clean <- analysis_data_hcv_clean %>%
  group_by(id) %>%
  mutate(id_seq = row_number())

analysis_data_hcv_bl <- analysis_data_hcv_clean %>%
  group_by(id) %>%
  mutate(id_seq = row_number(),
         SexWork6Mo = max(SexWork6Mo, na.rm = TRUE),
         SexWMen6Mo = max(SexWMen6Mo, na.rm = TRUE),
         TPrisJail6Mo = max(TPrisJail6Mo, na.rm = TRUE),
         MetBupPrg6M = max(MetBupPrg6M, na.rm = TRUE),
         Homeless = max(Homeless, na.rm = TRUE)) %>%
  ungroup() %>%
  subset(id_seq == 1)

# convert numeric to factor variables
analysis_data_hcv_bl$SEXBRTH <- factor(analysis_data_hcv_bl$SEXBRTH)
analysis_data_hcv_bl$Homeless <- factor(analysis_data_hcv_bl$Homeless)
analysis_data_hcv_bl$TPrisJail6Mo <- factor(analysis_data_hcv_bl$TPrisJail6Mo)
analysis_data_hcv_bl$MetBupPrg6M <- factor(analysis_data_hcv_bl$MetBupPrg6M)
analysis_data_hcv_bl$SexWork <- factor(analysis_data_hcv_bl$SexWork)
analysis_data_hcv_bl$SexWMen <- factor(analysis_data_hcv_bl$SexWMen)
analysis_data_hcv_bl$SexWork6Mo <- factor(analysis_data_hcv_bl$SexWork6Mo)
analysis_data_hcv_bl$SexWMen6Mo <- factor(analysis_data_hcv_bl$SexWMen6Mo)
analysis_data_hcv_bl$DGINJFQB <- factor(analysis_data_hcv_bl$DGINJFQB)
analysis_data_hcv_clean$days_risk <- as.numeric(analysis_data_hcv_clean$days_risk)

 # table
tab_bl_sw_sex_hcv <- tableby(SEXBRTH ~ AGE + SexWork6Mo + SexWMen6Mo + SexWork + SexWMen + AGE1INJ + Homeless + TPrisJail6Mo + MetBupPrg6M + DGINJFQB, data=analysis_data_hcv_bl)
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

# Function to calculate incidence rates and rate ratios
calculate_incidence_and_rate_ratio <- function(data, time_bin, group_label) {
  # selling sex work incidence rate
  total_days_hcv_sw <- sum(data$days_risk[data[[time_bin]] == 1])
  total_cases_sw <- sum(data$hcv_rslt[data[[time_bin]] == 1])
  incidence_rate_sw <- (total_cases_sw / total_days_hcv_sw) * 365.25 * 100

  # Calculate 95% confidence intervals for sex workers
  incidence_rate_sw_se <- sqrt(total_cases_sw) / total_days_hcv_sw * 365.25 * 100
  ci_lower_sw <- incidence_rate_sw - 1.96 * incidence_rate_sw_se
  ci_upper_sw <- incidence_rate_sw + 1.96 * incidence_rate_sw_se

  cat("Incidence rate of HCV per 100 person years among sex workers (", group_label, "):", incidence_rate_sw, "\n")
  cat("95% CI:", ci_lower_sw, "-", ci_upper_sw, "\n")

  # no sex work incidence rate
  total_days_hcv_nosw <- sum(data$days_risk[data[[time_bin]] == 0])
  total_cases_nosw <- sum(data$hcv_rslt[data[[time_bin]] == 0])
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
    group_by(!!sym(time_bin), Homeless, TPrisJail6Mo, DGINJFQB, AGE1INJ) %>%
    summarise(
      total_cases = sum(hcv_rslt),
      total_days = sum(days_risk)
    ) %>%
    mutate(
      rate = total_cases / total_days * 365.25 * 100
    )

  # Fit Poisson regression model controlling for Homeless and TPrisJail6Mo
  poisson_model1 <- glm(total_cases ~ get(time_bin) + Homeless + TPrisJail6Mo + offset(log(total_days)), 
                        family = poisson(link = "log"), 
                        data = summary_data)

  # Extract rate ratio and confidence intervals
  rate_ratio1 <- exp(coef(poisson_model1)[2])
  ci1 <- exp(confint(poisson_model1)[2, ])

  cat("Rate ratio of HCV (sex workers vs non-sex workers) controlling for Homeless and TPrisJail6Mo (", group_label, "):", rate_ratio1, "\n")
  cat("95% CI:", ci1[1], "-", ci1[2], "\n")

  # Fit Poisson regression model controlling for Homeless, TPrisJail6Mo, DGINJFQB, and AGE1INJ
  poisson_model2 <- glm(total_cases ~ get(time_bin) + Homeless + TPrisJail6Mo + DGINJFQB + AGE1INJ + offset(log(total_days)), 
                        family = poisson(link = "log"), 
                        data = summary_data)

  # Extract rate ratio and confidence intervals
  rate_ratio2 <- exp(coef(poisson_model2)[2])
  ci2 <- exp(confint(poisson_model2)[2, ])

  cat("Rate ratio of HCV (sex workers vs non-sex workers) controlling for Homeless, TPrisJail6Mo, DGINJFQB, and AGE1INJ (", group_label, "):", rate_ratio2, "\n")
  cat("95% CI:", ci2[1], "-", ci2[2], "\n")
}

# Calculate for recent exposure
analysis_data_hcv_clean <- analysis_data_hcv_clean %>%
  mutate(sw_time_bin_recent = ifelse(is.na(SexWork6Mo), 0, SexWork6Mo))

calculate_incidence_and_rate_ratio(analysis_data_hcv_clean, "sw_time_bin_recent", "recent exposure")

# Calculate for lifetime exposure
analysis_data_hcv_clean <- analysis_data_hcv_clean %>%
  mutate(sw_time_bin_lifetime = ifelse(is.na(SexWork), 0, SexWork))

calculate_incidence_and_rate_ratio(analysis_data_hcv_clean, "sw_time_bin_lifetime", "lifetime exposure")

# Function to calculate incidence rates and rate ratios
calculate_incidence_and_rate_ratio <- function(data, time_bin, group_label) {
  # selling sex work incidence rate
  total_days_hcv_sw <- sum(data$days_risk[data[[time_bin]] == 1])
  total_cases_sw <- sum(data$hcv_rslt[data[[time_bin]] == 1])
  incidence_rate_sw <- (total_cases_sw / total_days_hcv_sw) * 365.25 * 100

  # Calculate 95% confidence intervals for sex workers
  incidence_rate_sw_se <- sqrt(total_cases_sw) / total_days_hcv_sw * 365.25 * 100
  ci_lower_sw <- incidence_rate_sw - 1.96 * incidence_rate_sw_se
  ci_upper_sw <- incidence_rate_sw + 1.96 * incidence_rate_sw_se

  cat("Incidence rate of HCV per 100 person years among sex workers (", group_label, "):", incidence_rate_sw, "\n")
  cat("95% CI:", ci_lower_sw, "-", ci_upper_sw, "\n")

  # no sex work incidence rate
  total_days_hcv_nosw <- sum(data$days_risk[data[[time_bin]] == 0])
  total_cases_nosw <- sum(data$hcv_rslt[data[[time_bin]] == 0])
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
    group_by(!!sym(time_bin), Homeless, TPrisJail6Mo, DGINJFQB, AGE1INJ) %>%
    summarise(
      total_cases = sum(hcv_rslt),
      total_days = sum(days_risk)
    ) %>%
    mutate(
      rate = total_cases / total_days * 365.25 * 100
    )

  # Fit Poisson regression model controlling for Homeless and TPrisJail6Mo
  poisson_model1 <- glm(total_cases ~ get(time_bin) + Homeless + TPrisJail6Mo + offset(log(total_days)), 
                        family = poisson(link = "log"), 
                        data = summary_data)

  # Extract rate ratio and confidence intervals
  rate_ratio1 <- exp(coef(poisson_model1)[2])
  ci1 <- exp(confint(poisson_model1)[2, ])

  cat("Rate ratio of HCV (sex workers vs non-sex workers) controlling for Homeless and TPrisJail6Mo (", group_label, "):", rate_ratio1, "\n")
  cat("95% CI:", ci1[1], "-", ci1[2], "\n")

  # Fit Poisson regression model controlling for Homeless, TPrisJail6Mo, DGINJFQB, and AGE1INJ
  poisson_model2 <- glm(total_cases ~ get(time_bin) + Homeless + TPrisJail6Mo + DGINJFQB + AGE1INJ + offset(log(total_days)), 
                        family = poisson(link = "log"), 
                        data = summary_data)

  # Extract rate ratio and confidence intervals
  rate_ratio2 <- exp(coef(poisson_model2)[2])
  ci2 <- exp(confint(poisson_model2)[2, ])

  cat("Rate ratio of HCV (sex workers vs non-sex workers) controlling for Homeless, TPrisJail6Mo, DGINJFQB, and AGE1INJ (", group_label, "):", rate_ratio2, "\n")
  cat("95% CI:", ci2[1], "-", ci2[2], "\n")
}

# Calculate for recent exposure
analysis_data_hcv_clean <- analysis_data_hcv_clean %>%
  mutate(sw_time_bin_recent = ifelse(is.na(SexWork6Mo), 0, SexWork6Mo))

# Calculate for lifetime exposure
analysis_data_hcv_clean <- analysis_data_hcv_clean %>%
  mutate(sw_time_bin_lifetime = ifelse(is.na(SexWork), 0, SexWork))

# Stratify by sex and calculate incidence and rate ratios
analysis_data_hcv_men <- subset(analysis_data_hcv_clean, SEXBRTH == 1)
analysis_data_hcv_women <- subset(analysis_data_hcv_clean, SEXBRTH == 2)

# Recent exposure
calculate_incidence_and_rate_ratio(analysis_data_hcv_men, "sw_time_bin_recent", "men - recent exposure")
calculate_incidence_and_rate_ratio(analysis_data_hcv_women, "sw_time_bin_recent", "women - recent exposure")

# Lifetime exposure
calculate_incidence_and_rate_ratio(analysis_data_hcv_men, "sw_time_bin_lifetime", "men - lifetime exposure")
calculate_incidence_and_rate_ratio(analysis_data_hcv_women, "sw_time_bin_lifetime", "women - lifetime exposure")

### MSM HCV incidence rate calculations ###

# Calculate person-time for each group
analysis_data_hcv_clean <- analysis_data_hcv_clean %>%
  mutate(msm_time_bin_recent = ifelse(is.na(SexWMen6Mo), 0, SexWMen6Mo),
         msm_time_bin_lifetime = ifelse(is.na(SexWMen), 0, SexWMen))

# Function to calculate incidence rates and rate ratios
calculate_incidence_and_rate_ratio <- function(data, time_bin, group_label) {
  # MSM incidence rate
  total_days_hcv_msm <- sum(data$days_risk[data[[time_bin]] == 1])
  total_cases_msm <- sum(data$hcv_rslt[data[[time_bin]] == 1])
  incidence_rate_msm <- (total_cases_msm / total_days_hcv_msm) * 365.25 * 100

  # Calculate 95% confidence intervals for MSM
  incidence_rate_msm_se <- sqrt(total_cases_msm) / total_days_hcv_msm * 365.25 * 100
  ci_lower_msm <- incidence_rate_msm - 1.96 * incidence_rate_msm_se
  ci_upper_msm <- incidence_rate_msm + 1.96 * incidence_rate_msm_se

  cat("Incidence rate of HCV per 100 person years among MSM (", group_label, "):", incidence_rate_msm, "\n")
  cat("95% CI:", ci_lower_msm, "-", ci_upper_msm, "\n")

  # non-MSM incidence rate
  total_days_hcv_nonmsm <- sum(data$days_risk[data[[time_bin]] == 0])
  total_cases_nonmsm <- sum(data$hcv_rslt[data[[time_bin]] == 0])
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
    group_by(!!sym(time_bin), Homeless, TPrisJail6Mo, DGINJFQB, AGE1INJ) %>%
    summarise(
      total_cases = sum(hcv_rslt),
      total_days = sum(days_risk)
    ) %>%
    mutate(
      rate = total_cases / total_days * 365.25 * 100
    )

  # Fit Poisson regression model controlling for Homeless and TPrisJail6Mo
  poisson_model1 <- glm(total_cases ~ get(time_bin) + Homeless + TPrisJail6Mo + offset(log(total_days)), 
                        family = poisson(link = "log"), 
                        data = summary_data)

  # Extract rate ratio and confidence intervals
  rate_ratio1 <- exp(coef(poisson_model1)[2])
  ci1 <- exp(confint(poisson_model1)[2, ])

  cat("Rate ratio of HCV (MSM vs non-MSM) controlling for Homeless and TPrisJail6Mo (", group_label, "):", rate_ratio1, "\n")
  cat("95% CI:", ci1[1], "-", ci1[2], "\n")

  # Fit Poisson regression model controlling for Homeless, TPrisJail6Mo, DGINJFQB, and AGE1INJ
  poisson_model2 <- glm(total_cases ~ get(time_bin) + Homeless + TPrisJail6Mo + DGINJFQB + AGE1INJ + offset(log(total_days)), 
                        family = poisson(link = "log"), 
                        data = summary_data)

  # Extract rate ratio and confidence intervals
  rate_ratio2 <- exp(coef(poisson_model2)[2])
  ci2 <- exp(confint(poisson_model2)[2, ])

  cat("Rate ratio of HCV (MSM vs non-MSM) controlling for Homeless, TPrisJail6Mo, DGINJFQB, and AGE1INJ (", group_label, "):", rate_ratio2, "\n")
  cat("95% CI:", ci2[1], "-", ci2[2], "\n")
}

# Calculate for men - recent MSM exposure
analysis_data_hcv_men <- subset(analysis_data_hcv_clean, SEXBRTH == 1)
calculate_incidence_and_rate_ratio(analysis_data_hcv_men, "msm_time_bin_recent", "men - recent MSM exposure")

# Calculate for men - lifetime MSM exposure
calculate_incidence_and_rate_ratio(analysis_data_hcv_men, "msm_time_bin_lifetime", "men - lifetime MSM exposure")