# load packages
pacman::p_load(dplyr, arsenal, survival)

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
analysis_data_hcv_clean$days_risk <- as.numeric(analysis_data_hcv_clean$days_risk)

# male and female
tab_bl_sw_all_hcv <- tableby(~ SEXBRTH + AGE + SexWork + AGE1INJ + Homeless + TPrisJail6Mo + MetBupPrg6M + DGINJFQB, data=analysis_data_hcv_bl)
summary(tab_bl_sw_all_hcv, text=TRUE)
 # stratified
tab_bl_sw_sex_hcv <- tableby(SEXBRTH ~ AGE + SexWork + SexWMen + AGE1INJ + Homeless + TPrisJail6Mo + MetBupPrg6M + DGINJFQB, data=analysis_data_hcv_bl)
summary(tab_bl_sw_sex_hcv, text=TRUE)

## incidence rate calculations

# overall incidence rate
total_days_hcv <- sum(analysis_data_hcv_clean$days_risk)
total_cases <- sum(analysis_data_hcv_clean$hcv_rslt)
incidence_rate <- (total_cases / total_days_hcv) * 365.25 *100

cat("Incidence rate of HCV per 100 person years:", incidence_rate)

# selling sex work incidence rate
analysis_data_hcv_clean$sw_time_bin <- analysis_data_hcv_clean$SexWork6Mo
analysis_data_hcv_clean$sw_time_bin[is.na(analysis_data_hcv_clean$sw_time_bin)] <- 0
total_days_hcv_sw <- sum(analysis_data_hcv_clean$days_risk[analysis_data_hcv_clean$sw_time_bin == 1])
total_cases_sw <- sum(analysis_data_hcv_clean$hcv_rslt[analysis_data_hcv_clean$sw_time_bin == 1])
incidence_rate_sw <- (total_cases_sw / total_days_hcv_sw) * 365.25 *100

cat("Incidence rate of HCV per 100 person years among sex workers:", incidence_rate_sw)

# no sex work incidence rate
analysis_data_hcv_clean$sw_time_bin <- analysis_data_hcv_clean$SexWork6Mo
analysis_data_hcv_clean$sw_time_bin <- ifelse(is.na(analysis_data_hcv_clean$SexWork6Mo), 1, analysis_data_hcv_clean$SexWork6Mo)
total_days_hcv_nosw <- sum(analysis_data_hcv_clean$days_risk[analysis_data_hcv_clean$sw_time_bin == 0])
total_cases_nosw <- sum(analysis_data_hcv_clean$hcv_rslt[analysis_data_hcv_clean$sw_time_bin == 0])
incidence_rate_nosw <- (total_cases_nosw / total_days_hcv_nosw) * 365.25 *100

cat("Incidence rate of HCV per 100 person years among non-sex workers:", incidence_rate_nosw)

######## Recent HCV estimates #########

# unadjusted hazard ratio sw HCV recent - men and women
sw_model_hcv_crude_rec = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork6Mo, 
  data = analysis_data_hcv_clean
)

summary(sw_model_hcv_crude_rec)

# adjusted hazard ratio sw HCV recent - men and women
sw_model_hcv_adj1_rec = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork6Mo + Homeless + TPrisJail6Mo, 
  data = analysis_data_hcv_clean
)

summary(sw_model_hcv_adj1_rec)

sw_model_hcv_adj2_rec = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork6Mo + Homeless + TPrisJail6Mo + DGINJFQB + AGE1INJ, 
  data = analysis_data_hcv_clean
)

summary(sw_model_hcv_adj2_rec)

# unadjusted hazard ratio sw HCV recent - men
analysis_data_hcv_men <- subset(analysis_data_hcv_clean, SEXBRTH == 1)

sw_model_hcv_crude_men_rec = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork6Mo, 
  data = analysis_data_hcv_men
)

summary(sw_model_hcv_crude_men_rec)

# adjusted hazard ratio sw HCV recent - men 
sw_model_hcv_adj1_men_rec = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork6Mo + Homeless + TPrisJail6Mo, 
  data = analysis_data_hcv_men
)

summary(sw_model_hcv_adj1_men_rec)

sw_model_hcv_adj2_men_rec = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork6Mo + Homeless + TPrisJail6Mo + DGINJFQB + AGE1INJ, 
  data = analysis_data_hcv_men
)

summary(sw_model_hcv_adj2_men_rec)

# unadjusted hazard ratio sw HCV recent - women
analysis_data_hcv_women <- subset(analysis_data_hcv_clean, SEXBRTH == 2)

sw_model_hcv_crude_women_rec = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork6Mo, 
  data = analysis_data_hcv_women
)

summary(sw_model_hcv_crude_women_rec)

# adjusted hazard ratio sw HCV recent - women 
sw_model_hcv_adj1_women_rec = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork6Mo + Homeless + TPrisJail6Mo, 
  data = analysis_data_hcv_women
)

summary(sw_model_hcv_adj1_women_rec)

sw_model_hcv_adj2_women_rec = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork6Mo + Homeless + TPrisJail6Mo + DGINJFQB + AGE1INJ, 
  data = analysis_data_hcv_women
)

summary(sw_model_hcv_adj2_women_rec)

# unadjusted hazard ratio MSM HCV recent
msm_model_hcv_crude_men_rec = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWMen6Mo, 
  data = analysis_data_hcv_men
)

summary(msm_model_hcv_crude_men_rec)

# adjusted hazard ratio MSM HCV recent
msm_model_hcv_adj1_men_rec = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWMen6Mo + Homeless + TPrisJail6Mo, 
  data = analysis_data_hcv_men
)

summary(msm_model_hcv_adj1_men_rec)

msm_model_hcv_adj2_men_rec = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWMen6Mo + Homeless + TPrisJail6Mo + DGINJFQB + AGE1INJ, 
  data = analysis_data_hcv_men
)

summary(msm_model_hcv_adj2_men_rec)

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

