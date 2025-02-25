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

# male and female
tab_bl_sw_all_hcv <- tableby(~ SEXBRTH + AGE + SexWork + AGE1INJ + Homeless + TPrisJail6Mo + MetBupPrg6M + DGINJFQB, data=analysis_data_hcv_bl)
summary(tab_bl_sw_all_hcv, text=TRUE)
 # stratified
tab_bl_sw_sex_hcv <- tableby(SEXBRTH ~ AGE + SexWork + SexWMen + AGE1INJ + Homeless + TPrisJail6Mo + MetBupPrg6M + DGINJFQB, data=analysis_data_hcv_bl)
summary(tab_bl_sw_sex_hcv, text=TRUE)


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




# unadjusted hazard ratio sw HCV ever
sw_model_hcv_crude_ever = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWork, 
  data = analysis_data_hcv_clean
)

summary(sw_model_hcv_crude_ever)

msm_model_hcv_crude_rec = coxph(
  Surv(time = days_risk_start, time2 = days_risk_end, event = hcv_rslt) ~ SexWMen, 
  data = analysis_data_hcv_clean
)

summary(msm_model_hcv_crude_rec)













# baseline characteristics sex work (men and women)
analysis_data_hiv_long <- analysis_data_hiv_long %>%
  group_by(id) %>%
  mutate(id_seq = row_number())

analysis_data_hiv_bl <- analysis_data_hiv_long
analysis_data_hiv_bl <- subset(analysis_data_hiv_bl, id_seq == 1)

analysis_data_hcv_long <- analysis_data_hcv_long %>%
  group_by(id) %>%
  mutate(id_seq = row_number())

analysis_data_hcv_bl <- analysis_data_hcv_long
analysis_data_hcv_bl <- subset(analysis_data_hcv_bl, id_seq == 1)

# convert numeric to factor variables
analysis_data_hiv_bl$SEXBRTH <- factor(analysis_data_hiv_bl$SEXBRTH)
analysis_data_hiv_bl$Homeless <- factor(analysis_data_hiv_bl$Homeless)
analysis_data_hiv_bl$TPrisJail6Mo <- factor(analysis_data_hiv_bl$TPrisJail6Mo)
analysis_data_hiv_bl$MetBupPrg6M <- factor(analysis_data_hiv_bl$MetBupPrg6M)
analysis_data_hiv_bl$SexWork <- factor(analysis_data_hiv_bl$SexWork)
analysis_data_hiv_bl$SexWMen <- factor(analysis_data_hiv_bl$SexWMen)

# men and women
tab_bl_sw_all_hiv  <- tableby(~ SEXBRTH + AGE + SexWork + AGE1INJ + Homeless + TPrisJail6Mo + MetBupPrg6M + DGINJFQB, data=analysis_data_hiv_bl)
summary(tab_bl_sw_all_hiv, text=TRUE)



#stratified by sex
tab_bl_sw_sex_hiv <- tableby(SEXBRTH ~ AGE + SexWork + SexWMen + AGE1INJ + Homeless + TPrisJail6Mo + MetBupPrg6M + DGINJFQB, data=analysis_data_hiv_bl)
summary(tab_bl_sw_sex, text=TRUE)



# unadjusted hazard ratio sw HIV recent
sw_model_hiv_crude_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_sw_6m, 
  data = analysis_data_hiv_clean
)

summary(sw_model_hiv_crude_rec)

# unadjusted hazard ratio sw HIV ever
sw_model_hiv_crude_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWork, 
  data = analysis_data_hiv_clean
)

summary(sw_model_hiv_crude_ever)



# men only HIV recent
analysis_data_hiv_men <- subset(analysis_data_hiv_long, SEXBRTH == 1)

sw_model_hiv_crude_men_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_sw_6m, 
  data = analysis_data_hiv_men
)

summary(sw_model_hiv_crude_men_rec)

# men only HIV ever
sw_model_hiv_crude_men_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWork, 
  data = analysis_data_hiv_men
)

summary(sw_model_hiv_crude_men_ever)



# men only HCV ever
sw_model_hcv_crude_men_ever = coxph(
  Surv(time = hcv_start_risk, time2 = hcv_end_risk, event = hcv_tst_rslt) ~ SexWork, 
  data = analysis_data_hcv_men
)

summary(sw_model_hcv_crude_men_ever)

# women only HIV recent
analysis_data_hiv_women <- subset(analysis_data_hiv_long, SEXBRTH == 2)

sw_model_hiv_crude_women_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_sw_6m, 
  data = analysis_data_hiv_women
)

summary(sw_model_hiv_crude_women_rec)

# women only HIV ever
sw_model_hiv_crude_women_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWork, 
  data = analysis_data_hiv_women
)

summary(sw_model_hiv_crude_women_ever)

# women only HCV ever
analysis_data_hcv_women <- subset(analysis_data_hcv_long, SEXBRTH == 2)

sw_model_hcv_crude_women_ever = coxph(
  Surv(time = hcv_start_risk, time2 = hcv_end_risk, event = hcv_tst_rslt) ~ hcv_sw_6m, 
  data = analysis_data_hcv_women
)

summary(sw_model_hcv_crude_women_ever)

# unadjusted hazard ratio msm HIV recent
msm_model_hiv_crude_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_msm_6m, 
  data = analysis_data_hiv_men
)

summary(msm_model_hiv_crude_rec)

# unadjusted hazard ratio msm HIV ever
msm_model_hiv_crude_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWMen, 
  data = analysis_data_hiv_men
)

summary(msm_model_hiv_crude_ever)

# unadjusted hazard ratio msm HCV recent
msm_model_hcv_crude_rec = coxph(
  Surv(time = hcv_start_risk, time2 = hcv_end_risk, event = hcv_tst_rslt) ~ hcv_msm_6m, 
  data = analysis_data_hcv_men
)

summary(msm_model_hcv_crude_rec)

# unadjusted hazard ratio msm HCV ever
msm_model_hcv_crude_ever = coxph(
  Surv(time = hcv_start_risk, time2 = hcv_end_risk, event = hcv_tst_rslt) ~ SexWMen, 
  data = analysis_data_hcv_men
)

summary(msm_model_hcv_crude_ever)


