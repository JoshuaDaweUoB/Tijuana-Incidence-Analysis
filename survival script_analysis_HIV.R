# load packages
pacman::p_load(dplyr, arsenal, survival)

######## Baseline characteristics ########

# convert numeric to factor variables
analysis_data_hiv$SEXBRTH <- factor(analysis_data_hiv$SEXBRTH)
analysis_data_hiv$Homeless <- factor(analysis_data_hiv$Homeless)
analysis_data_hiv$TPrisJail6Mo <- factor(analysis_data_hiv$TPrisJail6Mo)
analysis_data_hiv$MetBupPrg6M <- factor(analysis_data_hiv$MetBupPrg6M)
analysis_data_hiv$SexWork <- factor(analysis_data_hiv$SexWork)
analysis_data_hiv$SexWMen <- factor(analysis_data_hiv$SexWMen)
analysis_data_hiv$SexWMen6Mo_1 <- factor(analysis_data_hiv$SexWMen6Mo_1)

# male and female
tab_bl_sw_all_hiv <- tableby(~ SEXBRTH + SexWork + SexWork6Mo_1 + AGE1INJ + Homeless + TPrisJail6Mo + MetBupPrg6M + DGINJFQB, data=analysis_data_hiv)
summary(tab_bl_sw_all_hiv, text=TRUE)

# stratified
tab_bl_sw_sex_hiv <- tableby(SEXBRTH ~ SexWork + SexWork6Mo_1 + SexWMen + SexWMen6Mo_1 + AGE1INJ + Homeless + TPrisJail6Mo + MetBupPrg6M + DGINJFQB, data=analysis_data_hiv)
summary(tab_bl_sw_sex_hiv, text=TRUE)

######## Recent HIV estimates #########

# unadjusted hazard ratio sw HIV recent - men and women
sw_model_hiv_crude_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_sw_6m, 
  data = analysis_data_hiv_long
)

summary(sw_model_hiv_crude_rec)

# adjusted hazard ratio sw HIV recent - men and women
sw_model_hiv_adj1_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_sw_6m + hiv_home_6m + hiv_pris_6m, 
  data = analysis_data_hiv_long
)

summary(sw_model_hiv_adj1_rec)

sw_model_hiv_adj2_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_sw_6m + hiv_home_6m + hiv_pris_6m + DGINJFQB + AGE1INJ, 
  data = analysis_data_hiv_long
)

summary(sw_model_hiv_adj2_rec)

# unadjusted hazard ratio sw HIV recent - men
analysis_data_hiv_men <- subset(analysis_data_hiv_long, SEXBRTH == 1)

sw_model_hiv_crude_men_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_sw_6m, 
  data = analysis_data_hiv_men
)

summary(sw_model_hiv_crude_men_rec)

# adjusted hazard ratio sw HIV recent - men 
sw_model_hiv_adj1_men_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_sw_6m + hiv_home_6m + hiv_pris_6m, 
  data = analysis_data_hiv_men
)

summary(sw_model_hiv_adj1_men_rec)

sw_model_hiv_adj2_men_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_sw_6m + hiv_home_6m + hiv_pris_6m + DGINJFQB + AGE1INJ, 
  data = analysis_data_hiv_men
)

summary(sw_model_hiv_adj2_men_rec)

# unadjusted hazard ratio sw HIV recent - women
analysis_data_hiv_women <- subset(analysis_data_hiv_long, SEXBRTH == 2)

sw_model_hiv_crude_women_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_sw_6m, 
  data = analysis_data_hiv_women
)

summary(sw_model_hiv_crude_women_rec)

# adjusted hazard ratio sw HIV recent - women 
sw_model_hiv_adj1_women_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_sw_6m + hiv_home_6m + hiv_pris_6m, 
  data = analysis_data_hiv_women
)

summary(sw_model_hiv_adj1_women_rec)

sw_model_hiv_adj2_women_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_sw_6m + hiv_home_6m + hiv_pris_6m + DGINJFQB + AGE1INJ, 
  data = analysis_data_hiv_women
)

summary(sw_model_hiv_adj2_women_rec)

# unadjusted hazard ratio MSM HIV recent
msm_model_hiv_crude_men_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_msm_6m, 
  data = analysis_data_hiv_men
)

summary(msm_model_hiv_crude_men_rec)

# adjusted hazard ratio MSM HIV recent
msm_model_hiv_adj1_men_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_msm_6m + hiv_home_6m + hiv_pris_6m, 
  data = analysis_data_hiv_men
)

summary(msm_model_hiv_adj1_men_rec)

msm_model_hiv_adj2_men_rec = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ hiv_msm_6m + hiv_home_6m + hiv_pris_6m + DGINJFQB + AGE1INJ, 
  data = analysis_data_hiv_men
)

summary(msm_model_hiv_adj2_men_rec)

######## Ever HIV estimates #########

# unadjusted hazard ratio sw HIV ever - men and women
sw_model_hiv_crude_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWork, 
  data = analysis_data_hiv_long
)

summary(sw_model_hiv_crude_ever)

# adjusted hazard ratio sw HIV ever - men and women
sw_model_hiv_adj1_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWork + hiv_home_6m + hiv_pris_6m, 
  data = analysis_data_hiv_long
)

summary(sw_model_hiv_adj1_ever)

sw_model_hiv_adj2_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWork + hiv_home_6m + hiv_pris_6m + DGINJFQB + AGE1INJ, 
  data = analysis_data_hiv_long
)

summary(sw_model_hiv_adj2_ever)

# unadjusted hazard ratio sw HIV ever - men
analysis_data_hiv_men <- subset(analysis_data_hiv_long, SEXBRTH == 1)

sw_model_hiv_crude_men_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWork, 
  data = analysis_data_hiv_men
)

summary(sw_model_hiv_crude_men_ever)

# adjusted hazard ratio sw HIV ever - men 
sw_model_hiv_adj1_men_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWork + hiv_home_6m + hiv_pris_6m, 
  data = analysis_data_hiv_men
)

summary(sw_model_hiv_adj1_men_ever)

sw_model_hcv_adj2_men_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWork + hiv_home_6m + hiv_pris_6m + DGINJFQB + AGE1INJ, 
  data = analysis_data_hiv_men
)

summary(sw_model_hiv_adj2_men_ever)

# unadjusted hazard ratio sw HIV ever - women
analysis_data_hiv_women <- subset(analysis_data_hiv_long, SEXBRTH == 2)

sw_model_hiv_crude_women_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWork, 
  data = analysis_data_hiv_women
)

summary(sw_model_hiv_crude_women_ever)

# adjusted hazard ratio sw HIV ever - women 
sw_model_hiv_adj1_women_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWork + hiv_home_6m + hiv_pris_6m, 
  data = analysis_data_hiv_women
)

summary(sw_model_hiv_adj1_women_ever)

sw_model_hiv_adj2_women_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWork + hiv_home_6m + hiv_pris_6m + DGINJFQB + AGE1INJ, 
  data = analysis_data_hiv_women
)

summary(sw_model_hiv_adj2_women_ever)

# unadjusted hazard ratio MSM HIV ever
msm_model_hiv_crude_men_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWMen, 
  data = analysis_data_hiv_men
)

summary(msm_model_hiv_crude_men_ever)

# adjusted hazard ratio MSM HIV ever
msm_model_hiv_adj1_men_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWMen + hiv_home_6m + hiv_pris_6m, 
  data = analysis_data_hiv_men
)

summary(msm_model_hiv_adj1_men_ever)

msm_model_hiv_adj2_men_ever = coxph(
  Surv(time = hiv_start_risk, time2 = hiv_end_risk, event = hiv_tst_rslt) ~ SexWMen + hiv_home_6m + hiv_pris_6m + DGINJFQB + AGE1INJ, 
  data = analysis_data_hiv_men
)

summary(msm_model_hiv_adj2_men_ever)