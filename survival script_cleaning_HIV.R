## load packages
pacman::p_load(dplyr, haven, tidyr, sas7bdat, writexl)

## set wd
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Publications/Sex work and risk of HIV and HCV/Emails to authors/Tijuana data/Data")

## open data for visits 1-6
visit_1 <- read_sas("lfvisit1data.sas7bdat")
visit_2 <- read_sas("lfvisit2data.sas7bdat")
visit_3 <- read_sas("lfvisit3data.sas7bdat")
visit_4 <- read_sas("lfvisit4data.sas7bdat")
visit_5 <- read_sas("lfvisit5data.sas7bdat")
visit_6 <- read_sas("lfvisit6data.sas7bdat")
# Note: data stored separately for each visit

# rename inconsistent homeless variable (capitalised in visits 5 & 6)
visit_5 <- visit_5 %>% 
  rename("Homeless" = "HOMELESS")
visit_6 <- visit_6 %>% 
  rename("Homeless" = "HOMELESS")

# manually fix incorrectly entered hiv test dates for id #319 in visits 2 and 3
visit_2$DateHiv[visit_2$id == "LF319"] <- as.Date("2022-01-07")
visit_2$DateHcv[visit_2$id == "LF319"] <- as.Date("2022-01-07")

visit_3$DateHiv[visit_3$id == "LF319"] <- as.Date("2022-02-01")
visit_3$DateHcv[visit_3$id == "LF319"] <- as.Date("2022-02-01")

## create visit 1-6 df with id, date HIV, result HIV, date HCV, result HCV
baseline = subset(visit_1, select = c(id, AGE, SEXBRTH, AGE1INJ, SexWork, TPrisJail, SexWMen, DGINJFQB, Homeless, TPrisJail6Mo, MetBupPrg6M))
visit_1.2 = subset(visit_1, select = c(id, DateHiv, DateHcv, HivConfStatus, HcvConfStatus, SexWork6Mo, SexWMen6Mo, Homeless, TPrisJail6Mo, MetBupPrg6M, DGINJFQB))
visit_2.2 = subset(visit_2, select = c(id, DateHiv, DateHcv, HivConfStatus, HcvConfStatus, SexWork6Mo, SexWMen6Mo, Homeless, TPrisJail6Mo, MetBupPrg6M, DGINJFQB))
visit_3.2 = subset(visit_3, select = c(id, DateHiv, DateHcv, HivConfStatus, HcvConfStatus, SexWork6Mo, SexWMen6Mo, Homeless, TPrisJail6Mo, MetBupPrg6M, DGINJFQB))
visit_4.2 = subset(visit_4, select = c(id, DateHiv, DateHcv, HivConfStatus, HcvConfStatus, SexWork6Mo, SexWMen6Mo, Homeless, TPrisJail6Mo, MetBupPrg6M, DGINJFQB))
visit_5.2 = subset(visit_5, select = c(id, DateHiv, DateHcv, HivConfStatus, HcvConfStatus, SexWork6Mo, SexWMen6Mo, Homeless, TPrisJail6Mo, MetBupPrg6M, DGINJFQB))
visit_6.2 = subset(visit_6, select = c(id, DateHiv, DateHcv, HivConfStatus, HcvConfStatus, SexWork6Mo, SexWMen6Mo, Homeless, TPrisJail6Mo, MetBupPrg6M, DGINJFQB))

# rename vars

# list of dataframe names
dataframe_names <- c("visit_1.2", "visit_2.2", "visit_3.2", "visit_4.2", "visit_5.2", "visit_6.2")

# loop through the dataframe names
for (i in seq_along(dataframe_names)) {
  current_df <- get(dataframe_names[i])  # Get the current dataframe
  
  # new column names
  new_col_names <- c("DateHiv", "DateHcv", "HivConfStatus", "HcvConfStatus", "SexWork6Mo", 
                     "SexWMen6Mo", "Homeless", "TPrisJail6Mo", "MetBupPrg6M", "DGINJFQB")
  new_col_names <- paste0(new_col_names, "_", i)
  
  # rename columns
  assign(dataframe_names[i], dplyr::rename(current_df, !!!setNames(names(current_df)[2:11], new_col_names)))
  
  message(paste("Columns renamed in", dataframe_names[i]))
}

## add visits 2-6 df as columns
analysis_data <- baseline %>%
  left_join(visit_1.2, by = "id") %>%
  left_join(visit_2.2, by = "id") %>%
  left_join(visit_3.2, by = "id") %>%
  left_join(visit_4.2, by = "id") %>%
  left_join(visit_5.2, by = "id") %>%
  left_join(visit_6.2, by = "id")

## replace sw_6m when equal to 1 in next column
for (i in 1:nrow(analysis_data)) {
  if (!is.na(analysis_data$SexWork6Mo_1[i]) && !is.na(analysis_data$SexWork6Mo_2[i]) && 
      analysis_data$SexWork6Mo_1[i] == "0" && analysis_data$SexWork6Mo_2[i] == "1") {
    analysis_data$SexWork6Mo_1[i] <- "1"
  }
  
  if (!is.na(analysis_data$SexWork6Mo_2[i]) && !is.na(analysis_data$SexWork6Mo_3[i]) && 
      analysis_data$SexWork6Mo_2[i] == "0" && analysis_data$SexWork6Mo_3[i] == "1") {
    analysis_data$SexWork6Mo_2[i] <- "1"
  }
  
  if (!is.na(analysis_data$SexWork6Mo_3[i]) && !is.na(analysis_data$SexWork6Mo_4[i]) && 
      analysis_data$SexWork6Mo_3[i] == "0" && analysis_data$SexWork6Mo_4[i] == "1") {
    analysis_data$SexWork6Mo_3[i] <- "1"
  }
  
  if (!is.na(analysis_data$SexWork6Mo_4[i]) && !is.na(analysis_data$SexWork6Mo_5[i]) && 
      analysis_data$SexWork6Mo_4[i] == "0" && analysis_data$SexWork6Mo_5[i] == "1") {
    analysis_data$SexWork6Mo_4[i] <- "1"
  }
  
  if (!is.na(analysis_data$SexWork6Mo_5[i]) && !is.na(analysis_data$SexWork6Mo_6[i]) && 
      analysis_data$SexWork6Mo_5[i] == "0" && analysis_data$SexWork6Mo_6[i] == "1") {
    analysis_data$SexWork6Mo_5[i] <- "1"
  }
  
}

# create hiv dataset
analysis_data_hiv <- subset(analysis_data, select = c(id, SEXBRTH, AGE1INJ, SexWork, TPrisJail, SexWMen, DGINJFQB, Homeless, TPrisJail6Mo, MetBupPrg6M, 
                                                      HivConfStatus_1, DateHiv_1, 
                                                      HivConfStatus_2, DateHiv_2, 
                                                      HivConfStatus_3, DateHiv_3, 
                                                      HivConfStatus_4, DateHiv_4, 
                                                      HivConfStatus_5, DateHiv_5, 
                                                      HivConfStatus_6, DateHiv_6, 
                                                      SexWork6Mo_1, SexWMen6Mo_1, TPrisJail6Mo_1, Homeless_1, MetBupPrg6M_1, DGINJFQB_1, 
                                                      SexWork6Mo_2, SexWMen6Mo_2, TPrisJail6Mo_2, Homeless_2, MetBupPrg6M_2, DGINJFQB_2, 
                                                      SexWork6Mo_3, SexWMen6Mo_3, TPrisJail6Mo_3, Homeless_3, MetBupPrg6M_3, DGINJFQB_3, 
                                                      SexWork6Mo_4, SexWMen6Mo_4, TPrisJail6Mo_4, Homeless_4, MetBupPrg6M_4, DGINJFQB_4, 
                                                      SexWork6Mo_5, SexWMen6Mo_5, TPrisJail6Mo_5, Homeless_5, MetBupPrg6M_5, DGINJFQB_5))
                                                      
# Replace dates with missing when test is missing
for (i in 2:6) {
  col_date_hiv <- paste0("DateHiv_", i)
  col_hiv_conf_status <- paste0("HivConfStatus_", i)
  
  analysis_data_hiv[[col_date_hiv]][is.na(analysis_data_hiv[[col_hiv_conf_status]])] <- NA
}

## find last recorded hcv and hiv test
analysis_data_hiv$hiv_last <- with(analysis_data_hiv, (pmax(DateHiv_6, DateHiv_5, DateHiv_4, DateHiv_3, DateHiv_2, na.rm = TRUE)))

# drop participants with no f/u test and/or positive at baseline
analysis_data_hiv <- analysis_data_hiv %>% drop_na(hiv_last)

# Calculate time between tests
for (i in 2:6) {
  col_date_hiv_dif <- paste0("DateHiv_dif_", i)
  col_date_hiv <- paste0("DateHiv_", i)
  col_prev_date_hiv <- paste0("DateHiv_", i - 1)
  
  analysis_data_hiv[[col_date_hiv_dif]] <- ifelse(is.na(analysis_data_hiv[[col_date_hiv]]), 
                                                  NA, 
                                                  analysis_data_hiv[[col_date_hiv]] - analysis_data_hiv[[col_prev_date_hiv]])
}

# calculate time at start of row
calculate_times <- function(data, prefix) {
  for (i in 1:5) {
    start_col <- paste0(prefix, "_start_", i)
    end_col <- paste0(prefix, "_end_", i)
    dif_col <- paste0(prefix, "_dif_", i + 1)
    
    if (i == 1) {
      data[[start_col]] <- 0
      data[[end_col]] <- data[[start_col]] + data[[dif_col]]
    } else {
      data[[start_col]] <- data[[paste0(prefix, "_start_", i - 1)]] + data[[paste0(prefix, "_dif_", i)]]
      data[[end_col]] <- data[[start_col]] + data[[dif_col]]
    }
  }
  return(data)
}
analysis_data_hiv <- calculate_times(analysis_data_hiv, "DateHiv")

# reshape test dates wide to long
analysis_data_hiv_date <- analysis_data_hiv %>%
  select(id, DateHiv_2, DateHiv_3, DateHiv_4, DateHiv_5, DateHiv_6) %>%
  pivot_longer(cols = -id, names_to = "date_hiv", values_to = "hiv_tst_dte") %>%
  select(id, hiv_tst_dte) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(test_seq = row_number())

# reshape test results wide to long
analysis_data_hiv_conf <- analysis_data_hiv %>%
  select(id, HivConfStatus_2, HivConfStatus_3, HivConfStatus_4, HivConfStatus_5, HivConfStatus_6) %>%
  pivot_longer(cols = -id, names_to = "conf_hiv", values_to = "hiv_tst_rslt") %>%
  select(id, hiv_tst_rslt) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(test_seq = row_number())

# reshape sex work past 6 months wide to long
analysis_data_hiv_sw <- analysis_data_hiv %>%
  select(id, SexWork6Mo_1, SexWork6Mo_2, SexWork6Mo_3, SexWork6Mo_4, SexWork6Mo_5) %>%
  pivot_longer(cols = -id, names_to = "sw_hiv", values_to = "hiv_sw_6m") %>%
  select(id, hiv_sw_6m) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(test_seq = row_number())

# reshape MSM past 6 months wide to long
analysis_data_hiv_msm <- analysis_data_hiv %>%
  select(id, SexWMen6Mo_1, SexWMen6Mo_2, SexWMen6Mo_3, SexWMen6Mo_4, SexWMen6Mo_5) %>%
  pivot_longer(cols = -id, names_to = "msm_hiv", values_to = "hiv_msm_6m") %>%
  select(id, hiv_msm_6m) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(test_seq = row_number())

# incarceration 
analysis_data_hiv_pris <- analysis_data_hiv %>%
  select(id, TPrisJail6Mo_1, TPrisJail6Mo_2, TPrisJail6Mo_3, TPrisJail6Mo_4, TPrisJail6Mo_5) %>%
  pivot_longer(cols = -id, names_to = "pris_hiv", values_to = "hiv_pris_6m") %>%
  select(id, hiv_pris_6m) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(test_seq = row_number())

# homeless
analysis_data_hiv_home <- analysis_data_hiv %>%
  select(id, Homeless_1, Homeless_2, Homeless_3, Homeless_4, Homeless_5) %>%
  pivot_longer(cols = -id, names_to = "home_hiv", values_to = "hiv_home_6m") %>%
  select(id, hiv_home_6m) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(test_seq = row_number())

# OAT
analysis_data_hiv_oat <- analysis_data_hiv %>%
  select(id, MetBupPrg6M_1, MetBupPrg6M_2, MetBupPrg6M_3, MetBupPrg6M_4, MetBupPrg6M_5) %>%
  pivot_longer(cols = -id, names_to = "oat_hiv", values_to = "hiv_oat_6m") %>%
  select(id, hiv_oat_6m) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(test_seq = row_number())

## injecting frequency
analysis_data_hiv_inj <- analysis_data_hiv %>%
  select(id, DGINJFQB_1, DGINJFQB_2, DGINJFQB_3, DGINJFQB_4, DGINJFQB_5) %>%
  pivot_longer(cols = -id, names_to = "inj_hiv", values_to = "hiv_inj_6m") %>%
  select(id, hiv_inj_6m) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(test_seq = row_number())

# reshape difference between tests (days) wide to long
analysis_data_hiv_start <- analysis_data_hiv %>%
  select(id, DateHiv_start_1, DateHiv_start_2, DateHiv_start_3, DateHiv_start_4, DateHiv_start_5) %>%
  pivot_longer(cols = -id, names_to = "start_hiv", values_to = "hiv_start_risk") %>%
  select(id, hiv_start_risk) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(test_seq = row_number())

# reshape difference between total days at risk wide to long
analysis_data_hiv_end <- analysis_data_hiv %>%
  select(id, DateHiv_end_1, DateHiv_end_2, DateHiv_end_3, DateHiv_end_4, DateHiv_end_5) %>%
  pivot_longer(cols = -id, names_to = "end_hiv", values_to = "hiv_end_risk") %>%
  select(id, hiv_end_risk) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(test_seq = row_number())

# merge long data frames
analysis_data_hiv_long <- analysis_data_hiv_date %>%
  left_join(analysis_data_hiv_conf, by = c("id", "test_seq")) %>%
  left_join(analysis_data_hiv_sw, by = c("id", "test_seq")) %>%
  left_join(analysis_data_hiv_msm, by = c("id", "test_seq")) %>%
  left_join(analysis_data_hiv_start, by = c("id", "test_seq")) %>%
  left_join(analysis_data_hiv_end, by = c("id", "test_seq")) %>%
  left_join(analysis_data_hiv_pris, by = c("id", "test_seq")) %>%
  left_join(analysis_data_hiv_home, by = c("id", "test_seq")) %>%
  left_join(analysis_data_hiv_oat, by = c("id", "test_seq")) %>%
  left_join(analysis_data_hiv_inj, by = c("id", "test_seq")) %>%
  left_join(baseline, by = c("id")) 

# calculate midpoint for cases
analysis_data_hiv_long <- analysis_data_hiv_long %>%
  mutate(midpoint_case = ifelse(hiv_tst_rslt == 0, NA, (hiv_end_risk - hiv_start_risk)/2))

# replace last test with midpoint for cases
analysis_data_hiv_long$hiv_end_risk[!is.na(analysis_data_hiv_long$midpoint_case)] <- 
  analysis_data_hiv_long$hiv_start_risk[!is.na(analysis_data_hiv_long$midpoint_case)] + 
  analysis_data_hiv_long$midpoint_case[!is.na(analysis_data_hiv_long$midpoint_case)]

# keep vars of interest
analysis_data_hiv_long <- subset(analysis_data_hiv_long, select = c(id, AGE, SEXBRTH, AGE1INJ, SexWork, TPrisJail, SexWMen, DGINJFQB, Homeless, TPrisJail6Mo, MetBupPrg6M, hiv_start_risk, hiv_end_risk, midpoint_case, hiv_sw_6m, hiv_msm_6m, hiv_pris_6m, hiv_home_6m, hiv_oat_6m, hiv_inj_6m, hiv_tst_rslt))

# drop rows after positive test
analysis_data_hiv_long <- analysis_data_hiv_long[!is.na(analysis_data_hiv_long$hiv_tst_rslt), ]

# drop if test is missing
analysis_data_hiv_long <- analysis_data_hiv_long[!is.na(analysis_data_hiv_long$hiv_end_risk), ]

# calculate difference between end and start dates
analysis_data_hiv_long <- analysis_data_hiv_long %>%
  mutate(days_risk = hiv_end_risk-hiv_start_risk)

# delete Homeless and TPrisJail6Mo columns
analysis_data_hiv_long <- analysis_data_hiv_long %>%
  select(-Homeless, -TPrisJail6Mo, -MetBupPrg6M, -DGINJFQB)

# rename columns for consistency 
analysis_data_hiv_long <- analysis_data_hiv_long %>%
  rename(SexWMen6Mo = hiv_msm_6m) %>%
  rename(SexWork6Mo = hiv_sw_6m) %>%
  rename(TPrisJail6Mo = hiv_pris_6m) %>%  
  rename(Homeless = hiv_home_6m) %>%
  rename(MetBupPrg6M = hiv_oat_6m) %>%    
  rename(DGINJFQB = hiv_inj_6m) %>%
  rename(hiv_rslt = hiv_tst_rslt)
  
# save data
write_xlsx(analysis_data_hiv_long,"HIV_data_clean.xlsx")


