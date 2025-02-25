## load packages
pacman::p_load(dplyr, haven, tidyr, sas7bdat, writexl, survival)

## set wd
# setwd("C:/Users/joshua.dawe/OneDrive - University of Bristol/Documents/Publications/Sex work and risk of HIV and HCV/Emails to authors/Tijuana data/Data")
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

# manually fix incorrectly entered test dates for id #319 and #159
visit_2$DateHcv[visit_2$id == "LF319"] <- as.Date("2022-01-07")
visit_3$DateHcv[visit_3$id == "LF319"] <- as.Date("2022-02-01")
visit_1$DateHcv[visit_1$id == "LF159"] <- as.Date("2021-04-26")

## subset visits 1-6 df with time varying columns
baseline = subset(visit_1, select = c(id, AGE, SEXBRTH, AGE1INJ, SexWork, TPrisJail, SexWMen))
visit_1 = subset(visit_1, select = c(id, DateHcv, HcvConfStatus, SexWork6Mo, SexWMen6Mo, Homeless, TPrisJail6Mo, MetBupPrg6M, DGINJFQB))
visit_2 = subset(visit_2, select = c(id, DateHcv, HcvConfStatus, SexWork6Mo, SexWMen6Mo, Homeless, TPrisJail6Mo, MetBupPrg6M, DGINJFQB))
visit_3 = subset(visit_3, select = c(id, DateHcv, HcvConfStatus, SexWork6Mo, SexWMen6Mo, Homeless, TPrisJail6Mo, MetBupPrg6M, DGINJFQB))
visit_4 = subset(visit_4, select = c(id, DateHcv, HcvConfStatus, SexWork6Mo, SexWMen6Mo, Homeless, TPrisJail6Mo, MetBupPrg6M, DGINJFQB))
visit_5 = subset(visit_5, select = c(id, DateHcv, HcvConfStatus, SexWork6Mo, SexWMen6Mo, Homeless, TPrisJail6Mo, MetBupPrg6M, DGINJFQB))
visit_6 = subset(visit_6, select = c(id, DateHcv, HcvConfStatus, SexWork6Mo, SexWMen6Mo, Homeless, TPrisJail6Mo, MetBupPrg6M, DGINJFQB))

# Combine the dataframes into a single longitudinal dataframe
analysis_data_hcv_long <- bind_rows(visit_1, visit_2, visit_3, visit_4, visit_5, visit_6) %>%
  arrange(id, DateHcv)

# Sort the combined dataframe by 'id' and 'date' columns
analysis_data_hcv_long <- analysis_data_hcv_long %>%
  arrange(id, DateHcv)

# Merge analysis_data_hcv_long with the baseline dataframe based on 'id'
analysis_data_hcv_long <- left_join(analysis_data_hcv_long, baseline, by = "id")

# find first negative hcv test
analysis_data_hcv_long <- analysis_data_hcv_long %>%
  filter(HcvConfStatus == 0) %>%
  group_by(id) %>%
  summarise(date_first_neg = min(DateHcv, na.rm = TRUE)) %>%
  left_join(analysis_data_hcv_long, ., by = "id")

# remove participants with no negative tests and positive tests before first negative
analysis_data_hcv_long <- analysis_data_hcv_long %>%
  filter(!(is.na(date_first_neg) | DateHcv < date_first_neg))

# find first positive hcv test
analysis_data_hcv_long <- analysis_data_hcv_long %>%
  filter(HcvConfStatus == 1) %>%
  group_by(id) %>%
  summarise(date_first_pos = min(DateHcv, na.rm = TRUE)) %>%
  left_join(analysis_data_hcv_long, ., by = "id")

# delete tests subsequent to first positive hcv test
analysis_data_hcv_long <- analysis_data_hcv_long %>%
  filter(is.na(date_first_pos) | DateHcv <= date_first_pos)

## testing data

# keep columns of interest
testing_df <- subset(analysis_data_hcv_long, select = c(id, DateHcv, HcvConfStatus)) 

# create lag of test (using lead function)
testing_df <- testing_df %>%
  arrange(id, DateHcv) %>%  
  group_by(id) %>%
  mutate(DateHcv_end = lead(DateHcv),
         hcv_rslt = lead(HcvConfStatus))  

# rename start date to avoid confusion
testing_df <- testing_df %>%
  rename(DateHcv_start = DateHcv)

# remove Nth row for each participant (because of lead)
testing_df <- testing_df %>%
  group_by(id) %>%
  filter(row_number() < n())

# calculate difference between end and start dates
testing_df <- testing_df %>%
  mutate(days_risk = DateHcv_end-DateHcv_start)

testing_df <- testing_df %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(days_risk_end = cumsum(as.numeric(days_risk, units = "days")))

testing_df <- testing_df %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(days_risk_start = lag(days_risk_end, default = 0))

testing_df <- subset(testing_df, select = c(id, days_risk, days_risk_start, days_risk_end, hcv_rslt)) 

write_xlsx(testing_df,"HCV_testing_data.xlsx")

## exposure data

# replace lifetime sex work with 1 if recent sex work is 1, also replacing subsequent rows with 1
analysis_data_hcv_long <- analysis_data_hcv_long %>%
  group_by(id) %>%
  arrange(DateHcv) %>%
  mutate(SexWork = ifelse(SexWork6Mo == 1 | lag(SexWork6Mo, default = 0) == 1, 1, SexWork)) %>%
  ungroup()

# keep columns of interest
exposure_df <- subset(analysis_data_hcv_long, select = c(id, SexWork6Mo, SexWMen6Mo, DGINJFQB, Homeless, TPrisJail6Mo, MetBupPrg6M, SexWMen, SexWork, AGE, SEXBRTH, AGE1INJ)) 

exposure_df <- exposure_df %>%
  arrange(id) %>%
  group_by(id) %>%
  slice(-n())  

# create analysis df

# sequence by id for merge
exposure_df <- exposure_df %>%
  arrange(id) %>%
  mutate(id_seq = row_number())

testing_df <- testing_df %>%
  arrange(id) %>%
  mutate(id_seq = row_number())

# merge
analysis_data_hcv_clean <- left_join(exposure_df, testing_df, by = c("id", "id_seq"))

# save data
write_xlsx(analysis_data_hcv_clean,"HCV_data_clean.xlsx")
