# QA 
visit_1_hcv <- visit_1 %>%
  filter(HcvConfStatus == 1) ## 68 positive at baseline
visit_1_hcv = subset(visit_1_hcv, select = c(id, DateHcv, HcvConfStatus, SexWork6Mo, SexWMen6Mo))
visit_1_hcv <- visit_1_hcv %>% 
  rename("HcvConfStatus_1" = "HcvConfStatus") %>%
  rename("DateHcv_1" = "DateHcv") %>%
  rename("SexWork6Mo_1" = "SexWork6Mo") %>%
  rename("SexWMen6Mo_1" = "SexWMen6Mo")


visit_2_hcv = subset(visit_2, select = c(id, DateHcv, HcvConfStatus, SexWork6Mo, SexWMen6Mo))
visit_2_hcv <- visit_2_hcv %>% 
  rename("HcvConfStatus_2" = "HcvConfStatus") %>%
  rename("DateHcv_2" = "DateHcv") %>%
  rename("SexWork6Mo_2" = "SexWork6Mo") %>%
  rename("SexWMen6Mo_2" = "SexWMen6Mo")

visit_3_hcv = subset(visit_3, select = c(id, DateHcv, HcvConfStatus, SexWork6Mo, SexWMen6Mo))
visit_3_hcv <- visit_3_hcv %>% 
  rename("HcvConfStatus_3" = "HcvConfStatus") %>%
  rename("DateHcv_3" = "DateHcv") %>%
  rename("SexWork6Mo_3" = "SexWork6Mo") %>%
  rename("SexWMen6Mo_3" = "SexWMen6Mo")

visit_4_hcv = subset(visit_4, select = c(id, DateHcv, HcvConfStatus, SexWork6Mo, SexWMen6Mo))
visit_4_hcv <- visit_4_hcv %>% 
  rename("HcvConfStatus_4" = "HcvConfStatus") %>%
  rename("DateHcv_4" = "DateHcv") %>%
  rename("SexWork6Mo_4" = "SexWork6Mo") %>%
  rename("SexWMen6Mo_4" = "SexWMen6Mo")

visit_5_hcv = subset(visit_5, select = c(id, DateHcv, HcvConfStatus, SexWork6Mo, SexWMen6Mo))
visit_5_hcv <- visit_5_hcv %>% 
  rename("HcvConfStatus_5" = "HcvConfStatus") %>%
  rename("DateHcv_5" = "DateHcv") %>%
  rename("SexWork6Mo_5" = "SexWork6Mo") %>%
  rename("SexWMen6Mo_5" = "SexWMen6Mo")

visit_6_hcv = subset(visit_6, select = c(id, DateHcv, HcvConfStatus, SexWork6Mo, SexWMen6Mo))
visit_6_hcv <- visit_6_hcv %>% 
  rename("HcvConfStatus_6" = "HcvConfStatus") %>%
  rename("DateHcv_6" = "DateHcv") %>%
  rename("SexWork6Mo_6" = "SexWork6Mo") %>%
  rename("SexWMen6Mo_6" = "SexWMen6Mo")

visit_hcv_all <- left_join(visit_1_hcv, visit_2_hcv, by = "id")
visit_hcv_all <- left_join(visit_hcv_all, visit_3_hcv, by = "id")
visit_hcv_all <- left_join(visit_hcv_all, visit_4_hcv, by = "id")
visit_hcv_all <- left_join(visit_hcv_all, visit_5_hcv, by = "id")
visit_hcv_all <- left_join(visit_hcv_all, visit_6_hcv, by = "id")

# cleared virus
visit_1_hcv_clear <- visit_hcv_all  %>% 
  filter(HcvConfStatus_2 == 0 | HcvConfStatus_3 == 0 | HcvConfStatus_4 == 0 | HcvConfStatus_5 == 0 | HcvConfStatus_6 == 0) ## 21 participants cleared
