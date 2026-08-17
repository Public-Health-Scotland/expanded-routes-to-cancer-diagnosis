#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# boss_smr06_cwt.R
# Calum Purdie
# 10/06/2026
# Testing linkage between BoSS, SMR06 and CWT
# Written/run on Posit Workbench
# R version 4.4.2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### 1 Housekeeping ----

# Load packages

library(odbc)
library(dplyr)
library(janitor)
library(tidyr)
library(haven)
library(lubridate)
library(stringr)
library(openxlsx)
library(tidylog)

# Define todays date

todays_date <- today()

# Define output path

output_path <- paste0("/PHI_conf/CancerGroup1/Topics/DiagnosisRoutes/", 
                      "20250505-ExpandedRoutesToDiagnosis/Output")

# Set up a DVPROD connection for CIP

dvprod_con <- dbConnect(odbc(),
                        dsn = "DVPROD",
                        uid = Sys.getenv("USER"),
                        pwd = .rs.askForPassword("What is your LDAP password?"))



### 2 Data Extraction ----

# Read in bowel extract and filter for records where a cancer was diagnosed for
# people invited in 2022

boss_extract <- readRDS(paste0("/PHI_conf/CancerGroup1/Topics/BowelScreening/", 
                             "Data/Programme/2025-11/combined_extract_all.rds")) |> 
  filter(cancer == "01" & year(invdate) == 2022 & !is.na(icd_10)) |> 
  select(chinum, invdate, screres, screresdate, screresdat, precolas, dateprecol,
         coloffered, datecoloff, colperf, datecolperf, colreason, barenctc, 
         barectalt, barctdat, icd_10)

# Read in CWT extract for anyone with a colorectal cancer referred through a
# national screening programme

cwt_extract <- tbl(dvprod_con, dbplyr::in_schema("cancer_intelligence", 
                                                 "cancer_waiting_time")) |> 
  filter(cancer_type_code == "02" & urgency_and_ref_source_code == "02") |> 
  select(pat_upi, urgency_and_ref_source_code, urgency_and_ref_source_desc,
         date_receipt_of_ref, date_decision_to_treat, date_first_trt,
         cancer_type_desc, cancer_type_code, urgent_with_suspicion_of_cancer) |>
  collect()

# Read in SMR06 extract for anyone with a colorectal cancer and method of 
# detection is screening

smr06_extract <- tbl(dvprod_con, dbplyr::in_schema("cancer_intelligence", 
                                                   "cancer_with_staging")) |> 
  filter(tumour_site_icd10_3char_code %in% c("C18", "C19", "C20") & 
         !is.na(pat_upi) & method_first_detection_code == "1") |> 
  select(pat_upi, tumour_id, date_incidence_scottish, date_incidence_encr, 
         encr_event_code, tumour_site_icd10_3char_code, 
         tumour_site_icd10_4char_code, method_first_detection_desc) |>
  collect()



### 3 BoSS and CWT Linkage ----

# Join BoSS and CWT
# Calculate time between screening and referral and also result date to referral
# Flag any matched records
# A match is defined as a screening or result date on or before referral where
# the result was within 365 days before referral
# Filter to select these matched records
# Group by CHI and take the record with the earliest decision to treat date

boss_cwt <- boss_extract |> 
  left_join(cwt_extract, by = c("chinum" = "pat_upi")) |> 
  mutate(time_result_to_referral = time_length(screresdate %--% date_receipt_of_ref, "days"), 
         time_pos_to_referral = time_length(screresdat %--% date_receipt_of_ref, "days")) |> 
  mutate(cwt_match_flag = case_when((time_result_to_referral >= 0 | time_pos_to_referral >= 0) & time_pos_to_referral < 365 ~ 1, 
                                    TRUE ~ 0)) |> 
  filter(cwt_match_flag == 1) |> 
  group_by(chinum) |> 
  filter(date_decision_to_treat == min(date_decision_to_treat) | is.na(date_decision_to_treat)) |> 
  ungroup()

boss_cwt |> count(time_result_to_referral == 0)
boss_cwt |> count(time_pos_to_referral == 0)

# 531 matches where screening result date matches referral date



### 4 BoSS and SMR06 Linkage ----

# Join BoSS and SMR06
# Calculate time between colonoscopy, result date, ct colonography and incidence
# Recode ICD10 codes to only keep alphanumeric characters
# Flag whether ICD10 codes match between BoSS and SMR06
# Flag any matched records
# A match is defined as a colonoscopy, result or ct date on or before incidence
# where the colonoscopy was within 365 days before referral
# Filter to select these matched records
# Group by CHI and prioritise records where four character ICD10 codes match
# Group by CHI and prioritise records where three character ICD10 codes match

boss_smr06 <- boss_extract |> 
  left_join(smr06_extract, by = c("chinum" = "pat_upi")) |> 
  mutate(time_col_to_inc = time_length(datecolperf %--% date_incidence_encr, "days"), 
         time_result_to_inc = time_length(screresdat %--% date_incidence_encr, "days"), 
         time_barct_to_inc = time_length(barctdat %--% date_incidence_encr, "days"), 
         icd_10 = gsub("[^[:alnum:]]", "", icd_10), 
         tumour_site_icd10_4char_code = gsub("[^[:alnum:]]", "", tumour_site_icd10_4char_code), 
         icd10_3_match = if_else(stringr::str_sub(icd_10, 1, 3) == tumour_site_icd10_3char_code, T, F), 
         icd10_4_match = if_else(stringr::str_sub(icd_10, 1, 4) == tumour_site_icd10_4char_code, T, F)) |> 
  mutate(smr06_match_flag = case_when((time_col_to_inc >= 0 | time_result_to_inc >= 0 | time_barct_to_inc >= 0) & (time_col_to_inc < 365 | is.na(time_col_to_inc)) ~ 1, 
                                      TRUE ~ 0)) |> 
  filter(smr06_match_flag == 1) |> 
  group_by(chinum) |> 
  filter((icd10_4_match == max(icd10_4_match)) | is.na(icd_10)) |> 
  ungroup() |> 
  group_by(chinum) |> 
  filter((icd10_3_match == max(icd10_3_match)) | is.na(icd_10)) |> 
  ungroup()

boss_smr06 |> count(icd10_3_match)
boss_smr06 |> count(icd10_4_match)



### 5 Combine BoSS, CWT and SMR06 ----

# Join SMR06 and CWT to BoSS
# Join by CHI and invite date

boss_linked <- boss_extract |> 
  left_join(boss_smr06 |> 
              select(chinum, invdate, 
                     date_incidence_scottish:smr06_match_flag)) |> 
  left_join(boss_cwt |> 
              select(chinum, invdate, urgency_and_ref_source_code:cwt_match_flag))

boss_linked |> count(smr06_match_flag, cwt_match_flag)

# 642 records in both
# 82 in SMR06 but not CWT
# 43 in CWT but not SMR06
# 34 in neither

# Review if the 84+35 have any CWT records
# Review if the 43+35 have any SMR06 records

# Filter for records with no CWT

chi_no_cwt <- boss_linked |> 
  filter(is.na(cwt_match_flag))

# Filter for records with no SMR06

chi_no_smr06 <- boss_linked |> 
  filter(is.na(smr06_match_flag))

# Filter for records with no CWT or SMR06

chi_no_cwt_smr06 <- boss_linked |> 
  filter(is.na(cwt_match_flag) & is.na(smr06_match_flag))



### 6 Review Records with No CWT or SMR06 ----

# Read in CWT data for relevant columns
# Filter for UPIs with no CWT data in main extract
# Join on BoSS data and calculate time between screening and result dates and
# referral date
# Derive the absolute value difference between result and referral date
# Flag any records where cancer type is colorectal
# Prioritise colorectal records over others and then take the closest record
# Add a flag to identify these as other CWT records

cwt_extract_missing <- tbl(dvprod_con, 
                           dbplyr::in_schema("cancer_intelligence", 
                                             "cancer_waiting_time")) |> 
  select(pat_upi, urgency_and_ref_source_code, urgency_and_ref_source_desc,
         date_receipt_of_ref, date_decision_to_treat, date_first_trt,
         cancer_type_desc, cancer_type_code, urgent_with_suspicion_of_cancer) |>
  collect() |> 
  filter(pat_upi %in% chi_no_cwt$chinum) |> 
  left_join(boss_extract |> select(pat_upi = chinum, screresdate, screresdat)) |> 
  mutate(time_pos_to_referral = time_length(screresdat %--% date_receipt_of_ref, "days")) |> 
  mutate(time_pos_to_referral_abs = abs(time_pos_to_referral), 
         cwt_bowel_flag = if_else(cancer_type_code == "02", 1, 0)) |> 
  group_by(pat_upi) |> 
  filter(cwt_bowel_flag == max(cwt_bowel_flag)) |> 
  filter(time_pos_to_referral_abs == min(time_pos_to_referral_abs) | 
         is.na(time_pos_to_referral)) |> 
  ungroup() |> 
  mutate(other_cwt_flag = 1)

cwt_extract_missing |> count(urgency_and_ref_source_desc)

cwt_extract_missing |>
  mutate(time_pos_to_referral_group = case_when(time_pos_to_referral < 0 ~ "CWT referral before BoSS referral", 
                                                 time_pos_to_referral >= 0 & time_pos_to_referral <= 90 ~ "Within 90 days", 
                                                 time_pos_to_referral > 90 & time_pos_to_referral <= 180 ~ "Within 180 days", 
                                                 time_pos_to_referral > 180 & time_pos_to_referral <= 365 ~ "Within one year", 
                                                 time_pos_to_referral > 365 ~ "Over one year")) |> 
  count(time_pos_to_referral_group)

cwt_extract_missing |> count(cancer_type_desc)


# Read in CWT data for relevant columns
# Filter for UPIs with no SMR06 data in main extract
# Join on BoSS data and calculate time between result and incidence date
# Derive the absolute value difference between result and incidence date
# Flag invasive and non-invasive colorectal tumours
# Prioritise colorectal records over non-colorectal
# Prioritise invasive colorectal records over non-invasive colorectal
# Add a flag to identify these as other SMR06 records

smr06_extract_missing <- tbl(dvprod_con, dbplyr::in_schema("cancer_intelligence", 
                                                           "cancer_with_staging")) |> 
  select(pat_upi, date_incidence_scottish, pat_gender_code, 
         tumour_site_icd10_3char_code, tumour_site_icd10_4char_code, 
         pat_age_at_incidence, pat_postcode, tumour_id, date_incidence_encr, 
         method_first_detection_code, method_first_detection_desc, 
         death_out_of_scotland_flag, death_certificate_only_code, 
         encr_event_code) |>
  collect() |> 
  filter(pat_upi %in% chi_no_smr06$chinum) |> 
  left_join(boss_extract |> select(pat_upi = chinum, screresdate, screresdat)) |> 
  mutate(date_incidence_combined = case_when(!is.na(date_incidence_encr) ~ date_incidence_encr, 
                                             is.na(date_incidence_encr) ~ date_incidence_scottish)) |> 
  mutate(time_screening_to_inc = time_length(screresdate %--% date_incidence_combined, "days")) |> 
  mutate(time_screening_to_inc_abs = abs(time_screening_to_inc)) |> 
  mutate(smr06_bowel_flag = case_when(tumour_site_icd10_3char_code %in% c("C18", "C19", "C20") ~ 1, 
                                tumour_site_icd10_4char_code %in% c("D010", "D011", "D012") ~ 2, 
                                TRUE ~ 3)) |> 
  mutate(icd10_flag = case_when(tumour_site_icd10_3char_code %in% c("C18", "C19", "C20") ~ 1, 
                                tumour_site_icd10_4char_code %in% c("D010", "D011", "D012") ~ 1, 
                                TRUE ~ 2)) |> 
  group_by(pat_upi) |> 
  filter(icd10_flag == min(icd10_flag)) |> 
  filter(time_screening_to_inc_abs == min(time_screening_to_inc_abs) | is.na(time_screening_to_inc_abs)) |> 
  filter(smr06_bowel_flag == min(smr06_bowel_flag)) |> 
  ungroup() |> 
  mutate(other_smr06_flag = 1)

smr06_extract_missing |> count(method_first_detection_desc)

smr06_extract_missing |> 
  mutate(time_screening_to_inc_group = case_when(time_screening_to_inc < 0 ~ "Incidence before screening", 
                                                 time_screening_to_inc >= 0 & time_screening_to_inc <= 90 ~ "Within 90 days", 
                                                 time_screening_to_inc > 90 & time_screening_to_inc <= 180 ~ "Within 180 days", 
                                                 time_screening_to_inc > 180 & time_screening_to_inc <= 365 ~ "Within one year", 
                                                 time_screening_to_inc > 365 ~ "Over one year")) |> 
  count(time_screening_to_inc_group)

smr06_extract_missing |> 
  count(tumour_site_icd10_3char_code)



### 7 Combine All Flags ----

# Join the other CWT and other SMR06 flags onto the linked data

boss_linked <- boss_linked |> 
  left_join(cwt_extract_missing |> select(chinum = pat_upi, other_cwt_flag)) |> 
  left_join(smr06_extract_missing |> select(chinum = pat_upi, other_smr06_flag))

boss_linked |> 
  count(smr06_match_flag, cwt_match_flag, other_cwt_flag, other_smr06_flag)

# Every record has some sort of match
# Only 34 records with no direct match
# Seems reasonable for linkage going forward but need to understand why these 34
# have no direct matches



### 8 Output ----

### 8.1 Set Up Workbook ----

# Create a workbook

wb <- createWorkbook()

# Define a header and title style for workbook

hs <- createStyle(fontColour = "#ffffff", fgFill = "#0078D4",
                  halign = "center", valign = "center", 
                  textDecoration = "bold", border = "TopBottomLeftRight")

title_style <- createStyle(fontSize = 14, textDecoration = "bold")


### 8.2 Add SMR06 Matches ----

addWorksheet(wb, "SMR06 Matches", gridLines = FALSE)

writeData(wb, sheet = "SMR06 Matches", startCol = 1, startRow = 1, 
          boss_smr06, 
          borders = "all", headerStyle = hs, colNames = TRUE)


### 8.3 Add CWT Matches ----

addWorksheet(wb, "CWT Matches", gridLines = FALSE)

writeData(wb, sheet = "CWT Matches", startCol = 1, startRow = 1, 
          boss_cwt, 
          borders = "all", headerStyle = hs, colNames = TRUE)


### 8.4 Add SMR06 Non-Matches ----

addWorksheet(wb, "SMR06 Non-Matches", gridLines = FALSE)

writeData(wb, sheet = "SMR06 Non-Matches", startCol = 1, startRow = 1, 
          smr06_extract_missing, 
          borders = "all", headerStyle = hs, colNames = TRUE)


### 8.5 Add CWT Non-Matches ----

addWorksheet(wb, "CWT Non-Matches", gridLines = FALSE)

writeData(wb, sheet = "CWT Non-Matches", startCol = 1, startRow = 1, 
          cwt_extract_missing, 
          borders = "all", headerStyle = hs, colNames = TRUE)


### 8.6 Output ----

saveWorkbook(wb, paste0(output_path, "/boss_cwt_smr06_linkage_cohorts.xlsx"), 
             overwrite = TRUE)
