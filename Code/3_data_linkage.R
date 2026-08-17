#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3_data_linkage.R
# Calum Purdie
# 29/07/2026
# Join data for routes to diagnosis
# Written/run on Posit Workbench
# R version 4.4.2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


### 1 Housekeeping ----

# Run housekeeping script to get dates

source(here::here("Code/1_housekeeping.R"))



### 2 Read in Saved Files ----

cancer_data <- readRDS(paste0(data_path, "/2026-07-30_cancer_data.rds"))

acute_data <- readRDS(paste0(data_path, "/2026-07-30_acute_data.rds")) |> 
  mutate(acute_id = row_number())

cwt_data <- readRDS(paste0(data_path, "/2026-07-30_cwt_data.rds")) |> 
  mutate(cwt_id = row_number())

outpatient_data <- readRDS(paste0(data_path, "/2026-07-30_outpatient_data.rds")) |> 
  mutate(outpatient_id = row_number())



### 3 Define Initial Endpoints ----

#### 3.1 Acute ----

# Join cancer and acute
# Calculate the time between admission date and incidence date
# Create an ID field for admission type, with emergency prioritised
# Arrange data and group by UPI and CIS marker
# Flag the first admission type for each group
# Categorise records into emergency, elective, transfer or non-emergency
# This uses the standard ICBP definition of an emergency within 30 days of
# incidence plus the PHS addition of diagnoses where the first admission was an
# emergency and the continuous inpatient stay carried into the 30 day window
# Create an emergency flag ID, with records prioritised as emergency > elective > transfer > non-emergency

cancer_acute_flags <- cancer_data |>  
  inner_join(acute_data) |>  
  mutate(time_adm_to_encr = time_length(date_admission %--% date_incidence_combined, "days")) |> 
  mutate(admission_flag_id = case_when(admission == "emergency" ~ 1, 
                                       admission == "elective" ~ 2)) |> 
  arrange(pat_upi, date_admission, admission_flag_id) |>  
  group_by(pat_upi, cis) |>  
  mutate(first_admission_type = first(admission)) |>  
  ungroup() |>  
  mutate(acute_admission_route = case_when(
    time_adm_to_encr >= 0 & time_adm_to_encr <= 30 & first_admission_type == "emergency" ~ "Emergency",
    time_adm_to_encr >= 0 & time_adm_to_encr <= 30 & admission == "emergency" ~ "Emergency", 
    time_adm_to_encr >= 0 & time_adm_to_encr <= 30 & admission == "elective" & admission_type_code == "11" ~ "Elective", 
    time_adm_to_encr >= 0 & time_adm_to_encr <= 30 & admission == "elective" & admission_type_code == "18" ~ "Transfer", 
    TRUE ~ "Non-Emergency")) |> 
  mutate(acute_admission_route_id = case_when(acute_admission_route == "Emergency" ~ 1, 
                                              acute_admission_route == "Elective" ~ 2, 
                                              acute_admission_route == "Transfer" ~ 3, 
                                              acute_admission_route == "Non-Emergency" ~ 4))

# Filter data to keep records where time between admission and incidence is within 183 days
# Flag if the admission and incidence date are the same day or within 30 days
# Add a count of tumour ID

cancer_acute_flags <- cancer_acute_flags |> 
  filter(time_adm_to_encr >= 0 & time_adm_to_encr <= 183) |> 
  mutate(acute_same_day = if_else(time_adm_to_encr == 0, 1, 0), 
         acute_within_30 = if_else(time_adm_to_encr >= 0 & time_adm_to_encr <= 30, 1, 0)) |> 
  add_count(tumour_id, name = "n_tumours")

# Group by tumour ID and calculate the minimum and maximum admission dates
# Create a flag to identify tumours which link to any emergency admission
# Calculate the number of unique admission types for each tumour
# Filter to exclude cases where a tumour links to multiple admissions, the admission is flagged as a transfer, 
# the admission date is not the first admission date and there is more than one type of admission
# These are transfers which happen after a previous admission and we are more interested in what happened at the original admission

cancer_acute_flags <- cancer_acute_flags |> 
  group_by(tumour_id) |> 
  mutate(acute_min_date = min(date_admission, na.rm = T), 
         acute_max_date = max(date_admission, na.rm = T)) |> 
  mutate(ip_em_flag = if_else(any(acute_admission_route == "Emergency"), 1, 0)) |> 
  mutate(unique_adm_types = n_distinct(acute_admission_route)) |>
  filter(!(n_tumours != 1 & acute_admission_route == "Transfer" & date_admission != acute_min_date & unique_adm_types != 1)) |>
  ungroup()

# Drop n_tumours and recount this variable now we have excluded some transfers
# Filter to exclude cases where a tumour links to multiple admissions, the admission is flagged as a transfer, 
# the admission date is the first admission date and there is more than one type of admission
# These are transfers which happen on the same day as the first admission and we prioritise the admission itself rather than the transfer

cancer_acute_flags <- cancer_acute_flags |> 
  select(-n_tumours) |> 
  add_count(tumour_id, name = "n_tumours") |> 
  group_by(tumour_id) |> 
  filter(!(n_tumours != 1 & acute_admission_route == "Transfer" & date_admission == acute_max_date & unique_adm_types != 1)) |>
  ungroup() |> 
  select(-n_tumours)

# Group by tumour ID
# Filter to take the prioritised admission type
# Filter to take admission closest to incidence date
# Filter to take take the prioritised admission if both on the same day
# Filter to take earliest discharge date
# Select columns and take distinct data
# Group by tumour ID and take the first row for each
# This deduplicates for cases where a patient has two of the same admission type on the same admission and discharge date
# Add an acute flag for use when linking to other datasets

cancer_acute_flags_filtered <- cancer_acute_flags |> 
  group_by(tumour_id) |> 
  filter(acute_admission_route_id == min(acute_admission_route_id)) |> 
  filter(time_adm_to_encr == min(time_adm_to_encr)) |> 
  filter(admission_flag_id == min(admission_flag_id)) |> 
  filter(date_discharge == min(date_discharge)) |> 
  ungroup() |> 
  select(tumour_id, acute_min_date, acute_max_date, time_adm_to_encr, 
         acute_same_day, acute_within_30, acute_id) |> 
  distinct() |> 
  group_by(tumour_id) |> 
  slice(1) |> 
  ungroup() |> 
  mutate(acute_flag = 1)


#### 3.2 Outpatient ----

# Join cancer and outpatient and take distinct rows
# Filter to only keep records where a patient attended
# Calculate time between clinic date and ENCR incidence date
# Flag any records referred by GP or A&E
# Filter to keep records where incidence date is within 183 days of clinic date
# Flag if the clinic date and incidence date are the same day or within 28 days

cancer_outpatient_flags <- cancer_data |> 
  inner_join(outpatient_data) |> 
  distinct() |> 
  mutate(time_cd_to_encr = time_length(clinic_date %--% date_incidence_combined, "days"), 
         gp_flag = case_when(referral_source_code == "1" ~ 1, 
                             TRUE ~ 0), 
         em_flag = case_when(referral_source_code == "A" ~ 1, 
                             TRUE ~ 0)) |> 
  filter(time_cd_to_encr >= 0 & time_cd_to_encr <= 183) |> 
  mutate(outpatient_same_day = case_when(time_cd_to_encr == 0 ~ 1, 
                                    TRUE ~ 0), 
         outpatient_within_28 = case_when(time_cd_to_encr >= 0 & time_cd_to_encr <= 28 ~ 1, 
                                     TRUE ~ 0))

# Group by tumour ID
# Calculate the minimum and maximum clinic dates
# Filter to prioritise emergency referrals
# Filter to prioritise GP referrals
# Filter to take the minimum time between clinic date and incidence
# Filter to take the minimum clinic type - 1 and 2 are consultant-led (1 = consultant, 2 = dentist)
# Others are non-consultant such as nurse-led clinics and AHPs (3 = nurse, unsure on 4)
# Select columns and take distinct data
# Group by tumour ID and take the first row for each
# This deduplicates for cases where a patient has two of the same clinic types on the same clinic date
# Add an outpatient flag for use when linking to other datasets

cancer_outpatient_flags_filtered <- cancer_outpatient_flags |> 
  group_by(tumour_id) |> 
  mutate(outpatient_min_date = min(clinic_date, na.rm = T), 
         outpatient_max_date = max(clinic_date, na.rm = T)) |>
  filter(em_flag == max(em_flag)) |> 
  filter(gp_flag == max(gp_flag)) |> 
  filter(time_cd_to_encr == min(time_cd_to_encr)) |>
  filter(clinic_type_code == min(clinic_type_code)) |>
  ungroup() |> 
  select(tumour_id, outpatient_min_date, outpatient_max_date, time_cd_to_encr, 
         outpatient_same_day, outpatient_within_28, outpatient_id) |> 
  distinct() |> 
  group_by(tumour_id) |> 
  slice(1) |> 
  ungroup() |> 
  mutate(outpatient_flag = 1)


#### 3.3 CWT ----

# Add a cancer type code to cancer to match to CWT
# Join cancer and CWT

cancer_cwt_flags <- cancer_data |> 
  mutate(cancer_cancer_type_code = case_when(tumour_site_icd10_3char_code == "C50" ~ "01", 
                                             tumour_site_icd10_3char_code %in% c("C18", "C19", "C20") ~ "02", 
                                             tumour_site_icd10_3char_code %in% c("C00", "C01", "C02", 
                                                                                 "C03", "C04", "C05", 
                                                                                 "C06", "C07", "C08", 
                                                                                 "C09", "C10", "C11", 
                                                                                 "C12", "C13", "C14", 
                                                                                 "C30", "C31", "C32") ~ "03",
                                             tumour_site_icd10_4char_code == "C760" ~ "03", 
                                             tumour_site_icd10_3char_code %in% c("C33", "C34") ~ "04",
                                             tumour_site_icd10_3char_code %in% c("C81", "C82", "C83", 
                                                                                 "C84", "C85") ~ "05",
                                             tumour_site_icd10_3char_code == "C43" ~ "06",
                                             tumour_site_icd10_3char_code %in% c("C48", "C56") ~ "07",
                                             tumour_site_icd10_3char_code %in% c("C22", "C23", "C24", 
                                                                                 "C25") ~ "08",
                                             tumour_site_icd10_3char_code %in% c("C15", "C16") ~ "09",
                                             tumour_site_icd10_4char_code == "C170" ~ "09", 
                                             tumour_site_icd10_3char_code == "C67" ~ "10",
                                             tumour_site_icd10_3char_code == "C61" ~ "11",
                                             tumour_site_icd10_3char_code %in% c("C60", "C62", "C63", 
                                                                                 "C64", "C65", "C66", 
                                                                                 "C68") ~ "12",
                                             tumour_site_icd10_3char_code == "C53" ~ "13",
                                             tumour_site_icd10_4char_code == "C541" ~ "14",
                                             tumour_site_icd10_3char_code %in% c("C70", "C71", "C72") ~ "15",
                                             tumour_site_icd10_4char_code == "C753" ~ "15",
                                             tumour_site_icd10_4char_code == "C900" ~ "17",
                                             tumour_site_icd10_3char_code == "C45" ~ "18")) |> 
  inner_join(cwt_data, by = c("pat_upi" = "pat_upi"))

# Calculate time between date referral was received and ENCR incidence date
# Calculate time between date decision to treat and ENCR incidence date
# Flag if the cancer type matches between cancer and CWT
# Take distinct rows
# Filter to keep records where date receipt of referral was within 183 days of
# incidence date, or where date receipt of referral is blank
# Flag if the referral and incidence date are the same day or within 28 days

cancer_cwt_flags <- cancer_cwt_flags |> 
  mutate(time_dror_to_encr = time_length(date_receipt_of_ref %--% date_incidence_combined, "days"), 
         # time_dtt_to_encr = time_length(date_decision_to_treat %--% date_incidence_combined, "days"), 
         type_match = if_else(cancer_cancer_type_code == cancer_type_code, 1, 0)) |> 
  filter((time_dror_to_encr >= 0 & time_dror_to_encr <= 183) | is.na(time_dror_to_encr)) |> 
  mutate(cwt_same_day = case_when(time_dror_to_encr == 0 ~ 1, 
                                  TRUE ~ 0), 
         cwt_within_28 = case_when(time_dror_to_encr >= 0 & time_dror_to_encr <= 28 ~ 1, 
                                   TRUE ~ 0))

# Group by tumour ID
# Calculate the minimum and maximum referral dates
# Filter to prioritise records where the cancer types match
# Filter to prioritise urgent referrals
# Filter to prioritise the referral closest to incidence date
# Filter to prioritise referrals on the same day or within 28 days
# Select columns and take distinct data
# Group by tumour ID and take the first row for each
# This deduplicates for cases where a patient has two of the same CWT referrals on the same referral date
# Add a CWT flag for use when linking to other datasets

cancer_cwt_flags_filtered <- cancer_cwt_flags |> 
  group_by(tumour_id) |> 
  mutate(cwt_min_date = min(date_receipt_of_ref, na.rm = T), 
         cwt_max_date = max(date_receipt_of_ref, na.rm = T)) |>
  filter(type_match == max(type_match) | is.na(type_match)) |>
  filter(urgency_and_ref_source_code == min(urgency_and_ref_source_code)) |> 
  filter(time_dror_to_encr == min(time_dror_to_encr) | is.na(time_dror_to_encr)) |> 
  filter(cwt_same_day == max(cwt_same_day)) |> 
  filter(cwt_within_28 == max(cwt_within_28)) |> 
  ungroup() |> 
  select(tumour_id, cwt_min_date, cwt_max_date, time_dror_to_encr, cwt_same_day, 
         cwt_within_28, cwt_id) |> 
  distinct() |> 
  group_by(tumour_id) |> 
  slice(1) |> 
  ungroup() |> 
  mutate(cwt_flag = 1)



### 4 Join Endpoints Together ----

# Join cancer, acute, outpatient and CWT together
# Flag any DCO or unknown endpoints for patients who don't appear in other data

cancer_endpoints <- cancer_data |> 
  left_join(cancer_acute_flags_filtered) |> 
  left_join(cancer_outpatient_flags_filtered) |> 
  left_join(cancer_cwt_flags_filtered) |> 
  mutate(endpoint = case_when(is.na(acute_flag) & 
                              is.na(outpatient_flag) & 
                              is.na(cwt_flag) & 
                              death_certificate_only_code == 1 ~ "death_certificate_only", 
                              is.na(acute_flag) & 
                              is.na(outpatient_flag) & 
                              is.na(cwt_flag) & 
                              death_certificate_only_code != 1 ~ "unknown"))

# Flag records where there is an inpatient admission on the same day as incidence
# and also either an outpatient attendance or referral on that date as inpatient
# Flag records where there is an inpatient admission within 30 days as inpatient
# Flag records where there is not an inpatient admission within 30 days but
# there is an outpatient attendance or referral within 28 days as outpatient

cancer_endpoints <- cancer_endpoints |> 
  mutate(endpoint = case_when(acute_same_day == 1 & (outpatient_same_day == 1 | cwt_same_day == 1) ~ "ip", 
                              TRUE ~ endpoint)) |> 
  mutate(endpoint = case_when(acute_within_30 == 1 ~ "ip", 
                              (acute_within_30 == 0 | is.na(acute_within_30)) & (outpatient_within_28 == 1 | cwt_within_28 == 1) ~ "op", 
                              TRUE ~ endpoint))

# Flag records which have an inpatient admission but no outpatient attendance or referral as inpatient
# Flag records which have an outpatient attendance but no inpatient admission or referral as outpatient
# Flag records which have a referral but no inpatient admission or outpatient attendance as outpatient

cancer_endpoints <- cancer_endpoints |> 
  mutate(endpoint = case_when(is.na(endpoint) & !is.na(acute_flag) & is.na(outpatient_flag) & is.na(cwt_flag) ~ "ip", 
                              is.na(endpoint) & is.na(acute_flag) & !is.na(outpatient_flag) & is.na(cwt_flag) ~ "op", 
                              is.na(endpoint) & is.na(acute_flag) & is.na(outpatient_flag) & !is.na(cwt_flag) ~ "op", 
                              TRUE ~ endpoint))

# Calculate the maximum date across acute, outpatient and CWT
# Flag records where the maximum date is an inpatient admission as inpatient
# Flag records where the maximum date is an outpatient attendance as outpatient
# Flag records where the maximum date is a referral as outpatient

cancer_endpoints <- cancer_endpoints |> 
  mutate(max_date = pmax(acute_max_date, outpatient_max_date, cwt_max_date, na.rm = T)) |> 
  mutate(endpoint = case_when(is.na(endpoint) & max_date == acute_max_date ~ "ip", 
                              is.na(endpoint) & max_date == outpatient_max_date ~ "op", 
                              is.na(endpoint) & max_date == cwt_max_date ~ "op", 
                              TRUE ~ endpoint))

# Count endpoints to ensure every tumour has been categorised

cancer_endpoints |> count(endpoint)



### 5 Inpatient Endpoint ----

# Filter to get inpatient endpoint flagged records
# Join acute, outpatient and CWT onto the inpatient endpoint records

ip_data <- cancer_endpoints |> 
  filter(endpoint == "ip") |> 
  left_join(cancer_acute_flags) |> 
  left_join(cancer_outpatient_flags) |> 
  left_join(cancer_cwt_flags)

# Flag elective inpatient admissions where there is an earlier outpatient attendance or referral as outpatient endpoint
# Flag non-emergency admissions where there is an earlier outpatient attendance or referral as outpatient endpoint

ip_data <- ip_data |> 
  mutate(endpoint = case_when(acute_admission_route == "Elective" & outpatient_min_date < acute_min_date ~ "op", 
                              acute_admission_route == "Elective" & cwt_min_date < acute_min_date ~ "op", 
                              acute_admission_route == "Non-Emergency" & outpatient_min_date < acute_min_date ~ "op", 
                              acute_admission_route == "Non-Emergency" & cwt_min_date < acute_min_date ~ "op", 
                              TRUE ~ endpoint))

# Categorise emergency admissions as the emergency route
# Categorise transfers as the elective route
# Flag elective inpatient admissions where there is no outpatient attendance or referral as the elective route
# Flag elective inpatient admissions where the outpatient attendance or referral are later as the elective route
# Flag non-emergency inpatient admissions where there is no outpatient attendance or referral as the elective route
# Flag non-emergency inpatient admissions where the outpatient attendance or referral are later as the elective route

ip_data <- ip_data |> 
  mutate(route = case_when(acute_admission_route == "Emergency" ~ "emergency", 
                           acute_admission_route == "Transfer" ~ "elective", 
                           acute_admission_route == "Elective" & is.na(outpatient_min_date) & is.na(cwt_min_date) ~ "elective", 
                           acute_admission_route == "Elective" & outpatient_min_date >= acute_min_date ~ "elective", 
                           acute_admission_route == "Elective" & cwt_min_date >= acute_min_date ~ "elective",
                           acute_admission_route == "Non-Emergency" & is.na(outpatient_min_date) & is.na(cwt_min_date) ~ "elective", 
                           acute_admission_route == "Non-Emergency" & outpatient_min_date >= acute_min_date ~ "elective", 
                           acute_admission_route == "Non-Emergency" & cwt_min_date >= acute_min_date ~ "elective"))

# Count endpoint and routes

ip_data |> count(endpoint)
ip_data |> count(route)
ip_data |> count(route, endpoint)



### 6 Outpatient Endpoint ----

# Filter for outpatient endpoint records
# Add on any records reclassified as outpatient through the inpatient processing
# Join on outpatient and CWT

op_data <- cancer_endpoints |> 
  filter(endpoint == "op") |> 
  bind_rows(ip_data |>
              filter(endpoint == "op") |>
              select(tumour_id:max_date)) |>
  left_join(cancer_outpatient_flags) |> 
  left_join(cancer_cwt_flags)

# Categorise emergency outpatient referrals as the emergency route
# Categorise GP outpatient referrals as the GP referral route
# Categorise GP CWT referrals as the GP referral route
# Categorise any other outpatient attendances as the other outpatient route

op_data <- op_data |> 
  mutate(route = case_when(em_flag == 1 ~ "emergency", 
                           gp_flag == 1 ~ "gp_referral", 
                           urgency_and_ref_source_code == "16" ~ "gp_referral", 
                           urgency_and_ref_source_code == "17" ~ "other_outpatient", 
                           !is.na(clinic_type_code) ~ "other_outpatient"))

# Count route

op_data |> count(route)

# Check any blank routes are urgent suspected cancers
# These will be populated later

op_data |> count(is.na(route), urgent_with_suspicion_of_cancer)



### 7 Assign Final Routes ----

# Filter to take any DCO or unknown endpoints and define their route
# Add on inpatient and outpatient endpoint data

final_endpoints <- cancer_endpoints |> 
  filter(endpoint %in% c("death_certificate_only", "unknown")) |> 
  mutate(route = endpoint) |> 
  bind_rows(ip_data |> 
              filter(endpoint == "ip") |> 
              select(tumour_id:max_date, route)) |> 
  bind_rows(op_data |> 
              select(tumour_id:max_date, route))

# Check all rows are included

nrow(cancer_data) - nrow(final_endpoints)

cancer_data |> select(tumour_id) |> n_distinct()
final_endpoints |> select(tumour_id) |> n_distinct()

# Join on acute and CWT flags to help derive final routes
# Calculate time between admission date and date decision to treat
# If a patient has both an emergency inpatient admission and urgent suspected
# cancer referral, class them as an emergency if the admission is within 30 days
# of the date decision to treat

final_endpoints <- final_endpoints |> 
  left_join(cancer_acute_flags |> 
              select(tumour_id, acute_id, acute_admission_route, date_admission)) |> 
  left_join(cancer_cwt_flags |> 
              select(tumour_id, cwt_id, urgent_with_suspicion_of_cancer, 
                     date_decision_to_treat)) |> 
  mutate(time_adm_to_dec = time_length(date_admission %--% date_decision_to_treat, 
                                       "days"), 
         em_flag_priority = case_when(acute_admission_route == "Emergency" & 
                                        urgent_with_suspicion_of_cancer == "YES" & 
                                        time_adm_to_dec >= 0 & 
                                        time_adm_to_dec <= 30 ~ 1, 
                                      TRUE ~ 0))

# Derive the final routes
# Categorise death certificate only records as that route
# Even if these diagnoses have an alternative route, we have chosen to keep them as DCO but maybe check these with the registry
# Categorise screen detected cancers as the screen detected route
# Categorise any records where the method of first detection is first outpatient clinic as the other outpatient route
# Select columns

final_endpoints <- final_endpoints |> 
  mutate(route = case_when(death_certificate_only_code == 1 ~ "death_certificate_only", 
                           method_first_detection_code == 1 & tumour_site_icd10_3char_code %in% c("C18", "C19", "C20", "C50", "C53") ~ "screen_detected", 
                           em_flag_priority == 1 ~ "emergency", 
                           urgent_with_suspicion_of_cancer == "YES" ~ "urgent_suspected_cancer", 
                           # Unsure on below as technically a 1st outpatient clinic because of the malignancy could be after diagnosis?
                           # route == "unknown" & incidence_year >= 2019 & incidence_year <= 2021 & encr_event == 3 ~ "other_outpatient",
                           # route == "unknown" & incidence_year >= 2022 & encr_event == 4 ~ "other_outpatient",
                           TRUE ~ route)) |> 
  select(tumour_id, pat_upi, date_incidence_combined, 
         date_incidence_combined_year, tumour_site_icd10_3char_code, 
         tumour_site_icd10_4char_code, pat_gender_code, pat_age_at_incidence, 
         method_first_detection_code, death_certificate_only_code, route)

final_endpoints |> count(route)
