#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2_data_extraction
# Calum Purdie
# 23/10/2025
# Extract and join data for routes to diagnosis
# Written/run on Posit Workbench
# R version 4.4.2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


### 1 Housekeeping ----

# Run housekeeping script to get dates

source(here::here("Code/1_housekeeping.R"))

# Set up a DVPROD connection for CIP

dvprod_con <- dbConnect(odbc(),
                        dsn = "DVPROD",
                        uid = Sys.getenv("USER"),
                        pwd = .rs.askForPassword("What is your LDAP password?"))



### 2 Cancer with Staging ----

# Extract data from cancer_with_staging view
# Exclude records with blank upi_number and where out_of_scotland is not blank
# Exclude C44 records and only keep invasive tumours
# Filter for relevant incidence dates
# Select columns and collect data

cancer_data <- tbl(dvprod_con, 
                  dbplyr::in_schema("cancer_intelligence", 
                                    "cancer_with_staging")) |> 
  filter(!is.na(pat_upi) & is.na(death_out_of_scotland_flag)) |> 
  filter(tumour_site_icd10_3char_code != "C44") |> 
  filter(substr(tumour_site_icd10_3char_code, 1, 1) == "C") |> 
  filter(date_incidence_combined_year >= start & 
         date_incidence_combined_year <= end) |> 
  select(tumour_id, pat_upi, date_incidence_combined, 
         tumour_site_icd10_3char_code, tumour_site_icd10_4char_code, 
         pat_gender_code, pat_age_at_incidence, date_incidence_combined_year, 
         method_first_detection_code, death_certificate_only_code) |>
  collect()

# Create 3 character cancer_site code based on icd10s_cancer_site
# Define incidence type, format date and calculate incidence year
# Remove men with breast or gynaecological cancer and women with prostate cancer

cancer_data <- cancer_data |> 
  mutate(tumour_type = case_when(tumour_site_icd10_3char_code %in% c("C00", "C01", "C02",
                                                                     "C03", "C04", "C05",
                                                                     "C06", "C07", "C08",
                                                                     "C09", "C10", "C11",
                                                                     "C12", "C13", "C14",
                                                                     "C30", "C31", "C32") ~ "Head and Neck",
                                  tumour_site_icd10_3char_code %in% c("C15") ~ "Oesophagus",
                                  tumour_site_icd10_3char_code %in% c("C16") ~ "Stomach",
                                  tumour_site_icd10_3char_code %in% c("C18", "C19", "C20") ~ "Colorectal",
                                  tumour_site_icd10_3char_code %in% c("C22") ~ "Liver",
                                  tumour_site_icd10_3char_code %in% c("C25") ~ "Pancreas",
                                  tumour_site_icd10_3char_code %in% c("C33", "C34") ~ "Trachea, Bronchus and Lung",
                                  tumour_site_icd10_3char_code %in% c("C43") ~ "Malignant Melanoma of the Skin",
                                  tumour_site_icd10_3char_code %in% c("C50") ~ "Breast",
                                  tumour_site_icd10_3char_code %in% c("C53") ~ "Cervix Uteri",
                                  tumour_site_icd10_3char_code %in% c("C54") ~ "Corpus Uteri",
                                  tumour_site_icd10_3char_code %in% c("C56") ~ "Ovary",
                                  tumour_site_icd10_3char_code %in% c("C61") ~ "Prostate",
                                  tumour_site_icd10_3char_code %in% c("C64", "C65") ~ "Kidney",
                                  tumour_site_icd10_3char_code %in% c("C67") ~ "Bladder",
                                  tumour_site_icd10_3char_code %in% c("C73") ~ "Thyroid", 
                                  tumour_site_icd10_3char_code %in% c("C70", "C71", "C72") |
                                  tumour_site_icd10_4char_code %in% c("C751", "C752", "C753") ~ "Brain and other CNS", 
                                  tumour_site_icd10_3char_code %in% c("C91", "C92", "C93", "C94", "C95") ~ "Leukaemias")) |> 
  mutate(keep_flag = case_when(tumour_type == "Breast" & pat_gender_code == "1" ~ 0, 
                               tumour_type == "Cervix Uteri" & pat_gender_code == "1" ~ 0, 
                               tumour_type == "Ovary" & pat_gender_code == "1" ~ 0, 
                               tumour_type == "Prostate" & pat_gender_code == "2" ~ 0, 
                               tumour_type == "Corpus Uteri" & pat_gender_code == "1" ~ 0, 
                               TRUE ~ 1)) |> 
  filter(keep_flag == 1)

gc()

saveRDS(cancer_data, paste0(data_path, "/", todays_date, "_cancer_data.rds"))



### 3 Acute ----

# Extract data from acute view
# Exclude blank UPIs and filter for relevant admission dates
# Select columns and collect data
# Filter for UPIs in cancer_data
# Define emergency and elective admissions

acute_data <- tbl(dvprod_con, 
                  dbplyr::in_schema("cancer_intelligence", "acute")) |> 
  filter(!is.na(pat_upi)) |> 
  filter(date_admission >= adm_start & date_admission <= adm_end) |> 
  select(pat_upi, date_admission, date_discharge, 
         treatment_nhs_board_code_current, admission_type_code, cis) |>
  collect() |> 
  filter(pat_upi %in% cancer_data$pat_upi) |> 
  mutate(admission = case_when(admission_type_code %in% c("20", "21", "22", "30", 
                                                          "31", "32", "33", "34", 
                                                          "35", "36", "38", 
                                                          "39") ~ "emergency", 
                               admission_type_code %in% c("10", "11", "12", "18", 
                                                          "19") ~ "elective"))

gc()

saveRDS(acute_data, paste0(data_path, "/", todays_date, "_acute_data.rds"))



### 4 Cancer Waiting Times ----

# Extract data from cancer_waiting_time view
# Exclude blank UPIs, select columns and collect data
# Filter for UPIs in cancer_data

cwt_data <- tbl(dvprod_con, dbplyr::in_schema("cancer_intelligence", 
                                              "cancer_waiting_time")) |> 
  filter(!is.na(pat_upi)) |> 
  select(pat_upi, urgency_and_ref_source_code, urgency_and_ref_source_desc,
         date_receipt_of_ref, date_decision_to_treat, date_first_trt,
         cancer_type_desc, cancer_type_code, urgent_with_suspicion_of_cancer) |>
  collect() |> 
  clean_names() |> 
  filter(pat_upi %in% cancer_data$pat_upi)

gc()

saveRDS(cwt_data, paste0(data_path, "/", todays_date, "_cwt_data.rds"))



### 5 Outpatient ----

# Extract data from outpatient view
# Exclude blank UPIs and filter for records where a patient attended
# Filter for relevant clinic dates
# Select columns and collect data
# Filter for UPIs in cancer_data

outpatient_data <- tbl(dvprod_con, 
                  dbplyr::in_schema("cancer_intelligence", "outpatient")) |> 
  filter(!is.na(pat_upi)) |> 
  filter(clinic_attendance_status_code == "1") |> 
  filter(clinic_date >= ref_start & clinic_date <= ref_end) |>
  select(pat_upi, date_referral_received, referral_source_code, 
         referral_type_code, clinic_date, clinic_type_code) |>
  collect() |> 
  filter(pat_upi %in% cancer_data$pat_upi)

gc()

saveRDS(outpatient_data, paste0(data_path, "/", todays_date, "_outpatient_data.rds"))
