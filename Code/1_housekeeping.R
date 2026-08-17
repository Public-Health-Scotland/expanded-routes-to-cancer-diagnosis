#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1_housekeeping
# Calum Purdie
# 03/11/2021
# Housekeeping script for use in later scripts
# Written/run on Posit Workbench
# R version 4.4.2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### 1 Housekeeping ----

# Load packages

library(odbc)
library(dplyr)
library(haven)
library(janitor)
library(tidyr)
library(lubridate)
library(ggplot2)
library(here)
library(purrr)
library(glue)
library(phsverse)
library(stringr)
library(openxlsx)
library(sjlabelled)
library(PHEindicatormethods)
library(forcats)
library(ggrepel)
library(tidylog)

# Define dates
# These are the start and end dates for admissions and cancer diagnoses
# Admission start date should be at least 30 days before cancer start date to
# allow for emergencies 30 days before admission to be identified

adm_start <- "2016-12-01"
adm_end <- "2022-12-31"
cancer_start <- "2018-01-01"
cancer_end <- "2022-12-31"
ref_start <- "2017-04-01"
ref_end <- "2023-08-30"
  
# Define years from start and end dates

start <- str_sub(cancer_start, 1, 4)
end <- str_sub(cancer_end, 1, 4)

# Stop scientific notation for small numbers

options(scipen = 999)

# Define data path

data_path <- paste0("/PHI_conf/CancerGroup1/Topics/DiagnosisRoutes/", 
                    "20250505-ExpandedRoutesToDiagnosis/Data")

# Define todays date

todays_date <- today()



### 2 Functions ----

# Define function for calculating age groups for standard populations

standard_pop_age_groups <- function(current_col){
  
  # Define an age group for each current_col value
  
  case_when(current_col >= 0 & current_col <= 4 ~ 0, 
            current_col >= 5 & current_col <= 9 ~ 1, 
            current_col >= 10 & current_col <= 14 ~ 2, 
            current_col >= 15 & current_col <= 19 ~ 3, 
            current_col >= 20 & current_col <= 24 ~ 4, 
            current_col >= 25 & current_col <= 29 ~ 5, 
            current_col >= 30 & current_col <= 34 ~ 6, 
            current_col >= 35 & current_col <= 39 ~ 7, 
            current_col >= 40 & current_col <= 44 ~ 8, 
            current_col >= 45 & current_col <= 49 ~ 9, 
            current_col >= 50 & current_col <= 54 ~ 10, 
            current_col >= 55 & current_col <= 59 ~ 11, 
            current_col >= 60 & current_col <= 64 ~ 12, 
            current_col >= 65 & current_col <= 69 ~ 13, 
            current_col >= 70 & current_col <= 74 ~ 14, 
            current_col >= 75 & current_col <= 79 ~ 15, 
            current_col >= 80 & current_col <= 84 ~ 16, 
            current_col >= 85 & current_col <= 89 ~ 17, 
            current_col >= 90 ~ 18)
  
}
