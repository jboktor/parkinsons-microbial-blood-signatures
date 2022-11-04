# Joe Boktor
# Caltech - Mazmanian Lab

source("notebooks/R_scripts/_load_packages.R")
source("notebooks/R_scripts/_plot-functions.R")
idvars <- c("participant_id", "visit_name", "visit_month")

# Define list filtering function
remove_idvars <- function(l) {
  l %>% subset(l %nin% c("participant_id", "visit_name", "visit_month"))
}

#_______________________________________________________________________________
#                              longitudinal data                           ####
#_______________________________________________________________________________

## Medical History  ----

medical_history <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/PD_Medical_History.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
medical_history_vars <-
  unique(colnames(medical_history)) %>% remove_idvars()

family_history <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/Family_History_PD.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
family_history_vars <-
  unique(colnames(family_history)) %>% remove_idvars()


genetic_status <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/Clinically_Reported_Genetic_Status.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-c(GUID, visit_name, visit_month)) %>%
  mutate_all(na_if, "") %>%
  mutate_all(na_if, "Unknown/Not collected as enrollment criterion")

genetic_status_vars <-
  unique(colnames(genetic_status)) %>% remove_idvars()

medical_history_metadata <-
  list(
    "medical_history" = medical_history_vars,
    "family_history" = family_history_vars,
    "genetic_status" = genetic_status_vars
  )

## Clinical Assessments ----
mds1 <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/MDS_UPDRS_Part_I.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
mds2 <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/MDS_UPDRS_Part_II.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
mds3 <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/MDS_UPDRS_Part_III.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
mds4 <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/MDS_UPDRS_Part_IV.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
mdsupdrs_vars <-
  unique(c(
    colnames(mds1),
    colnames(mds2),
    colnames(mds3),
    colnames(mds4)
  )) %>% remove_idvars()

updrs <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/UPDRS.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
updrs_vars <- unique(colnames(updrs)) %>% remove_idvars()

upsit <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/UPSIT.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
upsit_vars <- unique(colnames(upsit)) %>% remove_idvars()

mmse <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/MMSE.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
mmse_vars <- unique(colnames(mmse)) %>% remove_idvars()

moca <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/MOCA.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-c(
    GUID,
    code_education_12years_complete,
    education_12years_complete
  ))
moca_vars <- unique(colnames(moca)) %>% remove_idvars()

modswb_adl <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/Modified_Schwab___England_ADL.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
modswb_adl_vars <- unique(colnames(modswb_adl)) %>% remove_idvars()

rem_mayo <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/REM_Sleep_Behavior_Disorder_Questionnaire_Mayo.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
rem_mayo_vars <- unique(colnames(rem_mayo)) %>% remove_idvars()

rem_sk <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/REM_Sleep_Behavior_Disorder_Questionnaire_Stiasny_Kolster.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
rem_sk_vars <- unique(colnames(rem_sk)) %>% remove_idvars()

sleep_scale <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/Epworth_Sleepiness_Scale.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
sleep_scale_vars <-
  unique(colnames(sleep_scale)) %>% remove_idvars()


clincal_assessments_metadata <- list(
  "MDS-UPDRS" = mdsupdrs_vars,
  "UPDRS" = updrs_vars,
  # "MOCA" = moca_vars,
  "Modified_Schwab_ADL" = modswb_adl_vars,
  "REM_Mayo" = rem_mayo_vars,
  "REM_Stiasny_Kolster" = rem_sk_vars,
  "Epworth_Sleepiness_Scale" = sleep_scale_vars
)

## Brain Imaging ----

mri <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/MRI.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
mri_vars <- unique(colnames(mri)) %>% remove_idvars()

datscan <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/DaTSCAN_SBR.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
datscan_vars <- unique(colnames(datscan)) %>% remove_idvars()

dti <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/DTI.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
dti_vars <- unique(colnames(dti)) %>% remove_idvars()

brain_imaging_metadata <- list("MRI" = mri_vars,
                               # "DTI" = dti_vars,
                               "DATSCAN" = datscan_vars)

## Dietary/Behavioral questionnaires ----

smoking_alcohol <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/Smoking_and_alcohol_history.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
smoking_alcohol_vars <-
  unique(colnames(smoking_alcohol)) %>% remove_idvars()

caffeine <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/Caffeine_history.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
caffeine_vars <- unique(colnames(caffeine)) %>% remove_idvars()

pdq <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/PDQ_39.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-GUID)
pdq_vars <- unique(colnames(pdq)) %>% remove_idvars()

dietary_behavioral_metadata <- list(
  "smoking_and_alcohol" = smoking_alcohol_vars,
  "caffeine_intake" = caffeine_vars,
  "PDQ" = pdq_vars
)

longitudinal_data <-
  ## Clinical Assessments
  mds1 %>%
  full_join(mds2, by = idvars) %>%
  full_join(mds3, by = idvars) %>%
  full_join(mds4, by = idvars) %>%
  full_join(updrs, by = idvars) %>%
  full_join(upsit, by = idvars) %>%
  # full_join(moca, by = idvars) %>%     # needs additional wrangling - some duplicate values at M0
  full_join(modswb_adl, by = idvars) %>%
  full_join(rem_mayo, by = idvars) %>%
  full_join(rem_sk, by = idvars) %>%
  full_join(sleep_scale, by = idvars) %>%
  ## Brain Imaging
  full_join(mri, by = idvars) %>%
  full_join(datscan, by = idvars) %>%
  # full_join(dti, by = idvars) %>%         # needs additional wrangling to pivot into wide format
  ## Dietary/Behavioral questionnaires
  full_join(smoking_alcohol, by = idvars) %>%
  full_join(caffeine, by = idvars) %>%
  full_join(pdq, by = idvars) %>%
  ## Medical History
  full_join(medical_history, by = idvars) %>%
  full_join(family_history, by = idvars) %>%
  full_join(genetic_status, by = "participant_id")

records_wo_visit <- longitudinal_data %>%
  filter(visit_name == "LOG") %>%
  dplyr::select(-c("visit_name", "visit_month")) %>%
  mutate_all(na_if, "") %>%
  remove_empty_cols() %>%
  distinct()

static_start_data <- longitudinal_data %>%
  filter(visit_name == "M0") %>%
  mutate_all(na_if, "") %>% remove_empty_cols() %>%
  dplyr::select(-c("visit_name", "visit_month")) %>%
  distinct()
static_start_data2repair <- static_start_data %>%
  filter(participant_id %in% records_wo_visit$participant_id)

records_wo_visit2repair <- records_wo_visit %>%
  filter(participant_id %in% static_start_data2repair$participant_id)
records_wo_visit_extra <- records_wo_visit %>%
  filter(participant_id %nin% records_wo_visit2repair$participant_id)
repaired_static <-
  coalesce(static_start_data2repair, records_wo_visit2repair)

static_start_data_v2 <- static_start_data %>%
  filter(participant_id %nin% static_start_data2repair$participant_id) %>%
  full_join(repaired_static) %>%
  full_join(records_wo_visit_extra)

#_______________________________________________________________________________
#                          WGS Metadata (LOG & M0)                        ####
#_______________________________________________________________________________

sampdat <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/amp_pd_participants.csv",
           stringsAsFactors = F,
           header = T) %>%
  dplyr::select(participant_id, study)
case_control <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/amp_pd_case_control.csv",
           stringsAsFactors = F,
           header = TRUE)

# Data wrangling Demographics data

dem <-
  read.csv(file = "data/input/metadata/2021_v2-5release_0510/clinical/Demographics.csv",
           stringsAsFactors = F,
           header = TRUE) %>%
  dplyr::select(-c(visit_name, visit_month, GUID)) %>%
  mutate_all(na_if, "Unknown") %>%
  distinct()
# process rows with varying responses over time
dem_responded <- dem %>%
  filter(!is.na(education_level_years))
dem_dups <- dem_responded %>%
  get_dupes(participant_id) %>%
  filter(education_level_years == "Greater than 16 years") %>%
  dplyr::select(-dupe_count)
dem_true_na <- dem %>%
  filter(is.na(education_level_years)) %>%
  filter(participant_id %nin% dem_responded$participant_id)
dem_false_na <- dem %>%
  filter(is.na(education_level_years)) %>%
  filter(participant_id %in% dem_responded$participant_id)

dem_repaired <- dem_responded %>%
  filter(participant_id %nin% dem_dups$participant_id) %>%
  full_join(dem_dups) %>%
  filter(participant_id %in% dem_false_na$participant_id) %>%
  coalesce(dem_false_na)

dem_v2 <- dem_responded %>%
  filter(participant_id %nin% dem_dups$participant_id) %>%
  full_join(dem_dups) %>%
  full_join(dem_repaired) %>%
  full_join(dem_true_na)

demographic_vars <- list("demographic_vars" =
                           unique(c(
                             colnames(dem_v2),
                             colnames(sampdat),
                             colnames(case_control)
                           )))

saveRDS(dem_v2, file = glue("data/interim/metadata/{Sys.Date()}_demographics.rds"))

# Combining data ----

sample_info <- sampdat %>%
  full_join(case_control, by = "participant_id") %>%
  left_join(dem_v2, by = "participant_id") %>%
  full_join(static_start_data_v2, by = "participant_id") %>%
  mutate_all(na_if, "Unknown")

metadata_categories <- list(
  "medical_history" = medical_history_metadata,
  "clincal_assessments" = clincal_assessments_metadata,
  "brain_imaging" = brain_imaging_metadata,
  "dietary_behavioral" = dietary_behavioral_metadata,
  "demographics" = demographic_vars
)

saveRDS(sample_info, file = glue("data/interim/metadata/{Sys.Date()}_static_metdata.rds"))


# Longitudinal Metadata ----

sample_info_long <- sampdat %>%
  full_join(case_control, by = "participant_id") %>%
  left_join(dem_v2, by = "participant_id") %>%
  full_join(longitudinal_data, by = "participant_id") %>%
  mutate_all(na_if, "Unknown")

saveRDS(sample_info_long, file = glue("data/interim/metadata/{Sys.Date()}_longitudinal_metadata.rds"))

#_______________________________________________________________________________
# Metadata Categories ----
#_______________________________________________________________________________

dt_list <-
  map(metadata_categories, as.data.table) %>% suppressWarnings()
metadata_categories_df <-
  rbindlist(dt_list, fill = TRUE, idcol = T) %>%
  pivot_longer(!.id) %>%
  dplyr::rename(
    metadata_class = .id,
    metadata_subclass = name,
    metadata = value
  ) %>%
  drop_na(metadata) %>%
  distinct()


saveRDS(metadata_categories, file = glue("data/interim/metadata/{Sys.Date()}_metadata_categories.rds"))
saveRDS(metadata_categories_df, file = glue("data/interim/metadata/{Sys.Date()}_metadata_categories_dataframe.rds"))