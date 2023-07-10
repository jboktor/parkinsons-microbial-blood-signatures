# Joe Boktor
# Caltech - Mazmanian Lab

source("notebooks/R_scripts/_load_packages.R")
source("notebooks/R_scripts/_plot-functions.R")
idvars <- c("participant_id", "visit_name", "visit_month", "GUID")
mdir <- "data/input/metadata/2021_v2-5release_0510"
cdir <- glue("{mdir}/clinical")

# Utility functions ----

# Define list filtering function
remove_idvars <- function(l) {
  l %>% subset(l %nin% idvars)
}
read_format_csv <- function(path) {
  read.csv(file = path,
           stringsAsFactors = F,
           header = TRUE) %>%
    as_tibble()
}
collect_cols <- function(df) {
  unique(colnames(df)) %>% subset(. %nin% idvars)
}

#_______________________________________________________________________________
#                              longitudinal data                           ####
#_______________________________________________________________________________
## Medical History  ----
medical_history <-
  read_format_csv(glue("{cdir}/PD_Medical_History.csv")) %>% select(-GUID)
family_history <-
  read_format_csv(glue("{cdir}/Family_History_PD.csv"))  %>% select(-GUID)
genetic_status <-
  read_format_csv(glue("{cdir}/Clinically_Reported_Genetic_Status.csv")) %>% 
  select(-GUID)
medical_history_metadata <-
  list(
    "medical_history" = collect_cols(medical_history),
    "family_history" = collect_cols(family_history),
    "genetic_status" = collect_cols(genetic_status)
  )

## Clinical Assessments ----
mds1 <- read_format_csv(glue("{cdir}/MDS_UPDRS_Part_I.csv"))
mds2 <- read_format_csv(glue("{cdir}/MDS_UPDRS_Part_II.csv"))
mds3 <- read_format_csv(glue("{cdir}/MDS_UPDRS_Part_III.csv"))
mds4 <- read_format_csv(glue("{cdir}/MDS_UPDRS_Part_IV.csv"))
mdsupdrs_vars <-
  unique(c(
    collect_cols(mds1),
    collect_cols(mds2),
    collect_cols(mds3),
    collect_cols(mds4)
  ))

updrs <- read_format_csv(glue("{cdir}/UPDRS.csv"))
upsit <- read_format_csv(glue("{cdir}/UPSIT.csv"))
mmse <- read_format_csv(glue("{cdir}/MMSE.csv"))
rem_sk <-
  read_format_csv(glue(
    "{cdir}/REM_Sleep_Behavior_Disorder_Questionnaire_Stiasny_Kolster.csv"
  ))
sleep_scale <-
  read_format_csv(glue("{cdir}/Epworth_Sleepiness_Scale.csv"))
# some duplicate rows to clean in these dfs
rem_mayo <-
  read_format_csv(glue(
    "{cdir}/REM_Sleep_Behavior_Disorder_Questionnaire_Mayo.csv"
  )) %>%
  dplyr::group_by(participant_id, visit_name) %>%
  slice_max(order_by = msq01_act_out_dreams,
            n = 1,
            with_ties = FALSE) %>%
  ungroup()
moca <- read_format_csv(glue("{cdir}/MOCA.csv")) %>%
  dplyr::group_by(participant_id, visit_name) %>%
  slice_max(order_by = moca_abstraction_subscore,
            n = 1,
            with_ties = FALSE) %>%
  ungroup()
modswb_adl <-
  read_format_csv(glue("{cdir}/Modified_Schwab___England_ADL.csv")) %>%
  dplyr::group_by(participant_id, visit_name) %>%
  slice_max(order_by = mod_schwab_england_pct_adl_score,
            n = 1,
            with_ties = FALSE) %>%
  ungroup()

clincal_assessments_metadata <- list(
  "MDS-UPDRS" = mdsupdrs_vars,
  "UPDRS" = collect_cols(updrs),
  "UPSIT" = collect_cols(upsit),
  "MOCA" = collect_cols(moca),
  "Modified_Schwab_ADL" = collect_cols(modswb_adl),
  "REM_Mayo" = collect_cols(rem_mayo),
  "REM_Stiasny_Kolster" = collect_cols(rem_sk),
  "Epworth_Sleepiness_Scale" = collect_cols(sleep_scale)
)

## Brain Imaging ----
mri <- read_format_csv(glue("{cdir}/MRI.csv")) %>% select(-GUID)
datscan <-
  read_format_csv(glue("{cdir}/DaTSCAN_SBR.csv")) %>% select(-GUID)
# not longitudinal
dti <- read_format_csv(glue("{cdir}/DTI.csv")) %>%
  select(-c(GUID, visit_month)) %>%
  filter(!grepl("#", dti_measure)) %>%
  pivot_wider(names_from = "dti_measure", values_from = contains(c("roi", "ref")))

brain_imaging_metadata <- list(
  "MRI" = collect_cols(mri),
  "DTI" = collect_cols(datscan),
  "DATSCAN" = collect_cols(dti)
)

## Dietary/Behavioral questionnaires ----
smoking_alcohol <-
  read_format_csv(glue("{cdir}/Smoking_and_alcohol_history.csv")) %>% select(-GUID)
caffeine <-
  read_format_csv(glue("{cdir}/Caffeine_history.csv")) %>% select(-GUID)
pdq <- read_format_csv(glue("{cdir}/PDQ_39.csv")) %>%
  dplyr::group_by(participant_id, visit_name) %>%
  slice_max(order_by = pdq39_stigma_score,
            n = 1,
            with_ties = FALSE) %>%
  ungroup()

dietary_behavioral_metadata <- list(
  "smoking_and_alcohol" = collect_cols(smoking_alcohol),
  "caffeine_intake" = collect_cols(caffeine),
  "PDQ" = collect_cols(pdq)
)

testdupfunc <- function(df) {
  df %>%
    get_dupes(participant_id, visit_name, visit_month)
}

longitudinal_data <-
  ## Clinical Assessments
  mds1 %>%
  full_join(mds2) %>%
  full_join(mds3) %>%
  full_join(mds4) %>%
  full_join(updrs) %>%
  full_join(upsit) %>%
  full_join(moca) %>%
  full_join(modswb_adl) %>%
  full_join(rem_mayo) %>%
  full_join(rem_sk) %>%
  full_join(sleep_scale) %>%
  ## Brain Imaging
  full_join(mri) %>%
  full_join(datscan) %>%
  full_join(dti) %>%
  ## Dietary/Behavioral questionnaires
  full_join(smoking_alcohol) %>%
  full_join(caffeine) %>%
  full_join(pdq) %>%
  ## Medical History
  full_join(medical_history) %>%
  full_join(family_history) %>%
  full_join(genetic_status) %>%
  ## NOTE: removing secondary visit entries for specific visit
  filter(!grepl("#", visit_name))

# double check that there are no duplicates here
longitudinal_data %>% get_dupes(participant_id, visit_month)

records_wo_visit <- longitudinal_data %>%
  filter(visit_name == "LOG") %>%
  dplyr::select(-c("visit_name", "visit_month")) %>%
  remove_empty_cols() %>%
  distinct()
static_start_data <- longitudinal_data %>%
  filter(visit_name == "M0") %>%
  remove_empty_cols() %>%
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
  read_format_csv(glue("{mdir}/amp_pd_participants.csv")) %>%
  dplyr::select(participant_id, guid, study)
case_control <-
  read_format_csv(glue("{mdir}/amp_pd_case_control.csv"))
# Data wrangling Demographics data
ed_ranking <- list(
  "Greater than 16 years" = 5,
  "12-16 years" = 4,
  "Less than 12 years" = 3,
  "0 years" = 2,
  "Unknown" = 1
)
dem <- read_format_csv(glue("{cdir}/Demographics.csv")) %>%
  dplyr::select(-c(visit_name, visit_month, GUID)) %>%
  distinct() %>%
  mutate(
    education_level_years = case_when(
      education_level_years == "" ~ "Unknown",
      TRUE ~ education_level_years,
    )
  ) %>%
  mutate(ed_rank = education_level_years %>%
           purrr::map_dbl(~ unlist(ed_ranking[[.]]))) %>%
  group_by(participant_id) %>%
  slice_max(order_by = ed_rank, n = 1) %>%
  select(-ed_rank)

demographic_vars <- list("demographic_vars" =
                           unique(c(
                             colnames(dem),
                             colnames(sampdat),
                             colnames(case_control)
                           )))

core_meta <- sampdat %>%
  full_join(case_control) %>%
  full_join(dem)
saveRDS(core_meta,
        file = glue("data/interim/metadata/{Sys.Date()}_core-metadata.rds"))

# Combining data ----
sample_info <- core_meta %>%
  full_join(static_start_data_v2)

metadata_categories <- list(
  "medical_history" = medical_history_metadata,
  "clincal_assessments" = clincal_assessments_metadata,
  "brain_imaging" = brain_imaging_metadata,
  "dietary_behavioral" = dietary_behavioral_metadata,
  "demographics" = demographic_vars
)

# double checking there are no redundant participant IDs
sample_info$participant_id %>% length ==
  sample_info$participant_id %>% unique() %>% length

saveRDS(sample_info,
        file = glue("data/interim/metadata/{Sys.Date()}_static_metdata.rds"))

# Longitudinal Metadata ----
sample_info_long <- core_meta %>%
  full_join(longitudinal_data) # %>% mutate_all(na_if, "Unknown")

# # double checking there are no redundant participant ID x timepoints
sample_info_long %>% get_dupes(participant_id, visit_month) #%>% View

saveRDS(
  sample_info_long,
  file = glue(
    "data/interim/metadata/{Sys.Date()}_longitudinal_metadata.rds"
  )
)

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

saveRDS(
  metadata_categories,
  file = glue(
    "data/interim/metadata/{Sys.Date()}_metadata_categories.rds"
  )
)
saveRDS(
  metadata_categories_df,
  file = glue(
    "data/interim/metadata/{Sys.Date()}_metadata_categories_dataframe.rds"
  )
)
