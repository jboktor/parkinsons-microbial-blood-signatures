# Joe Boktor
# Caltech - Mazmanian Lab

home_dir <- "/central/groups/MazmanianLab/joeB"
pdmbs_dir <- paste0(home_dir, "/PDMBS")
ref_dir <- paste0(home_dir, "/Downloads/RefDBs")
wkdir <- paste0(pdmbs_dir, "/parkinsons-microbial-blood-signatures")
src_dir <- paste0(wkdir, "/notebooks")
source(paste0(src_dir, "/R_scripts/_load-core-pkgs.R"))
source(paste0(src_dir, "/R_scripts/_misc_functions.R"))
library(data.table)

idvars <- c("participant_id", "visit_name", "visit_month")
v4_metapath <- glue("{wkdir}/data/input/metadata/2023_v4release_1027")
v4_clincial <- glue("{v4_metapath}/clinical")
v4_genetics <- glue(
  "{wkdir}/data/input/metadata/",
  "genetic_status_tier2/2023_v4release_1027")
arcashla_interim <- glue("{wkdir}/data/interim/airr/arcasHLA")
qc_dir <- glue("{wkdir}/data/input/metadata/QC")
hla_dir <- glue("{wkdir}/data/interim/airr/arcasHLA/analysis_cleaned")
pca_dir <- glue("{hla_dir}/pca_data")

#_______________________________________________________________________________
#                              helper functions                           ####
#_______________________________________________________________________________
# Define list filtering function
remove_idvars <- function(l) {
  l %>% subset(l %nin% c("participant_id", "visit_name", "visit_month"))
}

get_colnames <- function(df) {
  unique(colnames(df)) %>% remove_idvars()
}

#' This function converges metadata from LOG, M0, and SC visits into a single
#' visit annotated as M0
#' unique values are compiled while conflicting column values are sampled
#' If there is a conflict, the M0 value is preferred
#' @param df a data frame with a column for visit_name
#' @return a data frame with a single row per participant
#' @export
converge_metadata <- function(df) {
  # filter for baseline, LOG, and SC visits
  df_filt <- df %>% 
    filter(visit_name %in% c("M0", "LOG") | grepl("SC", visit_name))
  
  set.seed(42)
  df_filt %>%
      # one row per participant
      dplyr::group_by(participant_id) %>%
      # put the “M0” row first      
      dplyr::arrange(desc(visit_name == "M0"), .by_group = TRUE) %>%
      # pick first non-NA (now from M0 if it exists)
      dplyr::summarise(dplyr::across(everything(), ~ first(na.omit(.x))), 
        .groups = "drop") %>%
      # sample 1 row per participant
      dplyr::mutate(visit_name = "M0", visit_month = 0) %>%
      distinct() %>%
      dplyr::group_by(participant_id, visit_name) %>%
      dplyr::slice_sample(n = 1) %>%
      glimpse()
}

test_df_format <- function(df) {
  sample_timept_cts <- df %>%
    dplyr::group_by(participant_id, visit_name, visit_month) %>% 
    dplyr::summarize(n = n(), .groups = "drop")  %>% 
    dplyr::arrange(desc(n))

  if (any(sample_timept_cts$n > 1)) {
    message("\n⚠️⚠️  Warning: Duplicate rows found in the data frame ⚠️⚠️\n")
    return(sample_timept_cts)
  } else {
    message("\n ✅ No duplicate rows found in the data frame\n")
  }
}

#_______________________________________________________________________________
#                              essential metadata
#_______________________________________________________________________________

cohort_meta <- read.csv(glue("{v4_metapath}/amp_pd_participants.csv"),
  stringsAsFactors = FALSE, header = TRUE
) %>%
  dplyr::select(participant_id, study) %>%
  glimpse()

case_control <- read.csv(glue("{v4_metapath}/amp_pd_case_control.csv"),
  stringsAsFactors = FALSE, header = TRUE
) %>%
  glimpse()
case_control_vars <- get_colnames(case_control)

# Data wrangling Demographics data ----
dem <- read.csv(glue("{v4_clincial}/Demographics.csv"),
           stringsAsFactors = FALSE,
           header = TRUE) %>%
  dplyr::select(-c(visit_name, visit_month, GUID)) %>%
  dplyr::mutate_if(is.character, ~na_if(., "Unknown")) %>%
  glimpse()

sample_inv <- bind_rows(
    read_csv(glue("{v4_metapath}/rnaseq_WB-RWTS-VHBS_sample_inventory.csv")),
    read_csv(glue("{v4_metapath}/rnaseq_WB-RWTS_sample_inventory.csv"))) %>% 
    dplyr::select(sample_id, participant_id, visit_month) %>% 
    glimpse()

amppd_meta_core <- read_csv(
    glue("{v4_metapath}/amp_pd_participants.csv")) %>% 
    dplyr::select(participant_id, study) %>% 
    right_join(sample_inv) %>%
    glimpse()
#_______________________________________________________________________________
#                              longitudinal data
#_______________________________________________________________________________

## Medical History  ----
medical_history_raw <-
  read.csv(file = glue("{v4_clincial}/PD_Medical_History.csv"),
           stringsAsFactors = FALSE,
           header = TRUE) %>%
  mutate_if(is.character, ~na_if(., "")) %>%
  dplyr::select(-c(GUID)) %>%
  glimpse() # 10,000+ rows

harmonized_med_baseline <- converge_metadata(medical_history_raw)
medical_history <- bind_rows(
  harmonized_med_baseline,
  medical_history_raw %>%
    filter(visit_name %nin% c("M0", "LOG") & !grepl("SC", visit_name))
  ) %>%
  glimpse()
medical_history_vars <- get_colnames(medical_history)
test_df_format(medical_history)

## Family History  ----
family_history <-
  read.csv(file = glue("{v4_clincial}/Family_History_PD.csv"),
           stringsAsFactors = FALSE,
           header = TRUE) %>%
  dplyr::select(-GUID)
family_history_vars <- get_colnames(family_history)

## Genetic Status  ----
genetic_status <-
  read.csv(
  glue("{v4_genetics}/clinical/Clinically_Reported_Genetic_Status.csv"),
           stringsAsFactors = FALSE,
           header = TRUE) %>%
  dplyr::select(-c(GUID, visit_name, visit_month)) %>%
  mutate_all(na_if, "") %>%
  mutate_all(na_if, "Unknown/Not collected as enrollment criterion") %>%
  distinct() %>% 
  glimpse()
genetic_status_vars <- get_colnames(genetic_status)

## Mutation Status  ----
mutation_status <- read.csv(glue(
  "{v4_genetics}/amp_pd_participant_mutations.csv")
  )
mutation_status_vars <- get_colnames(mutation_status)

## APOE Status  ----
apoe_status <- read.csv(glue(
  "{v4_genetics}/amp_pd_participant_apoe_mutations.csv")
  )
apoe_status_vars <- get_colnames(apoe_status)

## Clinical Assessments ----
mds1 <- read.csv(glue("{v4_clincial}/MDS_UPDRS_Part_I.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  dplyr::select(-GUID)
mds2 <- read.csv(glue("{v4_clincial}/MDS_UPDRS_Part_II.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  dplyr::select(-GUID)
mds3 <- read.csv(glue("{v4_clincial}/MDS_UPDRS_Part_III.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  dplyr::select(-GUID)
mds4 <- read.csv(glue("{v4_clincial}/MDS_UPDRS_Part_IV.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  dplyr::select(-GUID)

mdsupdrs_vars <-
  unique(c(
    colnames(mds1),
    colnames(mds2),
    colnames(mds3),
    colnames(mds4)
  )) %>% remove_idvars()

updrs <- read.csv(glue("{v4_clincial}/UPDRS.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  dplyr::select(-c(GUID, updrs1_on_medication, updrs2_on_medication)) %>% 
  left_join(sample_inv) %>% 
  filter(!is.na(updrs1_ment_behav_mood_score) |
    !is.na(updrs3_motor_examination_score)) %>%
  dplyr::select(-sample_id) %>%
  glimpse()
test_df_format(updrs)
updrs_vars <- get_colnames(updrs)

upsit <- read.csv(glue("{v4_clincial}/UPSIT.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  dplyr::select(-GUID)
upsit_vars <- get_colnames(upsit)

# filtering for duplicate entries for 14 samples in the MMSE dataset
set.seed(42)
mmse <- read.csv(glue("{v4_clincial}/MMSE.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  dplyr::select(-GUID) %>% 
  mutate_if(is.character, ~na_if(., "")) %>%
  mutate(col_entries = rowSums(!is.na(.))) %>%
  group_by(participant_id, visit_name, visit_month) %>% 
  slice_max(col_entries, n = 1, with_ties = FALSE) %>%
  select(-col_entries) %>% 
  ungroup() %>% 
  glimpse()
test_df_format(mmse)
mmse_vars <- get_colnames(mmse)

set.seed(42) 
moca <- read.csv(glue("{v4_clincial}/MOCA.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  dplyr::select(-c(GUID)) %>% 
  mutate_if(is.character, ~na_if(., "")) %>%
  mutate(col_entries = rowSums(!is.na(.))) %>%
  group_by(participant_id, visit_name, visit_month) %>% 
  slice_max(col_entries, n = 1, with_ties = FALSE) %>%
  select(-col_entries) %>% 
  ungroup() %>% 
  glimpse()
test_df_format(moca)
moca_vars <- get_colnames(moca)

set.seed(42) 
modswb_adl <- read.csv(glue("{v4_clincial}/Modified_Schwab___England_ADL.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  dplyr::select(-GUID, -mod_schwab_england_on_off_med) %>% 
  filter(!is.na(mod_schwab_england_pct_adl_score)) %>%
  mutate_if(is.character, ~na_if(., "")) %>%
  mutate(col_entries = rowSums(!is.na(.))) %>%
  group_by(participant_id, visit_name, visit_month) %>% 
  slice_max(col_entries, n = 1, with_ties = FALSE) %>%
  select(-col_entries) %>% 
  ungroup() %>% 
  glimpse()
test_df_format(modswb_adl)
modswb_adl_vars <- get_colnames(modswb_adl)

#' one sample (SU-33113) has two entries at M24 - 
#' here we sample one of the two entries (either more info or random)
set.seed(42) 
rem_mayo <- read.csv(glue("{v4_clincial}/",
    "REM_Sleep_Behavior_Disorder_Questionnaire_Mayo.csv"),
           stringsAsFactors = FALSE,
           header = TRUE) %>%
  dplyr::select(-GUID) %>% 
  mutate_if(is.character, ~na_if(., "")) %>%
  mutate(col_entries = rowSums(!is.na(.))) %>%
  group_by(participant_id, visit_name, visit_month) %>% 
  slice_max(col_entries, n = 1, with_ties = FALSE) %>%
  select(-col_entries) %>% 
  ungroup() %>% 
  glimpse()
test_df_format(rem_mayo)
rem_mayo_vars <- get_colnames(rem_mayo)

rem_sk <- read.csv(glue("{v4_clincial}/",
    "REM_Sleep_Behavior_Disorder_Questionnaire_Stiasny_Kolster.csv"),
           stringsAsFactors = FALSE,
           header = TRUE) %>%
  dplyr::select(-GUID)
rem_sk_vars <- get_colnames(rem_sk)

sleep_scale <- read.csv(glue("{v4_clincial}/Epworth_Sleepiness_Scale.csv"),
           stringsAsFactors = FALSE,
           header = TRUE) %>%
  dplyr::select(-GUID)
sleep_scale_vars <- get_colnames(sleep_scale)

## Brain Imaging ----
mri <- read.csv(glue("{v4_clincial}/MRI.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  dplyr::select(-GUID)
mri_vars <- get_colnames(mri)

datscan <- read.csv(glue("{v4_clincial}/DaTSCAN_SBR.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  dplyr::select(-GUID)
datscan_vars <- get_colnames(datscan)

datscan_vis <- read.csv(glue("{v4_clincial}/DaTSCAN_visual_interpretation.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  dplyr::select(-c(GUID, scan_months_after_baseline, 
    visit_name, visit_month)) %>%
  distinct() %>%
  glimpse()
datscan_vis_vars <- get_colnames(datscan_vis)


dti_long <- read.csv(glue("{v4_clincial}/DTI.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  filter(dti_measure %in% 
    c(paste0("Eigenvalue", 1:3), "Fractional Anisotropy")) %>% 
  dplyr::select(-GUID, dti_brain_tissue)
dti_metrics <- colnames(dti_long) %>% keep(grepl("right|left", .))
dti <- dti_long %>%
  pivot_wider(
    id_cols = c(participant_id, visit_name, visit_month),
    names_from = dti_measure,
    values_from = all_of(dti_metrics),
    names_glue = "{.value}_{dti_measure}"
  ) %>%
  janitor::clean_names() %>% 
  # mutate(visit_month = 0) %>% 
  dplyr::select(-c(visit_name, visit_month)) %>% 
  glimpse()
dti_vars <- get_colnames(dti)

## Analytical Metrics ----
abeta_ptau <- read.csv(
  glue("{v4_clincial}/Biospecimen_analyses_CSF_abeta_tau_ptau.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
      mutate(new_col = 
        glue("{test_name}_{test_units}") %>% gsub("/mL", "_per_ml", .)) %>%
  dplyr::select(-GUID, -sample_type) %>%
  pivot_wider(
    id_cols = c(participant_id, visit_name, visit_month),
    names_from = new_col,
    values_from = test_value
  ) %>%
  janitor::clean_names() %>% 
  glimpse()
abeta_ptau_vars <- get_colnames(abeta_ptau)

gcase <- read.csv(
  glue("{v4_clincial}/Biospecimen_analyses_CSF_beta_glucocerebrosidase.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  mutate(new_col = "gcase_avg_pmol_min_ml") %>%
  filter(test_units != "%CV") %>% 
  dplyr::select(-GUID, -sample_type) %>%
  pivot_wider(
    id_cols = c(participant_id, visit_name, visit_month),
    names_from = new_col,
    values_from = test_value
  ) %>%
  janitor::clean_names() %>%
  dplyr::select(-visit_name) %>% 
  glimpse()
gcase_vars <- get_colnames(gcase)

#' SAA Results ----
#' Data cleaning rules - 
#' if there are duplicate entries for a given sample x timepoint
#' prefer the entry the 24 hour assay over the 24 hour assay
#' prefer a positive result over a negative result or inclusive 
#' prefer an inclusive result over a negative result
#' If all else equal - look at SAA Status
#' Type1 > Type2 > Undetermined

saa_result <- read.csv(
  glue("{wkdir}/data/input/metadata/",
    "SAA_Biospecimen_Analysis_Results_26Jun2025.csv"),
      stringsAsFactors = FALSE, header = TRUE) %>%
  filter(SAAMethod == "Amprion-24h alpha-synuclein-SAA") %>%
  dplyr::select(PATNO, CLINICAL_EVENT, SAA_Status, SAA_Type) %>%
  glimpse()

ppmi_events <- tribble(
  ~CLINICAL_EVENT, ~visit_name,
  "BL",  "M0",    # Baseline
  "V01", "M3",
  "V02", "M6",
  "V03", "M9",
  "V04", "M12",
  "V05", "M18",
  "V06", "M24",
  "V07", "M30",
  "V08", "M36",
  "V09", "M42",
  "V10", "M48",
  "V11", "M54",
  "V12", "M60",
  "V16", "M108"
)

saa_result <- saa_result %>%
  mutate(
    # Adjust these column names as needed
    result_priority = case_when(
      tolower(SAA_Status) == "Positive" ~ 1,
      tolower(SAA_Status) == "Inclusive" ~ 2,
      tolower(SAA_Status) == "Negative" ~ 3,
      TRUE ~ NA_real_
    ),
    status_priority = case_when(
      SAA_Type == "Type1" ~ 1,
      SAA_Type == "Type2" ~ 2,
      SAA_Type == "Undetermined" ~ 3,
      TRUE ~ NA_real_
    )
  )

# Now, for each sample x timepoint, keep the "best" row
saa_result_clean <- saa_result %>%
  group_by(PATNO, CLINICAL_EVENT) %>%
  arrange(result_priority, status_priority) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  select(-c(result_priority, status_priority)) %>% 
  mutate(participant_id = paste0("PP-", PATNO)) %>% 
  # map month names to visit_name
  left_join(ppmi_events, by = "CLINICAL_EVENT") %>%
  select(participant_id, visit_name, SAA_Status, SAA_Type) %>% 
  drop_na(SAA_Status, visit_name) %>%
  glimpse()
saa_vars <- get_colnames(saa_result_clean)

## Dietary/Behavioral questionnaires ----
#' here we sample one general entry per participant - disregarding the timeline -
#' the entry with greater info is selected

set.seed(42)
smoking_alcohol <-
  read.csv(file = glue("{v4_clincial}/Smoking_and_alcohol_history.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  dplyr::select(-GUID, -visit_name, -visit_month) %>% 
  mutate_if(is.character, ~na_if(., "")) %>%
  mutate(col_entries = rowSums(!is.na(.))) %>%
  group_by(participant_id) %>%
  slice_max(col_entries, n = 1, with_ties = FALSE) %>%
  select(-col_entries) %>% 
  ungroup() %>% 
  glimpse()
smoking_alcohol_vars <- get_colnames(smoking_alcohol)

caffeine <-
  read.csv(file = glue("{v4_clincial}/Caffeine_history.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  dplyr::select(-GUID) %>% 
  glimpse()
caffeine_vars <- get_colnames(caffeine)

#' duplicate entries for SU-32242, SU-33113, SU-33114, and SU-33368
#' Randomly sample one of the two entries for each of these samples
set.seed(42)
pdq <- read.csv(glue("{v4_clincial}/PDQ_39.csv"),
           stringsAsFactors = FALSE, header = TRUE) %>%
  dplyr::select(-GUID) %>% 
  mutate_if(is.character, ~na_if(., "")) %>%
  mutate(col_entries = rowSums(!is.na(.))) %>%
  group_by(participant_id, visit_name, visit_month) %>% 
  slice_max(col_entries, n = 1, with_ties = FALSE) %>%
  select(-col_entries) %>%
  ungroup() %>%
  glimpse()
test_df_format(pdq)
pdq_vars <- get_colnames(pdq)

#_______________________________________________________________________________
## Manually Curated Metadata ----
#_______________________________________________________________________________
#_______________________________________________________________________________
# Read stats and estimated T-cell reads per sample ----

main_cohort_stats <- read_tsv(glue(
    "{qc_dir}/rnaseq-WB-RWTS/multiqc_data/",
    "mqc_picard_aligned_reads_1.txt")) %>% 
    clean_names()

hbs_cohort_stats <- read_tsv(glue(
    "{qc_dir}/rnaseq-WB-RWTS-VHBS/multiqc_data/",
      "picard_alignment_summary_Aligned_Reads.txt")) %>% 
    clean_names()

picard_stats <- bind_rows(
    main_cohort_stats,
    hbs_cohort_stats) %>% 
    dplyr::rename(sample_id = sample) %>% 
    left_join(amppd_meta_core) %>%
    glimpse()

# Selecting HBS samples with the largest estimated T cell Reads
decon_meta <- readRDS(
    glue(
      "{wkdir}/data/interim/celltype_deconvolution/",
      "deconvoluted_cell_types_metadata_2025-06-23.rds")
)
decon_meta %>% glimpse()
picard_stats %>% glimpse()

read_stats <- picard_stats %>%
  dplyr::select(
    sample_id, participant_id, aligned_reads,
    unaligned_reads, study, visit_month) %>% 
  left_join(decon_meta) %>%
  mutate_at(vars(contains("abis_")), ~ . / 100) %>% 
  mutate(
    est_tcell_reads_abis = abis_t_cell * aligned_reads,
    est_tcell_reads_epic = epic_t_cell * aligned_reads
  ) %>%
  glimpse()

read_stats_vars <- get_colnames(read_stats) %>%
  discard(~ . %in% c("sample_id", "study"))

#_______________________________________________________________________________
# Inflammation gene program scores ----
deconv_outdir <- glue("{wkdir}/data/interim/celltype_deconvolution")
inf_scores <- readRDS(
    glue("{deconv_outdir}/",
        "inflammation_gene_program_scores_2025-07-06.rds"
    )
)
inf_scores_vars <- get_colnames(inf_scores) %>% discard(~ . %in% c("sample_id"))

#_______________________________________________________________________________
# HLA genotypes ----
#' Annotations with the highest number of reads are used for 
#' each sample x HLA gene combination. Some small inconsistences 
#' are found in the genotype table below -- reads

mhc1_alleles <- c("A", "B", "C")
mhc2_alleles <- c("DPA1", "DPB1", "DQA1", "DQB1", "DRB1")
major_hla_alleles <- c(mhc1_alleles, mhc2_alleles)
paste0(rep(major_hla_alleles, each = 2), 1:2)

# hla_genotypes_cleaned <- readRDS(
#   glue("{arcashla_interim}/genotypes_cleaned_2025-07-06.rds")
#   ) %>%
#   dplyr::select(participant_id, 
#     all_of(paste0(rep(major_hla_alleles, each = 2), 1:2))) %>% 
#   glimpse()

# hla_genotypes_cleaned_annot <- hla_genotypes_cleaned %>%
#   mutate(
#     DRB1_grouped = case_when(
#       # Homozygous cases
#       grepl("DRB1\\*15:01", DRB11) & grepl("DRB1\\*15:01", DRB12) ~ 
#         "DRB1*15:01/DRB1*15:01",
#       grepl("DRB1\\*04:01", DRB11) & grepl("DRB1\\*04:01", DRB12) ~ 
#         "DRB1*04:01/DRB1*04:01",
#       grepl("DRB1\\*04:04", DRB11) & grepl("DRB1\\*04:04", DRB12) ~ 
#         "DRB1*04:04/DRB1*04:04",

#       # Heterozygous combinations of the three alleles
#       (grepl("DRB1\\*15:01", DRB11) & grepl("DRB1\\*04:01", DRB12)) |
#         (grepl("DRB1\\*04:01", DRB11) & grepl("DRB1\\*15:01", DRB12)) ~ 
#           "DRB1*15:01/DRB1*04:01",
#       (grepl("DRB1\\*15:01", DRB11) & grepl("DRB1\\*04:04", DRB12)) |
#         (grepl("DRB1\\*04:04", DRB11) & grepl("DRB1\\*15:01", DRB12)) ~ 
#           "DRB1*15:01/DRB1*04:04",
#       (grepl("DRB1\\*04:01", DRB11) & grepl("DRB1\\*04:04", DRB12)) |
#         (grepl("DRB1\\*04:04", DRB11) & grepl("DRB1\\*04:01", DRB12)) ~ 
#           "DRB1*04:01/DRB1*04:04",

#       # One allele of interest with another allele
#       grepl("DRB1\\*15:01", DRB11) | grepl("DRB1\\*15:01", DRB12) ~ 
#         "DRB1*15:01/DRB1X",
#       grepl("DRB1\\*04:01", DRB11) | grepl("DRB1\\*04:01", DRB12) ~ 
#         "DRB1*04:01/DRB1X",
#       grepl("DRB1\\*04:04", DRB11) | grepl("DRB1\\*04:04", DRB12) ~ 
#         "DRB1*04:04/DRB1X",
#       # None of the alleles of interest
#       TRUE ~ "DRB1X/DRB1X"
#     ),
#     DRB1_risk = case_when(
#       grepl("DRB1\\*15:01", DRB1_grouped) & grepl("DRB1\\*04", DRB1_grouped) ~ 
#         "mixed",
#       grepl("DRB1\\*15:01", DRB1_grouped) ~ "risk",
#       grepl("DRB1\\*04", DRB1_grouped) ~ "protective",
#       TRUE ~ "none"
#     )
#   ) %>%
#   glimpse()
# hla_genotypes_cleaned_vars <- get_colnames(hla_genotypes_cleaned_annot)

# Loading in HLA risk score
hla_dir <- glue("{wkdir}/data/interim/airr/arcasHLA/analysis_cleaned")
hla_score <- readRDS(glue("{hla_dir}/hla_risk_score_2025-07-21.rds"))
hla_score_vars <- get_colnames(hla_score)

# Loading in MHC PCA data
pca_mhc_all_df <- readRDS(glue("{pca_dir}/df_pca_MHC_all_2025-07-06.rds")) %>% 
    dplyr::select(1:10) %>% 
    rownames_to_column(var = "participant_id") %>% 
    glimpse()
pca_mhc_1_df <- readRDS(glue("{pca_dir}/df_pca_MHCI_2025-07-06.rds")) %>% 
    dplyr::select(1:10) %>% 
    rownames_to_column(var = "participant_id") %>% 
    glimpse()
pca_mhc_2_df <- readRDS(glue("{pca_dir}/df_pca_MHCII_2025-07-06.rds")) %>% 
    dplyr::select(1:10) %>% 
    rownames_to_column(var = "participant_id") %>% 
    glimpse()
pca_mhc_all_vars <- get_colnames(pca_mhc_all_df)
pca_mhc_1_vars <- get_colnames(pca_mhc_1_df)
pca_mhc_2_vars <- get_colnames(pca_mhc_2_df)

# Loading biallelic HLA genotype count matrix
hla_alleles_matrix <- readRDS(
  glue("{hla_dir}/MHC_all_cleaned_unfiltered_2025-07-21.rds")) %>%
  rownames_to_column(var = "participant_id") %>%
  glimpse()
hla_allele_matrix_vars <- get_colnames(hla_alleles_matrix)


#_______________________________________________________________________________

metadata <- sample_inv %>%
  full_join(cohort_meta) %>% 
  full_join(case_control) %>%
  full_join(dem) %>%
  ## Clinical Assessments
  full_join(mds1) %>%
  full_join(mds2) %>%
  full_join(mds3) %>%
  full_join(mds4) %>%
  left_join(updrs) %>%
  # prodromal assessments
  full_join(upsit) %>%
  full_join(rem_mayo) %>%
  full_join(rem_sk) %>%
  full_join(sleep_scale) %>%
  # cognitive assessments
  full_join(mmse) %>%
  full_join(moca) %>%
  full_join(modswb_adl) %>%
  ## Dietary/Behavioral questionnaires
  full_join(smoking_alcohol) %>%
  full_join(caffeine) %>%
  full_join(pdq) %>%
  ## Medical History
  full_join(medical_history) %>%
  full_join(family_history) %>%
  full_join(genetic_status) %>%
  full_join(mutation_status) %>%
  full_join(apoe_status) %>%
  ## Brain Imaging
  full_join(mri) %>%
  full_join(datscan) %>%
  full_join(datscan_vis) %>%
  full_join(dti) %>%
  # other analytical metrics
  full_join(abeta_ptau) %>%
  full_join(gcase) %>%
  full_join(saa_result_clean) %>%
  full_join(read_stats) %>%
  # Bulk RNA-Seq Inflammatory signatures
  full_join(inf_scores) %>%
  # HLA genotypes
  left_join(hla_alleles_matrix) %>%
  full_join(hla_score) %>%
  # HLA PCA
  full_join(pca_mhc_all_df) %>%
  full_join(pca_mhc_1_df) %>%
  full_join(pca_mhc_2_df) %>%
  # Final filter for samples that have RNA-Seq data
  mutate_if(is.character, ~na_if(., "")) %>%
  right_join(sample_inv) %>%
  mutate(years_since_diagnosis =  age_at_baseline - age_at_diagnosis) %>%
  glimpse()

test_df_format(metadata)

# # viewing the number of non-na entries per variable
# nrow(metadata) - colSums(is.na(metadata))

#_______________________________________________________________________________
# Metadata Categories ----
#_______________________________________________________________________________

demographic_meta <- list(
  "Demographic" =
  unique(c(
      colnames(dem),
      colnames(cohort_meta),
      "sample_id", "visit_month", "visit_name"
    )))

clinical_history_meta <- list(
  "Family History" = family_history_vars,
  "Medical History" = 
    c(medical_history_vars, case_control_vars, "years_since_diagnosis")
)

lifestyle_factors_meta <- list(
  "Smoking and Alcohol" = smoking_alcohol_vars,
  "Caffeine Intake" = caffeine_vars
)

motor_severity_meta <- list(
  "MDS-UPDRS" = mdsupdrs_vars,
  "UPDRS" = updrs_vars
)

daily_function_qol_meta <- list(
  "PDQ" = pdq_vars,
  "Modified_Schwab_ADL" = modswb_adl_vars
)

cognition_meta <- list(
  "MMSE" = mmse_vars,
  "MOCA" = moca_vars
)

prodromal_meta <- list(
  "UPSIT" = upsit_vars,
  "REM Mayo" = rem_mayo_vars,
  "REM Stiasny Kolster" = rem_sk_vars,
  "Epworth Sleepiness Scale" = sleep_scale_vars
)

biomarkers_neuro_meta <- list(
  "MRI" = mri_vars,
  "DTI" = dti_vars,
  "DATSCAN" = c(datscan_vars, datscan_vis_vars)
  )

biomarkers_csf_meta <- list(
  "A-Beta" = abeta_ptau_vars,
  "GCase" = gcase_vars,
  "SAA" = saa_vars
)

biomarkers_genetic_meta <- list(
  "Genetic Status" = genetic_status_vars,
  "Mutation Status" = mutation_status_vars,
  "ApoE Status" = apoe_status_vars,
  "HLA I&II PCA" = pca_mhc_all_vars,
  "HLA I PCA" = pca_mhc_1_vars,
  "HLA II PCA" = pca_mhc_2_vars,
  "HLA Alleles" = hla_allele_matrix_vars,
  "HLA Risk Score" = hla_score_vars
)

biomarkers_omics_meta <- list(
  "RNA-Seq Reads" = read_stats_vars,
  "RNA-Seq Inf. Signatures" = inf_scores_vars
)

# Aggregating a metadata into a list of lists
metadata_categories <- list(
  "Demographics" = demographic_meta,
  "Clinical History" = clinical_history_meta,
  "Lifestyle Factors" = lifestyle_factors_meta,
  "Motor Severity" = motor_severity_meta,
  "Daily Function/QoL" = daily_function_qol_meta,
  "Cognition" = cognition_meta,
  "Prodromal Assessments" = prodromal_meta,
  "Biomarkers Neuroimaging" = biomarkers_neuro_meta,
  "Biomarkers CSF" = biomarkers_csf_meta,
  "Biomarkers Genetics" = biomarkers_genetic_meta,
  "Biomarkers Omics" = biomarkers_omics_meta
)

# Joining metadata categories into a single dataframe
dt_list <- purrr::map(metadata_categories, as.data.table) %>% suppressWarnings()
metadata_categories_df <-
  rbindlist(dt_list, fill = TRUE, idcol = TRUE) %>%
  pivot_longer(!.id) %>%
  dplyr::rename(
    metadata_class = .id,
    metadata_subclass = name,
    metadata = value
  ) %>%
  drop_na(metadata) %>%
  distinct()

metadata_categories_df$metadata_class %>% table()
metadata_categories_df$metadata_subclass %>% table()

# Adding a color scheme to the metadata categories
len_categories <- length(unique(metadata_categories_df$metadata_class))
class_colors <- pal_simpsons("springfield")(len_categories)
names(class_colors) <- unique(metadata_categories_df$metadata_class)
class_colors


metadata_categories_df %<>%
  mutate(class_color = class_colors[metadata_class]) %>% 
  glimpse()

# Should be  - and is empty
colnames(metadata) %>% discard(. %in% metadata_categories_df$metadata)

# Saving metadata 
saveRDS(metadata,
  glue("{wkdir}/data/interim/metadata/",
    "compiled_metadata_v4_release_{Sys.Date()}.rds")
  )
# Saving metadata categories
saveRDS(metadata_categories,
  glue("{wkdir}/data/interim/metadata/",
    "metadata_categories_{Sys.Date()}.rds")
    )
saveRDS(metadata_categories_df,
  glue("{wkdir}/data/interim/metadata/",
    "metadata_categories_dataframe_{Sys.Date()}.rds")
    )



# GRAVE YARD TOSS ME when safe

#     mutate(raw_total_reads_binned = cut(raw_total_reads,
#         breaks = quantile(raw_total_reads, probs = seq(0, 1, 0.1), na.rm = TRUE),
#         include.lowest = TRUE,
#         labels = paste0(1:10)
#     ) %>% as.numeric()
#     ) %>%
#     mutate(
#         tcell_reads_x_cell = raw_total_reads * x_cell_t_cell,
#         tcell_reads_epic = raw_total_reads * epic_t_cell,
#         tcell_reads_abis = raw_total_reads * abis_t_cell
#     ) %>%
#     glimpse()


# records_wo_visit <- longitudinal_data %>%
#   filter(visit_name == "LOG") %>%
#   dplyr::select(-c("visit_name", "visit_month")) %>%
#   dplyr::mutate_if(is.character, ~na_if(., "")) %>%
#   janitor::remove_empty("cols") %>%
#   distinct() %>% glimpse()

# static_start_data <- longitudinal_data %>%
#   filter(visit_name == "M0") %>%
#   dplyr::mutate_if(is.character, ~na_if(., "")) %>%
#   janitor::remove_empty("cols") %>%
#   dplyr::select(-c("visit_name", "visit_month")) %>%
#   distinct() %>% glimpse()

# static_start_data2repair <- static_start_data %>%
#   filter(participant_id %in% records_wo_visit$participant_id) #%>% glimpse()

# records_wo_visit2repair <- records_wo_visit %>%
#   filter(participant_id %in% static_start_data2repair$participant_id)
# records_wo_visit_extra <- records_wo_visit %>%
#   filter(participant_id %nin% records_wo_visit2repair$participant_id)
# repaired_static <-
#   coalesce(static_start_data2repair, records_wo_visit2repair)

# static_start_data_v2 <- static_start_data %>%
#   filter(participant_id %nin% static_start_data2repair$participant_id) %>%
#   full_join(repaired_static) %>%
#   full_join(records_wo_visit_extra) %>% 
#   glimpse()

# dir.create("data/interim/metadata/", showWarnings=FALSE, recursive = TRUE)
# saveRDS(dem_v2, file = glue("data/interim/metadata/{Sys.Date()}_demographics.rds"))
# Combining data ----

# sample_info <- cohort_meta %>%
#   full_join(case_control, by = "participant_id") %>%
#   left_join(dem, by = "participant_id") %>%
#   full_join(static_start_data_v2, by = "participant_id") %>%
#   mutate_all(na_if, "Unknown")

# saveRDS(sample_info, file = glue("data/interim/metadata/{Sys.Date()}_static_metdata.rds"))



