# Joe Boktor
# Caltech - Mazmanian Lab
# Dec 2021

pdmbs_dir <- "/resnick/groups/MazmanianLab/jboktor/PDMBS"
wgs_wkdir <- paste0(pdmbs_dir, "/workflow/WGS")
wkdir <- paste0(pdmbs_dir, "/pdairr")
source(paste0(wkdir, "/notebooks/R_scripts/_load-core-pkgs.R"))
source(paste0(wkdir, "/notebooks/R_scripts/_misc_functions.R"))
library(phyloseq)
library(microbiome)
library(gt)
library(gtsummary)
library(webshot2)

figures_dir <- glue("{wkdir}/figures/sample_summary")
ps <- readRDS(
  glue(
    "{wkdir}/data/processed/phyloseq_objects/raw/",
    "2023-07-14_WGS_RefSeqPlusPF_All_phyloseq.rds"
  )
)
rnaseq_samples <- sample_names(ps)

# for (seq_method in c("WGS", "RNASEQ")) {
sample_metadata <- readRDS(
  glue("{wkdir}/data/interim/metadata/2023-07-14_phyloseq-metadata.rds")
)

sample_metadata %<>% purrr::map(
  ~ select(.,
  sample_id, participant_id, study, diagnosis_at_baseline,
  case_control_other_at_baseline, age_at_baseline,
  case_control_other_latest,
  diagnosis_latest,
  ethnicity, sex, race
  ))
# core_meta <- sample_metadata[["WGS"]] %>% select(-sample_id)


seq_method <- "RNASEQ"
sample_metadata[[seq_method]] %>%
  select(study, case_control_other_latest, sex, age_at_baseline) %>%
  tbl_summary() %>%
  bold_labels() %>%
  as_gt() %>%
  gt::gtsave(
    filename = glue(
      "{figures_dir}/{seq_method}_summary-stats-table.png"
    )
  )
sample_metadata[[seq_method]] %>%
  select(diagnosis_latest, case_control_other_latest) %>%
  filter(case_control_other_latest == "Case") %>%
  tbl_summary(by = case_control_other_latest,
              statistic =  all_categorical() ~ "{n} / {N} ({p}%)") %>%
  as_gt() %>%
  gt::gtsave(
    filename =
    glue("{figures_dir}/{seq_method}_group-status_Case.png")
)

sample_metadata[[seq_method]] %>%
  select(diagnosis_latest, case_control_other_latest) %>%
  filter(case_control_other_latest == "Control") %>%
  tbl_summary(by = case_control_other_latest,
              statistic =  all_categorical() ~ "{n} / {N} ({p}%)") %>%
  as_gt() %>%
    gt::gtsave(
      filename =
        glue("{figures_dir}/{seq_method}_group-status_Control.png")
    )


sample_metadata[[seq_method]] %>%
  select(diagnosis_latest, case_control_other_latest) %>%
  filter(case_control_other_latest == "Other") %>%
  tbl_summary(by = case_control_other_latest,
              statistic =  all_categorical() ~ "{n} / {N} ({p}%)") %>%
  as_gt() %>%
    gt::gtsave(
      filename =
        glue("{figures_dir}/{seq_method}_group-status_Other.png")
    )


# sample_metadata[[seq_method]] %>%
#   select(diagnosis_latest, case_control_other_latest) %>%
#   # filter(case_control_other_latest == "Other") %>%
#   tbl_summary(by = case_control_other_latest,
#               statistic =  all_categorical() ~ "{n} / {N} ({p}%)") %>%
#   as_gt()
