# Joe Boktor
# Caltech - Mazmanian Lab
# Dec 2021

source("src/_load_packages.R")
base::load("data/Phyloseq_Objects/Phyloseq_all_outliers_removed.RData") #phyloseq_objs
dat.obj <- refDB_phyloseq_orm[["UHGG"]]$Species

metadat <- meta(dat.obj) %>% 
  select(study, case_control_other_latest, sex, age_at_baseline)

metadat %>%
  tbl_summary() %>%
  bold_labels() %>% 
  as_gt() %>%
  gt::gtsave(filename = "figures/misc/summary-stats-table.png")



meta(dat.obj) %>% 
  select(diagnosis_latest, case_control_other_latest) %>%
  filter(case_control_other_latest == "Case") %>% 
  tbl_summary(by = case_control_other_latest, 
              statistic =  all_categorical() ~ "{n} / {N} ({p}%)") %>% 
  as_gt() %>%
  gt::gtsave(filename = "figures/misc/summary-stats-table_group-status_Case.png")


meta(dat.obj) %>% 
  select(diagnosis_latest, case_control_other_latest) %>%
  filter(case_control_other_latest == "Control") %>% 
  tbl_summary(by = case_control_other_latest, 
              statistic =  all_categorical() ~ "{n} / {N} ({p}%)") %>% 
  as_gt() %>%
  gt::gtsave(filename = "figures/misc/summary-stats-table_group-status_Control.png")


meta(dat.obj) %>% 
  select(diagnosis_latest, case_control_other_latest) %>%
  filter(case_control_other_latest == "Other") %>% 
  tbl_summary(by = case_control_other_latest, 
              statistic =  all_categorical() ~ "{n} / {N} ({p}%)") %>% 
  as_gt() %>%
  gt::gtsave(filename = "figures/misc/summary-stats-table_group-status_Other.png")

