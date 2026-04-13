# Load the required library
library(readxl)
library(tidyverse)

metadata <- readRDS(
  glue("{wkdir}/data/interim/metadata/compiled_metadata_v4_release_2025-07-21.rds")) %>% 
  dplyr::filter(case_control_other_latest != "Other") %>% 
  dplyr::filter(study == "PPMI") %>% 
  glimpse()
pilot_pids <- metadata$participant_id %>% unique()

# Read the Excel file
pbmc_inventory <- read_excel("/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/input/temp/PBMC_Inventory_ForResearcher_9-15-2025.xlsx") %>% 
    janitor::clean_names() %>% 
    mutate(participant_id = glue::glue("PP-{id}")) %>% 
    glimpse()


# (pbmc_inventory$participant_id %in% pilot_pids) %>% table()
# FALSE  TRUE 
#   274   530 
# > 

pbmc_inventory$cohort_subgroup %>% table()


pbmc_inventory_filtered <- pbmc_inventory %>% 
    filter(participant_id %in% pilot_pids) %>% 
    glimpse()

pbmc_inventory_filtered$cohort_subgroup %>% table()

# View(pbmc_inventory_filtered)

# HOw many aliquotes in total:cohort_subgroup
sum(pbmc_inventory_filtered$grand_total) #899 aliquotes

pbmc_inventory_filtered$grand_total


# View the data
View(pbmc_inventory)
head(pbmc_inventory)
str(pbmc_inventory)


# Looking for PD and helathy control samples with longitudinal samples availabile 





pbmc_inventory_filtered %>% glimpse()
# pbmc_inventory_filtered$cohort_subgroup %>% table()
# .
#    HC:Healthy Control                PD:GBA              PD:LRRK2 
#                    54                    38                    82 
#        PD:LRRK2 + GBA               PD:PRKN               PD:SNCA 
#                     3                     3                    15 
#        PD:Sporadic PD         Prodromal:GBA   Prodromal:GBA + RBD 
#                   164                    60                     1 
#    Prodromal:Hyposmia       Prodromal:LRRK2 Prodromal:LRRK2 + GBA 
#                     8                    79                     8 
#         Prodromal:RBD        Prodromal:SNCA            SWEDD:PRKN 
#                     8                     4                     1 
#           SWEDD:SWEDD 
#                     2 
# > 


# 54 healthy control samples, (all of them)
# b 


pbmc_inventory_filtered_noprodromal <- pbmc_inventory_filtered %>% 
  filter(!grepl("Prodromal|SWEDD", cohort_subgroup)) %>%
  glimpse()

pbmc_inventory_filtered_noprodromal$cohort_subgroup %>% table()


# 38 GBA samples, (all of them)
# 82 LRRK2 samples, (all of them)
# 164 Sporadic samples, (all of them)
# 1 GBA + RBD samples, (all of them)

# pbmc_inventory_filtered_noprodromal %>% 


pbmc_inventory_filtered_noprodromal %>% 
  filter(grand_total > 1) %>%
  pull(cohort_subgroup) %>% table()

# include two longitudinal samples from each of the following cohorts:
# longituindal samples from 8 healthy and 
# 10 GBA, 10 LRRK2, 10 Sporadic


pbmc_inventory_filtered_noprodromal %>% 
  filter(grand_total > 1) %>% 
  group_by(cohort_subgroup) %>% 
  slice_min(order_by = date_collected, n = 2) %>% 
  glimpse()


# select the earliest timepoints available 

pbmc_inventory_filtered_noprodromal %>% glimpse()
pbmc_inventory_filtered_noprodromal$cohort_subgroup %>% table()
set.seed(123)
pbmc_inventory_filtered_noprodromal_longitudinal <- 
  pbmc_inventory_filtered_noprodromal %>% 
  filter(grand_total > 1)

# randomly select 10 samples from each cohort with longitudinal samples for 2 samples
longitudinal_id_list <- list()
for (cg in unique(pbmc_inventory_filtered_noprodromal_longitudinal$cohort_subgroup)) {
  id_list <- pbmc_inventory_filtered_noprodromal_longitudinal %>% 
    filter(cohort_subgroup == cg) %>% 
    pull(id) 
  if (length(id_list) >= 10) {
    set.seed(123)
    longitudinal_id_list[[cg]] <- id_list %>% sample(10) %>% unique()
  } else {
    longitudinal_id_list[[cg]] <- id_list
  }
}


#' For each of the id's in this list- we will select two samples resulting in a total
#' 96 longitudinal samples (leaving 360-96 = 264 unique donors not is this set to select from)
long_sids <- longitudinal_id_list %>% unlist()
# long_sids %>% length() * 2

remaining_healthy_sids <- pbmc_inventory_filtered_noprodromal %>% 
  filter(id %nin% long_sids) %>% 
  filter(cohort_subgroup == "HC:Healthy Control") %>% 
  pull(id) %>% unique()

# 264-length(remaining_healthy_sids)
remaining_pd_sids <- pbmc_inventory_filtered_noprodromal %>% 
  filter(id %nin% long_sids) %>% 
  filter(cohort_subgroup != "HC:Healthy Control") %>% 
  pull(id) %>% unique()  %>% 
  sample(264-length(remaining_healthy_sids))


# for single timepoint samples, select the earliest timepoint available
single_timepoint_samples <- pbmc_inventory_filtered_noprodromal %>% 
  dplyr::rename("0" = "bl") %>% 
  filter(id %in% remaining_healthy_sids| id %in% remaining_pd_sids) %>% 
  pivot_longer(!c(id, participant_id, cohort_subgroup, grand_total), names_to = "timepoint_name") %>% 
  mutate(timepoint = gsub("v0", "", timepoint_name) %>% gsub("v", "", .) %>% as.numeric()) %>% 
  drop_na(value) %>% 
  filter(value > 0) %>% 
  group_by(id) %>% 
  slice_min(order_by = timepoint, n = 1) %>% 
  glimpse()

longit_timepoint_samples <- pbmc_inventory_filtered_noprodromal_longitudinal %>% 
  dplyr::rename("0" = "bl") %>% 
  filter(id %in% unname(long_sids)) %>% 
  pivot_longer(!c(id, participant_id, cohort_subgroup, grand_total), names_to = "timepoint_name") %>% 
  mutate(timepoint = gsub("v0", "", timepoint_name) %>% gsub("v", "", .) %>% as.numeric()) %>% 
  drop_na(value) %>% 
  filter(value > 0) %>% 
  group_by(id) %>% 
  sample_n(size = 2) %>% 
  glimpse()

full_sample_set <- bind_rows(single_timepoint_samples, longit_timepoint_samples) %>%
  mutate(timepoint_name = if_else(timepoint_name == "0", "Baseline", timepoint_name)) %>% 
  glimpse()

sample_info_table_for_grant <- full_sample_set %>% 
  group_by(cohort_subgroup, timepoint_name) %>% 
  tally()

readr::write_csv(sample_info_table_for_grant, 
  "/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/input/temp/sample_info_table_for_grant.csv"
  )


# longit_timepoint_samples %>% 
#   select(participant_id, cohort_subgroup) %>% 
#   distinct() %>% 
#   pull(cohort_subgroup) %>% table()

single_timepoint_samples %>% 
  select(participant_id, cohort_subgroup) %>% 
  distinct() %>% 
  pull(cohort_subgroup) %>% table()
