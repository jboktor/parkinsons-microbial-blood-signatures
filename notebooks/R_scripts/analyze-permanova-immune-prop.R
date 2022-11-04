# # Joe Boktor
# Caltech - Mazmanian Lab

rm(list = ls())
source("src/_load_packages.R")
source("src/_plot-functions.R")
metadata_categories <- readRDS("data/Metadata/metadata_categories.rds")
metadata_categories_df <- readRDS("data/Metadata/metadata_categories_df.rds")
base::load("data/Phyloseq_Objects/Phyloseq_all_outliers_removed.RData")

# load immune cell data
res_xcell <- readRDS(file = "data/Immune_Cell_Proportions/xCell.rds")
# load metadata
longitudinal_metadata <- readRDS(file = "data/Metadata/longitudinal_metadata.rds")
rna_sample_inv <-
  read.csv(file = "input_files/2021_v2-5release_0510/rna_sample_inventory.csv",
           stringsAsFactors = F, header = TRUE)
# format metadata
longitudinal_metadata %<>%
  mutate(visit_name = gsub("#2", "", visit_name)) %>%
  left_join(rna_sample_inv, by = c("participant_id", "visit_month")) %>% 
  mutate_all(na_if,"")




RefDB <- "UHGG"
level <- "Genus"
ps_obj <- refDB_phyloseq_orm[[RefDB]][[level]]
# Define metadata matrix
metadat <- meta(ps_obj)











#_______________________________________________________________________________
#####                       PERMANOVA Analysis                             ##### 
#_______________________________________________________________________________

permanova_analysis <- 
  phyloseq_permanova(ps_object = ps_obj, 
                     nperm = 99999,
                     unlist(metadata_categories,
                            recursive = T, use.names = F))


permanova_df <- permanova_analysis %>% 
  left_join(metadata_categories_df) %>% 
  dplyr::mutate(R2_perc = R2 * 100) %>%
  mutate(FDR = p.adjust(p_value, method = 'BH')) %>%
  ungroup()


# dir.create(file.path("data/Analyses/community_composition/"), showWarnings = FALSE)
# save(permanova_df, file = paste0("data/Analyses/community_composition/PERMANOVA_", 
#                                  RefDB, "_", level, ".RData"))
# openxlsx::write.xlsx(permanova_df,
#                      file = paste0('data/Analyses/community_composition/PERMANOVA_', 
#                                    RefDB, '_', level, '.xlsx'), overwrite = T)