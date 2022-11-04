# Joe Boktor
# Caltech - Mazmanian Lab

rm(list = ls())
source("src/_load_packages.R")
source("src/_plot-functions.R")
source("src/_analysis_funcs.R")
metadata_categories <- readRDS("data/Metadata/metadata_categories.rds")
metadata_categories_df <- readRDS("data/Metadata/metadata_categories_df.rds")
base::load("data/Phyloseq_Objects/Phyloseq_all_outliers_removed.RData")

RefDB <- "UHGG"
level <- "Genus"
ps_obj <- refDB_phyloseq_orm[[RefDB]][[level]]

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

dir.create(file.path("data/Analyses/community_composition/"), showWarnings = FALSE)
save(permanova_df, file = paste0("data/Analyses/community_composition/PERMANOVA_", 
                                         RefDB, "_", level, ".RData"))
openxlsx::write.xlsx(permanova_df,
                     file = paste0('data/Analyses/community_composition/PERMANOVA_', 
                                   RefDB, '_', level, '.xlsx'), overwrite = T)


