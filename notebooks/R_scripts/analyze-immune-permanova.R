# Joe Boktor
# Caltech - Mazmanian Lab

source("src/_load_packages.R")
source("src/_misc_functions.R")
source("src/_analysis_funcs.R")
source("src/_plot-functions.R")
library(ggpointdensity)
library(umap)
library(stringr)

# load Data
res_quantiseq <- readRDS(file = "data/Immune_Cell_Proportions/quanTIseq.rds")
res_epic <- readRDS(file = "data/Immune_Cell_Proportions/EPIC.rds")
res_xcell <- readRDS(file = "data/Immune_Cell_Proportions/xCell.rds")
res_mcp_counter <- readRDS(file = "data/Immune_Cell_Proportions/MCP-counter.rds")
# load Metadata
longitudinal_metadata <- readRDS(file = "data/Metadata/longitudinal_metadata.rds")
rna_sample_inv <-
  read.csv(file = "input_files/2021_v2-5release_0510/rna_sample_inventory.csv",
           stringsAsFactors = F, header = TRUE)
rna_sample_inv %>% head()

# format metadata
longitudinal_metadata %<>%
  mutate(visit_name = gsub("#2", "", visit_name)) %>%
  left_join(rna_sample_inv, by = c("participant_id", "visit_month")) %>% 
  mutate_all(na_if,"")
longitudinal_metadata %>% dim()

longitudinal_metadata %<>%
  mutate(
    age_at_baseline_factor =
      cut(
        longitudinal_metadata$age_at_baseline,
        breaks = quantile(longitudinal_metadata$age_at_baseline, na.rm = T),
        labels = c(1, 2, 3, 4),
        include.lowest = F
      )
  )

#_______________________________________________________________________________
# Select  
res_input <- res_xcell
#_______________________________________________________________________________


immune_long <-
  res_input %>%
  pivot_longer(!cell_type, names_to = "sample_id", values_to = "score") %>% 
  right_join(longitudinal_metadata) %>% 
  filter(case_control_other_latest %in% c("Case", "Control"))

immune_avg_df <-
  immune_long %>%
  drop_na(score) %>%
  group_by(participant_id, cell_type) %>% 
  summarize(score_avg = mean(score)) %>% 
  pivot_wider(names_from = "cell_type", values_from = "score_avg") %>% 
  column_to_rownames(var = "participant_id")
# saveRDS(immune_avg_df, "data/Immune_Cell_Proportions/EPIC_averages.rds") # CHANGE DEPENDING ON RES INPUT

basic_meta <- 
  longitudinal_metadata %>% 
  filter(participant_id %in% rownames(immune_avg_df)) %>% 
  dplyr::select(participant_id, study, diagnosis_latest, case_control_other_latest,
                age_at_baseline, sex, ethnicity, race) %>% distinct()

# check all row names align between metadata and immune cell abund. matrix
all(basic_meta$participant_id == rownames(immune_avg_df))



immune.pca <- prcomp(immune_avg_df, center = TRUE,scale. = TRUE)

immune_pca_df <- immune.pca$x %>% as.data.frame() %>% 
  rownames_to_column("participant_id") %>% 
  left_join(basic_meta)

immune_pca_df %>% 
  ggplot(aes(PC1, PC2, color = age_at_baseline)) +
  geom_point() +
  scale_color_viridis() +
  theme_minimal() +
  theme(legend.position = "top")

immune_pca_df %>% 
  ggplot(aes(PC1, PC2, fill = case_control_other_latest)) +
  geom_point(shape=21, stroke = 0.2) +
  scale_fill_manual(values = cco_pal) +
  theme_minimal() +
  theme(legend.position = "top")

immune_pca_df %>% 
  ggplot(aes(PC1, PC2, fill = sex)) +
  geom_point(shape=21, stroke = 0.2) +
  scale_fill_calc() +
  theme_minimal() +
  theme(legend.position = "top")

# UMAP Analysis
immune.umap <- umap(immune_avg_df)

immune_umap_df <- immune.umap$layout %>% 
  as.data.frame() %>% 
  rownames_to_column("participant_id") %>% 
  left_join(basic_meta)

immune_umap_df %>% 
  ggplot(aes(V1, V2, color = age_at_baseline)) +
  geom_point() +
  scale_color_viridis() +
  theme_minimal() +
  theme(legend.position = "top")

immune_umap_df %>% 
  ggplot(aes(V1, V2, color = case_control_other_latest)) +
  geom_point() +
  scale_color_manual(values = cco_pal) +
  theme_minimal() +
  theme(legend.position = "top")


# ggbiplot(immune.pca, ellipse=TRUE, obs.scale = 1, var.scale = 1,
#          groups=basic_meta$case_control_other_latest) +
#   theme_minimal() +
#   theme(legend.position = "top")
# screeplot(immune.pca, type = "line", main = "Scree plot")
# plot(immune.pca)

base::load("data/Phyloseq_Objects/Phyloseq_all_outliers_removed.RData")
RefDB <- "UHGG"
level <- "Genus"
ps_obj <- refDB_phyloseq_orm[[RefDB]][[level]]

# prep immune cell prop data for permanova
immune_avg_df %<>% janitor::clean_names()
# Quantize average abundance 
immune_avg_df_quant <- list()
for (celltype in colnames(immune_avg_df)) {
  print(celltype)
  quantile_breaks <- unique(quantile(immune_avg_df[[celltype]], na.rm = T))
  
  immune_avg_df_quant[[celltype]] <- 
    cut(immune_avg_df[[celltype]],
    breaks = quantile_breaks,
    labels = c(1:(length(quantile_breaks)-1)),
    include.lowest = T)

}

immune_avg_df_quant %<>% as.data.frame()
rownames(immune_avg_df_quant) <- rownames(immune_avg_df)
# create a metadata phyloseq obj for immune cell averages
immune_data_phy <- phyloseq::sample_data(immune_avg_df_quant)
# merge immune cell averages with phyloseq metadata
ps_obj_new <- merge_phyloseq(ps_obj, immune_data_phy)


#_______________________________________________________________________________
#####        PERMANOVA Analysis --Average Immune Cell Proportions        ##### 
#_______________________________________________________________________________

permanova_analysis_immune <- 
  phyloseq_permanova(ps_object = ps_obj_new, 
                     nperm = 1000,
                     colnames(immune_avg_df_quant))


permanova_analysis_immune_df <- permanova_analysis_immune %>% 
  dplyr::mutate(R2_perc = R2 * 100) %>%
  mutate(FDR = p.adjust(p_value, method = 'BH')) %>%
  ungroup()

permanova_analysis_immune_df %>% 
  ggplot(aes(R2_perc, FDR)) +
  geom_point()

permanova_analysis_immune_df %>% 
  ggplot(aes(x= R2_perc, y = reorder(metadata, R2))) +
  geom_segment(aes(x=0, xend=R2_perc, y=reorder(metadata, R2), yend=reorder(metadata, R2)),
               color="grey") +
  geom_point(aes(fill = FDR), shape = 21, size = 4) +
  scale_fill_viridis_c(option = "F", direction = -1) + 
  labs(x = expression(R^"2"~"(%)"), y = "", fill = "P-adj") +
  clean_theme()


# save(permanova_df, file = paste0("data/Analyses/community_composition/PERMANOVA_IMMUNE", 
#                                  RefDB, "_", level, ".RData"))
# openxlsx::write.xlsx(permanova_df,
#                      file = paste0('data/Analyses/community_composition/PERMANOVA_IMMUNE', 
#                                    RefDB, '_', level, '.xlsx'), overwrite = T)



