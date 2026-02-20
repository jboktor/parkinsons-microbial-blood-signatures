# Joe Boktor
# Caltech - Mazmanian Lab
pdmbs_dir <- "/central/groups/MazmanianLab/joeB/PDMBS"
wgs_wkdir <- paste0(pdmbs_dir, "/workflow/WGS")
wkdir <- paste0(pdmbs_dir, "/parkinsons-microbial-blood-signatures")
source(paste0(wkdir, "/notebooks/R_scripts/_load-core-pkgs.R"))
source(paste0(wkdir, "/notebooks/R_scripts/_misc_functions.R"))
source(paste0(wkdir, "/notebooks/R_scripts/_plot-functions.R"))
source(paste0(wkdir, "/notebooks/R_scripts/_analysis_funcs.R"))
library(phyloseq)
library(microbiome)
library(ggsci)
library(ggpubr)
library(ggdist)
library(gghalves)

refDB <- "UHGG"
level <- "Species"
seq_method <- "WGS"
decon_status <- "decontaminated"
refDBs <- c("UHGG", "RefSeqPlusPF")
taxa_levels <- c("Species", "Genus")

# Read in ps
ps <- readRDS(
  glue(
    "{wkdir}/data/processed/phyloseq_objects/{decon_status}/",
    "2023-07-17_{seq_method}_{refDB}_{level}_phyloseq.rds"
  )
) %>%
  prune_samples(sample_sums(.) > 1, .)
env <- meta(ps)

alpha_figures_dir <- glue("{wkdir}/figures/alpha-diversity")
alpha_data_dir <- glue("{wkdir}/data/processed/alpha-diversity")
dir.create(alpha_figures_dir, showWarnings = FALSE)
dir.create(alpha_data_dir, showWarnings = FALSE, recursive = TRUE)

alpha_stats_df <- tibble()
ps_objs <- list()
for (refDB in refDBs) {
  for (level in taxa_levels) {
    ps <- readRDS(
      glue(
        "{wkdir}/data/processed/phyloseq_objects/{decon_status}/",
        "2023-07-17_{seq_method}_{refDB}_{level}_phyloseq.rds"
      )
    )
    ## Calculate Alpha Diversity Metrics and add cols to df
    env <- meta(ps)
    alpha_stats <- microbiome::alpha(abundances(ps), index = "observed") %>%
      rownames_to_column(var = "participant_id")
    stats_df <- env %>%
      left_join(alpha_stats) %>%
      dplyr::select(case_control_other_latest, colnames(alpha_stats)) %>%
      pivot_longer(!c(case_control_other_latest, participant_id),
        names_to = "alpha_metric"
      ) %>%
      mutate(DB = refDB, rank = level)
    alpha_stats_df %<>% bind_rows(stats_df)
  }
}

alpha_stats_df %<>%
  mutate(grouping = glue("{DB}_{rank}")) %>%
  dplyr::left_join(env) %>%
  mutate(case_control_other_latest = factor(
    case_control_other_latest,
    levels = c("Control", "Case", "Other"),
  ))

saveRDS(alpha_stats_df,
  glue("{alpha_data_dir}/{Sys.Date()}_{seq_method}_alpha_stats_df.rds")
)



alpha_stats_df <- readRDS(
  glue("{alpha_data_dir}/2023-07-18_{seq_method}_alpha_stats_df.rds")
)
grp_comparisons <- list(
  c("Control", "Case"),
  c("Case", "Other"),
  c("Control", "Other")
)
p_alpha_obs_spec <- alpha_stats_df %>%
  filter(rank == "Species") %>%
  ggplot(aes(case_control_other_latest, value))  +
    theme_bw() +
    scale_color_brewer(palette = "Greys") +
    stat_interval(
      .width = c(.25, 0.5, .75, 1),
      height = 5, show.legend = T
    ) +
    stat_halfeye(aes(fill = case_control_other_latest),
      .width = 0,
      scale = 0.5, size = 0.2, alpha = 0.4, height = 0.4,
      justification = -.25,
      point_alpha = 1, point_fill = "white",
      point_color = "white", shape = 23
    ) +
    geom_half_point(
      side = "l", range_scale = .3, alpha = 0.4,
      shape = 21, color = "black", size = 0.5
    ) +
    facet_wrap(~grouping, nrow = 1, scales = "free") +
    stat_compare_means(comparisons = grp_comparisons) +
    stat_compare_means() +
    # stat_compare_means(method = "wilcox.test", label.y.npc = 0.9) +
    labs(x = NULL, y = "Unique Taxa per sample", color = "interval") +
      guides(size = "none", fill = "none") +
      scale_fill_d3() +
      scale_y_log10() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.border = element_blank()
      )
# p_alpha_obs_spec

ggsave(
  filename = glue(
    "{alpha_figures_dir}/",
    "{Sys.Date()}_{seq_method}_Species_{decon_status}",
    "_alpha-diversity_Observed.png"
  ),
  p_alpha_obs_spec,
  width = 9, height = 8
)

p_alpha_obs_genus <- alpha_stats_df %>%
  filter(rank == "Genus") %>% 
  ggplot(aes(case_control_other_latest, value))  +
    theme_bw() +
    scale_color_brewer(palette = "Greys") +
    stat_interval(
      .width = c(.25, 0.5, .75, 1),
      height = 5, show.legend = T
    ) +
    stat_halfeye(aes(fill = case_control_other_latest),
      .width = 0,
      scale = 0.5, size = 0.2, alpha = 0.4, height = 0.4,
      justification = -.25,
      point_alpha = 1, point_fill = "white",
      point_color = "white", shape = 23
    ) +
    geom_half_point(
      side = "l", range_scale = .3, alpha = 0.4,
      shape = 21, color = "black", size = 0.5
    ) +
    facet_wrap(~grouping, nrow = 1, scales = "free") +
    stat_compare_means(comparisons = grp_comparisons) +
    stat_compare_means() +
    labs(x = NULL, y = "Unique Taxa per sample", color = "interval") +
      guides(size = "none", fill = "none") +
      scale_fill_d3() +
      # scale_y_log10() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.border = element_blank()
      )
# p_alpha_obs_genus

ggsave(
  filename = glue(
    "{alpha_figures_dir}/",
    "{Sys.Date()}_{seq_method}_Genus_{decon_status}",
    "_alpha-diversity_Observed.png"
  ),
  p_alpha_obs_genus,
  width = 9, height = 8
)

meta_df <- readRDS(glue(
  "{wkdir}/data/interim/metadata/2023-06-06_metadata_categories_dataframe.rds"
))

alpha_stats_meta_df <- alpha_stats_df %>%
  dplyr::left_join(env) %>%
  dplyr::select(
    case_control_other_latest, participant_id,
    alpha_metric, value, grouping, study,
    any_of(meta_df$metadata)
  )

# alpha_stats_meta_df %>%
#   select(where(is.numeric)) %>%
#   colnames()
numeric_outcomes <- alpha_stats_meta_df %>%
  select(any_of(meta_df$metadata)) %>%
  select(where(is.numeric)) %>%
  colnames()
categorical_outcomes <- alpha_stats_meta_df %>%
  select(any_of(meta_df$metadata)) %>%
  select(!where(is.numeric)) %>%
  colnames()


# nested_stats <- 
alpha_stats_df %>%
  dplyr::left_join(env) %>%
  dplyr::select(
    case_control_other_latest, participant_id,
    alpha_metric, value, grouping, study,
    all_of(numeric_outcomes)
  ) %>%
  pivot_longer(
    !c(
      case_control_other_latest, participant_id,
      alpha_metric, value, grouping, study
    ),
    names_to = "metadata", values_to = "metadata_value"
  ) %>%
    head(10) %>%
      View
  
  group_by(grouping,) %>%
  nest()

library(broom)
library(lme4)
library(easystats)

nested_models <- nested_stats %>%
  mutate(
    models = map(
      data, ~ lme4::glmer(value ~ case_control_other_latest + (1 | study),
        family = "gaussian", #"poisson",
        na.action = na.omit,
        data = .
      )
    )
  )

nested_models_gaus <- nested_stats %>%
  mutate(
    models = map(
      data, ~ lme4::lmer(value ~ case_control_other_latest + (1 | study),
        na.action = na.omit,
        data = .
      )
    )
  )

nested_models  %>% glimpse
nested_models$models[[1]] %>% summary
nested_models$models[[1]] %>% check_model()




meta_df$metadata_class %>% unique
meta_df$metadata_subclass %>% unique

meta_df$metadata %>%
  purrr::map()

alpha_stats_df %>%
  mutate(grouping = glue("{DB}_{rank}")) %>%
  dplyr::left_join(env) %>%
  ggplot(aes(value, mds_updrs_part_i_summary_score)) +
  geom_point(aes(color = case_control_other_latest)) +
  geom_smooth(aes(color = case_control_other_latest), method = "loess") +
  facet_wrap(~grouping, ncol = 2, scales = "free") +
  theme_bw()



env %>% glimpse





# ## Calculate Alpha Diversity Metrics and add cols to df
# env <- meta(ps)
# env$diagnosis_latest <- factor(env$diagnosis_latest, levels=c("No PD Nor Other Neurological Disorder", "Idiopathic PD"))
# env$Observed <- microbiome::alpha(abundances(ps), 'observed')$observed
# env$Shannon <- microbiome::alpha(abundances(ps), 'shannon')$diversity_shannon
# env$Evenness <- evenness(abundances(ps), 'simpson')$simpson
# 
# 
# env %>% 
#   ggplot(aes(diagnosis_latest, Observed))  +
#   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
#   geom_point(aes(fill = diagnosis_latest), shape = 21, 
#              position = position_jitterdodge(), stroke = 0.1) +
#   clean_theme() +
#   scale_fill_npg()
# env %>% 
#   ggplot(aes(diagnosis_latest, Shannon))  +
#   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
#   geom_point(aes(color = diagnosis_latest), position = position_jitterdodge()) +
#   clean_theme() +
#   scale_color_npg()
# env %>% 
#   ggplot(aes(diagnosis_latest, Evenness))  +
#   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
#   geom_point(aes(color = diagnosis_latest), position = position_jitterdodge()) +
#   clean_theme() +
#   scale_color_npg()
# 
# env %>% 
#   ggplot(aes(Observed, Shannon, color = diagnosis_latest))  +
#   geom_point(position = position_jitterdodge()) +
#   geom_smooth(method = "loess", se = F) +
#   clean_theme() +
#   scale_color_npg()

#_______________________________________________________________________________
#####                    Limma - VOOM Normalization                         ##### 
#_______________________________________________________________________________

library(limma)
library(DEFormats)
library(DESeq2)
library(edgeR)

for (refDB in refDBs) {
  for (level in taxa_levels) {
    ps <- readRDS(
      glue(
        "{wkdir}/data/processed/phyloseq_objects/{decon_status}/",
        "2023-07-17_{seq_method}_{refDB}_{level}_phyloseq.rds"
      )
    ) %>%
      prune_samples(sample_sums(.) > 1, .)
    voom_matrix <-
      model.matrix(~ 0 + case_control_other_latest,
        data = meta(ps)
      )
    dge <- ps %>% 
      phyloseq::phyloseq_to_deseq2(~ case_control_other_latest) %>% 
      DEFormats::as.DGEList()
    dge.norm <- dge %>% 
      calcNormFactors(method = "TMM") %>% 
      voom(design = voom_matrix, plot = TRUE, 
      save.plot = TRUE, normalize.method="none")
    # print(dim(t(dge.norm$E)))
    # Principal Components Analysis
    pca_df <- prcomp(t(dge.norm$E))
    
    pca_plot_study <- 
      pca_df %>% 
      ggbiplot::ggbiplot(choices = c(1,2), #obs.scale = 1, var.scale = 1,
                        groups =  meta(ps)$study,
                        ellipse = TRUE, 
                        alpha = 0.4,
                        circle = TRUE, var.axes=F,
                        size = 0.01
                        ) +
      scale_color_nejm() +
      clean_theme() +
      theme(legend.direction = "horizontal", legend.position = "top")
    ggsave(
      glue(
        "{wkdir}/figures/community_composition/",
        "{Sys.Date()}_{seq_method}_{refDB}_{level}_",
        "PCA-voom_normalized_StudyColored.png"
      ),
      pca_plot_study,
      width = 12, height = 9
    )
  }
}

library(plotly)
library(ggfortify)


refDB <- "UHGG"
level <- "Species"
seq_method <- "WGS"
decon_status <- "decontaminated"
refDBs <- c("UHGG", "RefSeqPlusPF")
taxa_levels <- c("Species", "Genus")

# Read in ps
ps <- readRDS(
  glue(
    "{wkdir}/data/processed/phyloseq_objects/{decon_status}/",
    "2023-07-17_{seq_method}_{refDB}_{level}_phyloseq.rds"
  )
) %>%
  prune_samples(sample_sums(.) > 1, .) %>%
  microbiome::transform("clr")

env <- meta(ps)
counts_df <- ps %>%
  abundances() %>%
  t() %>% 
  as.data.frame()
# counts_df %>% glimpse
full_df <-
  bind_cols(
    counts_df,
    env,
  )
# full_df %>% glimpse

pca_res <- prcomp(
  counts_df,
  center = TRUE #, scale. = TRUE
)

pca_df <- pca_res$x[, 1:3] %>%
  as.data.frame() %>%
  bind_cols(env)

# pca_df %>%
#   plot_ly(
#     x = ~PC1,
#     y = ~PC2,
#     z = ~PC3,
#     text = ~ paste("ID:", participant_id),
#     mode = "markers", marker = list(size = 6)
#   ) %>% add_markers(
#     color = ~ study,
#     colors = pal_npg(palette = c("nrc"), alpha = 1)(9)
#   ) %>%
#   layout(scene = list(
#     xaxis = list(title = glue("PC1 ({round(pca_res$sdev[1], 2)})")),
#     yaxis = list(title = "Polarizability Parameter"),
#     zaxis = list(title = "Relative Mutability")
#   ))




# p <-
autoplot(
  pca_res,
  data = full_df,
  colour = "study"
)


pca_df <- prcomp(t(counts))
dge.norm$E


require(FactoMineR)
require(factoextra)
require(MASS)
require(reshape2)

# Run the PCA
pca1 <- PCA(
  counts_df,
  quali.sup = c(8:10), graph = FALSE
)
plot.PCA(pca1)




env <- meta(ps)
pca_df_meta <- cbind(env, pca_df$x[,1:3])

var_explained <- pca_df$sdev^2/sum(pca_df$sdev^2)
var_explained[1:5]

pca_df_meta %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = study), alpha = 0.7, size = 1) +
  theme_bw(base_size = 32) +
  labs(
    x = paste0("PC1: ", round(var_explained[1] * 100, 1), "%"),
    y = paste0("PC2: ", round(var_explained[2] * 100, 1), "%")
  ) +
      guides(
      color =
        guide_legend(override.aes = list(size = 4))
    ) +
  scale_color_nejm() +
  # coord_fixed() +
  theme(legend.position = "top")

ggsave(
  filename = glue("{wkdir}/figures/community_composition/",
  "PCA-voom_normalized_StudyColored.png"),
       width = 12, height = 11, dpi = 600)







# pca_plot_case <- 
#   pca_df %>% 
#   ggbiplot::ggbiplot(choices = c(1,2), obs.scale = 1, var.scale = 1,
#                      groups =  meta(ps)$case_control_other_latest,
#                      ellipse = TRUE,
#                      alpha = 0.4,
#                      circle = TRUE,var.axes=FALSE) +
#   scale_color_aaas() +
#   clean_theme() +
#   theme(legend.direction = "horizontal", legend.position = "top")
# pca_plot_case
# ggsave(pca_plot_case, filename = "figures/community_composition/PCA-voom_normalized_CaseColored.png",
#        width = 6, height = 4)


# pca_plot_sex <- 
#   pca_df %>% 
#   ggbiplot::ggbiplot(choices = c(1,2), obs.scale = 1, var.scale = 1,
#                      groups =  meta(ps)$sex,
#                      ellipse = TRUE,
#                      alpha = 0.4,
#                      circle = TRUE,var.axes=FALSE) +
#   scale_color_aaas() +
#   clean_theme() +
#   theme(legend.direction = "horizontal", legend.position = "top")
# pca_plot_sex
# ggsave(pca_plot_sex, filename = "figures/community_composition/PCA-voom_normalized_SexColored.png",
#        width = 6, height = 4)






#___________________________________________________________________________
# PCA without normalization


# Principal Components Analysis
pca_df <- ps %>% 
  microbiome::transform("clr") %>% 
  microbiome::abundances() %>% 
  t() %>% 
  prcomp()

pca_plot_study <- 
  pca_df %>% 
  ggbiplot::ggbiplot(choices = c(1,2), #obs.scale = 1, var.scale = 1,
                     groups =  meta(ps)$study,
                     ellipse = TRUE, 
                     alpha = 0.4,
                     circle = TRUE, var.axes=F) +
  # scale_color_aaas() +
  clean_theme() +
  theme(legend.direction = "horizontal", legend.position = "top")
pca_plot_study



distance.metric <- "euclidean"
obj_dist <- ps
taxa_n <- obj_dist %>% microbiome::alpha("Observed")
obj_dist %<>% microbiome::transform(transform = "clr")
iDist <- phyloseq::distance(obj_dist, method = distance.metric)
iMDS <- phyloseq::ordinate(obj_dist, "MDS", distance = iDist)
p <-
  plot_ordination(
    obj_dist,
    iMDS,
    # color = "line",
    # shape = "organism",
    axes = c(1, 2, 3)
  )
pcoa_data <- p$data %>% cbind(taxa_n)
pcoa_data %>% glimpse



# ggsave(pca_plot_study, filename = "figures/community_composition/PCA-voom_normalized_StudyColored.png",
#        width = 6, height = 4)


#_______________________________________________________________________________
#                          SNM Correction        ----
#_______________________________________________________________________________

# library(snm)

# rawdat <- obj.plot %>% abundances()
rawdat <- dge.norm$E
adj.var <- model.matrix(~ study, data=metadat)
bio.var.input <- model.matrix(~ sex, data = metadat)

snm.obj <-
  snm(
    raw.dat = rawdat,
    bio.var = bio.var.input,
    adj.var = adj.var,
    rm.adj = TRUE,
    verbose = TRUE,
    diagnose = TRUE
  )

snm.dge.norm <- t(snm.obj$norm.dat)
dim(snm.dge.norm)



#_______________________________________________________________________________
#                         PERMNOVA + PCA Vis        ----
#_______________________________________________________________________________

# base::load("data/Analyses/community_composition/PERMANOVA_UHGG_Species.RData")
permanova_pca_vis(permanova_df = species_permanova_df,
                  pca_df = prcomp(t(dge.norm$E)), 
                  fileName = "figures/community_composition/PCA_PERMANOVA_hits_UHGG_Voom-Norm_")

permanova_pca_vis(permanova_df = species_permanova_df,
                  pca_df = prcomp(t(snm.obj$norm.dat)), 
                  fileName = "figures/community_composition/PCA_PERMANOVA_hits_UHGG_Voom-SNM-Norm_")




#_______________________________________________________________________________
#              Principal Variance Component Analysis (PVCA)       ----
#_______________________________________________________________________________

# library(golubEsets)
# data(Golub_Merge)

# Create ExpressionSet Object from phyloseq Object
dat.exp <- ExpressionSet(assayData=otu_table(ps),
                         phenoData=AnnotatedDataFrame(data.frame(sample_data(ps))))

tst <- data.frame(sample_data(ps))

library(pvca)
pct_threshold <- 0.6
batch.factors <- c("sex", "study", "case_control_other_latest", "age_at_baseline",
                   "ethnicity", "race", "education_level_years")

pvcaObj <- pvcaBatchAssess (dat.exp, batch.factors, pct_threshold)

pvca_df <- do.call(rbind, Map(data.frame, metadata=pvcaObj$label, prop.var=pvcaObj$dat))
  
pvca_df %>% 
  filter(metadata != "resid") %>% 
  filter(!grepl(":", metadata)) %>% 
  ggplot(aes(x = prop.var, y = reorder(metadata, prop.var))) +
  geom_col(width = 0.4) +
  geom_text(aes(x = prop.var + .003, label = round(prop.var, digits = 3)),
            size = 3) +
  clean_theme()


bp <- barplot(pvcaObj$dat, xlab = "Effects",
                ylab = "Weighted average proportion variance",
                ylim= c(0,1.1),col = c("blue"), las=2,
                main="PVCA estimation bar chart")

axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 0.5, las=2)
values = pvcaObj$dat
new_values = round(values , 3)
text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.8)




#_______________________________________________________________________________
#                           Taxa Barplots       ----
#_______________________________________________________________________________

library(fantaxtic)

# base::load("data/Phyloseq_Objects/UHGG/Genus_counts.RData")

# dat.dataframe <- dat.genus %>% psmelt()
# dat.agr = aggregate(Abundance~study+diagnosis_latest+OTU, data=dat.dataframe, FUN=mean)
# ggplot(dat.agr, aes(x=diagnosis_latest, y=Abundance, fill=OTU)) + 
#   geom_bar(stat="identity") + 
#   facet_grid(~study, scale="free")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
prune.dat_top20 <- transform_sample_counts(dat.genus, function(OTU) OTU/sum(OTU))
prune.dat_top20 <- prune_taxa(top20, prune.dat_top20)
genus_plot <- plot_bar(prune.dat_top20, x="participant_id", fill="OTU") + 
  facet_wrap(~diagnosis_latest, scales="free_x")
print(genus_plot)


# Plot all Samples
barplt1 <- 
  fantaxtic_bar(
    ps,
    color_by = "Genus",
    label_by = "Genus",
    other_label = "Other",
    facet_by = "diagnosis",
    grid_by = "cohort",
    facet_cols = 3,
    order_alg = "hclust",
    # base_color = "#5b9bd5", 
    palette = barcols
    # color_levels = barcol_ID
  ) +
  labs(y = "Relative Abundance") +
  theme(axis.text.x = element_blank())
barplt1
# ggsave(barplt1, filename = "data/Community_Composition/Stacked_Barplots/Top30_Genera_Cohort_Facet.png",
#        width = 12, height = 6)


# obj_dist <- 
#   microbiome::transform(ps, "compositional")
#   # microbiome::transform(ps, "clr")
# 
# iDist <- phyloseq::distance(obj_dist, method="bray")
# # iDist <- phyloseq::distance(obj_dist, method="euclidean")
# 
# dist_label <- "Aitchison"
# cat("Processing", dist_label, "Distance:", z[cnt], "\n")
# iMDS  <- phyloseq::ordinate(obj_dist, "MDS", distance=iDist)
# p <- plot_ordination(obj_dist, iMDS, color="description", axes = c(1, 2))
# 
# df12 = p$data
# p <- ggplot(df12, aes(Axis.1, Axis.2, fill = diagnosis_latest, color=diagnosis_latest))
# p <- p + geom_point(shape=21, size=3, alpha=0.7)
# ord <- p + 
#   theme_bw() + 
#   labs(fill="Donor Group") +
#   xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1])*100, digits = 2), "%)")) +
#   ylab(paste0("PCoA 2 (", round((iMDS$values$Relative_eig[2])*100, digits = 2), "%)")) +
#   labs(fill="Donor Group") +
#   scale_color_aaas() +
#   scale_fill_aaas() +
#   # scale_fill_manual(values = cols.pdpchc.dark) +
#   # scale_color_manual(values = cols.pdpchc.rim) +
#   theme(plot.title = element_text(hjust = 0.5), 
#         panel.grid = element_blank())
# ord
# ggsave(ord, 
#        filename = "figures/community_composition/beta-diversity_BrayCurtis.png", 
#        width = 7, height = 4)

