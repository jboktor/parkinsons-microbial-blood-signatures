# Joe Boktor
# Caltech - Mazmanian Lab

rm(list = ls())
source("src/_load_packages.R")
source("src/_plot-functions.R")
base::load("data/Phyloseq_Objects/WoL/Species_counts.RData")
base::load("data/Metadata/metadata_categories.RData") #metadata_categories
base::load("data/Metadata/metadata_categories_df.RData") #metadata_categories_df

# Define metadata matrix
metadat <- meta(dat.species)
alpha_diversity <- dat.species %>% abundances() %>% 
  microbiome::alpha('shannon')
metadat$diversity_shannon <- alpha_diversity$diversity_shannon


# ## Calculate Alpha Diversity Metrics and add cols to df
# env <- meta(dat.species)
# alpha_stats <- microbiome::alpha(abundances(dat.species)) %>% 
#   rownames_to_column(var = "participant_id")
# 
# alpha_df <- env %>% left_join(alpha_stats) %>% 
#   dplyr::select(case_control_other_latest, colnames(alpha_stats)) %>% 
#   pivot_longer(!c(case_control_other_latest, participant_id), 
#                names_to = "alpha_metric")
# 
# alpha_matrix <- alpha_df %>% 
#   filter(case_control_other_latest != "Other") %>% 
#   ggplot(aes(case_control_other_latest, value))  +
#   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), lwd = 0.3) +
#   geom_point(aes(fill = case_control_other_latest), shape = 21, 
#              position = position_jitterdodge(jitter.width = 0.3), stroke = 0.1, alpha =0.3) +
#   facet_wrap(~alpha_metric, scales = "free") +
#   stat_compare_means(method = "wilcox.test", label.y.npc = 0.9) +
#   clean_theme() +
#   scale_fill_aaas()
# alpha_matrix
# ggsave(alpha_matrix, 
#        filename = "figures/community_composition/alpha_matrix.png", 
#        width = 13, height = 13)


# ## Calculate Alpha Diversity Metrics and add cols to df
# env <- meta(dat.species)
# env$diagnosis_latest <- factor(env$diagnosis_latest, levels=c("No PD Nor Other Neurological Disorder", "Idiopathic PD"))
# env$Observed <- microbiome::alpha(abundances(dat.species), 'observed')$observed
# env$Shannon <- microbiome::alpha(abundances(dat.species), 'shannon')$diversity_shannon
# env$Evenness <- evenness(abundances(dat.species), 'simpson')$simpson
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

voom_matrix <-
  model.matrix( ~ 0 + case_control_other_latest + sex,
                data = meta(dat.species))

dge <- 
  dat.species %>% 
  phyloseq_to_deseq2(~ case_control_other_latest) %>% 
  as.DGEList()
dge_highQC <- filterByExpr(dge, voom_matrix)
dge.filtered <- dge[dge_highQC,,keep.lib.sizes=FALSE] 
dge.norm <- dge.filtered %>% 
  calcNormFactors(method = "TMM") %>% 
  voom(design = voom_matrix, plot = TRUE, save.plot = TRUE, normalize.method="none")
print(dim(t(dge.norm$E)))



# Principal Components Analysis
pca_df <- prcomp(t(dge.norm$E))

pca_plot_study <- 
  pca_df %>% 
  ggbiplot::ggbiplot(choices = c(1,2), #obs.scale = 1, var.scale = 1,
                     groups =  meta(dat.species)$study,
                     ellipse = TRUE, 
                     alpha = 0.4,
                     circle = TRUE, var.axes=F) +
  scale_color_aaas() +
  clean_theme() +
  theme(legend.direction = "horizontal", legend.position = "top")
pca_plot_study
ggsave(pca_plot_study, filename = "figures/community_composition/PCA-voom_normalized_StudyColored.png",
       width = 6, height = 4)


pca_plot_case <- 
  pca_df %>% 
  ggbiplot::ggbiplot(choices = c(1,2), obs.scale = 1, var.scale = 1,
                     groups =  meta(dat.species)$case_control_other_latest,
                     ellipse = TRUE,
                     alpha = 0.4,
                     circle = TRUE,var.axes=FALSE) +
  scale_color_aaas() +
  clean_theme() +
  theme(legend.direction = "horizontal", legend.position = "top")
pca_plot_case
ggsave(pca_plot_case, filename = "figures/community_composition/PCA-voom_normalized_CaseColored.png",
       width = 6, height = 4)


pca_plot_sex <- 
  pca_df %>% 
  ggbiplot::ggbiplot(choices = c(1,2), obs.scale = 1, var.scale = 1,
                     groups =  meta(dat.species)$sex,
                     ellipse = TRUE,
                     alpha = 0.4,
                     circle = TRUE,var.axes=FALSE) +
  scale_color_aaas() +
  clean_theme() +
  theme(legend.direction = "horizontal", legend.position = "top")
pca_plot_sex
ggsave(pca_plot_sex, filename = "figures/community_composition/PCA-voom_normalized_SexColored.png",
       width = 6, height = 4)






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
dat.exp <- ExpressionSet(assayData=otu_table(dat.species),
                         phenoData=AnnotatedDataFrame(data.frame(sample_data(dat.species))))

tst <- data.frame(sample_data(dat.species))

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

# library(fantaxtic)

base::load("data/Phyloseq_Objects/UHGG/Genus_counts.RData")
prune.dat <- prune_taxa(taxa_sums(dat.genus) > 10, dat.genus) #remove less than 2%

# dat.dataframe <- dat.genus %>% psmelt()
# dat.agr = aggregate(Abundance~study+diagnosis_latest+OTU, data=dat.dataframe, FUN=mean)
# ggplot(dat.agr, aes(x=diagnosis_latest, y=Abundance, fill=OTU)) + 
#   geom_bar(stat="identity") + 
#   facet_grid(~study, scale="free")

top20 <- names(sort(taxa_sums(dat.genus), decreasing=TRUE))[1:20]
prune.dat_top20 <- transform_sample_counts(dat.genus, function(OTU) OTU/sum(OTU))
prune.dat_top20 <- prune_taxa(top20, prune.dat_top20)
genus_plot <- plot_bar(prune.dat_top20, x="participant_id", fill="OTU") + 
  facet_wrap(~diagnosis_latest, scales="free_x")
print(genus_plot)


# Plot all Samples
barplt1 <- 
  fantaxtic_bar(
    dat.genus,
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
#   microbiome::transform(dat.species, "compositional")
#   # microbiome::transform(dat.species, "clr")
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

