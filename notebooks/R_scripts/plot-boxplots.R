#  plot - abundance
source("src/_load_packages.R")
source("src/_misc_functions.R")
source("src/_plot-functions.R")

base::load("data/Phyloseq_Objects/Phyloseq_all_outliers_removed.RData") #phyloseq_objs
dat.obj <- refDB_phyloseq_orm[["RefSeqPlusPF"]]$Genus



#_______________________________________________________________________________
#  BoxPlot function
#_______________________________________________________________________________


boxplots_maaslin2 <- function(refDB, level, feature){
  
  # # TROUBLE
  # refDB = "UHGG"
  # level = "Genus"
  # feature = "g__MGYG.HGUT.01660"
  
  dat.obj <- refDB_phyloseq_orm[[refDB]][[level]] %>% 
    subset_samples(case_control_other_latest != "Other")
  model_matrix <-
    model.matrix( ~ 0 + case_control_other_latest + sex,
                  data = meta(dat.obj))
  
  dge <-
    dat.obj %>%
    phyloseq_to_deseq2(~ case_control_other_latest) %>%
    as.DGEList()
  # dge_highQC <- filterByExpr(dge, model_matrix)
  # dge.filtered <- dge[dge_highQC,,keep.lib.sizes=FALSE]
  dge.norm <- dge %>% #dge.filtered %>%
    calcNormFactors(method = "TMM") %>%
    voom(design = model_matrix, plot = F, save.plot = F, normalize.method="none")
  print(dim(t(dge.norm$E)))
  
  df_plot <- t(dge.norm$E) %>%
    as.data.frame() %>%
    rename_all(~(
      gsub("-", ".", .) %>% gsub(" ", ".", .)
    )) %>% 
    rownames_to_column(var = "participant_id") %>%
    left_join(meta(dat.obj))
  
  plot <-
    df_plot %>%
    ggplot(aes(x = case_control_other_latest,
               y = !!sym(feature),
               color = case_control_other_latest)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.8) +
    geom_boxplot(alpha = 0.4, outlier.alpha = 0, color = "black") +  clean_theme() +
    scale_color_aaas() +
    theme(axis.title.x = element_blank(), 
          legend.position = "none")
  return(plot)
}










refDB <- "RefSeqPlusPF"
level <- "Species"
dat.obj <- refDB_phyloseq_orm[[refDB]][[level]] %>% 
  subset_samples(case_control_other_latest != "Other")
#_______________________________________________________________________________

model_matrix <-
  model.matrix( ~ 0 + case_control_other_latest + sex,
                data = meta(dat.obj))

dge <-
  dat.obj %>%
  phyloseq_to_deseq2(~ case_control_other_latest) %>%
  as.DGEList()
# dge_highQC <- filterByExpr(dge, model_matrix)
# dge.filtered <- dge[dge_highQC,,keep.lib.sizes=FALSE]
dge.norm <- dge %>% #dge.filtered %>%
  calcNormFactors(method = "TMM") %>%
  voom(design = model_matrix, plot = TRUE, save.plot = TRUE, normalize.method="none")
print(dim(t(dge.norm$E)))

df_plot <- t(dge.norm$E) %>%
  as.data.frame() %>%
  rownames_to_column(var = "participant_id") %>%
  left_join(meta(dat.obj))

df_plot_raw <-
  dat.obj %>%
  abundances() %>% t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "participant_id") %>%
  left_join(meta(dat.obj))
#_______________________________________________________________________________


# Group by Sex
F_nuc <-
  # df_plot_raw %>%
  df_plot %>%
  ggplot(aes(x = case_control_other_latest,
             y = `s__Fusobacterium nucleatum_E`,
             color = sex)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.8) +
  geom_boxplot(aes(position = sex),
               alpha = 0.4, outlier.alpha = 0, color = "black") +  clean_theme() +
  facet_grid(cols = vars(study), scales = "free_x", space = "free") +
  scale_color_aaas() +
  labs(x = NULL)
F_nuc
ggsave(F_nuc, filename = "figures/misc/s__Fusobacterium nucleatum_E - UHGG.png",
       width = 7, height = 4, dpi = 600)

C_dif <-
  # df_plot_raw %>%
  df_plot %>%
  ggplot(aes(x = case_control_other_latest,
             y = `s__Clostridioides difficile`,
             color = sex)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.8) +
  geom_boxplot(aes(position = sex),
    alpha = 0.4, outlier.alpha = 0, color = "black") +
  clean_theme() +
  facet_grid(cols = vars(study), scales = "free_x", space = "free") +
  scale_color_aaas() +
  labs(x = NULL)
C_dif
ggsave(C_dif, filename = "figures/misc/s__Clostridioides difficile - UHGG.png",
       width = 7, height = 4, dpi = 600)


HS <- 
  df_plot_raw %>%
  # df_plot %>%
  ggplot(aes(x = case_control_other_latest,
             y = `Homo sapiens`,
             color = sex)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.8) +
  geom_boxplot(aes(position = sex),
               alpha = 0.4, outlier.alpha = 0, color = "black") +
  clean_theme() +
  facet_grid(cols = vars(study), scales = "free_x", space = "free") +
  scale_color_aaas() +
  scale_y_log10() +
  labs(x = NULL)
HS
ggsave(HS, filename = "figures/misc/Homo sapiens - RefSeqPlusPF.png",
       width = 9, height =3.5, dpi = 600)

TSPG <- 
  df_plot %>%
  # df_plot_raw %>%
  ggplot(aes(x = case_control_other_latest,
             y = `Toxoplasma gondii`,
             color = sex)) +
  geom_point(position = position_jitterdodge(), alpha = 0.8) +
  geom_boxplot(aes(position = sex),
               alpha = 0.4, outlier.alpha = 0, color = "black") +  # scale_y_log10() +
  facet_grid(cols = vars(study), scales = "free_x", space = "free") +
  clean_theme() +
  scale_color_aaas()
TSPG
ggsave(TSPG, filename = "figures/misc/Toxoplasma gondii - RefSeqPlusPF.png",
       width = 7, height = 4, dpi = 600)


# # df_plot$`s__Fusobacterium nucleatum_E`
# df_plot_raw$`Saccharomyces cerevisiae`
# df_plot_raw$`s__Clostridioides difficile`
# df_plot$`Fusobacterium nucleatum`
# df_plot_raw$
