# Joe Boktor
# Caltech - Mazmanian Lab

source("src/_load_packages.R")
source("src/_plot-functions.R")
source("src/plot-PCA.R")
# base::load("data/Phyloseq_Objects/UHGG/Species_counts.RData")
base::load("data/Phyloseq_Objects/Phyloseq_all_outliers_removed.RData")

refDB <- "UHGG"
level <- "Genus"
base::load(paste0("data/Analyses/community_composition/PERMANOVA_", refDB, "_", level, ".RData"))
# Temp fix; delete later
permanova_df <- species_permanova_df
ps_object <- refDB_phyloseq_orm[[refDB]][[level]]
metadat <- meta(ps_object)
#_______________________________________________________________________________
#####                    Limma - VOOM Normalization                         ##### 
#_______________________________________________________________________________

voom_matrix <-
  model.matrix( ~ 0 + case_control_other_latest + sex,
                data = meta(ps_object))

dge <- 
  ps_object %>% 
  phyloseq_to_deseq2(~ case_control_other_latest) %>% 
  as.DGEList()
dge_highQC <- filterByExpr(dge, voom_matrix)
dge.filtered <- dge[dge_highQC,,keep.lib.sizes=FALSE] 
dge.norm <- dge.filtered %>% 
  calcNormFactors(method = "TMM") %>% 
  voom(design = voom_matrix, plot = TRUE, save.plot = TRUE, normalize.method="none")
print(dim(t(dge.norm$E)))


#_______________________________________________________________________________
#                          SNM Correction        ----
#_______________________________________________________________________________

rawdat <- dge.norm$E
adj.var <- model.matrix(~ study, data=metadat)
bio.var.input <- model.matrix(~ sex + age_at_baseline, data = metadat)

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
#                         PERMNAOVA + PCA Visualization        ----
#_______________________________________________________________________________

# base::load("data/Analyses/community_composition/PERMANOVA_UHGG_Species.RData")
permanova_pca_vis(permanova_df = permanova_df,
                  pca_df = prcomp(t(dge.norm$E), center = T, scale. = T), 
                  fileName = paste0(
                    "figures/community_composition/PCA_PERMANOVA_hits_", refDB, 
                    "_", level, "_Voom-Norm_"))

permanova_pca_vis(permanova_df = permanova_df,
                  pca_df = prcomp(t(snm.obj$norm.dat), center = T, scale. = T), 
                  fileName = paste0(
                    "figures/community_composition/PCA_PERMANOVA_hits_", refDB, 
                    "_", level, "_Voom-SNM-Norm_"))
