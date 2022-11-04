
source("src/_load_packages.R")
library(edgeR)
library(limma)
library(DEFormats)
library(DESeq2)
library(apeglm)
library(Maaslin2)
library(VennDiagram)
base::load("data/Phyloseq_Objects/Phyloseq_all_outliers_removed.RData")

refDB <- "UHGG"
level <- "Species"
dat.obj <- refDB_phyloseq_orm[[ref_DB]][[level]] %>% 
  subset_samples(case_control_other_latest != "Other")
# filter out low prevalence features

dat.analysis <- dat.obj %>%
  core(detection = 5, prevalence = 0.1)
sample_data(dat.analysis)$case_control_other_latest <-
  factor(sample_data(dat.analysis)$case_control_other_latest,
         levels = c( "Control", "Case")) # Set control first as reference level
table(sample_data(dat.analysis)$case_control_other_latest)
dat.analysis.dseq2 <- phyloseq_to_deseq2(dat.analysis, ~ case_control_other_latest)

#_______________________________________________________________________________
# LIMMA / VOOM

dge <- as.DGEList(dat.analysis.dseq2)
#store TMM norm factor
dge <- calcNormFactors(dge, method = "TMM")
head(dge$samples$norm.factors)

# tst <- cpm(dge, log=TRUE, prior.count=3)

#construct model matrix
mm <- model.matrix(~ group, dge$samples)
head(mm)
table(mm[, 2])
#obtain Voom weights
y <- voom(dge, mm, plot = T)
#fit lm with limma
fit <- lmFit(y, mm)
fit <- eBayes(fit)
head(coef(fit))
limma_res_df <- data.frame(topTable(fit, coef = "groupCase", number = Inf))    #extract results

fdr_limma <- limma_res_df %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  rownames_to_column(var = "features")
dim(fdr_limma)

ggplot(fdr_limma, aes(x = features, y = logFC, color = features)) +
  geom_point(size = 4) +
  labs(y = "\nLog2 Fold-Change for PD vs. Controls", x = "") +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype="dotted")

#_______________________________________________________________________________
# DESeq2

dds <- DESeq(dat.analysis.dseq2, test = "Wald", fitType = "local", sfType = "poscounts")
plotDispEsts(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
plotMA(dds)
deseq_res_df <- data.frame(res) %>%
  rownames_to_column(var = "features") %>%
  dplyr::arrange(padj)

fdr_deseq <- deseq_res_df %>%
  dplyr::filter(padj < 0.05)

dim(fdr_deseq)
head(fdr_deseq)

ggplot(fdr_deseq, aes(
  x =  reorder(features, log2FoldChange),
  y = log2FoldChange,
  color = features
)) +
  geom_point() +
  labs(y = "\nLog2 Fold-Change for PD vs. Control", x = "") +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype="dotted")



#_______________________________________________________________________________
# CORNCOB

corn_da <- differentialTest(formula = ~ case_control_other_latest,
                            phi.formula = ~ case_control_other_latest,
                            formula_null = ~ 1,
                            phi.formula_null = ~ case_control_other_latest,
                            data = dat.analysis,
                            test = "Wald", boot = TRUE, B = 100,
                            fdr_cutoff = 0.05)

fdr_corncob <- corn_da$significant_taxa
dim(data.frame(fdr_corncob))
head(sort(corn_da$p_fdr))


corn_dv <- differentialTest(formula = ~ case_control_other_latest,
                            phi.formula = ~ case_control_other_latest,
                            formula_null = ~ case_control_other_latest,
                            phi.formula_null = ~ 1,
                            data = dat.analysis,
                            test = "Wald", boot = TRUE, B = 1000,
                            fdr_cutoff = 0.05)

fdr_corncob_dv <- corn_dv$significant_taxa
dim(data.frame(fdr_corncob_dv))
head(sort(corn_dv$p_fdr))


#_______________________________________________________________________________
# Maaslin2 ----

input_data_df <- data.frame(otu_table(dat.analysis))
input_metadata_df <- data.frame(sample_data(dat.analysis))
rownames(input_metadata_df) <- gsub("-", ".", input_metadata_df$participant_id)

mas_1 <- Maaslin2(
  input_data = input_data_df,
  input_metadata = input_metadata_df,
  output = "data/Analyses/differential_abundance/testrun/",
  min_abundance = 0.0,
  min_prevalence = 0.0,
  normalization = "TMM",
  transform = "NONE",
  analysis_method = "NEGBIN",
  max_significance = 0.05,
  fixed_effects = c("case_control_other_latest", "sex"),
  correction = "BH",
  standardize = FALSE,
  cores = 8)

mas_res_df <- mas_1$results
fdr_mas <- mas_res_df %>%
  dplyr::filter(qval < 0.05,
                metadata == "case_control_other_latest")
dim(fdr_mas)
head(fdr_mas)


#_______________________________________________________________________________
# ANCOM-BC ----

ancom_da <- ancombc(phyloseq = dat.analysis, formula = "case_control_other_latest + sex",
                    p_adj_method = "fdr", zero_cut = 0.90, lib_cut = 1000,
                    group = "case_control_other_latest", struc_zero = F,
                    neg_lb = F, tol = 1e-5, max_iter = 100, conserve = F,
                    alpha = 0.05, global = FALSE)

# tst <- unlist(ancom_da[["res"]][["beta"]], )


ancom_res_df <- data.frame(
  features = row.names(ancom_da$res$beta),
  beta = unlist(ancom_da$res$beta),
  se = unlist(ancom_da$res$se),
  W = unlist(ancom_da$res$W),
  p_val = unlist(ancom_da$res$p_val),
  q_val = unlist(ancom_da$res$q_val),
  diff_abn = unlist(ancom_da$res$diff_abn)) %>%
  rownames_to_column(var = "covariate")

fdr_ancom <- ancom_res_df %>%
  dplyr::filter(q_val < 0.05,
                grepl("case_control_other_latest", covariate))
dim(fdr_ancom)
head(fdr_ancom)

#_______________________________________________________________________________
# Overlap

Reduce(intersect, list(fdr_limma$features, fdr_deseq$features, fdr_corncob, fdr_mas$feature, fdr_ancom$features))   #67 features


da_venn <- venn.diagram(
  x = list(fdr_limma$features, fdr_deseq$features, 
           #fdr_corncob, 
           fdr_mas$feature, fdr_ancom$features),
  category.names = c("Limma-Voom" , "DESeq2", 
                     #"corncob", 
                     "MaAsLin2", "ANCOM-BC"),
  filename = NULL,
  fill = c("#8DD3C7", "#FFFFB3", 
           #"#BEBADA", 
           "#FB8072", "#80B1D3"),
  margin = 0.1)

grid::grid.newpage()
grid::grid.draw(da_venn)
