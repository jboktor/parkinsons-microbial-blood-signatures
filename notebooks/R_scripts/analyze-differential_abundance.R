# Joe Boktor
# Caltech - Mazmanian Lab

source("src/_load_packages.R")
# base::load("data/Phyloseq_Objects/UHGG/Species_counts.RData")
data_output_dir <- "data/Analyses/differential_abundance/"
dir.create(file.path(data_output_dir), showWarnings = FALSE)
base::load("data/Phyloseq_Objects/Phyloseq_all_outliers_removed.RData")


#_______________________________________________________________________________

# Loop through phyloseq objects for differential abundance testing
refDB_list <- c( "UHGG", "WoL", "RefSeqPlusPF" )#,
  # "UHGG", "WoL", "RefSeqPlusPF") # names(refDB_phyloseq_orm)
level_list <- c("Species") # names(refDB_phyloseq_orm[[ref_DB]])

for (ref_DB in refDB_list){
  for (level in level_list){
    cat("Analyzing: ", ref_DB, level, "\n")
    dir.create(file.path(paste0(data_output_dir, ref_DB)), showWarnings = FALSE)
    
    dat.obj <- refDB_phyloseq_orm[[ref_DB]][[level]] %>% 
      subset_samples(case_control_other_latest != "Other")
    # input_data_df <- data.frame(otu_table(dat.obj))
    # input_metadata_df <- data.frame(sample_data(dat.obj))
    # rownames(input_metadata_df) <- gsub("-", ".", input_metadata_df$participant_id)
    
    # mas_1 <- Maaslin2(
    #   input_data = input_data_df,
    #   input_metadata = input_metadata_df,
    #   output = paste0("data/Analyses/differential_abundance/", ref_DB, "/Maaslin2_output_", level),
    #   min_abundance = 0.0,
    #   min_prevalence = 0.05,
    #   normalization = "TMM",
    #   transform = "NONE",
    #   analysis_method = "NEGBIN",
    #   max_significance = 0.05,
    #   fixed_effects = c("case_control_other_latest", "sex"),
    #   random_effects = "study",
    #   correction = "BH",
    #   standardize = FALSE,
    #   cores = 8)
    
    # out = ancombc(phyloseq = dat.obj, formula = "case_control_other_latest + sex",
    #               p_adj_method = "BH", zero_cut = 0.90, lib_cut = 10,
    #               group = "diagnosis_latest", struc_zero = TRUE, neg_lb = FALSE,
    #               tol = 1e-5, max_iter = 100, conserve = TRUE,
    #               alpha = 0.05, global = TRUE)
    
    
    # Adding a covariate within
    # (formula) - test for differential abundance
    # (phi.formula) - tests for differential dispersion
    # To test a model for only differential abundance while controlling for dispersion
    # add the covariate to (formula, & (phi.formula & phi.formula.null) for variability control)
    # To test exclusively for differential variability (abundance should be controlled for)
    # This is done by adding covariate to (phi.formula (variability), & (formula & formula.null) for abundance control)

      # Testing for differential abundance across diagnosis,
      # controlling effect of diagnosis on dispersion
      model_diagnosis_disp <-
        differentialTest(
          formula = ~ case_control_other,
          phi.formula = ~ 1,
          formula_null = ~ 1,
          phi.formula_null = ~ 1,
          data = dat.obj,
          test = "Wald",
          boot = FALSE,
          fdr_cutoff = 0.05
        )
      summary(model_diagnosis_disp)
      plot(model_diagnosis_disp, B = 1000, total = TRUE)
      save(model_diagnosis_disp, file = paste0(data_output_dir, ref_DB,
                                               "/", level, "_corncob_case_control_other", 
                                               Sys.Date(),".RData"))

      model_diagnosisSex_disp <-
        differentialTest(
          formula = ~ sex,
          phi.formula = ~ 1,
          formula_null = ~ 1,
          phi.formula_null = ~ 1,
          data = dat.obj,
          test = "Wald",
          boot = FALSE,
          fdr_cutoff = 0.05
        )
      summary(model_diagnosisSex_disp)
      plot(model_diagnosisSex_disp, B = 1000, total = TRUE)
      save(model_diagnosisSex_disp, file = paste0(data_output_dir, ref_DB,
                                                  "/", level, "_corncob_sex.RData"))

  }
}


#______________________________________________________________________________

# 
# # # Run ancombc function
# # out = ancombc(phyloseq = dat.species, formula = "diagnosis_latest",
# #               p_adj_method = "BH", zero_cut = 0.90, lib_cut = 10,
# #               group = "diagnosis_latest", struc_zero = TRUE, neg_lb = FALSE,
# #               tol = 1e-5, max_iter = 100, conserve = TRUE,
# #               alpha = 0.05, global = TRUE)
# # 
# # res = out$res$diff_abn
# # res_global = out$res_global
# 
# 
# corncob analysis

# set.seed(2021)
# # Testing for differential abundance across diagnosis, not controlling for any other variables
# model_diagnosis <- differentialTest(
#   formula = ~ diagnosis_latest,
#   phi.formula = ~ 1,
#   formula_null = ~ 1,
#   phi.formula_null = ~ 1,
#   data = dat.species,
#   test = "Wald",
#   boot = FALSE,
#   fdr_cutoff = 0.05
# )
# summary(model_diagnosis)
# plot(model_diagnosis, B = 1000, color = "diagnosis_latest", total = TRUE)


# Adding a covariate within
# (phi.formula) - tests for differential dispersion
# (phi.formula & phi.formula.null) - controls for dispersion
# (formula_null) - controls for abundance

# Testing for differential abundance across diagnosis,
# controlling effect of diagnosis on dispersion
# model_diagnosis_disp <-
#   differentialTest(
#     formula = ~ diagnosis_latest,
#     phi.formula = ~ diagnosis_latest,
#     formula_null = ~ 1,
#     phi.formula_null = ~ diagnosis_latest,
#     data = dat.species,
#     test = "Wald",
#     boot = FALSE,
#     fdr_cutoff = 0.05
#   )
# summary(model_diagnosis_disp)
# plot(model_diagnosis_disp, B = 1000, color = "sex", total = TRUE)
# save(model_diagnosis_disp, file = paste0(data_output_dir,  "UHGG",
#                                          "/corncob_diagnosis_latest_DispControl_", 
#                                          Sys.Date(),".RData"))
# 
# 
# # Jointly testing for differential abundance across diagnosis & sex, 
# # controlling effect of diagnosis & sex on dispersion
# model_diagnosisSex_disp <-
#   differentialTest(
#     formula = ~ diagnosis_latest + sex,
#     phi.formula = ~ diagnosis_latest + sex,
#     formula_null = ~ 1,
#     phi.formula_null = ~ diagnosis_latest + sex,
#     data = dat.species,
#     test = "Wald",
#     boot = FALSE,
#     fdr_cutoff = 0.05
#   )
# summary(model_diagnosisSex_disp)
# plot(model_diagnosisSex_disp, B = 1000, color = "diagnosis_latest", total = TRUE)
# 
# 
# # Jointly testing for differential abundance and differential variability across diagnosis 
# model_diagnosis_abudvar <-
#   differentialTest(
#     formula = ~ diagnosis_latest,
#     phi.formula = ~ diagnosis_latest,
#     formula_null = ~ 1,
#     phi.formula_null = ~ 1,
#     data = dat.species,
#     test = "Wald",
#     boot = FALSE,
#     fdr_cutoff = 0.05
#   )
# summary(model_diagnosis_abudvar)
# plot(model_diagnosis_abudvar, B = 1000, color = "diagnosis_latest", total = TRUE)
# 
# corn_da_qvals <- corn_da$p_fdr %>% as.data.frame()
# fdr_corncob <- corn_da$significant_taxa
# dim(data.frame(fdr_corncob))
# 
# plot(corn_da, B = 1000, color = "sex", total = TRUE)
# summary(corn_da)
# 


