# Machine Learning - Model Comparison Analysis

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
source("src/ml_models.R")
source("src/ml_plots.R")
load("files/low_quality_samples.RData")


tst <- tax_table(phyloseq_objs[["Species"]]) %>% 
  as.data.frame()

# load("files/Phyloseq_Merged_ML.RData") #phyloseq_objs
load("files/Phyloseq_Merged_ML_Rarefied.RData") #phyloseq_objs_rare
csv_prefix <- "data/Machine_Learning_Analysis/model_stats/"
levlist <- c(
  # "Species", "Genus",  "Pathways.slim", 
  "KOs.slim", "PFAMs.slim", "eggNOGs.slim")
ml_models <- c(
  "randomforest"
  # "xgboost"
  )

obj_PD <- subset_samples(dat.species, diagnosis_latest == "Idiopathic PD")
obj_Controls <- subset_samples(dat.species, diagnosis_latest == "No PD Nor Other Neurological Disorder")


source("src/ml_models.R")
firstml <-
  ml_summary(
    obj_A = obj_PD,
    obj_B = obj_Controls,
    label_A = "PD",
    label_B = "Controls",
    analysis = "LOSO",
    model_type = "lasso",
    featSelection =  NULL
  )


# 
# # ML Model loop including feature selection
# for (level in levlist){
#   
#   print_line()
#   print(level)
#   
#   obj <- subset_samples(phyloseq_objs_rare[[level]], donor_id %ni% low_qc[[1]])
#   print(obj)
#   
#   obj_noShanghai <- subset_samples(obj, cohort != "Shanghai")
#   obj_Shanghai <- subset_samples(obj, cohort == "Shanghai")
#   obj_noRush <- subset_samples(obj, cohort != "Rush")
#   obj_Rush <- subset_samples(obj, cohort == "Rush")
#   obj_noTBC <- subset_samples(obj, cohort != "TBC")
#   obj_TBC <- subset_samples(obj, cohort == "TBC")
#   obj_noBonn <- subset_samples(obj, cohort != "Bonn")
#   obj_Bonn <- subset_samples(obj, cohort == "Bonn")
#   
#   for(model_eng in ml_models){
#     for (N in c(25, 100)){
#       
#       df_loso <-
#         ml_loso(
#           Shanghai = obj_Shanghai,
#           noShanghai = obj_noShanghai,
#           TBC = obj_TBC,
#           noTBC = obj_noTBC,
#           Rush = obj_Rush,
#           noRush = obj_noRush,
#           Bonn = obj_Bonn,
#           noBonn = obj_noBonn,
#           model_type = model_eng,
#           featSelection = F,
#           data_type = level,
#           nfeats =  N
#         )
#       
#       write.csv(df_loso, 
#                 file = paste0(csv_prefix, "LOSO_", level, "_", model_eng, "_", N, 
#                               "Rarefied.csv"))
#       df_loso.plot <- cohort_comparison_bars(df_loso)
#       print(df_loso.plot)
#       print_line()
#       ggsave(df_loso.plot,
#              filename = paste("data/Machine_Learning_Analysis/LOSO", level, 
#                               model_eng, "barplot", N, "Rarefied.png", sep = "_"),
#              width = 3.5, height = 4)
#       
#       
#       df_s2s <- 
#         ml_s2s(
#           Shanghai = obj_Shanghai,
#           TBC = obj_TBC,
#           Rush = obj_Rush,
#           Bonn = obj_Bonn,
#           model_type = model_eng,
#           featSelection = F, 
#           data_type = level,
#           nfeats = N
#         )
#       write.csv(df_s2s, 
#                 file = paste0(csv_prefix, "S2S_", level, "_", model_eng, 
#                               "_", N, ".csv"))
#       df_s2s.plot <- ml_heatmap_summary(df_s2s)
#       print(df_s2s.plot)
#       print_line()
#       ggsave(df_s2s.plot, 
#              filename = paste("data/Machine_Learning_Analysis/S2S", level, 
#                               model_eng, "heatmap ", N, "Rarefied.png", 
#                               sep = "_"),
#              width = 4, height = 3)
#     }
#   }
# }
# 

  