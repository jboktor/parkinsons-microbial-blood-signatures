# Joe Boktor
# Caltech - Mazmanian lab

source("src/_load_packages.R")
source("src/_plot-functions.R")
library(tidyverse)
library(tidymodels)
library(skimr)
library(finetune)

base::load("data/Phyloseq_Objects/Phyloseq_all_outliers_removed.RData") #phyloseq_objs
refDB <- "UHGG"
level <- "Species"
dat.obj <- refDB_phyloseq_orm[[refDB]][[level]] %>% 
  subset_samples(case_control_other_latest != "Other")


metadat <- meta(dat.obj)
model_matrix <-
  model.matrix( ~ 0 + case_control_other_latest + sex + study,
                data = metadat)

dge <- 
  dat.obj %>% 
  phyloseq_to_deseq2(~ case_control_other_latest + sex + study) %>% 
  as.DGEList()
# dge_highQC <- filterByExpr(dge, model_matrix)
# dge.filtered <- dge[dge_highQC,,keep.lib.sizes=FALSE]
dge.norm <- dge %>%  #dge.filtered %>% 
  calcNormFactors(method = "TMM") %>% 
  voom(design = model_matrix, plot = TRUE, save.plot = TRUE, normalize.method="none")
print(dim(t(dge.norm$E)))


# rawdat <- dge.norm$E
# adj.var <- model.matrix(~ study, data=metadat)
# bio.var.input <- model.matrix(~ sex + age_at_baseline, data = metadat)
# snm.obj <-
#   snm(
#     raw.dat = rawdat,
#     bio.var = bio.var.input,
#     adj.var = adj.var,
#     rm.adj = TRUE,
#     verbose = TRUE,
#     diagnose = TRUE
#   )


#______________________________________________________________________________
## Build Models ----


set.seed(42)
metadata_df <- metadat %>% 
  dplyr::select(participant_id, case_control_other_latest, sex, study)
df_plot <- 
  t(dge.norm$E) %>%
  # t(snm.obj$norm.dat) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "participant_id") %>%
  left_join(metadata_df) %>%
  initial_split(strata = case_control_other_latest)
# df_plot <- 
#   dat.obj %>% 
#   abundances() %>% t() %>% 
#   as.data.frame() %>% 
#   rownames_to_column(var = "participant_id") %>% 
#   left_join(metadata_df) %>%
#   initial_split(strata = case_control_other_latest)

ml_train <- training(df_plot)
ml_test <- testing(df_plot)
ml_metrics <- metric_set(accuracy, roc_auc, mn_log_loss)


set.seed(42)
train_folds <- vfold_cv(ml_train, v = 10, strata = case_control_other_latest)
train_folds

ml_recipe <- recipe(case_control_other_latest ~., data = ml_train) %>%
  update_role(participant_id, new_role = "participant_id") %>%
  update_role(sex, new_role = "sex") %>%
  update_role(study, new_role = "study") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_predictors())
ml_recipe

#______________________________________________________________________________
#                             LASSO Model ----
#______________________________________________________________________________

{
  

# Define tunable Lasso model
tune_spec <- logistic_reg(penalty = tune(), mixture = 1) %>%
  set_engine("glmnet")

wf <- workflow() %>% 
  add_recipe(ml_recipe) %>% 
  add_model(tune_spec)


#_____________________________________________________
#                TUNE LASSO MODEL  
#_____________________________________________________
# Create a grid of penalty values to test
lambda_grid <- grid_regular(penalty(), levels = 50)
# ctrl <- control_grid(save_pred = TRUE, verbose = TRUE)
doParallel::registerDoParallel()

set.seed(42)
lasso_grid <- tune_grid(wf, resamples = train_folds, grid = lambda_grid)

## Evaluate model results
show_best(lasso_grid)
# select optimal penalty by filtering largest rocauc
best_aucroc <- select_best(lasso_grid, "roc_auc")
# visualize model metrics of grid 
model_performance <- 
  lasso_grid %>% 
  collect_metrics() %>% 
  ggplot(aes(penalty, mean, color = .metric)) +
  geom_errorbar(aes(ymin = mean - std_err,
                    ymax = mean + std_err),
                alpha = 0.5) +
  geom_line(size = 1.25, show.legend = F) +
  facet_wrap(~.metric, scales = "free", nrow = 2) +
  theme_bw() + 
  scale_x_log10() +
  scale_color_viridis_d(option = "cividis", begin = .9, end = 0) +
  theme(legend.position = "none")
print(model_performance)

ggsave(model_performance, filename = paste0("figures/machine_learning/lasso_vanilla/training_perf.png"),
       dpi = 600, width = 5, height =3)

# Finalize and fit workflow with tuned parameters
train_lasso <-
  wf %>%
  finalize_workflow(best_aucroc) %>%
  fit(ml_train)

# Predictions on test data
hold_out_set <-
  wf %>%
  finalize_workflow(best_aucroc) %>%
  last_fit(df_plot)
hold_out_set

conf_matrix <- 
  hold_out_set %>% 
  collect_predictions() %>% 
  conf_mat(case_control_other_latest, .pred_class) %>% 
  autoplot(type = "heatmap")
print(conf_matrix)

aurocplot <- 
  hold_out_set %>% 
  collect_predictions() %>% 
  roc_curve(case_control_other_latest, .pred_Case) %>%
  autoplot()
aurocplot

aupr.plot <- 
  hold_out_set %>% 
  collect_predictions() %>% 
  pr_curve(case_control_other_latest, .pred_Case) %>%
  autoplot()
aupr.plot

ggsave(conf_matrix, filename = paste0("figures/machine_learning/lasso_vanilla/",
                                      refDB, "_", level, "_confusion.png"),
       dpi = 600, width = 3.25, height =3)
ggsave(aurocplot, filename = paste0("figures/machine_learning/lasso_vanilla/",
                                    refDB, "_", level, "_AUROC.png"),
       dpi = 600, width = 4.5, height =4.5)
ggsave(aupr.plot, filename = paste0("figures/machine_learning/lasso_vanilla/",
                                    refDB, "_", level, "_AUPR.png"),
       dpi = 600, width = 4.5, height =4.5)

pred_accuracy <- hold_out_set %>% 
  collect_predictions() %>% 
  metrics(case_control_other_latest, .pred_class) %>% 
  filter(.metric == "accuracy")
pred_auroc <- hold_out_set %>% 
  collect_predictions() %>% 
  roc_auc(case_control_other_latest, .pred_Case)
pred_prauc <- hold_out_set %>% 
  collect_predictions() %>% 
  pr_auc(case_control_other_latest, .pred_Case)
pred_ppv <- hold_out_set %>% 
  collect_predictions() %>% 
  ppv(case_control_other_latest, .pred_class)
pred_npv <- hold_out_set %>% 
  collect_predictions() %>% 
  npv(case_control_other_latest, .pred_class)
pred_sensitivity <- hold_out_set %>% 
  collect_predictions() %>% 
  sensitivity(case_control_other_latest, .pred_class)
pred_specificity <- hold_out_set %>% 
  collect_predictions() %>% 
  specificity(case_control_other_latest, .pred_class)



model_stats <- bind_rows(pred_accuracy, pred_auroc, pred_prauc, 
                         pred_ppv, pred_npv, pred_sensitivity, pred_specificity)

model_stats_summary <- model_stats %>% 
  ggplot(aes(x = .estimate, y = .metric)) +
  geom_col(width = 0.5) +
  scale_y_discrete(
    labels = c(
      "accuracy" = "Accuracy",
      "roc_auc" = "AUROC",
      "pr_auc" = "AUPR",
      "sens" = "Sensitivity",
      "spec" = "Specificity",
      "ppv" = "PPV",
      "npv" = "NPV")) +
  labs(x = "", y = "") +
  geom_text(aes(label = round(.estimate, digits = 3), x = .estimate - 0.075, y = .metric), 
            color = "white") +
  theme_bw() 
model_stats_summary
ggsave(model_stats_summary, filename = paste0("figures/machine_learning/lasso_vanilla/",
                                              refDB, "_", level, "_model_stats_summary.png"),
       dpi = 600, width = 4.5, height =3)



}
#______________________________________________________________________________
#                             XGBoost Model ----
#______________________________________________________________________________


# Define tunable XGBoost model
stopping_spec <-
  boost_tree(
    trees = 1000,
    mtry = tune(),
    learn_rate = tune(),
    stop_iter = tune()
  ) %>%
  set_engine("xgboost", validation = 0.2) %>%
  set_mode("classification")

stopping_grid <-
  grid_latin_hypercube(
    mtry(range = c(5L, 20L)), ## depends on number of columns in data
    learn_rate(range = c(-5, -1)), ## keep pretty big
    stop_iter(range = c(10L, 50L)), ## bigger than default
    size = 50
  )

wf <- workflow() %>% 
  add_recipe(ml_recipe) %>% 
  add_model(stopping_spec)

#_____________________________________________________
#                TUNE XGBoost MODEL  
#_____________________________________________________

doParallel::registerDoParallel()
set.seed(42)
XGBoost_grid <- tune_grid(
  wf,
  train_folds,
  grid = stopping_grid,
  metrics = ml_metrics
)
## Evaluate model results
show_best(XGBoost_grid)
autoplot(XGBoost_grid) # + theme_light(base_family = "IBMPlexSans")
# select optimal penalty by filtering largest rocauc
best_aucroc <- select_best(XGBoost_grid, "roc_auc")

# ggsave(model_performance, filename = paste0("figures/machine_learning/XGBoost_vanilla/training_perf.png"),
#        dpi = 600, width = 5, height =3)

# Finalize and fit workflow with tuned parameters
train_XGBoost <-
  wf %>%
  finalize_workflow(best_aucroc) %>%
  fit(ml_train)

# Predictions on test data
hold_out_set <-
  wf %>%
  finalize_workflow(best_aucroc) %>%
  last_fit(df_plot)
hold_out_set

# Save Models
xgboost_models <-
  list(
    "10_fold_CV" = XGBoost_grid,
    "trained_model" = train_XGBoost,
    "testset_pred" = hold_out_set
  )
# saveRDS(xgboost_models, file = paste0("data/ML_models/XGBoost_", refDB, "_", level, ".rds"))


pred_accuracy <- hold_out_set %>% 
  collect_predictions() %>% 
  metrics(case_control_other_latest, .pred_class) %>% 
  filter(.metric == "accuracy")
pred_auroc <- hold_out_set %>% 
  collect_predictions() %>% 
  roc_auc(case_control_other_latest, .pred_Case)
pred_prauc <- hold_out_set %>% 
  collect_predictions() %>% 
  pr_auc(case_control_other_latest, .pred_Case)
pred_ppv <- hold_out_set %>% 
  collect_predictions() %>% 
  ppv(case_control_other_latest, .pred_class)
pred_npv <- hold_out_set %>% 
  collect_predictions() %>% 
  npv(case_control_other_latest, .pred_class)
pred_sensitivity <- hold_out_set %>% 
  collect_predictions() %>% 
  sensitivity(case_control_other_latest, .pred_class)
pred_specificity <- hold_out_set %>% 
  collect_predictions() %>% 
  specificity(case_control_other_latest, .pred_class)



model_stats <- bind_rows(pred_accuracy, pred_auroc, pred_prauc, 
                         pred_ppv, pred_npv, pred_sensitivity, pred_specificity)

model_stats_summary <- model_stats %>% 
  ggplot(aes(x = .estimate, y = .metric)) +
  geom_col(width = 0.5) +
  scale_y_discrete(
    labels = c(
      "accuracy" = "Accuracy",
      "roc_auc" = "AUROC",
      "pr_auc" = "AUPR",
      "sens" = "Sensitivity",
      "spec" = "Specificity",
      "ppv" = "PPV",
      "npv" = "NPV")) +
  labs(x = "", y = "") +
  geom_text(aes(label = round(.estimate, digits = 3), x = .estimate - 0.075, y = .metric), 
            color = "white") +
  theme_bw() 
model_stats_summary

ggsave(model_stats_summary, filename = paste0("figures/machine_learning/XGBoost_vanilla/",
                                              refDB, "_", level, "_model_stats_summary.png"),
       dpi = 600, width = 4.5, height =3)






conf_matrix <- 
  hold_out_set %>% 
  collect_predictions() %>% 
  conf_mat(case_control_other_latest, .pred_class) %>% 
  autoplot() +
  scale_fill_gradient(low="#D6EAF8",high = "#2E86C1")
# autoplot(type = "heatmap")
print(conf_matrix)

aurocplot <- 
  hold_out_set %>% 
  collect_predictions() %>% 
  roc_curve(case_control_other_latest, .pred_Case) %>%
  autoplot() +
  annotate("text", x = 0.6, y = 0.1, 
           label = paste0("AUROC: ", round(pred_auroc$.estimate, digits = 2))) +
  theme_bw() +
  coord_fixed() +
  theme(panel.grid = element_blank())
aurocplot

aupr.plot <- 
  hold_out_set %>% 
  collect_predictions() %>% 
  pr_curve(case_control_other_latest, .pred_Case) %>%
  autoplot() +
  annotate("text", x = 0.6, y = 0.1, 
           label = paste0("AUPR: ", round(pred_prauc$.estimate, digits = 2))) +
  theme_bw() +
  coord_fixed() +
  theme(panel.grid = element_blank())
aupr.plot

# Variables of importance
library(vip)
## use this fitted workflow `extract_workflow(stopping_fit)` to predict on new data
extract_workflow(hold_out_set) %>%
  extract_fit_parsnip() %>%
  vip(num_features = 15, geom = "point")

figOutputDir <- "figures/machine_learning/"
ggsave(conf_matrix, filename = paste0(figOutputDir, "XGBoost_", refDB, "_", level, "_confusion.svg"),
       dpi = 600, width = 3.25, height =3)
ggsave(aurocplot, filename = paste0(figOutputDir, "XGBoost_", refDB, "_", level, "_AUROC.svg"),
       dpi = 600, width = 2.5, height =2.5)
ggsave(aupr.plot, filename = paste0(figOutputDir, "XGBoost_", refDB, "_", level, "_AUPR.svg"),
       dpi = 600, width = 2.5, height = 2.5)

