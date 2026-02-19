# Joe Boktor
# Caltech - Mazmanian lab

# source("src/_load_packages.R")
# source("src/_plot-functions.R")
# library(tidymodels)
# library(parsnip)
# # library(skimr)
# library(finetune)

# ______________________________________________________________________________
## Build Models ----

library(tidyverse)
library(parsnip)
library(rsample)
library(yardstick)
library(recipes)
library(workflows)
library(dials)
library(tune)

ps_trim_norm <- readRDS(
  glue(
    "{wkdir}/data/processed/",
    "phyloseq_objects/decontaminated/",
    "2024-02-27_WGS_KrakenUnique_phyloseq_humanreadnorm.rds"
  )
)

# load normalized phyloseq object and filter for cohorts of interest
ps_model <- ps_trim_norm %>%
  # subset_samples(study %in% c("BioFIND", "HBS", "PDBP", "PPMI")) %>%
  subset_samples(case_control_other_latest != "Other") %>%
  core(detection = 0, prevalence = 1/100)

# collect metdata
metadata_df <- microbiome::meta(ps_model) %>%
  dplyr::select(participant_id, case_control_other_latest, sex, study) %>%
  glimpse


set.seed(42)
df_plot <- abundances(ps_model) %>%
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "participant_id") %>%
  left_join(metadata_df) %>%
  rsample::initial_split(strata = case_control_other_latest)

ml_train <- rsample::training(df_plot)
ml_test <- rsample::testing(df_plot)
ml_metrics <- yardstick::metric_set(accuracy, roc_auc, mn_log_loss)


# Set up 10-fold cross-validation recipe
set.seed(42)
train_folds <- vfold_cv(ml_train, v = 10, strata = case_control_other_latest)
ml_recipe <- recipes::recipe(case_control_other_latest ~., data = ml_train) %>%
  recipes::update_role(participant_id, new_role = "participant_id") %>%
  recipes::update_role(sex, new_role = "sex") %>%
  recipes::update_role(study, new_role = "study") %>%
  recipes::step_nzv(all_predictors()) %>%
  recipes::step_normalize(all_predictors())

# Define tunable Lasso model
tune_spec <- parsnip::logistic_reg(penalty = tune(), mixture = 1) %>%
  parsnip::set_engine("glmnet")

# defining workflow
wf <- workflows::workflow() %>%
  workflows::add_recipe(ml_recipe) %>%
  workflows::add_model(tune_spec)

# Create a grid of penalty values to test (LASSO MODEL)
lambda_grid <- dials::grid_regular(penalty(), levels = 100)

# TUNE MODEL
doParallel::registerDoParallel()
lasso_grid <- tune_grid(wf, resamples = train_folds, grid = lambda_grid)

## Evaluate model results
tune::show_best(lasso_grid, metric = "roc_auc")
# select optimal penalty by filtering largest rocauc
best_aucroc <- select_best(lasso_grid, "roc_auc")
# visualize model metrics of grid 

# model_performance <- 
#   lasso_grid %>% 
#   collect_metrics() %>% 
#   ggplot(aes(penalty, mean, color = .metric)) +
#   geom_errorbar(aes(ymin = mean - std_err,
#                     ymax = mean + std_err),
#                 alpha = 0.5) +
#   geom_line(linewidth = 1.25, show.legend = F) +
#   facet_wrap(~.metric, scales = "free", nrow = 2) +
#   theme_bw() + 
#   scale_x_log10() +
#   scale_color_viridis_d(option = "cividis", begin = .9, end = 0) +
#   theme(legend.position = "none")
# model_performance

# ggsave(model_performance,
#   filename = glue("{wkdir}/figures/ml/lasso_vanilla_training_perf.png"),
#   width = 5, height = 3
# )

# Finalize and fit workflow with tuned parameters
train_lasso <- finalize_workflow(wf, best_aucroc) %>% fit(ml_train)

# Predictions on test data
hold_out_set <- finalize_workflow(wf, best_aucroc) %>% last_fit(df_plot)

conf_matrix <- collect_predictions(hold_out_set) %>% 
  conf_mat(case_control_other_latest, .pred_class) %>% 
  autoplot(type = "heatmap")
print(conf_matrix)

aurocplot <- collect_predictions(hold_out_set) %>% 
  roc_curve(case_control_other_latest, .pred_Case) %>%
  autoplot()
aurocplot

aupr_plot <- collect_predictions(hold_out_set) %>% 
  pr_curve(case_control_other_latest, .pred_Case) %>%
  autoplot()
aupr_plot

# ggsave(conf_matrix, filename = paste0("figures/machine_learning/lasso_vanilla/",
#                                       refDB, "_", level, "_confusion.png"),
#        dpi = 600, width = 3.25, height =3)
# ggsave(aurocplot, filename = paste0("figures/machine_learning/lasso_vanilla/",
#                                     refDB, "_", level, "_AUROC.png"),
#        dpi = 600, width = 4.5, height =4.5)
# ggsave(aupr.plot, filename = paste0("figures/machine_learning/lasso_vanilla/",
#                                     refDB, "_", level, "_AUPR.png"),
#        dpi = 600, width = 4.5, height =4.5)


collect_metrics <- function(hold_out_set, group) {
  preds <- collect_predictions(hold_out_set)
  pred_accuracy <- metrics(preds, group, .pred_class) %>%
    filter(.metric == "accuracy")
  pred_auroc <- roc_auc(preds, group, .pred_Case)
  pred_prauc <- pr_auc(preds, group, .pred_Case)
  pred_ppv <- ppv(preds, group, .pred_class)
  pred_npv <- npv(preds, group, .pred_class)
  pred_sensitivity <- sensitivity(preds, group, .pred_class)
  pred_specificity <- specificity(preds, group, .pred_class)
  model_stats <- dplyr::bind_rows(
    pred_accuracy, pred_auroc, pred_prauc,
    pred_ppv, pred_npv, pred_sensitivity, pred_specificity
  )
  return(model_stats)
}

model_stats <- collect_metrics(hold_out_set, "case_control_other_latest")




model_stats_summary <- model_stats %>% 
  ggplot(aes(x = .estimate, y = .metric)) +
  geom_col(width = 0.5) +
  scale_y_discrete(
    labels = c(
      "accuracy" = "Accuracy",
      "roc_auc" = "AUROC",
      "pr_auc" = "AUPR",
      "Sensitivity" = "Sensitivity",
      "Specificity" = "Specificity",
      "PPV" = "PPV",
      "NPV" = "NPV")) +
  labs(x = "", y = "") +
  geom_text(aes(label = round(.estimate, digits = 3), x = .estimate - 0.075, y = .metric), 
            color = "white") +
  theme_bw() 
model_stats_summary

ggsave(model_stats_summary, filename = paste0("figures/machine_learning/lasso_vanilla/",
                                              refDB, "_", level, "_model_stats_summary.png"),
       dpi = 600, width = 4.5, height =3)


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

