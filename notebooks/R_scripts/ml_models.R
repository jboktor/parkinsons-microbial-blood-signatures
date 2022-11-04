# Machine-Learning-Models

# Packages for xgboost
# library(glue)
library(ModelMetrics)
library(OpenMPController) # for Kaggle backend
library(xgboost)
library(parallel)
library(doParallel)
library(mlbench)
library(caret)
library(MLeval)

library(tidyverse)
library(parsnip)
library(rsample)
library(yardstick)
library(recipes)
library(workflows)
library(dials)
library(tune)
library(janitor)
library(remotes)
library(glmnet)
# lightgbm and catboost for parsnip
# remotes::install_github("curso-r/treesnip", dependencies = T) 
# library(treesnip)







#_______________________________________________________________________________
#                           ML Functions
#_______________________________________________________________________________

prep_ml_input <- function(obj){
  all_abund <- abundances(obj) %>%
    t() %>% as.data.frame() %>% 
    rownames_to_column("participant_id")
  all_meta <- meta(obj) %>% 
    dplyr::select(participant_id, diagnosis_latest)
  ml_input <- inner_join(all_meta, all_abund) %>% 
    as_tibble() %>% 
    dplyr::select(-participant_id)
  return(ml_input)
}
#______________________________________

# prep_mRMR_input <- function(obj){
#   
#   all_abund <- abundances(obj) %>%
#     t() %>% as.data.frame() %>%
#     rownames_to_column("donor_id")
#   all_meta <- process_meta(obj, cohort = "Merged_ML") %>%
#     dplyr::select(donor_id, PD)
#   mRMR_input <- left_join(all_meta, all_abund) %>%
#     as.data.frame() %>%
#     column_to_rownames(var = "donor_id") %>% 
#     mutate(PD = factor(PD, levels = c("No", "Yes"))) %>%
#     mutate(PD = as.numeric(PD)) 
#   
#   return(mRMR_input)
# }

#______________________________________
#             mRMRe wrapper 
#______________________________________

# feature_selection <- function(df, n_features, n_algorithm = 1){
#   
#   mdat <- mRMR.data(data = data.frame(df))
#   mRMR <- mRMR.ensemble(data = mdat, target_indices = 1, 
#                         feature_count = n_features, 
#                         method = "bootstrap", bootstrap_count = 3,
#                         solution_count = n_algorithm)
#   selections <- df[, solutions(mRMR)$`1`] %>% colnames()
#   return(selections)
# }

#______________________________________
#            AUROC filter method 
#______________________________________

feature_selection_auroc <- function(df, n_features, subgroup, data_type){
  
  feats <- df %>% 
    dplyr::filter(data_level == data_type) %>%
    dplyr::filter(study == subgroup) %>%
    dplyr::mutate(auroc_center = auroc - 0.5) %>% 
    slice_max(order_by = abs(auroc_center), n = n_features) %>% 
    dplyr::select(feature)
  return(feats$feature)
}

#_______________________________________________________________________________
#                       Tidymodels LASSO 
#_______________________________________________________________________________

lasso_model <- function(train, test){
  
  # # TROUBLE
  # obj_PD <- subset_samples(dat.species, diagnosis_latest == "Idiopathic PD")
  # obj_Controls <- subset_samples(dat.species, diagnosis_latest == "No PD Nor Other Neurological Disorder")
  # ml_input_A <- prep_ml_input(obj_PD)
  # ml_input_B <- prep_ml_input(obj_Controls)
  # train = ml_input_A
  # test = ml_input_B
  ml_input <- prep_ml_input(dat.species)
  

  # combined <- bind_rows(train, test)
  # ind <- list(
  #   analysis = seq(nrow(train)),
  #   assessment = nrow(train) + seq(nrow(test)))
  set.seed(2021)
  splits <- initial_split(ml_input,  strata = diagnosis_latest)
  ml_train <- training(splits)
  ml_test <- testing(splits)
  
  # # Generate 10-fold CV model with 5 repetitions, stratified by PD status
  set.seed(42)
  boot_splits <- rsample::bootstraps(ml_train, times = 25, strata = diagnosis_latest)
  
  # Define tunable Lasso model
  tune_spec <- logistic_reg(penalty = tune(), mixture = 1) %>%
    set_engine("glmnet")
  
  ml_recipe <- recipe(diagnosis_latest ~ ., data = ml_train) %>% 
    # update_role(participant_id, new_role = "participant_id") %>% 
    update_role(diagnosis_latest, new_role = "diagnosis_latest") %>%
    step_nzv(all_numeric(), -all_outcomes()) %>%
    step_normalize(all_numeric(), -all_outcomes())
  
  filter_obj <- prep(ml_recipe, training = ml_train)
  filtered_te <- bake(filter_obj, ml_test)
  any(names(filtered_te) == "sparse")
  
  wf <- workflow() %>% 
    add_recipe(ml_recipe) %>% 
    add_model(tune_spec)
  
  #_____________________________________________________
  #                TUNE LASSO MODEL  
  #_____________________________________________________
  # Create a grid of penalty values to test
  lambda_grid <- grid_regular(penalty(), levels = 25)
  # ctrl <- control_grid(save_pred = TRUE, verbose = TRUE)
  doParallel::registerDoParallel()
  
  set.seed(42)
  lasso_grid <- tune_grid(wf, resamples = boot_splits, grid = lambda_grid)
  
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
  
  # Extract metrics from optimal models (average of N bootstraps)
  train_metrics <- 
    lasso_grid %>% 
    collect_metrics() %>% 
    filter(penalty == best_aucroc$penalty)
  
  # Finalize and fit workflow with tuned parameters
  train_lasso <-
    wf %>%
    finalize_workflow(best_aucroc) %>%
    fit(ml_train)
  
  # Predictions on test data
  test_lasso <- 
    last_fit(train_lasso, split = splits) 
  test_metrics <- test_lasso %>%
    collect_metrics()
  
  conf_matrix <- 
    test_lasso %>% 
    collect_predictions() %>% 
    conf_mat(PD, .pred_class) %>% 
    autoplot()
  print(conf_matrix)
  
  output <- list(
    "bootstrap_tune_model" = lasso_grid,
    "train_lasso" = train_lasso,
    "test_lasso" = test_lasso,
    "train_metrics" = train_metrics,
    "test_metrics" = test_metrics
  )
  return(output)
}


#_______________________________________________________________________________
#                   Tidymodels RANDOM FOREST 
#_______________________________________________________________________________

rf_model <- function(train, test){
  
  # # TROUBLE
  # train = ml_input_A
  # test = ml_input_B
  
  combined <- bind_rows(train, test)
  ind <- list(
    analysis = seq(nrow(train)),
    assessment = nrow(train) + seq(nrow(test)))
  splits <- make_splits(ind, combined)
  ml_train <- training(splits)
  ml_test <- testing(splits)

  # Generate Grid of bootstrapped testing values for training 
  set.seed(42)
  boot_splits <- rsample::bootstraps(ml_train, times = 25, strata = PDxcohort)

  # Define tunable RF model
  ranger_spec <-
    rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
    set_mode("classification") %>%
    set_engine("ranger") 
    
  ml_recipe <- recipe(formula = PD ~ ., data = ml_train) %>% 
    update_role(donor_id, new_role = "donor_id") %>% 
    update_role(cohort, new_role = "cohort") %>% 
    update_role(paired, new_role = "paired") %>%
    update_role(PDxcohort, new_role = "PDxcohort") %>%
    step_zv(all_numeric(), -all_outcomes()) %>% 
    step_normalize(all_numeric(), -all_outcomes())
  
  wf <- workflow() %>% 
    add_recipe(ml_recipe) %>% 
    add_model(ranger_spec)
  
  #_____________________________________________________
  #                TUNE RF MODEL  
  #_____________________________________________________

  doParallel::registerDoParallel()
  set.seed(42)
  ranger_tune <- tune_grid(wf, resamples = boot_splits, grid = 25)
  
  # select optimal penalty by filtering largest rocauc
  best_aucroc <- select_best(ranger_tune, "roc_auc")

  # visualize model metrics of grid 
  print(autoplot(ranger_tune))
  
  # Extract metrics from optimal models (average of N bootstraps)
  train_metrics <- 
    ranger_tune %>% 
    collect_metrics() %>% 
    filter(mtry == best_aucroc$mtry,
           min_n == best_aucroc$min_n)
  
  # Finalize and fit workflow with tuned parameters
  train_rf <-
    wf %>%
    finalize_workflow(best_aucroc) %>%
    fit(ml_train)
  
  # Predictions on test data
  test_rf <- 
    last_fit(train_rf, split = splits) 
  test_metrics <- test_rf %>%
    collect_metrics()
  
  conf_matrix <- 
    test_rf %>% 
    collect_predictions() %>% 
    conf_mat(PD, .pred_class) %>% 
    autoplot()
  print(conf_matrix)
  
  output <- list(
    "bootstrap_tune_model" = ranger_tune,
    "train_rf" = train_rf,
    "test_rf" = test_rf,
    "train_metrics" = train_metrics,
    "test_metrics" = test_metrics
  )
  return(output)
}


#_______________________________________________________________________________
#                       Tidymodels XGBoost 
#_______________________________________________________________________________


xgb_model <- function(train, test){
  
  # # TROUBLE
  # obj_A <- obj_Shanghai
  # obj_B <- obj_TBC
  # ml_input_A <- prep_ml_input(obj_A)
  # ml_input_B <- prep_ml_input(obj_B)
  # train = ml_input_A
  # test = ml_input_B
  
  combined <- bind_rows(train, test)
  ind <- list(
    analysis = seq(nrow(train)),
    assessment = nrow(train) + seq(nrow(test)))
  splits <- make_splits(ind, combined)
  ml_train <- training(splits)
  ml_test <- testing(splits)
  
  set.seed(42)
  # Generate Grid of bootstrapped testing values for training 
  boot_splits <- rsample::bootstraps(ml_train, times = 25, strata = PDxcohort)
  
  # Define tunable model
  xgboost_spec <- 
    boost_tree(
      trees = 1000,
      min_n = tune(),
      tree_depth = tune(),
      learn_rate = tune(),
      loss_reduction = tune(),
      sample_size = tune(),
      mtry = tune(),
    ) %>%
    set_mode("classification") %>% 
    set_engine("xgboost") 
  
  ml_recipe <- recipe(formula = PD ~ ., data = ml_train) %>% 
    update_role(donor_id, new_role = "donor_id") %>% 
    update_role(cohort, new_role = "cohort") %>% 
    update_role(paired, new_role = "paired") %>%
    update_role(PDxcohort, new_role = "PDxcohort") %>%
    step_zv(all_numeric(), -all_outcomes()) %>% 
    step_normalize(all_numeric(), -all_outcomes())
  
  wf <- workflow() %>% 
    add_recipe(ml_recipe) %>% 
    add_model(xgboost_spec)
  
  #_____________________________________________________
  #                TUNE MODEL  
  #_____________________________________________________
  
  xgb_grid <-
    grid_max_entropy(
      # trees = 1000,
      min_n(),
      learn_rate(),
      loss_reduction(),
      sample_size = sample_prop(),
      tree_depth(),
      finalize(mtry(), ml_train),
      size = 25
    )
  
  doParallel::registerDoParallel()
  set.seed(42)
  xgb_tune <- tune_grid(wf, resamples = boot_splits, grid = xgb_grid)
  
  # select optimal penalty by filtering largest rocauc
  best_aucroc <- select_best(xgb_tune, "roc_auc")

  # visualize model metrics of grid 
  xgb_tune.plot <- xgb_tune %>% 
    collect_metrics() %>% 
    filter(.metric == "roc_auc") %>% 
    dplyr::select(mean, mtry:sample_size) %>% 
    pivot_longer(mtry:sample_size, names_to = "parameter", values_to = "value") %>% 
    ggplot(aes(value, mean, color = parameter)) +
    geom_point(show.legend = F) +
    facet_wrap(~ parameter, scales = "free_x", nrow = 2)
  print(xgb_tune.plot)
  
  # Extract metrics from optimal models (average of N bootstraps)
  train_metrics <- 
    xgb_tune %>% 
    collect_metrics() %>% 
    filter(mtry == best_aucroc$mtry,
           min_n == best_aucroc$min_n)
  
  # Finalize trained workflow - use on hold out test sets
  train_xgb <-
    wf %>%
    finalize_workflow(best_aucroc) %>%
    fit(ml_train)
  
  # Predictions on test data
  test_xgb <- 
    last_fit(train_xgb, split = splits) 
  test_metrics <- test_xgb %>%
    collect_metrics()
  
  conf_matrix <- 
    test_xgb %>% 
    collect_predictions() %>% 
    conf_mat(PD, .pred_class) %>% 
    # roc_curve(PD, .pred_No) %>% 
    autoplot()
  print(conf_matrix)
  
  output <- list(
    "bootstrap_tune_model" = xgb_tune,
    "train_xgb" = train_xgb,
    "test_xgb" = test_xgb,
    "train_metrics" = train_metrics,
    "test_metrics" = test_metrics
  )
  return(output)
}


#_______________________________________________________________________________
#                       Tidymodels LightGBM 
#_______________________________________________________________________________


lgbm_model <- function(train, test){
  
  # # TROUBLE
  # obj_A <- obj_Shanghai
  # obj_B <- obj_TBC
  # ml_input_A <- prep_ml_input(obj_A)
  # ml_input_B <- prep_ml_input(obj_B)
  # train = ml_input_A
  # test = ml_input_B

  combined <- bind_rows(train, test)
  ind <- list(
    analysis = seq(nrow(train)),
    assessment = nrow(train) + seq(nrow(test)))
  splits <- make_splits(ind, combined)
  ml_train <- training(splits)
  ml_test <- testing(splits)
  
  set.seed(42)
  # Generate Grid of bootstrapped testing values for training 
  boot_splits <- rsample::bootstraps(ml_train, times = 25, strata = PDxcohort)
  
  # Define tunable model
  lgbm_spec <- 
    boost_tree(
      trees = tune(),
      min_n = tune(),
      tree_depth = tune(),
      learn_rate = tune(),
      loss_reduction = tune(),
      sample_size = tune(),
      mtry = tune(),
    ) %>%
    set_mode("classification") %>% 
    set_engine("lightgbm") 
  
  ml_recipe <- recipe(formula = PD ~ ., data = ml_train) %>% 
    update_role(donor_id, new_role = "donor_id") %>% 
    update_role(cohort, new_role = "cohort") %>% 
    update_role(paired, new_role = "paired") %>%
    update_role(PDxcohort, new_role = "PDxcohort") %>%
    step_zv(all_numeric(), -all_outcomes()) %>% 
    step_normalize(all_numeric(), -all_outcomes())
  
  wf <- workflow() %>% 
    add_recipe(ml_recipe) %>% 
    add_model(lgbm_spec)
  
  #_____________________________________________________
  #                TUNE MODEL  
  #_____________________________________________________
  
  lgbm_grid <-
    grid_max_entropy(
      trees(),
      min_n(),
      learn_rate(),
      loss_reduction(),
      sample_size = sample_prop(),
      tree_depth(),
      finalize(mtry(), ml_train),
      size = 30
    )
  
  doParallel::registerDoParallel()
  set.seed(42)
  lgbm_tune <- tune_grid(wf, resamples = boot_splits, grid = lgbm_grid)

  # select optimal penalty by filtering largest rocauc
  best_aucroc <- select_best(lgbm_tune, "roc_auc")
  
  # visualize model metrics of grid 
  lgbm_tune.plot <- lgbm_tune %>% 
    collect_metrics() %>% 
    filter(.metric == "roc_auc") %>% 
    dplyr::select(mean, mtry:sample_size) %>% 
    pivot_longer(mtry:sample_size, names_to = "parameter", values_to = "value") %>% 
    ggplot(aes(value, mean, color = parameter)) +
    geom_point(show.legend = F) +
    facet_wrap(~ parameter, scales = "free_x", nrow = 2)
  print(lgbm_tune.plot)
  
  # Extract metrics from optimal models (average of N bootstraps)
  train_metrics <- 
    lgbm_tune %>% 
    collect_metrics() %>% 
    filter(mtry == best_aucroc$mtry,
           min_n == best_aucroc$min_n)
  
  # Finalize trained workflow - use on hold out test sets
  train_lgbm <-
    wf %>%
    finalize_workflow(best_aucroc) %>%
    fit(ml_train)
  
  # Predictions on test data
  test_lgbm <- 
    last_fit(train_lgbm, split = splits) 
  test_metrics <- test_lgbm %>%
    collect_metrics()
  
  conf_matrix <- 
    test_lgbm %>% 
    collect_predictions() %>% 
    conf_mat(PD, .pred_class) %>%
    # roc_curve(PD, .pred_No) %>%
    autoplot()
  print(conf_matrix)
  
  output <- list(
    "bootstrap_tune_model" = lgbm_tune,
    "train_lgbm" = train_lgbm,
    "test_lgbm" = test_lgbm,
    "train_metrics" = train_metrics,
    "test_metrics" = test_metrics
  )
  return(output)
}



#_______________________________________________________________________________
#                      Feature Selection Function
#_______________________________________________________________________________

mRMR_selection <- function(input, feats){
  cat("Running mRMR feature selection .. \n")
  trimmed_input <- input %>% 
    dplyr::select(donor_id, PD, cohort, paired, matches(feats))
  cat("Feature N input:", length(feats) ,"Features selected : ", ncol(trimmed_input)-4, "\n")
  return(trimmed_input)
}


#_______________________________________________________________________________
#                      ML Cohort Cross Testing
#_______________________________________________________________________________

ml_summary <- function(obj_A,
                       obj_B,
                       label_A,
                       label_B,
                       analysis,
                       featSelection = NULL,
                       model_type = "lasso") {

  
  # TROUBLE
  obj_A = obj_PD
  obj_B = obj_Controls
  label_A = "PD"
  label_B = "Controls"
  analysis = "LOSO"
  model_type = "lasso"
  featSelection =  NULL
  
  # ml_input_A <- prep_ml_input(obj_A)
  # ml_input_B <- prep_ml_input(obj_B)
  
  if (model_type == "lasso"){
    A_vs_B <- lasso_model(train = ml_input_A, test = ml_input_B)
  } else if (model_type == "randomforest"){
    A_vs_B <- rf_model(train = ml_input_A, test = ml_input_B)
  } else if (model_type == "xgboost"){
    A_vs_B <- xgb_model(train = ml_input_A, test = ml_input_B)
  } else if (model_type == "lightgbm"){
    A_vs_B <- lgbm_model(train = ml_input_A, test = ml_input_B)
  } 
  
  summary_A <-
    bind_rows(
      A_vs_B$train_metrics %>%
        mutate(train = label_A, test = label_A, model = "Bootstrap validation"),
      A_vs_B$test_metrics %>%
        mutate(train = label_A, test = label_B, model = "Prediction")) %>% 
    mutate(analysis = analysis)

  return(summary_A)
}



#_______________________________________________________________________________
#                      LASSO Cohort Cross Testing TRIMMED FUNCTION
#_______________________________________________________________________________

lasso_cohort_summary_trim <- function(ml_input_A, ml_input_B, 
                                 label_A, label_B, analysis, featSelection = NULL){

  A_vs_B <- lasso_model(train = ml_input_A,
                        test = ml_input_B)
  summary_A <-
    bind_rows(
      A_vs_B$train_metrics %>%
        mutate(train = label_A, test = label_A, model = "10-fold CV"),
      A_vs_B$test_metrics %>%
        mutate(train = label_A, test = label_B, model = "Prediction")) %>% 
    mutate(analysis = analysis)
  
  return(summary_A)
}


#_______________________________________________________________________________
#                      Study to Study Transfer
#_______________________________________________________________________________

ml_s2s <- function(Shanghai,
                   TBC,
                   Rush,
                   Bonn,
                   model_type = "lasso",
                   featSelection = NULL,
                   nfeats = 100, 
                   data_type = data_type) {
  
  # ## TROUBLE
  # Shanghai = obj_Shanghai
  # TBC = obj_TBC
  # Rush = obj_Rush
  # Bonn = obj_Bonn
  # model_type = "lasso"
  # data_type = "KOs.slim"
  # nfeats = 100
  # featSelection = T
  
  #____________________________________________________________________________
  #                        SWAP FOR RAREFIED DATA ------
  # load("data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs.RData")
  load("data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs_Rarefied.RData")
  #____________________________________________________________________________
  
  feats <- vector("list", length = 4)

  if(!is.null(featSelection)){
    featselectStart <- Sys.time()
    featSelect_Shanghai <- feature_selection_auroc(df = aucs, n_features = nfeats,
                                       subgroup = "Shanghai", data_type = data_type)[1:nfeats]
    featSelect_TBC <- feature_selection_auroc(df = aucs, n_features = nfeats,
                                        subgroup = "TBC", data_type = data_type)[1:nfeats]
    featSelect_Rush <- feature_selection_auroc(df = aucs, n_features = nfeats,
                                         subgroup = "Rush", data_type = data_type)[1:nfeats]
    featSelect_Bonn <- feature_selection_auroc(df = aucs, n_features = nfeats,
                                         subgroup = "Bonn", data_type = data_type)[1:nfeats]
    # MRMR methods
      # feature_selection(df = prep_mRMR_input(Bonn), n_features = nfeats, n_algorithm = 1)
    
    feats[[1]] <- featSelect_Shanghai
    feats[[2]] <- featSelect_TBC
    feats[[3]] <- featSelect_Rush
    feats[[4]] <- featSelect_Bonn
    names(feats) <- c("featSelect_Shanghai", "featSelect_TBC", "featSelect_Rush", "featSelect_Bonn")
    
    featselectEnd <- Sys.time()
    cat("\nfeature selection complete : ", featselectEnd - featselectStart,
        attr(featselectEnd - featselectStart, "units"), "\n")
  }
  rm(aucs)
  
  modelsStart <- Sys.time()
  # Study to Study Transfer
  s2s_Shanghai_x_TBC <-
    ml_summary(
      obj_A = Shanghai,
      obj_B = TBC,
      label_A = "Shanghai",
      label_B = "TBC",
      analysis = "S2S",
      model_type = model_type,
      featSelection = feats[["featSelect_Shanghai"]]
    )
  s2s_Shanghai_x_Rush <-
    ml_summary(
      obj_A = Shanghai,
      obj_B = Rush,
      label_A = "Shanghai",
      label_B = "Rush",
      analysis = "S2S",
      model_type = model_type,
      featSelection = feats[["featSelect_Shanghai"]]
    )
  s2s_Shanghai_x_Bonn <-
    ml_summary(
      obj_A = Shanghai,
      obj_B = Bonn,
      label_A = "Shanghai",
      label_B = "Bonn",
      analysis = "S2S",
      model_type = model_type,
      featSelection = feats[["featSelect_Shanghai"]]
      )
  
  s2s_TBC_x_Shanghai <-
    ml_summary(
      obj_A = TBC,
      obj_B = Shanghai,
      label_A = "TBC",
      label_B = "Shanghai",
      analysis = "S2S",
      model_type = model_type,
      featSelection = feats[["featSelect_TBC"]]
    )
  s2s_TBC_x_Rush <-
    ml_summary(
      obj_A = TBC,
      obj_B = Rush,
      label_A = "TBC",
      label_B = "Rush",
      analysis = "S2S",
      model_type = model_type,
      featSelection = feats[["featSelect_TBC"]]
    )
  s2s_TBC_x_Bonn <-
    ml_summary(
      obj_A = TBC,
      obj_B = Rush,
      label_A = "TBC",
      label_B = "Bonn",
      analysis = "S2S",
      model_type = model_type,
      featSelection = feats[["featSelect_TBC"]]
    )
  
  s2s_Rush_x_Shanghai <-
    ml_summary(
      obj_A = Rush,
      obj_B = Shanghai,
      label_A = "Rush",
      label_B = "Shanghai",
      analysis = "S2S",
      model_type = model_type,
      featSelection =  feats[["featSelect_Rush"]]
    )
  s2s_Rush_x_TBC <-
    ml_summary(
      obj_A = Rush,
      obj_B = TBC,
      label_A = "Rush",
      label_B = "TBC",
      analysis = "S2S",
      model_type = model_type,
      featSelection =  feats[["featSelect_Rush"]]
    )
  s2s_Rush_x_Bonn <-
    ml_summary(
      obj_A = Rush,
      obj_B = TBC,
      label_A = "Rush",
      label_B = "Bonn",
      analysis = "S2S",
      model_type = model_type,
      featSelection =  feats[["featSelect_Rush"]]
    )
  
  s2s_Bonn_x_TBC <-
    ml_summary(
      obj_A = Bonn,
      obj_B = TBC,
      label_A = "Bonn",
      label_B = "TBC",
      analysis = "S2S",
      model_type = model_type,
      featSelection = feats[["featSelect_Bonn"]]
    )
  s2s_Bonn_x_Rush <-
    ml_summary(
      obj_A = Bonn,
      obj_B = Rush,
      label_A = "Bonn",
      label_B = "Rush",
      analysis = "S2S",
      model_type = model_type,
      featSelection = feats[["featSelect_Bonn"]]
    )
  s2s_Bonn_x_Shanghai <-
    ml_summary(
      obj_A = Bonn,
      obj_B = Rush,
      label_A = "Bonn",
      label_B = "Shanghai",
      analysis = "S2S",
      model_type = model_type,
      featSelection = feats[["featSelect_Bonn"]]
    )
  
  modelsEnd <- Sys.time()
  cat("\n", model_type, " Study to Study Transfer Analysis complete : ", 
      modelsEnd - modelsStart, attr(modelsEnd - modelsStart, "units"), "\n")
  
  s2s_summary <- bind_rows(
    s2s_Shanghai_x_TBC, s2s_Shanghai_x_Rush, s2s_Shanghai_x_Bonn,
    s2s_TBC_x_Shanghai, s2s_TBC_x_Rush, s2s_TBC_x_Bonn,
    s2s_Rush_x_Shanghai, s2s_Rush_x_TBC, s2s_Rush_x_Bonn,
    s2s_Bonn_x_TBC, s2s_Bonn_x_Rush, s2s_Bonn_x_Shanghai
  )
  return(s2s_summary)
}

#_______________________________________________________________________________
#                    LOSO - Leave One Study Out
#_______________________________________________________________________________


ml_loso <- function(Shanghai,
                   noShanghai,
                   TBC,
                   noTBC,
                   Rush,
                   noRush,
                   Bonn,
                   noBonn,
                   model_type = "lasso",
                   featSelection = NULL,
                   nfeats = 100, 
                   data_type){
  
  
  # # TROUBLE
  # Shanghai = obj_Shanghai
  # noShanghai = obj_noShanghai
  # TBC = obj_TBC
  # noTBC = obj_noTBC
  # Rush = obj_Rush
  # noRush = obj_noRush
  # Bonn = obj_Bonn
  # noBonn = obj_noBonn
  # model_type = "lasso"
  # featSelection = TRUE
  # nfeats = 100
  
  #____________________________________________________________________________
  #                        SWAP FOR RAREFIED DATA ------
  # load("data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs_LODO.RData")
  load("data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs_LODO_Rarefied.RData")
  #____________________________________________________________________________
  
  feats <- vector("list", length = 4)
  
  if(!is.null(featSelection)){
    
    featselectStart <- Sys.time()
    featSelect_noShanghai <- feature_selection_auroc(df = aucs, n_features = nfeats,
                                                      subgroup = "no_Shanghai", data_type = data_type)[1:nfeats]
    featSelect_noTBC <- feature_selection_auroc(df = aucs, n_features = nfeats,
                                                 subgroup = "no_TBC", data_type = data_type)[1:nfeats]
    featSelect_noRush <- feature_selection_auroc(df = aucs, n_features = nfeats,
                                               subgroup = "no_Rush", data_type = data_type)[1:nfeats]
    featSelect_noBonn <- feature_selection_auroc(df = aucs, n_features = nfeats,
                                               subgroup = "no_Bonn", data_type = data_type)[1:nfeats]
    # mRMR_noShanghai <-
    #   feature_selection(df = prep_mRMR_input(noShanghai), n_features = nfeats, n_algorithm = 1)

    feats[[1]] <- featSelect_noShanghai
    feats[[2]] <- featSelect_noTBC
    feats[[3]] <- featSelect_noRush
    feats[[4]] <- featSelect_noBonn
    names(feats) <- c("featSelect_noShanghai", "featSelect_noTBC", 
                      "featSelect_noRush", "featSelect_noBonn")
    
    featselectEnd <- Sys.time()
    cat("\nfeature selection complete : ", 
        featselectEnd - featselectStart,
        attr(featselectEnd - featselectStart, "units"), "\n")
  }
  print(summary(feats))
  rm(aucs)
  modelsStart <- Sys.time()
  
  loso_shanghai <-
    ml_summary(
      obj_A = noShanghai,
      obj_B = Shanghai,
      label_A = "LOSO_Shanghai",
      label_B = "Shanghai",
      analysis = "LOSO",
      model_type = model_type,
      featSelection =  feats[["featSelect_noShanghai"]]
    )
  loso_TBC <-
    ml_summary(
      obj_A = noTBC,
      obj_B = TBC,
      label_A = "LOSO_TBC",
      label_B = "TBC",
      analysis = "LOSO",
      model_type = model_type,
      featSelection =  feats[["featSelect_noTBC"]]
    )
  loso_Rush <-
    ml_summary(
      obj_A = noRush,
      obj_B = Rush,
      label_A = "LOSO_Rush",
      label_B = "Rush",
      analysis = "LOSO",
      model_type = model_type,
      featSelection =  feats[["featSelect_noRush"]]
    )
  loso_Bonn <-
    ml_summary(
      obj_A = noBonn,
      obj_B = Bonn,
      label_A = "LOSO_Bonn",
      label_B = "Bonn",
      analysis = "LOSO",
      model_type = model_type,
      featSelection =  feats[["featSelect_noBonn"]]
    )
  
  modelsEnd <- Sys.time()
  cat("\n ", model_type, " LOSO Analysis complete : ", 
      modelsEnd - modelsStart, attr(modelsEnd - modelsStart, "units"), "\n")
  
  loso_summary <- bind_rows(loso_shanghai, loso_TBC, loso_Rush, loso_Bonn)
  return(loso_summary)
}


#_______________________________________________________________________________
#                 LASSO Cohort x Donor Group Cross Testing
#_______________________________________________________________________________

lasso_cohort_x_group_summary <- function(obj_tbc_all, obj_rush_all, featSelection = NULL){
  
  ml_input_tbc <- prep_ml_input(obj_tbc_all)
  ml_input_rush <- prep_ml_input(obj_rush_all)
  
  if(!is.null(featSelection)){
    cat("Selecting mRMR features .. \n")
    ml_input_tbc <- ml_input_tbc %>% 
      dplyr::select(donor_id, PD, cohort, paired, matches(featSelection))
    ml_input_rush <- ml_input_rush %>% 
      dplyr::select(donor_id, PD, cohort, paired, matches(featSelection))
  }
  
  ml_input_tbc_house <- obj_tbc_all %>% 
    subset_samples(paired != "No") %>% prep_ml_input()
  ml_input_tbc_healthy <- obj_tbc_all %>% 
    subset_samples(donor_group != "HC") %>% prep_ml_input() 
  ml_input_rush_house <- obj_rush_all %>% 
    subset_samples(paired != "No") %>% prep_ml_input()
  ml_input_rush_healthy <- obj_rush_all %>% 
    subset_samples(donor_group != "HC") %>% prep_ml_input()
  
  
  rush_PC_vs_rush_HC <- lasso_model(train = ml_input_rush_healthy,
                                    test = ml_input_rush_house)
  rush_PC_vs_tbc_PC <- lasso_model(train = ml_input_rush_healthy,
                                   test = ml_input_tbc_healthy)
  rush_PC_vs_tbc_HC <- lasso_model(train = ml_input_rush_healthy,
                                   test = ml_input_tbc_house)
  
  rush_HC_vs_rush_PC <- lasso_model(train = ml_input_rush_house,
                                    test = ml_input_rush_healthy)
  rush_HC_vs_tbc_PC <- lasso_model(train = ml_input_rush_house,
                                   test = ml_input_tbc_healthy)
  rush_HC_vs_tbc_HC <- lasso_model(train = ml_input_rush_house,
                                   test = ml_input_tbc_house)
  
  tbc_PC_vs_tbc_HC <- lasso_model(train = ml_input_tbc_healthy,
                                  test = ml_input_tbc_house)
  tbc_PC_vs_rush_PC <- lasso_model(train = ml_input_tbc_healthy,
                                   test = ml_input_rush_healthy)
  tbc_PC_vs_rush_HC <- lasso_model(train = ml_input_tbc_healthy,
                                   test = ml_input_rush_house)
  
  tbc_HC_vs_tbc_PC <- lasso_model(train = ml_input_tbc_house,
                                  test = ml_input_tbc_healthy)
  tbc_HC_vs_rush_PC <- lasso_model(train = ml_input_tbc_house,
                                   test = ml_input_rush_healthy)
  tbc_HC_vs_rush_HC <- lasso_model(train = ml_input_tbc_house,
                                   test = ml_input_rush_house)
  
  
  summaryB <-
    bind_rows(
      # RUSH PC vs all
      rush_PC_vs_rush_HC$train_metrics %>%
        mutate(train = "Rush_PC", test = "Rush_PC"),
      rush_PC_vs_rush_HC$test_metrics %>%
        mutate(train = "Rush_PC", test = "Rush_HC"),
      rush_PC_vs_tbc_PC$test_metrics %>%
        mutate(train = "Rush_PC", test = "TBC_PC"),
      rush_PC_vs_tbc_HC$test_metrics %>%
        mutate(train = "Rush_PC", test = "TBC_HC"),
      # RUSH HC vs all
      rush_HC_vs_rush_PC$train_metrics %>%
        mutate(train = "Rush_HC", test = "Rush_HC"),
      rush_HC_vs_rush_PC$test_metrics %>%
        mutate(train = "Rush_HC", test = "Rush_PC"),
      rush_HC_vs_tbc_PC$test_metrics %>%
        mutate(train = "Rush_HC", test = "TBC_PC"),
      rush_HC_vs_tbc_HC$test_metrics %>%
        mutate(train = "Rush_HC", test = "TBC_HC"),
      # TBC PC vs all
      tbc_PC_vs_tbc_HC$train_metrics %>%
        mutate(train = "TBC_PC", test = "TBC_PC"),
      tbc_PC_vs_tbc_HC$test_metrics %>%
        mutate(train = "TBC_PC", test = "TBC_HC"),
      tbc_PC_vs_rush_PC$test_metrics %>%
        mutate(train = "TBC_PC", test = "Rush_PC"),
      tbc_PC_vs_rush_HC$test_metrics %>%
        mutate(train = "TBC_PC", test = "Rush_HC"),
      # TBC PC vs all
      tbc_HC_vs_tbc_PC$train_metrics %>%
        mutate(train = "TBC_HC", test = "TBC_HC"),
      tbc_HC_vs_tbc_PC$test_metrics %>%
        mutate(train = "TBC_HC", test = "TBC_PC"),
      tbc_HC_vs_rush_PC$test_metrics %>%
        mutate(train = "TBC_HC", test = "Rush_PC"),
      tbc_HC_vs_rush_HC$test_metrics %>%
        mutate(train = "TBC_HC", test = "Rush_HC")
    ) %>% 
    mutate(train_cohort = substr(train, 1, nchar(train)-3)) %>% 
    mutate(test_cohort = substr(test, 1, nchar(test)-3)) %>% 
    mutate(train_group = substr(train, nchar(train)-1, nchar(train))) %>% 
    mutate(test_group = substr(test, nchar(test)-1, nchar(test)))
  
  return(summaryB)
}



































# 
# 
# lasso_model <- function(train, test){
#   
#   # TROUBLE
#   ml_input_tbc <- prep_ml_input(obj_tbc_all)
#   ml_input_rush <- prep_ml_input(obj_rush_all)
#   train = ml_input_tbc
#   test = ml_input_rush
# 
#   
# 
#   combined <- bind_rows(train, test)
#   ind <- list(
#     analysis = seq(nrow(train)),
#     assessment = nrow(train) + seq(nrow(test)))
#   splits <- make_splits(ind, combined)
#   ml_train <- training(splits)
#   ml_test <- testing(splits)
#   
#   set.seed(42)
#   # Generate 10-fold CV model with 5 repetitions, stratified by PD status 
#   cv_splits <- rsample::vfold_cv(ml_train, v = 10, repeats = 10, strata = PD)
#   # Define tunable Lasso model
#   tune_spec <- logistic_reg(penalty = tune(), mixture = 1) %>%
#     set_engine("glmnet")
#   
#   ml_recipe <- recipe(PD ~ ., data = ml_train) %>% 
#     update_role(donor_id, new_role = "donor_id") %>% 
#     update_role(cohort, new_role = "cohort") %>% 
#     update_role(paired, new_role = "paired") %>%
#     step_zv(all_numeric(), -all_outcomes()) %>% 
#     step_normalize(all_numeric(), -all_outcomes())
#   
#   wf <- workflow() %>% 
#     add_recipe(ml_recipe) %>% 
#     add_model(tune_spec)
#   
#   #_____________________________________________________
#   #                TUNE LASSO MODEL  
#   #_____________________________________________________
#   # Create a grid of penalty values to test
#   lambda_grid <- grid_regular(penalty(), levels = 50)
#   ctrl <- control_grid(save_pred = TRUE, verbose = TRUE)
#   doParallel::registerDoParallel()
#   
#   set.seed(42)
#   lasso_grid <-
#     tune_grid(wf,
#               resamples = cv_splits,
#               grid = lambda_grid,
#               control = ctrl)
#   
#   # select optimal penalty by filtering largest rocauc
#   best_aucroc <- select_best(lasso_grid, "roc_auc")
#   # visualize model metrics of grid 
#   model_performance <- 
#     lasso_grid %>% 
#     collect_metrics() %>% 
#     ggplot(aes(penalty, mean, color = .metric)) +
#     geom_errorbar(aes(ymin = mean - std_err,
#                       ymax = mean + std_err),
#                   alpha = 0.5) +
#     geom_line(size = 1.25, show.legend = F) +
#     facet_wrap(~.metric, scales = "free", nrow = 2) +
#     theme_bw() + 
#     scale_x_log10() +
#     scale_color_viridis_d(option = "cividis", begin = .9, end = 0) +
#     theme(legend.position = "none")
#   print(model_performance)
#   
#   train_lasso2 <-
#     wf %>%
#     finalize_workflow(best_aucroc) %>% 
#     fit_resamples(cv_splits, control = control_resamples(save_pred = TRUE))
# 
#   train_lasso2 %>% 
#     collect_metrics() %>% 
#     filter(penalty == best_aucroc$penalty)
#   
#   
#   # Filtering ROCAUC valued from optimal 5x10fold CV model
#   # train_metrics <- 
#   #   lasso_grid %>% 
#   #   collect_metrics() %>% 
#   #   filter(penalty == best_aucroc$penalty)
#   
#   # Finalize trained workflow - use on hold out test sets
#   train_lasso <-
#     wf %>%
#     finalize_workflow(best_aucroc) %>%
#     fit(ml_train)
#   
#   # Predictions on test data
#   test_lasso <- 
#     last_fit(train_lasso2, split = splits) 
#   test_metrics <- test_lasso %>%
#     collect_metrics()
#   
#   output <- list("train_lasso" = train_lasso, 
#                  "test_lasso" = test_lasso, 
#                  "train_metrics" = train_metrics,
#                  "test_metrics" = test_metrics)
#   return(output)
# }
#' 
#' 
#' 
#' #_______________________________________________________________________________--------
#' ###                       RIDGE, LASSO, ENET LOGISTIC REGRESSION                     ###
#' #_______________________________________________________________________________--------
#' 
#' 
#' ridge.lasso.enet.regression.model <- function(obj, comparison = "PDvPC", model.type){
#'   
#'   #' Function to train a Ridge, lASSO, or ElasticNet Regression Model
#'   #' Input: Phyloseq Obj, comparison of interest, and model type
#'   #' Returns: a list including 
#'   #' 1) The fitted Model 
#'   #' 2) A dataframe of values with the optimal parameters (for ROC plot)
#'   #' 3) MLeval AUC-ROC values 
#'   #' 4) Other MLeval analysis parameters
#'   
#'   # Initalize variables
#'   model <- NULL
#'   model.input <- NULL
#'   output.list <- vector(mode="list", length=4)
#'   names(output.list) <- c("fitted_model", "optimal.df", "AUCROC", "MLevaldata")
#'   
#'   if (model.type == "ridge"){
#'     tune.grid = expand.grid(alpha = 0, lambda=seq(0, 1, 0.1))
#'   } else if (model.type == "lasso"){
#'     tune.grid = expand.grid(alpha = 1, lambda=seq(0, 1, 0.1))
#'   } else if (model.type == "enet"){
#'     tune.grid = expand.grid(alpha = 0.5, lambda=seq(0, 1, 0.1))
#'   }
#'   
#'   # Select samples and prep abundance DF 
#'   if (comparison == "PDvPC"){
#'     dat_pdpc = subset_samples(obj, donor_group !="HC")
#'     d <- dat_pdpc %>%
#'       microbiome::transform("compositional") %>% 
#'       microbiome::abundances() %>% 
#'       t() %>% 
#'       as.data.frame() 
#'     d <- asin(sqrt(d))
#'     pdPC.df <- group_col_from_ids(d, id= rownames(d))
#'     rownames(pdPC.df) <- rownames(d)
#'     pdPC.df$group <- factor(pdPC.df$group, levels = c("PC", "PD"))
#'     model.input <- pdPC.df
#'     
#'   } else if (comparison == "PDvHC"){
#'     dat_pdhc = subset_samples(obj, Paired !="No")
#'     d <- dat_pdhc %>%
#'       microbiome::transform("compositional") %>% 
#'       microbiome::abundances() %>% 
#'       t() %>% 
#'       as.data.frame() 
#'     d <- asin(sqrt(d))
#'     pdHC.df <- group_col_from_ids(d, id= rownames(d))
#'     rownames(pdHC.df) <- rownames(d)
#'     pdHC.df$group <- factor(pdHC.df$group, levels = c("HC", "PD"))
#'     model.input <- pdHC.df
#'   }
#'   
#'   # Model Parameters
#'   numbers <- 10
#'   repeats <- 5  
#'   set.seed(42)
#'   seed <- 42
#'   rcvSeeds <- setSeeds(method = "repeatedcv", numbers = numbers, repeats = repeats, seed = seed)
#'   
#'   mytrainControl <- 
#'     trainControl(method='repeatedcv',
#'                  number=numbers, 
#'                  repeats=repeats,
#'                  search='grid',
#'                  seeds = rcvSeeds,
#'                  savePredictions = TRUE, 
#'                  classProbs = TRUE, 
#'                  verboseIter = TRUE)
#'   
#'   # Run 10-fold CV Model with 5 Repitions
#'   model <-train(group ~.,
#'                 data=model.input, 
#'                 method='glmnet', 
#'                 metric='Accuracy', 
#'                 tuneGrid=tune.grid, 
#'                 trControl=mytrainControl,
#'                 verboseIter = T)
#'   cat("Model Summary")
#'   print(model)
#'   cat("\n\n")
#'   print(plot(model))
#'   
#'   # Select model with optimal Lambda
#'   selectedIndices <- model$pred$lambda == model$bestTune$lambda
#'   df.output <- model$pred[selectedIndices, ]
#'   
#'   # MLeval
#'   mleval <- evalm(model)
#'   mleval$roc
#'   output.list$fitted_model <- model
#'   output.list$optimal.df <- df.output
#'   output.list$AUCROC <- mleval$stdres$`Group 1`["AUC-ROC", "Score"]
#'   output.list$MLevaldata <- mleval$stdres
#'   
#'   
#'   return(output.list)
#'   
#' }
#' 
#' #_______________________________________________________________________________--------
#' ###                              RANDOM FOREST MODELS                                ###
#' #_______________________________________________________________________________--------
#' 
#' 
#' 
#' random.forest.model <- function(obj, comparison = "PDvPC"){
#'   
#'   #' Function to train a Random Forrest Model
#'   #' Input: Phyloseq Obj and comparison of interest 
#'   #' Returns: a list including 
#'   #' 1) The fitted Model 
#'   #' 2) A dataframe of values with the optimal parameters (for ROC plot)
#'   #' 3) MLeval AUC-ROC values 
#'   #' 4) Other MLeval analysis parameters
#'   
#'   intervalStart <- Sys.time()
#'   
#'   # Initalize variables
#'   model <- NULL
#'   output.list <- vector(mode="list", length=4)
#'   names(output.list) <- c("fitted_model", "optimal.df", "AUCROC", "MLevaldata")
#'   
#'   # Select samples and prep abundance DF 
#'   if (comparison == "PDvPC"){
#'     dat_pdpc = subset_samples(obj, donor_group !="HC")
#'     d <- dat_pdpc %>%
#'       microbiome::transform("compositional") %>% 
#'       microbiome::abundances() %>% 
#'       t() %>% 
#'       as.data.frame()
#'     d <- asin(sqrt(d))
#'     pdPC.df <- group_col_from_ids(d, id= rownames(d))
#'     rownames(pdPC.df) <- rownames(d)
#'     pdPC.df$group <- factor(pdPC.df$group, levels = c("PC", "PD"))
#'     model.input <- pdPC.df
#'     
#'   } else if (comparison == "PDvHC"){
#'     dat_pdhc = subset_samples(obj, Paired !="No")
#'     d <- dat_pdhc %>%
#'       microbiome::transform("compositional") %>% 
#'       microbiome::abundances() %>% 
#'       t() %>% 
#'       as.data.frame() 
#'     d <- asin(sqrt(d))
#'     pdHC.df <- group_col_from_ids(d, id= rownames(d))
#'     rownames(pdHC.df) <- rownames(d)
#'     pdHC.df$group <- factor(pdHC.df$group, levels = c("HC", "PD"))
#'     model.input <- pdHC.df
#'   }
#'   
#'   # Model Parameters
#'   numbers <- 10
#'   repeats <- 5
#'   tune.grid = expand.grid(.mtry =   seq(1, 2 * as.integer(sqrt(ncol(model.input) - 1)), by=2)) 
#'   set.seed(42)
#'   seed <- 42
#'   # Repeated cross validation
#'   rcvSeeds <- setSeeds(method = "repeatedcv", numbers = numbers, repeats = repeats, 
#'                        tunes = length(tune.grid$.mtry), seed = seed)
#'   # c('B + 1' = length(rcvSeeds), M = length(rcvSeeds[[1]]))
#'   # rcvSeeds[c(1, length(rcvSeeds))]
#'   mytrainControl <- 
#'     trainControl(method='repeatedcv',
#'                  number=numbers, 
#'                  repeats=repeats,
#'                  search='grid',
#'                  savePredictions = TRUE, 
#'                  classProbs = TRUE, 
#'                  verboseIter = TRUE,
#'                  allowParallel = TRUE,
#'                  seeds = rcvSeeds)
#'   
#'   # Run model using multiple threads 
#'   cluster <- parallel::makeCluster(detectCores() - 1, setup_strategy = "sequential")
#'   registerDoParallel(cluster)
#'   
#'   # Run 10-fold CV Model with 5 Repitions
#'   system.time(
#'     model <-train(group ~.,
#'                   data=model.input, 
#'                   method='rf', 
#'                   metric='Accuracy', 
#'                   tuneGrid=tune.grid, 
#'                   ntree = 1000,
#'                   trControl=mytrainControl,
#'                   importance = TRUE,
#'                   verboseIter = T)
#'   )
#'   
#'   stopCluster(cluster)
#'   unregister()
#'   
#'   cat("Model Summary")
#'   print(model)
#'   cat("\n\n")
#'   print(plot(model))
#'   
#'   # Select model with optimal __mtry__ 
#'   selectedIndices <- model$pred$mtry == model$bestTune$mtry
#'   df.output <- model$pred[selectedIndices, ]
#'   
#'   # MLeval
#'   mleval <- evalm(model)
#'   mleval$roc
#'   output.list$fitted_model <- model
#'   output.list$optimal.df <- df.output
#'   output.list$AUCROC <- mleval$stdres$`Group 1`["AUC-ROC", "Score"]
#'   output.list$MLevaldata <- mleval$stdres
#'   
#'   intervalEnd <- Sys.time()
#'   cat("Random Forest analysis took",
#'       intervalEnd - intervalStart,attr(intervalEnd - intervalStart,"units"))
#'   cat("\n\n")
#'   return(output.list)
#'   
#' }
#' 
#' 
#' 
#' 
#' #____________________________________________________________________________
#' ###                           XBOOST MODELS                           ###
#' #____________________________________________________________________________
#' 
#' 
#' xgboost.model <- function(obj, comparison = "PDvPC"){
#'   
#'   
#'   # ########    DELETE LATER  4 - TESTING ######## 
#'   # obj = dat
#'   # comparison = "PDvPC"
#'   # ########  ########  ########  ########  ######## 
#'   # 
#'   
#'   #' Function to train a Random Forrest Model
#'   #' Input: Phyloseq Obj and comparison of interest 
#'   #' Returns: a list including 
#'   #' 1) A dataframe of values with the optimal parameters (for ROC plot)
#'   #' 2) MLeval individual model ROC plot
#'   #' 3) MLeval AUC-ROC values 
#'   #' 4) Other MLeval analysis parameters
#'   
#'   intervalStart <- Sys.time()
#'   
#'   # Initalize variables
#'   model <- NULL
#'   output.list <- vector(mode="list", length=3)
#'   names(output.list) <- c("optimal.df", "AUCROC", "MLevaldata")
#'   
#'   set.seed(42)
#'   
#'   # Select samples and prep abundance DF 
#'   if (comparison == "PDvPC"){
#'     dat_pdpc = subset_samples(obj, donor_group !="HC")
#'     d <- dat_pdpc %>%
#'       microbiome::transform("compositional") %>% 
#'       microbiome::abundances() %>% 
#'       t() %>% 
#'       as.data.frame() 
#'     d <- asin(sqrt(d))
#'     pdPC.df <- group_col_from_ids(d, id= rownames(d))
#'     rownames(pdPC.df) <- rownames(d)
#'     pdPC.df$group <- factor(pdPC.df$group, levels = c("PC", "PD"))
#'     model.input <- pdPC.df
#'     
#'   } else if (comparison == "PDvHC"){
#'     dat_pdhc = subset_samples(obj, Paired !="No")
#'     d <- dat_pdhc %>%
#'       microbiome::transform("compositional") %>% 
#'       microbiome::abundances() %>% 
#'       t() %>% 
#'       as.data.frame() 
#'     d <- asin(sqrt(d))
#'     pdHC.df <- group_col_from_ids(d, id= rownames(d))
#'     rownames(pdHC.df) <- rownames(d)
#'     pdHC.df$group <- factor(pdHC.df$group, levels = c("HC", "PD"))
#'     model.input <- pdHC.df
#'   }
#'   
#'   # Set-up Parallel Processing Threads
#'   # omp_set_num_threads(3)
#'   # intervalStart <- Sys.time()
#'   
#'   #----------------------------------- 
#'   # Hyperparameter Tuning 
#'   #Informed by: https://www.kaggle.com/pelkoja/visual-xgboost-tuning-with-caret/report
#'   
#'   #_____________________________________________________ 
#'   # 1) nrounds & eta 
#'   #_____________________________________________________ 
#'   
#'   tune_grid <- expand.grid(
#'     nrounds = seq(from = 200, to = 1000, by = 50),
#'     eta = c(0.025, 0.05, 0.1, 0.3),
#'     max_depth = c(2, 3, 4, 5, 6),
#'     gamma = 0,
#'     colsample_bytree = 1,
#'     min_child_weight = 1,
#'     subsample = 1
#'   )
#'   
#'   tune_control <-
#'     trainControl(method='repeatedcv',
#'                  number=5,
#'                  repeats=3,
#'                  search='grid',
#'                  # savePredictions = TRUE,
#'                  # classProbs = TRUE,
#'                  # allowParallel = TRUE,
#'                  verboseIter = TRUE)
#'   
#'   xgb_tune <- caret::train(
#'     group ~.,
#'     data=model.input,
#'     trControl = tune_control,
#'     tuneGrid = tune_grid,
#'     method = "xgbTree",
#'     verbose = TRUE
#'   )
#'   
#'   cat("Tuning Step 1: Maximum Depth, learning rate, and nrounds baseline \n COMPLETE \n")
#'   print(tuneplot(xgb_tune))
#'   print(xgb_tune$bestTune)
#'   
#'   #_____________________________________________________ 
#'   # 2) Maximum Depth and Minimum Child Weight
#'   #_____________________________________________________ 
#'   
#'   if (xgb_tune$bestTune$max_depth == 2) {
#'     mxdpth <- 
#'       c(xgb_tune$bestTune$max_depth:4)
#'   }  else {
#'     mxdpth <-  
#'       c((xgb_tune$bestTune$max_depth - 1):(xgb_tune$bestTune$max_depth + 1))
#'   }
#'   
#'   tune_grid2 <- expand.grid(
#'     nrounds = seq(from = 50, to = 1000, by = 50),
#'     eta = xgb_tune$bestTune$eta,
#'     max_depth = mxdpth,
#'     gamma = 0,
#'     colsample_bytree = 1,
#'     min_child_weight = c(1, 2, 3),
#'     subsample = 1
#'   )
#'   
#'   xgb_tune2 <- caret::train(
#'     group ~.,
#'     data=model.input,
#'     trControl = tune_control,
#'     tuneGrid = tune_grid2,
#'     method = "xgbTree",
#'     verbose = TRUE
#'   )
#'   cat("Tuning Step 2: Maximum Depth and Minimum Child Weight \n COMPLETE \n")
#'   print(tuneplot(xgb_tune2))
#'   print(xgb_tune2$bestTune)
#'   
#'   #_____________________________________________________ 
#'   # 3)  Column and Row Sampling
#'   #_____________________________________________________ 
#'   
#'   tune_grid3 <- expand.grid(
#'     nrounds = seq(from = 50, to = 1000, by = 50),
#'     eta = xgb_tune$bestTune$eta,
#'     max_depth = xgb_tune2$bestTune$max_depth,
#'     gamma = 0,
#'     colsample_bytree = c(0.4, 0.6, 0.8, 1.0),
#'     min_child_weight = xgb_tune2$bestTune$min_child_weight,
#'     subsample = c(0.5, 0.75, 1.0)
#'   )
#'   
#'   xgb_tune3 <- caret::train(
#'     group ~.,
#'     data=model.input,
#'     trControl = tune_control,
#'     tuneGrid = tune_grid3,
#'     method = "xgbTree",
#'     verbose = TRUE
#'   )
#'   
#'   cat("Tuning Step 3: Column and Row Sampling \n COMPLETE \n")
#'   print(tuneplot(xgb_tune3, probs = .95))
#'   print(xgb_tune3$bestTune)
#'   
#'   #_____________________________________________________ 
#'   # 4)  Gamma
#'   #_____________________________________________________ 
#'   
#'   tune_grid4 <- expand.grid(
#'     nrounds = seq(from = 50, to = 1000, by = 50),
#'     eta = xgb_tune$bestTune$eta,
#'     max_depth = xgb_tune2$bestTune$max_depth,
#'     gamma = c(0, 0.05, 0.1, 0.5, 0.7, 0.9, 1.0),
#'     colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
#'     min_child_weight = xgb_tune2$bestTune$min_child_weight,
#'     subsample = xgb_tune3$bestTune$subsample
#'   )
#'   
#'   xgb_tune4 <- caret::train(
#'     group ~.,
#'     data=model.input,
#'     trControl = tune_control,
#'     tuneGrid = tune_grid4,
#'     method = "xgbTree",
#'     verbose = TRUE
#'   )
#'   
#'   cat("Tuning Step 4: Gamma \n COMPLETE \n")
#'   print(tuneplot(xgb_tune4))
#'   print(xgb_tune4$bestTune)
#'   
#'   #_____________________________________________________ 
#'   # 5)  Reducing the Learning Rate
#'   #_____________________________________________________ 
#'   
#'   tune_grid5 <- expand.grid(
#'     nrounds = seq(from = 100, to = 1000, by = 100),
#'     eta = c(0.01, 0.015, 0.025, 0.05, 0.1),
#'     max_depth = xgb_tune2$bestTune$max_depth,
#'     gamma = xgb_tune4$bestTune$gamma,
#'     colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
#'     min_child_weight = xgb_tune2$bestTune$min_child_weight,
#'     subsample = xgb_tune3$bestTune$subsample
#'   )
#'   
#'   xgb_tune5 <- caret::train(
#'     group ~.,
#'     data=model.input,
#'     trControl = tune_control,
#'     tuneGrid = tune_grid5,
#'     method = "xgbTree",
#'     verbose = TRUE
#'   )
#'   
#'   cat("Tuning Step 5: Final Learning Rate \n COMPLETE \n")
#'   print(tuneplot(xgb_tune5))
#'   print(xgb_tune5$bestTune)
#'   
#'   
#'   #_____________________________________________________ 
#'   # 5)  Fitting the Model
#'   #_____________________________________________________ 
#'   
#'   final_grid <- expand.grid(
#'     nrounds = xgb_tune5$bestTune$nrounds,
#'     eta = xgb_tune5$bestTune$eta,
#'     max_depth = xgb_tune5$bestTune$max_depth,
#'     gamma = xgb_tune5$bestTune$gamma,
#'     colsample_bytree = xgb_tune5$bestTune$colsample_bytree,
#'     min_child_weight = xgb_tune5$bestTune$min_child_weight,
#'     subsample = xgb_tune5$bestTune$subsample
#'   )
#'   
#'   mytrainControl <-
#'     trainControl(method='repeatedcv',
#'                  number=5,
#'                  repeats=5,
#'                  search='grid',
#'                  savePredictions = TRUE,
#'                  classProbs = TRUE,
#'                  allowParallel = TRUE,
#'                  verboseIter = TRUE)
#'   
#'   
#'   # Run 5-fold CV Model with 5 Repetition
#'   model <-train(group ~.,
#'                 data=model.input, 
#'                 method='xgbTree', 
#'                 metric='Accuracy', 
#'                 tuneGrid=final_grid, 
#'                 trControl=mytrainControl,
#'                 verbose = TRUE)
#'   
#'   
#'   cat("Model Summary")
#'   print(model)
#'   cat("\n\n")
#'   # print(plot(model))
#'   
#'   intervalEnd <- Sys.time()
#'   cat("XGBoost model tuning completed in: ",
#'       intervalEnd - intervalStart, attr(intervalEnd - intervalStart, "units"))
#'   cat("\n\n")
#'   
#'   # Select model with optimal hyper-parameters
#'   df.output <- model$pred 
#'   
#'   # MLeval
#'   mleval <- evalm(model)
#'   mleval$roc
#'   output.list$optimal.df <- df.output
#'   # output.list$AUCROC <- mleval$roc
#'   output.list$AUCROC <- mleval$stdres$`Group 1`["AUC-ROC", "Score"]
#'   output.list$MLevaldata <- mleval$stdres
#'   
#'   
#'   return(output.list)
#'   
#' }
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' #_______________________________________________________________________________--------
#' #_____________________________________________________
#' #_____________________________________________________
#' #_____________________________________________________
#' 
#' 
#' 
#' #____________________________________________________________________________
#' ###         RIDGE, LASSO, ENET LOGISTIC REGRESSION  --  For Curated Datasets            
#' #____________________________________________________________________________
#' 
#' 
#' ridge.lasso.enet.regression.model.DS <- function(model.input, model.type){
#'   
#'   #' Function to train a Rigde, lASSO, or ElasticNet Regression Model
#'   #' Input: Phyloseq Obj, comparison of interest, and model type
#'   #' Returns: a list including 
#'   #' 1) The fitted Model 
#'   #' 2) A dataframe of values with the optimal parameters (for ROC plot)
#'   #' 3) MLeval AUC-ROC values 
#'   #' 4) Other MLeval analysis parameters
#'   
#'   # Initialize variables
#'   model <- NULL
#'   output.list <- vector(mode="list", length=4)
#'   names(output.list) <- c("fitted_model", "optimal.df", "AUCROC", "MLevaldata")
#'   
#'   if (model.type == "ridge"){
#'     tune.grid = expand.grid(alpha = 0, lambda=seq(0, 1, 0.1))
#'   } else if (model.type == "lasso"){
#'     tune.grid = expand.grid(alpha = 1, lambda=seq(0, 1, 0.1))
#'   } else if (model.type == "enet"){
#'     tune.grid = expand.grid(alpha = 0.5, lambda=seq(0, 1, 0.1))
#'   }
#'   
#'   # AST Transformation on merged data
#'   groupcol <- factor(model.input$group)
#'   temp <- dplyr::select(model.input, -group)
#'   model.input <- asin(sqrt(temp))
#'   model.input$group <- groupcol
#'   
#'   # Model Parameters
#'   numbers <- 10
#'   repeats <- 5  
#'   set.seed(42)
#'   seed <- 42
#'   rcvSeeds <- setSeeds(method = "repeatedcv", numbers = numbers, repeats = repeats, seed = seed)
#'   
#'   mytrainControl <- 
#'     trainControl(method='repeatedcv',
#'                  number=numbers, 
#'                  repeats=repeats,
#'                  search='grid',
#'                  seeds = rcvSeeds,
#'                  savePredictions = TRUE, 
#'                  classProbs = TRUE, 
#'                  verboseIter = TRUE)
#'   
#'   # Run 10-fold CV Model with 5 Repitions
#'   model <-train(group ~.,
#'                 data=model.input, 
#'                 method='glmnet', 
#'                 metric='Accuracy', 
#'                 tuneGrid=tune.grid, 
#'                 trControl=mytrainControl,
#'                 verboseIter = T)
#'   
#'   cat("Model Summary")
#'   print(model)
#'   cat("\n\n")
#'   print(plot(model))
#'   
#'   # Select model with optimal Lambda
#'   selectedIndices <- model$pred$lambda == model$bestTune$lambda
#'   df.output <- model$pred[selectedIndices, ]
#'   
#'   # MLeval
#'   mleval <- evalm(model)
#'   mleval$roc
#'   output.list$fitted_model <- model
#'   output.list$optimal.df <- df.output
#'   output.list$AUCROC <- mleval$stdres$`Group 1`["AUC-ROC", "Score"]
#'   output.list$MLevaldata <- mleval$stdres
#'   
#'   return(output.list)
#'   
#' }
#' 
#' 
#' #____________________________________________________________________________
#' ###   RIDGE, LASSO, ENET LOGISTIC REGRESSION  --  For Curated Datasets & PD COMPARISON        
#' #____________________________________________________________________________
#' 
#' 
#' ridge.lasso.enet.regression.model.DSxPD <- function(disease.model.input, model.type, obj = dat){
#'   
#'   # # TROUBLESHOOTING
#'   # disease.model.input = VincentC_2016.model.input;
#'   # obj = dat;
#'   # model.type = "enet"
#'   # 
#'   #' Function to train a Rigde, lASSO, or ElasticNet Regression Model
#'   #' Input: Phyloseq Obj, comparison of interest, and model type
#'   #' Returns: a list including 
#'   #' 1) The fitted Model 
#'   #' 2) A dataframe of values with the optimal parameters (for ROC plot)
#'   #' 3) MLeval AUC-ROC values 
#'   #' 4) Other MLeval analysis parameters
#'   
#'   # Initalize variables
#'   model <- NULL
#'   output.list <- vector(mode="list", length=4)
#'   names(output.list) <- c("fitted_model", "optimal.df", "AUCROC", "MLevaldata")
#'   
#'   if (model.type == "ridge"){
#'     tune.grid = expand.grid(alpha = 0, lambda=seq(0, 1, 0.1))
#'   } else if (model.type == "lasso"){
#'     tune.grid = expand.grid(alpha = 1, lambda=seq(0, 1, 0.1))
#'   } else if (model.type == "enet"){
#'     tune.grid = expand.grid(alpha = 0.5, lambda=seq(0, 1, 0.1))
#'   }
#'   
#'   # Trim Control Samples from Disease Dataset 
#'   disease.model.input <- filter(disease.model.input, group != "control")
#'   # Load Native PD Dataset
#'   dat_pd = subset_samples(dat, donor_group == "PD")
#'   d <- dat_pd %>%
#'     microbiome::transform("compositional") %>% 
#'     microbiome::abundances() %>% 
#'     t() %>% 
#'     as.data.frame() 
#'   pd.model.input <- group_col_from_ids(d, id= rownames(d))
#'   rownames(pd.model.input) <- rownames(d)
#'   
#'   
#'   # Merge Disease and PD input datasets
#'   model.input <- full_join(pd.model.input, disease.model.input)
#'   groupcol <- factor(model.input$group)
#'   
#'   # option 1)   Replace all NAs with 0s 
#'   # model.input[is.na(model.input)] = 0
#'   # option 2)  Trim all features that aren't shared (Detected in both groups)
#'   model.input <- model.input[ ,colSums(is.na(model.input)) == 0]
#'   
#'   # AST Transformation on merged data
#'   temp <- dplyr::select(model.input, -group)
#'   model.input <- asin(sqrt(temp))
#'   model.input$group <- groupcol
#'   
#'   
#'   
#'   # Model Parameters
#'   numbers <- 10
#'   repeats <- 5  
#'   set.seed(42)
#'   seed <- 42
#'   rcvSeeds <- setSeeds(method = "repeatedcv", 
#'                        numbers = numbers, repeats = repeats, seed = seed)
#'   
#'   mytrainControl <- 
#'     trainControl(method='repeatedcv',
#'                  number=numbers, 
#'                  repeats=repeats,
#'                  search='grid',
#'                  seeds = rcvSeeds,
#'                  savePredictions = TRUE, 
#'                  classProbs = TRUE, 
#'                  verboseIter = TRUE)
#'   
#'   # Run 10-fold CV Model with 5 Repitions
#'   model <-train(group ~.,
#'                 data=model.input, 
#'                 method='glmnet', 
#'                 metric='Accuracy', 
#'                 tuneGrid=tune.grid, 
#'                 trControl=mytrainControl,
#'                 verboseIter = T)
#'   
#'   cat("Model Summary")
#'   print(model)
#'   cat("\n\n")
#'   print(plot(model))
#'   
#'   # Select model with optimal Lambda
#'   selectedIndices <- model$pred$lambda == model$bestTune$lambda
#'   df.output <- model$pred[selectedIndices, ]
#'   
#'   # MLeval
#'   mleval <- evalm(model)
#'   mleval$roc
#'   output.list$fitted_model <- model
#'   output.list$optimal.df <- df.output
#'   output.list$AUCROC <- mleval$stdres$`Group 1`["AUC-ROC", "Score"]
#'   output.list$MLevaldata <- mleval$stdres
#'   
#'   return(output.list)
#'   
#' }


