# R functions for running association screening with MCPAS and VDJDB
require(glue)
require(future)
require(tidyverse)
require(furrr)
require(nlme)

wkdir <- "/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr"
metrics_dir <- glue("{wkdir}/data/interim/airr/immunarch/metrics_tcrb")
mcpas_stats_dir <- glue("{metrics_dir}/association_screening_McPAS")
vdjdb_stats_dir <- glue("{metrics_dir}/association_screening_VDJdb")
refdir <- glue("{wkdir}/data/input/reference_datasets")

# ------------------------------------------------------------
# Loading data tables ----
# VDJdb ----
vdjdb <- read_tsv(
  glue("{refdir}/SearchTable-2025-06-11 00_47_14.856.tsv")) %>% 
  filter(Gene == "TRB") %>% 
  filter(Species == "HomoSapiens") %>% 
  janitor::clean_names() %>%
  filter(epitope_species != "synthetic") %>% 
  mutate(epitope_species = case_when(
    epitope_species == "Mtb" ~ "M.tuberculosis",
    epitope_species == "Hcv" ~ "Hepatitis C virus (HCV)",
    epitope_species == "HPV-16" ~ "HPV",
    epitope_species == "HCoV-HKU1" ~ "SARS-CoV",
    epitope_species == "SARS-CoV-2" ~ "SARS-CoV",
    TRUE ~ epitope_species
  ))
vdjdb_path_vars <- table(vdjdb$epitope_species)
testable_vdjdb_vars <- names(vdjdb_path_vars[vdjdb_path_vars > 30])

# McPAS ----
mcpas <- read_csv(
  glue("{refdir}/McPAS-TCR.csv")) %>% 
  filter(!is.na(CDR3.beta.aa)) %>% 
  filter(Species == "Human") %>% 
  mutate(Pathology = case_when(
    Pathology == "Hepatitis C virus" ~ "Hepatitis C virus (HCV)",
    Pathology == "M.Tuberculosis" ~ "M. tuberculosis",
    TRUE ~ Pathology
  ))
mcpas_path_vars <- table(mcpas$Pathology)
testable_mcpas_vars <- names(mcpas_path_vars[mcpas_path_vars > 30])

mcpas_map <- readRDS(
  glue("{metrics_dir}/McPAS_mapping_2025-07-15.rds")
)
vdjdb_map <- readRDS(
  glue("{metrics_dir}/VDJdb_mapping_2025-07-15.rds")
)

clonotype_volume_df <- readRDS(
    glue("{metrics_dir}/clonotype_volume_df_2025-07-18.rds")
    )
metadata <- readRDS(
  glue("{wkdir}/data/interim/metadata/",
    "compiled_metadata_v4_release_2025-07-21.rds")) %>%
  dplyr::filter(case_control_other_latest != "Other")

# ------------------------------------------------------------
# Functions ----
# Utility function to screen explained variance of metadata variables
test_metadata_variable <- function(data, response, metadata_col,
                                   random_effect = "participant_id") {
  # filter out non-NA data
  data <- data %>% 
    filter(!is.na(.data[[metadata_col]])) %>% 
    filter(!is.na(match_fraction))
  
  message(glue("Testing {metadata_col} with {nrow(data)} samples",
    " and {unique(data$participant_id) %>% length()} participants"))
  # Build formulas
  response <- as.name(response)
  metadata_col <- as.name(metadata_col)
  full_formula <- stats::as.formula(
    paste(deparse(response), "~", deparse(metadata_col))
  )
  null_formula <- stats::as.formula(
    paste(deparse(response), "~ 1")
  )

  # Fit models
  full_model <- nlme::lme(
    fixed = full_formula,
    random = stats::as.formula(paste("~1 |", random_effect)),
    data = data,
    na.action = na.omit,
    method = "ML" # Use ML for model comparison
  )

  null_model <- nlme::lme(
    fixed = null_formula,
    random = stats::as.formula(paste("~1 |", random_effect)),
    data = data,
    na.action = na.omit,
    method = "ML"
  )

  # Likelihood ratio test
  model_comparison <- stats::anova(null_model, full_model)
  
  # Extract metrics and include the full anova model comparison
  res <- list(
    LRT = model_comparison$L.Ratio[2],
    p_value = model_comparison$`p-value`[2],
    null_AIC = AIC(null_model),
    full_AIC = AIC(full_model),
    null_BIC = BIC(null_model),
    full_BIC = BIC(full_model),
    null_logLik = logLik(null_model)[1],
    full_logLik = logLik(full_model)[1],
    null_df = attr(logLik(null_model), "df"),
    full_df = attr(logLik(full_model), "df"),
    LRT_model = model_comparison,
    n_samples = unique(data$sample_id) %>% length(),
    n_participants = unique(data$participant_id) %>% length()
  )
  return(res)
}

# SLURM functions ----
run_mcpas_association_screening <- function(pathology, testable_meta,
                                            wkdir, nthreads) {
  require(glue)
  source(glue("{wkdir}/notebooks/R_scripts/_airr_mcpas_vdjdb_funcs.R"))
  options(future.globals.maxSize = 2 * 1024^3)  # 2 GB
  future::plan(future::multisession, workers = nthreads)
  
  outfile <- glue("{mcpas_stats_dir}/{pathology}_res_df_2025-07-17.rds")
  message(glue("Testing: {pathology}"))
  p_cdr3b <- dplyr::filter(mcpas, Pathology == pathology) %>%
    dplyr::pull(CDR3.beta.aa)

  formatted_res_df <- mcpas_map %>%
    dplyr::filter(CDR3.aa %in% p_cdr3b) %>%
    dplyr::select(-Samples) %>%
    tidyr::pivot_longer(!CDR3.aa,
      names_to = "sample_id",
      values_to = "mcpas_match"
    ) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarize(mcpas_total = sum(mcpas_match)) %>%
    dplyr::left_join(clonotype_volume_df) %>%
    dplyr::mutate(match_fraction = log((mcpas_total / tcrb_volume) + 1e-4)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(metadata)

  # Parallel LRT testing with furrr
  res_df <- furrr::future_map_dfr(testable_meta, function(m) {
    message(glue("Testing: {m}"))
    tryCatch(
      {
        df_input <- formatted_res_df %>% dplyr::filter(!is.na(.data[[m]]))
        res <- test_metadata_variable(df_input, "match_fraction", m)
        tibble(
          variable = m,
          LRT = res$LRT,
          p_value = res$p_value,
          null_AIC = res$null_AIC,
          full_AIC = res$full_AIC,
          null_BIC = res$null_BIC,
          full_BIC = res$full_BIC,
          null_logLik = res$null_logLik,
          full_logLik = res$full_logLik,
          null_df = res$null_df,
          full_df = res$full_df,
          LRT_model = res$LRT_model,
          n_samples = res$n_samples,
          n_participants = res$n_participants
        )
      },
      error = function(e) {
        warning(glue("Error testing {m}: {e$message}"))
        tibble(variable = m, LRT = NA, p_value = NA)
      }
    )
  })
  saveRDS(res_df, outfile)
}

run_vdjdb_association_screening <- function(epitope_group, testable_meta,
                                            wkdir, nthreads) {
  require(glue)
  source(glue("{wkdir}/notebooks/R_scripts/_airr_mcpas_vdjdb_funcs.R"))
  options(future.globals.maxSize = 2 * 1024^3)  # 2 GB
  future::plan(future::multisession, workers = nthreads)
  
  outfile <- glue("{vdjdb_stats_dir}/{epitope_group}_res_df_2025-07-17.rds")
  message(glue("Testing: {epitope_group}"))
  p_cdr3b <- dplyr::filter(vdjdb, epitope_species == epitope_group) %>%
    dplyr::pull(cdr3)

  formatted_res_df <- vdjdb_map %>%
    dplyr::filter(CDR3.aa %in% p_cdr3b) %>%
    dplyr::select(-Samples) %>%
    tidyr::pivot_longer(!CDR3.aa,
      names_to = "sample_id",
      values_to = "vdjdb_match"
    ) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarize(vdjdb_total = sum(vdjdb_match)) %>%
    dplyr::left_join(clonotype_volume_df) %>%
    dplyr::mutate(match_fraction = log((vdjdb_total / tcrb_volume) + 1e-4)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(metadata)

  # Parallel LRT testing with furrr
  res_df <- furrr::future_map_dfr(testable_meta, function(m) {
    message(glue("Testing: {m}"))
    tryCatch(
      {
        df_input <- formatted_res_df %>% dplyr::filter(!is.na(.data[[m]]))
        res <- test_metadata_variable(df_input, "match_fraction", m)
        tibble::tibble(variable = m,
          LRT = res$LRT,
          p_value = res$p_value,
          null_AIC = res$null_AIC,
          full_AIC = res$full_AIC,
          null_BIC = res$null_BIC,
          full_BIC = res$full_BIC,
          null_logLik = res$null_logLik,
          full_logLik = res$full_logLik,
          null_df = res$null_df,
          full_df = res$full_df,
          LRT_model = res$LRT_model,
          n_samples = res$n_samples,
          n_participants = res$n_participants
          )
      },
      error = function(e) {
        warning(glue("Error testing {m}: {e$message}"))
        tibble::tibble(variable = m, LRT = NA, p_value = NA)
      }
    )
  })
  saveRDS(res_df, outfile)
}

pull_epitope_matches <- function(epitope_group, db) {
  if (db == "McPAS") {
    cdr3s <- mcpas %>% 
      filter(Pathology == epitope_group) %>% 
      pull(CDR3.beta.aa)
    df <- mcpas_map
  }
  else if (db == "VDJdb") {
    cdr3s <- vdjdb %>% 
      filter(epitope_species == epitope_group) %>% 
      pull(cdr3)
    df <- vdjdb_map
  }
  else {
    stop("Invalid database")
  }
  final_df <- df %>%
    filter(CDR3.aa %in% cdr3s) %>% 
    dplyr::select(-Samples) %>% 
    pivot_longer(!CDR3.aa, names_to = "sample_id", 
      values_to = "match") %>% 
    group_by(sample_id) %>% 
    dplyr::summarize(match_total = sum(match)) %>% 
    left_join(clonotype_volume_df) %>% 
    mutate(match_fraction = log((match_total / tcrb_volume) + 1e-4)) %>% 
    ungroup() %>%
    left_join(metadata)
  return(final_df)
}






# run_mcpas_association_screening <- function(pathology, testable_meta, wkdir, nthreads) {
#   require(glue)
#   require(future)
#   require(tidyverse)
#   require(furrr)
#   source(glue("{wkdir}/notebooks/R_scripts/_airr_mcpas_vdjdb_funcs.R"))
#   options(future.globals.maxSize = 2 * 1024^3)  # 2 GB
#   future::plan(future::multisession, workers = nthreads)
#   metrics_dir <- glue("{wkdir}/data/interim/airr/immunarch/metrics_tcrb")
#   mcpas_stats_dir <- glue("{metrics_dir}/association_screening_McPAS")
#   refdir <- glue("{wkdir}/data/input/reference_datasets")
#   outfile <- glue("{mcpas_stats_dir}/{pathology}_res_df_2025-07-17.rds")
  
#   # if (file.exists(outfile)) next
#   message(glue("Testing: {pathology}"))
  
#   clonotype_volume_df <- readRDS(
#       glue("{metrics_dir}/clonotype_volume_df_2025-07-18.rds"))
#   metadata <- readRDS(
#     glue("{wkdir}/data/interim/metadata/compiled_metadata_v4_release_2025-07-17.rds")) %>% 
#     dplyr::filter(case_control_other_latest != "Other")

#   mcpas_map <- readRDS(glue("{metrics_dir}/McPAS_mapping_2025-07-15.rds"))
#   mcpas <- readr::read_csv(
#     glue("{refdir}/McPAS-TCR.csv")) %>% 
#     dplyr::filter(!is.na(CDR3.beta.aa)) %>% 
#     dplyr::filter(Species == "Human") %>% 
#     dplyr::mutate(Pathology = case_when(
#       Pathology == "Hepatitis C virus" ~ "Hepatitis C virus (HCV)",
#       Pathology == "M.Tuberculosis" ~ "M. tuberculosis",
#       TRUE ~ Pathology
#     ))
  
#   p_cdr3b <- dplyr::filter(mcpas, Pathology == pathology) %>% 
#     dplyr::pull(CDR3.beta.aa)

#   formatted_res_df <- mcpas_map %>%
#     dplyr::filter(CDR3.aa %in% p_cdr3b) %>% 
#     dplyr::select(-Samples) %>%
#     tidyr::pivot_longer(!CDR3.aa, names_to = "sample_id", 
#       values_to = "mcpas_match") %>% 
#     dplyr::group_by(sample_id) %>% 
#     dplyr::summarize(mcpas_total = sum(mcpas_match)) %>% 
#     dplyr::left_join(clonotype_volume_df) %>% 
#     dplyr::mutate(match_fraction = log((mcpas_total/tcrb_volume) + 1e-4)) %>% 
#     dplyr::ungroup() %>%
#     dplyr::left_join(metadata)
  
#     # Parallel LRT testing with furrr
#     res_df <- furrr::future_map_dfr(testable_meta, function(m) {
#       message(glue("Testing: {m}"))  
#       tryCatch({
#         df_input <- formatted_res_df %>% dplyr::filter(!is.na(.data[[m]]))
#         res <- test_metadata_variable(df_input, "match_fraction", m)
#         tibble(variable = m, LRT = res$L.Ratio[2], p_value = res$`p-value`[2])
#       }, error = function(e) {
#         warning(glue("Error testing {m}: {e$message}"))
#         tibble(variable = m, LRT = NA, p_value = NA)
#       })
#     })
#     saveRDS(res_df, outfile)
# }