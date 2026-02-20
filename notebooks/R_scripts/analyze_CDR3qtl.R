# Joe Boktor
# 2025-08-02

library(tidyverse)
library(glue)
library(tictoc)
library(furrr)
library(janitor)
library(MVLM)
library(broom)

wkdir <- "/central/groups/MazmanianLab/joeB/PDMBS/parkinsons-microbial-blood-signatures"
c3_pheno_dir <- glue("{wkdir}/data/interim/cdr3_qtl/cdr3_phenotypes")
gdata_dir <- glue("{wkdir}/data/interim/cdr3_qtl/cdr3_hla_genotypes")

# Loading formatted HLA-site data
full_aa_hla_sites <- readRDS(
  glue("{gdata_dir}/all_hla_AA_sites_2025-07-31.rds")) %>% 
  pillar::glimpse()

# phenotype_type <- "cdr3_p1_strat_filtered_cv1k"
# hsite <- "DRB1_13"
# cdr3length <- 13
# cdr3_pos <- "P109"


# This function aggregates the CDR3 phenotype data, genotype data, and metadata
# into a single data frame for analysis.
# phenotype_type: the type of CDR3 phenotype to load
# hsite: the HLA site to analyze
# cdr3_pos: the CDR3 position to analyze
# cdr3length: the CDR3 length to filter on
# Returns a list containing the aggregated data frame (M), genotype data (gdata),
# CDR3 phenotype data (pdata), and covariate data (covar_df).
# CDR3 phenotype data is normalized by inverse normal transformation
generate_model_inputs <- function(phenotype_type, hsite, cdr3_pos, cdr3length) {
  message("=== Starting CDR3 QTL Analysis ===")
  message(glue("Phenotype type: {phenotype_type}"))
  message(glue("HLA site: {hsite}"))
  message(glue("CDR3 position: {cdr3_pos}"))
  message(glue("CDR3 length filter: {if(is.na(cdr3length)) 'None' else cdr3length}"))
  
  # Loading CDR3 phenotypes based on phenotype_type parameter
  message("\n--- Loading CDR3 phenotype data ---")
  
  switch(phenotype_type,
    "cdr3_p1_strat_filtered_cv1k" = {
      message("Loading stratified CDR3 phenotype data (model 1)...")
      pheno_data <- readRDS(glue("{c3_pheno_dir}/",
        "cdr3_pheno_m1_cdr3len_stratified_filtered_cv1000_2025-08-07.rds")
      )
    },
    "cdr3_p1_unstrat_filtered_cv1k" = {
      message("Loading unstratified CDR3 phenotype data (model 1)...")
      pheno_data <- readRDS(glue("{c3_pheno_dir}/",
        "cdr3_pheno_m1_cdr3len_UNstratified_filtered_cv1000_2025-08-07.rds")
      )
    },
    "cdr3_p2_strat_filtered_cv1k" = {
      message("Loading stratified CDR3 phenotype data (model 2)...")
      pheno_data <- readRDS(glue("{c3_pheno_dir}/",
        "cdr3_pheno_m2_cdr3len_stratified_filtered_cv1000_2025-08-07.rds")
      )
    },
    "cdr3_p2_unstrat_filtered_cv1k" = {
      message("Loading unstratified CDR3 phenotype data (model 2)...")
      pheno_data <- readRDS(glue("{c3_pheno_dir}/",
        "cdr3_pheno_m2_cdr3len_UNstratified_filtered_cv1000_2025-08-07.rds")
      )
    },
    stop(glue("Invalid phenotype_type. Choose from: ",
        "cdr3_p1_strat_filtered_cv1k, cdr3_p1_unstrat_filtered_cv1k,",
        " cdr3_p2_strat_filtered_cv1k, cdr3_p2_unstrat_filtered_cv1k"))
  )

  message(glue("Loaded phenotype data with {nrow(pheno_data)} rows",
    " and {ncol(pheno_data)} columns"))

  message("\n--- Loading metadata and covariates ---")
  metadata <- readRDS(
    glue("{wkdir}/data/interim/metadata/",
      "compiled_metadata_v4_release_2025-07-21.rds")) %>% 
    dplyr::filter(case_control_other_latest != "Other")

  covar_df <- metadata %>% 
    dplyr::select(participant_id, pd_hla_risk_score,
      MHC_all_PC1, MHC_all_PC2, MHC_all_PC3) %>%
    dplyr::distinct() %>%
    pillar::glimpse()
  message(glue("Covariate data: {nrow(covar_df)} unique participants"))

  message("\n--- Processing genotype data ---")
  # first prepping genotype data into a wide format
  gdata <- full_aa_hla_sites %>% 
    dplyr::filter(hla_site == hsite) %>% 
    dplyr::filter(nchar(HLA_AA_site) == 1) %>%
    dplyr::select(-c(fdigit, aatype)) %>% 
    dplyr::group_by(participant_id, HLA_AA_site) %>% 
    dplyr::summarize(allele_dosage = sum(allele_dosage)) %>% 
    tidyr::pivot_wider(names_from = HLA_AA_site, 
        values_from = allele_dosage) %>%
    dplyr::rename_at(vars(-participant_id), ~ paste0(., "_allele")) %>% 
    pillar::glimpse()
  message(glue("Genotype data for {hsite}: {nrow(gdata)} participants,",
    " {ncol(gdata)-1} HLA alleles"))

  message("\n--- Processing CDR3 phenotype data ---")
  pdata <- pheno_data %>% 
    dplyr::filter(imgt_pos == cdr3_pos) %>% 
    dplyr::filter(if (is.na(cdr3length)) TRUE else cdr3len == cdr3length) %>%
    dplyr::select(participant_id, AA = aa, ratio) %>%
    pillar::glimpse()
  message(glue("Filtered phenotype data: {nrow(pdata)} rows,",
    " {length(unique(pdata$AA))} unique amino acids"))

  message("\n--- Performing inverse normal transformation ---")
  pall <- data.frame()
  pnamelist <- unique(pdata$AA)
  message(glue("Processing {length(pnamelist)} amino acids..."))

  for (ptarget in pnamelist) {
    pdata2 <- subset(pdata, AA == ptarget)
    pdata2 <- pdata2[, c("participant_id", "ratio")]
    x <- pdata2$ratio
    # inverse normal normalization
    pdata2$normrate <- qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
    df <- pdata2[, c("participant_id", "normrate")]
    df$AA <- ptarget
    pall <- rbind(pall, df)
  }
  message(glue("Transformation complete: {nrow(pall)} total observations"))

  message("\n--- Creating wide format phenotype matrix ---")
  cdr3_freq_df <- pall %>%
    tidyr::pivot_wider(names_from = AA, values_from = normrate) %>% 
    pillar::glimpse()
  message(glue("Phenotype matrix: {nrow(cdr3_freq_df)} participants,",
    " {ncol(cdr3_freq_df)-1} amino acid features"))

  message("\n--- Joining all data sources ---")
  # Joining data
  M <- cdr3_freq_df %>% 
    dplyr::left_join(gdata, by = "participant_id") %>% 
    dplyr::left_join(covar_df, by = "participant_id") %>% 
    pillar::glimpse()
  
  message(glue("Final dataset: {nrow(M)} participants,",
    " {ncol(M)} total features"))
  message(glue("Features breakdown: {ncol(cdr3_freq_df)-1} amino acids +",
    " {ncol(gdata)-1} HLA alleles + {ncol(covar_df)-1} covariates"))
  message("\n=== Analysis setup complete ===")
  
  return(list(M = M, gdata = gdata, pdata = pall, covar_df = covar_df))
}

generate_model_inputs_cdr3_riskscore <- function(phenotype_type, cdr3_pos, cdr3length) {
  message("=== Starting CDR3 QTL Analysis ===")
  message(glue("Phenotype type: {phenotype_type}"))
  message(glue("CDR3 position: {cdr3_pos}"))
  message(glue("CDR3 length filter: {if(is.na(cdr3length)) 'None' else cdr3length}"))
  
  # Loading CDR3 phenotypes based on phenotype_type parameter
  message("\n--- Loading CDR3 phenotype data ---")
  
  switch(phenotype_type,
    "cdr3_p1_strat_filtered_cv1k" = {
      message("Loading stratified CDR3 phenotype data (model 1)...")
      pheno_data <- readRDS(glue("{c3_pheno_dir}/",
        "cdr3_pheno_m1_cdr3len_stratified_filtered_cv1000_2025-08-07.rds")
      )
    },
    "cdr3_p1_unstrat_filtered_cv1k" = {
      message("Loading unstratified CDR3 phenotype data (model 1)...")
      pheno_data <- readRDS(glue("{c3_pheno_dir}/",
        "cdr3_pheno_m1_cdr3len_UNstratified_filtered_cv1000_2025-08-07.rds")
      )
    },
    "cdr3_p2_strat_filtered_cv1k" = {
      message("Loading stratified CDR3 phenotype data (model 2)...")
      pheno_data <- readRDS(glue("{c3_pheno_dir}/",
        "cdr3_pheno_m2_cdr3len_stratified_filtered_cv1000_2025-08-07.rds")
      )
    },
    "cdr3_p2_unstrat_filtered_cv1k" = {
      message("Loading unstratified CDR3 phenotype data (model 2)...")
      pheno_data <- readRDS(glue("{c3_pheno_dir}/",
        "cdr3_pheno_m2_cdr3len_UNstratified_filtered_cv1000_2025-08-07.rds")
      )
    },
    stop(glue("Invalid phenotype_type. Choose from: ",
        "cdr3_p1_strat_filtered_cv1k, cdr3_p1_unstrat_filtered_cv1k,",
        " cdr3_p2_strat_filtered_cv1k, cdr3_p2_unstrat_filtered_cv1k"))
  )

  message(glue("Loaded phenotype data with {nrow(pheno_data)} rows",
    " and {ncol(pheno_data)} columns"))

  message("\n--- Loading metadata and covariates ---")
  metadata <- readRDS(
    glue("{wkdir}/data/interim/metadata/",
      "compiled_metadata_v4_release_2025-07-21.rds")) %>% 
    dplyr::filter(case_control_other_latest != "Other")

  covar_df <- metadata %>% 
    dplyr::select(participant_id, pd_hla_risk_score,
      MHC_all_PC1, MHC_all_PC2, MHC_all_PC3) %>%
    dplyr::distinct() %>%
    pillar::glimpse()
  message(glue("Covariate data: {nrow(covar_df)} unique participants"))

  message("\n--- Processing CDR3 phenotype data ---")
  pdata <- pheno_data %>% 
    dplyr::filter(imgt_pos == cdr3_pos) %>% 
    dplyr::filter(if (is.na(cdr3length)) TRUE else cdr3len == cdr3length) %>%
    dplyr::select(participant_id, AA = aa, ratio) %>%
    pillar::glimpse()
  message(glue("Filtered phenotype data: {nrow(pdata)} rows,",
    " {length(unique(pdata$AA))} unique amino acids"))

  message("\n--- Performing inverse normal transformation ---")
  pall <- data.frame()
  pnamelist <- unique(pdata$AA)
  message(glue("Processing {length(pnamelist)} amino acids..."))

  for (ptarget in pnamelist) {
    pdata2 <- subset(pdata, AA == ptarget)
    pdata2 <- pdata2[, c("participant_id", "ratio")]
    x <- pdata2$ratio
    # inverse normal normalization
    pdata2$normrate <- qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
    df <- pdata2[, c("participant_id", "normrate")]
    df$AA <- ptarget
    pall <- rbind(pall, df)
  }
  message(glue("Transformation complete: {nrow(pall)} total observations"))

  message("\n--- Creating wide format phenotype matrix ---")
  cdr3_freq_df <- pall %>%
    tidyr::pivot_wider(names_from = AA, values_from = normrate) %>% 
    pillar::glimpse()
  message(glue("Phenotype matrix: {nrow(cdr3_freq_df)} participants,",
    " {ncol(cdr3_freq_df)-1} amino acid features"))

  message("\n--- Joining all data sources ---")
  # Joining data
  M <- cdr3_freq_df %>% 
    dplyr::left_join(covar_df, by = "participant_id") %>% 
    pillar::glimpse()
  
  message(glue("Final dataset: {nrow(M)} participants,",
    " {ncol(M)} total features"))
  message(glue("Features breakdown: {ncol(cdr3_freq_df)-1} amino acids +",
    " {ncol(covar_df)-1} covariates"))
  message("\n=== Analysis setup complete ===")
  
  return(list(M = M, pdata = pall, covar_df = covar_df))
}

# Function to build the formula for the MVLM model
# amino_acids: list of amino acids to include in the model
# hla_alleles: list of HLA alleles to include in the model
# covariates: list of covariates to include in the model
build_mvlm_formula <- function(amino_acids, hla_alleles, covariates) {
  # Build the response variables (amino acids)
  response_vars <- paste(amino_acids, collapse = ", ")
  
  # Build the predictor variables
  predictors <- c(hla_alleles, covariates)
  predictor_vars <- paste(predictors, collapse = " + ")
  
  # Construct the formula
  formula_str <- glue("cbind({response_vars}) ~ {predictor_vars}")
  
  message(glue("Building QTL model formula:"))
  message(glue("Response variables: {response_vars}"))
  message(glue("Predictor variables: {predictor_vars}"))
  
  return(as.formula(formula_str))
}

run_lm_site_pos_pair <- function(pdata, gdata, covar_df, hla_aa_allele, cdr3_aa) {
  require(glue)
  require(dplyr)
  
  # filtering the cdr3_aa data
  pdata2 <- subset(pdata, AA==cdr3_aa) %>% 
    dplyr::select(participant_id, normrate) %>% 
    glimpse()
  
  # preparing genotype data
  gdata2 <- gdata %>% dplyr::select(participant_id, {{hla_aa_allele}})
  
  # Joining data
  M_data <- pdata2 %>% 
    dplyr::left_join(gdata2, by = "participant_id") %>% 
    dplyr::left_join(covar_df, by = "participant_id")

  model_formula <- as.formula(
      glue("normrate ~ {hla_aa_allele} + MHC_all_PC1+MHC_all_PC2+MHC_all_PC3")
      )
  mod <- lm(model_formula, data = M_data)
  return(mod)
}

run_lm_site_pos_pair_riskscore <- function(pdata, covar_df, cdr3_aa) {
  require(glue)
  require(dplyr)
  
  # filtering the cdr3_aa data
  pdata2 <- subset(pdata, AA==cdr3_aa) %>% 
    dplyr::select(participant_id, normrate) %>% 
    glimpse()
  
  # Joining data
  M_data <- pdata2 %>% 
    dplyr::left_join(covar_df, by = "participant_id")

  model_formula <- as.formula(
      glue("normrate ~ pd_hla_risk_score + MHC_all_PC1+MHC_all_PC2+MHC_all_PC3")
      )
  mod <- lm(model_formula, data = M_data)
  return(mod)
}


run_qtl_stats <- function(results_dir, phenotype_type, hsite, cdr3_pos, cdr3length, wkdir) {
    require(tidyverse)
    require(purrr)
    
    source(paste0(wkdir, "/notebooks/R_scripts/analyze_CDR3qtl.R"))
    run_descriptor <- glue("{phenotype_type}__HLAsites-{hsite}__CDR3pos-{cdr3_pos}__CDR3len-{cdr3length}")
    dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

    dlist <- generate_model_inputs(
        phenotype_type = phenotype_type,
        hsite = hsite,
        cdr3_pos = cdr3_pos,
        cdr3length = cdr3length
    )
    M <- dlist$M
    gdata <- dlist$gdata
    pdata <- dlist$pdata
    covar_df <- dlist$covar_df

    #' get the names of HLA alleles AA's excluding the top most 
    #' abundant as a reference allele
    if (length(colnames(gdata[,-1])) > 1) {
        hla_aa_alleles_minustop <- gdata[,-1] %>% colSums() %>% 
            sort(decreasing = TRUE) %>% names() %>% .[-1]
    } else {
        hla_aa_alleles_minustop <- colnames(gdata[,-1])
    }

    full_formula <- build_mvlm_formula(
        unique(pdata$AA), 
        hla_aa_alleles_minustop,
        c("MHC_all_PC1", "MHC_all_PC2", "MHC_all_PC3")
    )
    null_formula <- build_mvlm_formula(
        unique(pdata$AA), 
        c(), 
        c("MHC_all_PC1", "MHC_all_PC2", "MHC_all_PC3")
    )
    # MANOVA p-value (Pillai's trace for model fit)
    mod0_mlm <- lm(null_formula, data = M)
    mod1_mlm <- lm(full_formula, data = M)
    mlm_res <- anova(mod1_mlm, mod0_mlm)


    # MVLM explained variance
    mod0_mvlm <- mvlm(null_formula, data = M)
    mod1_mvlm <- mvlm(full_formula, data = M)
    var_exp_null_mvlm <- mod0_mvlm$pseudo.rsq["Omnibus Effect",1]
    var_exp_full_mvlm <- mod1_mvlm$pseudo.rsq["Omnibus Effect",1]
    corrected_var_exp_mvlm <- var_exp_full_mvlm - var_exp_null_mvlm

    mlm_mvlm_res <- list(
        "manova_null" = mod0_mlm,
        "manova_full" = mod1_mlm,
        "manova_res" = mlm_res,
        "mvlm_null" = mod0_mvlm,
        "mvlm_full" = mod1_mvlm,
        "mvlm_var_exp_null" = var_exp_null_mvlm,
        "mvlm_var_exp_full" = var_exp_full_mvlm,
        "mvlm_corrected_var_exp" = corrected_var_exp_mvlm
    )

    # Save the results
    saveRDS(mlm_mvlm_res, glue("{results_dir}/mlm_mvlm_{run_descriptor}.rds"))

    lm_map_df <- expand.grid(
        "hla_aa_allele" = colnames(gdata[,-1]),
        "cdr3_aa" = unique(pdata$AA))

    # running the LM for each combination
    mod_list <- 1:nrow(lm_map_df) %>%
    purrr::set_names(
        glue("HLA-{lm_map_df$hla_aa_allele[.]}",
            " CDR3-{lm_map_df$cdr3_aa[.]}")
            ) %>% 
    purrr::map(
        ~ run_lm_site_pos_pair(
            pdata = pdata,
            gdata = gdata,
            covar_df = covar_df,
            hla_aa_allele = lm_map_df$hla_aa_allele[.], 
            cdr3_aa = lm_map_df$cdr3_aa[.]
        )
    )

    lm_results_df <- mod_list %>% 
        purrr::map(broom::tidy) %>% 
        bind_rows(.id = "model_name") %>% 
        mutate(cdr3_aa = strex::str_after_last(model_name, "CDR3-")) %>% 
        glimpse()

    lm_fileout <- glue("{results_dir}/lm_model_results__{run_descriptor}")
    write_csv(lm_results_df, glue("{lm_fileout}.csv"))
    saveRDS(lm_results_df, glue("{lm_fileout}.rds")) # may want to change this to save model list instead
}

calculate_cdr3_riskscore_betacoefs <- function(results_dir, phenotype_type, cdr3_pos, cdr3length) {
    require(tidyverse)
    require(purrr)
    run_descriptor <- glue(
      "{phenotype_type}",
      "__CDR3pos-{cdr3_pos}",
      "__CDR3len-{cdr3length}"
      )
    dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

    dlist <- generate_model_inputs_cdr3_riskscore(
        phenotype_type = phenotype_type,
        cdr3_pos = cdr3_pos,
        cdr3length = cdr3length
    )

    pdata <- dlist$pdata
    covar_df <- dlist$covar_df

    # running the LM for each combination
    mod_list <- unique(pdata$AA) %>%
      purrr::set_names(glue("PD-CDR3-Risk-Score_CDR3-{.}")) %>% 
      purrr::map(
          ~ run_lm_site_pos_pair_riskscore(
              pdata = pdata, 
              covar_df = covar_df,
              cdr3_aa = .
          )
      )

    lm_results_df <- mod_list %>% 
        purrr::map(broom::tidy) %>% 
        bind_rows(.id = "model_name") %>% 
        mutate(cdr3_aa = strex::str_after_last(model_name, "CDR3-")) %>% 
        glimpse()

    lm_fileout <- glue("{results_dir}/lm_model_results__{run_descriptor}")
    lm_modlist <- glue("{results_dir}/lm_model_list__{run_descriptor}")

    write_csv(lm_results_df, glue("{lm_fileout}.csv"))
    saveRDS(mod_list, glue("{lm_modlist}.rds"))
}
