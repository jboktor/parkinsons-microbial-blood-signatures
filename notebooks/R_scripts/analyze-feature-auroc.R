
source("src/_load_packages.R")
# refDBlist <- c("RefSeqPlusPF", "UHGG", "WoL")
base::load("data/Phyloseq_Objects/Phyloseq_all_outliers_removed.RData")

aucs <- tibble()
start_time <- Sys.time()

for (ref_DB in names(refDB_phyloseq_orm)){
  for (level in names(refDB_phyloseq_orm[[ref_DB]])) {
    
    #_____________________________________________________________________________
    #####                    Limma - VOOM Normalization                     ##### 
    #_____________________________________________________________________________
    ps_object <- refDB_phyloseq_orm[[ref_DB]][[level]]
    metadat <- meta(ps_object)
    print(ps_object)
    
    voom_matrix <-
      model.matrix( ~ 0 + case_control_other_latest + sex,
                    data = metadat)
    
    dge <- 
      ps_object %>% 
      phyloseq_to_deseq2(~ case_control_other_latest + sex) %>% 
      as.DGEList()
    # dge_highQC <- filterByExpr(dge, voom_matrix)
    # dge.filtered <- dge[dge_highQC,,keep.lib.sizes=FALSE] 
    
    dge.norm <- dge %>%  #dge.filtered %>% 
      calcNormFactors(method = "TMM") %>% 
      voom(design = voom_matrix, block = metadat$study, plot = TRUE, 
           save.plot = TRUE, normalize.method="none")
    print(dim(t(dge.norm$E)))
    otu_table(ps_object) <- otu_table(dge.norm$E, taxa_are_rows = T)
    
    dat.pd <-
      subset_samples(ps_object,
                     case_control_other_latest == "Case") %>%
      abundances()
    dat.cont <-
      subset_samples(ps_object,
                     case_control_other_latest == "Control") %>%
      abundances()
    
    require(foreach)
    require(doParallel)
    require(progress)
    
    iterations <- length(rownames(dat.pd))
    pb <- progress_bar$new(
      format = "  calculating auroc [:bar] :current/:total :percent eta: :eta",
      total = iterations,
      width = 60)
    
    # allows progress bar to print within foreach
    progress_report <- function() {
      pb$tick()
    }
    opts <- list(progress = progress_report)

    cores = detectCores()
    cl <- makeCluster(cores[1] - 2)
    registerDoParallel(cl)

    roc2add <-
      foreach(
        feat = rownames(dat.pd),
        .combine = 'rbind',
        .packages = c('magrittr', 'pROC'),
        .options.snow = opts
      ) %dopar% {
        
        x <- dat.pd[feat, ] %>% t()
        y <- dat.cont[feat, ] %>% t()
        
        # AUROC
        rocdata <-
          c(roc(
            controls = y,
            cases = x,
            direction = '<',
            ci = TRUE,
            auc = TRUE
          )$ci)
        data.frame(
          "feature" = feat,
          "data_level" = level,
          "ref_DB" = ref_DB,
          "ci_lower" = rocdata[1],
          "auroc" = rocdata[2],
          "ci_upper" = rocdata[3]
        )
      }
    
    stopCluster(cl)
    aucs <- rbind(aucs, roc2add)
  }
}


end_time <- Sys.time()
cat("AUROCs calculated in : ",
    end_time - start_time, attr(end_time - start_time, "units"), "\n")

dir.create(file.path("data/feature_AUROC"), showWarnings = FALSE)
save(aucs, file = paste0("data/Analyses/feature_AUROC/feature_AUROCs.RData"))
openxlsx::write.xlsx(aucs, 
                     file = 'data/Analyses/feature_AUROC/feature_AUROCs.xlsx', 
                     overwrite = T)




for (ref_DB in names(refDB_phyloseq_orm)){
  for (level in names(refDB_phyloseq_orm[[ref_DB]])) {
    
    #_____________________________________________________________________________
    #####                    Limma - VOOM Normalization                     ##### 
    #_____________________________________________________________________________
    ps_object <- refDB_phyloseq_orm[[ref_DB]][[level]]
    metadat <- meta(ps_object)
    print(ps_object)
    
    voom_matrix <-
      model.matrix( ~ 0 + case_control_other_latest + sex,
                    data = metadat)
    
    dge <- 
      ps_object %>% 
      phyloseq_to_deseq2(~ case_control_other_latest + sex) %>% 
      as.DGEList()
    # dge_highQC <- filterByExpr(dge, voom_matrix)
    # dge.filtered <- dge[dge_highQC,,keep.lib.sizes=FALSE] 
    
    dge.norm <- dge %>%  #dge.filtered %>% 
      calcNormFactors(method = "TMM") %>% 
      voom(design = voom_matrix, block = metadat$study, plot = TRUE, 
           save.plot = TRUE, normalize.method="none")
    print(dim(t(dge.norm$E)))
    otu_table(ps_object) <- otu_table(dge.norm$E, taxa_are_rows = T)
    
    dat.male <-
      subset_samples(ps_object,
                     sex == "Male") %>%
      abundances()
    dat.female <-
      subset_samples(ps_object,
                     sex == "Female") %>%
      abundances()
    
    require(foreach)
    require(doParallel)
    require(progress)
    
    iterations <- length(rownames(dat.male))
    pb <- progress_bar$new(
      format = "  calculating auroc [:bar] :current/:total :percent eta: :eta",
      total = iterations,
      width = 60)
    
    # allows progress bar to print within foreach
    progress_report <- function() {
      pb$tick()
    }
    opts <- list(progress = progress_report)
    
    cores = detectCores()
    cl <- makeCluster(cores[1] - 2)
    registerDoParallel(cl)
    
    roc2add <-
      foreach(
        feat = rownames(dat.male),
        .combine = 'rbind',
        .packages = c('magrittr', 'pROC'),
        .options.snow = opts
      ) %dopar% {

        x <- dat.male[feat, ] %>% t()
        y <- dat.female[feat, ] %>% t()
        
        # AUROC
        rocdata <-
          c(roc(
            controls = y,
            cases = x,
            direction = '<',
            ci = TRUE,
            auc = TRUE
          )$ci)
        data.frame(
          "feature" = feat,
          "data_level" = level,
          "ref_DB" = ref_DB,
          "ci_lower" = rocdata[1],
          "auroc" = rocdata[2],
          "ci_upper" = rocdata[3]
        )
      }
    
    stopCluster(cl)
    aucs <- rbind(aucs, roc2add)
  }
}


end_time <- Sys.time()
cat("AUROCs calculated in : ",
    end_time - start_time, attr(end_time - start_time, "units"), "\n")

dir.create(file.path("data/feature_AUROC"), showWarnings = FALSE)
save(aucs, file = paste0("data/Analyses/feature_AUROC/feature_AUROCs_sex.RData"))
openxlsx::write.xlsx(aucs, 
                     file = 'data/Analyses/feature_AUROC/feature_AUROCs_sex.xlsx', 
                     overwrite = T)


