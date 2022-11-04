# Joe Boktor
# Caltech - Mazmanian Lab
# Dec 2021

## TO DO - create log file for function
### Possibly use only one level to inform which samples to remove for all others

rm(list = ls())
source("src/_load_packages.R")
source("src/_plot-functions.R")
refDB_phyloseq <- readRDS("data/Phyloseq_Objects/Phyloseq_all.rds")

refDB_phyloseq_orm <- list()
for (refDB in c("UHGG", "WoL")){
  for (level in c("Genus", "Species")) {
    
    # # # TROUBLE
    # refDB <- "WoL"
    # level <- "Genus"

    cat("\nProcessing: ", refDB, level, "-------\n\n")
    ps_object <- refDB_phyloseq[[refDB]][[level]]
    print(ps_object)
    
    mdf <- meta(ps_object)
    counts <- ps_object %>% 
      abundances() %>% 
      as.data.frame() %>% 
      colSums()
    
    counts_df <- as.data.frame(counts) %>% 
      rownames_to_column("participant_id") %>% 
      left_join(mdf)
    
    # Saving figures in pdf
    fileName <- paste0("figures/quality_control/", refDB, "_", level)
    fileOut <- paste0(fileName, "_", Sys.Date(), ".pdf")
    grDevices::cairo_pdf(file = fileOut, onefile = TRUE, width =20, height = 14)
    
    count_eCDF <- 
      counts_df %>% 
      ggplot(aes(x=counts, color = study)) + 
      stat_ecdf(geom = "point", pad = FALSE, size = 1) + 
      scale_x_log10() +
      scale_color_aaas() + 
      clean_theme() +
      labs(x = "Mapped Reads", y = "eCDF", color = "") +
      theme(legend.background = element_blank(),
            legend.text  = element_text(size = 8))
    cowplot::plot_grid(count_eCDF, scale = 0.7)
  
    
    #_______________________________________________________________________________
    #####                    Limma - VOOM Normalization                         ##### 
    #_______________________________________________________________________________
    
    voom_matrix <-
      model.matrix( ~ 0 + case_control_other_latest + sex + study,
                    data = meta(ps_object))
    
    dge <- 
      ps_object %>% 
      phyloseq_to_deseq2(~ 0 + case_control_other_latest + sex + study) %>% 
      as.DGEList()
    dge_highQC <- filterByExpr(dge, voom_matrix)
    dge.filtered <- dge[dge_highQC, keep.lib.sizes=FALSE] 
    dge.norm <- dge.filtered %>% 
      calcNormFactors(method = "TMM") %>% 
      voom(design = voom_matrix, plot = TRUE, save.plot = TRUE, normalize.method="none")
    
    # Principal Components Analysis
    dat.pca <- prcomp(t(dge.norm$E), center = T, scale. = T)
    df_plot_pca <- dat.pca$x %>%
      as.data.frame() %>%
      rownames_to_column(var = "participant_id") %>%
      left_join(counts_df)
    
    seq_depth_reg <- 
      df_plot_pca %>%
        ggplot(aes(x = PC1, y = counts)) +
        geom_point() +
        geom_smooth(method = "lm", se = F) +
        
        df_plot_pca %>%
        ggplot(aes(x = PC2, y = counts)) +
        geom_point() +
        geom_smooth(method = "lm", se = F)
    print(seq_depth_reg)
    
    
    #_______________________________________________________________________________
    
    dat.pca.input <- 
      dat.pca$x %>% 
      as.data.frame() %>% 
      dplyr::select(1:4)
    
    # Basic Outlier test on 1O Principle Components
    #    defined as > 6 standard deviations from mean
    six.std.dev <- apply(dat.pca.input, 2, function(x) which( (abs(x - mean(x)) / sd(x)) > 6 )) %>%
      Reduce(union, .) 
    # Mahalanobis distance estimation of outliers ---
    mahalanobis_distance <- bigutilsr::dist_ogk(as.matrix(dat.pca.input))
    # Local Outlier Factor (LOF) ----
    llof <- LOF(dat.pca.input)  # log(LOF) by default
    
    
    
    # Outlier detection threshold (upper) based on Tukey's rule
    # adjusting for skewness and multiple testing
    basic.out <- (1:nrow(dat.pca.input) %in% six.std.dev) 
    cat("Outliers detected via Std Dev Threshold: ", sum(basic.out), "\n")
    maha.out <- (mahalanobis_distance > tukey_mc_up(mahalanobis_distance))
    cat("Outliers detected via Mahalanobis Distance: ", sum(maha.out), "\n")
    lof.out <- (llof > tukey_mc_up(llof, coef = 1.5))
    cat("Outliers detected via Local Outlier Factor Analysis: ", sum(lof.out), "\n")
    total_outliers <- as.logical(maha.out + lof.out + basic.out)
    cat("Total outliers detected: ", sum(total_outliers), "\n")
    
    #_______________________________________________________________________________
    #                  
    mahal_llof <- 
      qplot(mahalanobis_distance, llof) +
      geom_vline(xintercept = tukey_mc_up(mahalanobis_distance, coef = 1.5), color = "red") +
      geom_hline(yintercept = tukey_mc_up(llof, coef = 1.5),  color = "red")
    print(cowplot::plot_grid(mahal_llof, scale = 0.7))
    
    
    outlier_summary1 <- 
      patchwork::wrap_plots(
        qplot(dat.pca.input[, 1], dat.pca.input[, 2], 
              color = mahalanobis_distance, size = I(3)) +
          scale_color_viridis_c(trans = "log") +
          labs(x = "PC1", y = "PC2", color = "Mahalanobis\nDistance", alpha = "Outlier") + 
          clean_theme(),
        qplot(dat.pca.input[, 1], dat.pca.input[, 2], 
              color = total_outliers, size = I(3), alpha = total_outliers) +
          scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
          scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.2)) +
          labs(x = "PC1", y = "PC2", color = "Outlier", alpha = "Outlier") +
          clean_theme(),
        qplot(dat.pca.input[, 3], dat.pca.input[, 4],
              color = mahalanobis_distance, size = I(3)) +
          scale_color_viridis_c(trans = "log") +
          labs(x = "PC3", y = "PC4", color = "Mahalanobis\nDistance", alpha = "Outlier") + 
          clean_theme(),
        qplot(dat.pca.input[, 3], dat.pca.input[, 4],
              color = total_outliers, size = I(3), alpha = total_outliers) +
          scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
          scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.2)) +
          labs(x = "PC3", y = "PC4", color = "Outlier", alpha = "Outlier") + 
          clean_theme() #,
        # scale = 0.5
      )
    print(outlier_summary1)
    
    features <- t(dge.norm$E)[!total_outliers, ] %>%
      as.data.frame()
    filtered.pca <- prcomp(features)
    
    # Visualize data post-outlier removal
    outlier_summary2 <- 
      cowplot::plot_grid(
        
        dat.pca %>%
          ggbiplot::ggbiplot(
            choices = c(1, 2),
            alpha = 0.4,
            scale = 1,
            circle = F,
            var.axes = FALSE
          ) +
          labs(title = "Voom Normalized") +
          clean_theme(),
        
        filtered.pca %>%
          ggbiplot::ggbiplot(
            choices = c(1, 2),
            alpha = 0.4,
            scale = 1,
            circle = F,
            var.axes = FALSE
          ) +
          scale_color_aaas() +
          labs(title = "Voom Normalized + Outliers Removed") +
          clean_theme(), 
        ncol = 2, scale = 0.8, align = "hv"
        )
    print(outlier_summary2)
    
    dev.off()
    #_______________________________________________________________________________
    #                   Saving cleaned phyloseq objects
    #_______________________________________________________________________________
    
    # Outliers
    names(total_outliers) <- names(maha.out)
    
    dat.cleaned <- ps_object %>%
      subset_samples(!total_outliers)
    
    physeq_meta <- dat.cleaned %>%
      meta() %>% sample_data()
    
    features <- dat.cleaned %>%
      otu_table(taxa_are_rows = T)
    
    refDB_phyloseq_orm[[refDB]][[level]] <- phyloseq(features, physeq_meta)
  }
}

saveRDS(refDB_phyloseq_orm,
     file = "data/Phyloseq_Objects/Phyloseq_all_outliers_removed.rds")


