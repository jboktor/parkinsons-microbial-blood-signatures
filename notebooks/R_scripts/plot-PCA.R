# Joe Boktor
# Caltech - Mazmanian Lab

#______________________________________________________________________________
#                    Function to visualize top PERMANOVA hits in PCA
#______________________________________________________________________________

permanova_pca_vis <- function(permanova_df, pca_df, fileName) {
  
  
  top_permanova_hits <- species_permanova_df %>%
    filter(FDR <= 0.05) %>%
    filter(R2_perc > 1)
  
  fileOut <- paste0(fileName, Sys.Date(), ".pdf")
  grDevices::cairo_pdf(file = fileOut, onefile = TRUE, width = 8, height = 6)
  
  for (metahit in top_permanova_hits$metadata){
    
    if ( is.numeric(metadat[[metahit]]) ){
      pca_plot <- pca_df %>%
        ggbiplot::ggbiplot(choices = c(1,2), #obs.scale = 1, var.scale = 1,
                           groups =  metadat[[metahit]],
                           scale = 0,
                           ellipse = F,
                           alpha = 0.6,
                           circle = TRUE,var.axes=FALSE) +
        scale_color_viridis() +
        labs(title = metahit, color = "") +
        clean_theme()
      print(pca_plot)
    } else {
      
      pca_plot <- pca_df %>%
        ggbiplot::ggbiplot(choices = c(1,2),
                           groups =  metadat[[metahit]],
                           # scale = 0, 
                           obs.scale = 1, var.scale = 1,
                           ellipse = TRUE,
                           alpha = 0.6,
                           circle = TRUE,var.axes=FALSE) +
        scale_color_aaas() +
        labs(title = metahit, color = "") +
        clean_theme()
      print(pca_plot)
    }
  }
  dev.off()
}


