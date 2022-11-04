# Joe Boktor
# Caltech - Mazmanian Lab

#------------------------------------------------------------------------------
#                       Aesthetic Variables
#------------------------------------------------------------------------------

cco_pal <- c("Case" = "red", "Control" = "grey60")

#------------------------------------------------------------------------------

clean_theme <- function() {
  th <- ggplot2::theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5),
      # plot.margin = unit(c(1, 1, 1, 2), "cm")
    )
  return(th)
}


legend_plus_fill <- function(legsize) {
  th <- guides(fill = guide_legend(override.aes = list(size=legsize)))
  return(th)
}
legend_plus_color <- function(legsize) {
  th <- guides(color = guide_legend(override.aes = list(size=legsize)))
  return(th)
}

#_______________________________________________________________________________
# Save legend and plot as separate files to fix width:
save_me_cleanly <- function(ggobj, filename, plot_w, plot_h, leg_w, leg_h, 
                            res = 600, filetype = ".png"){
  
  require(cowplot)
  require(ggpubr)
  
  my_legend <- get_legend(ggobj) %>% as_ggplot()
  ggobj_out <- ggobj + theme(legend.position = "none")
  plot_filename <- paste0(filename, filetype)
  legend_filename <- paste0(filename, "__Legend", filetype)
  
  ggsave(ggobj_out, filename = plot_filename, width = plot_w, height = plot_h, dpi = res)
  ggsave(my_legend, filename = legend_filename, width = leg_w, height = leg_h, dpi = res)
}

#_______________________________________________________________________________

                         
color_loop_generator <- function(names_col){
  accent_pal <- 
    c('#7fc97f','#beaed4','#fdc086','#ffff99',
      '#386cb0','#f0027f','#bf5b17','#666666')
  color_output <- c()
  pal_length <- length(accent_pal)
  feats <- unique(names_col)
  n_cols <- length(feats)
  color_output <- c(rep(accent_pal, floor(n_cols/pal_length)), 
                    accent_pal[1:(n_cols %% pal_length)])
  names(color_output) <- feats
  
  return(color_output)
}

