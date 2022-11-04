# Aitchisons Distance Analysis

####### Load PhyloSeq objects  ####### 
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")

beta_diversity_summary <- function(cohort, dist = "Aitchisons"){
  
  load_betadiv_colors() 
  
  x <- c(dat.species)
  z <- c("Species")
  
  # --------------------------------------------------- 
  # Aitchisons PCoA/Ridgeline/Violin Plots Loop 
  # --------------------------------------------------- 
  
  cnt <- 1
  for (i in x){
    load("files/low_quality_samples.RData")
    i <- i %>% 
      subset_samples(donor_id %ni% low_qc[[1]])
    
    if (dist == "Aitchisons"){
      obj_dist <- 
        microbiome::transform(i, "clr")
      iDist <- phyloseq::distance(obj_dist, method="euclidean")
      dist_label <- "Aitchison"
    } else if (dist == "bray"){
      obj_dist <- microbiome::transform(i, "compositional") 
      iDist <- phyloseq::distance(obj_dist, method="bray")
      dist_label <- "Bray-Curtis"
    } else if (dist == "jaccard"){
      obj_dist <- microbiome::transform(i, "compositional") 
      iDist <- phyloseq::distance(obj_dist, method="jaccard")
      dist_label <- "Jaccard"
    }
    
    # Calculate ordination
    cat("Processing", dist_label, "Distance:", z[cnt], "\n")
    iMDS  <- phyloseq::ordinate(obj_dist, "MDS", distance=iDist)
    
    # --------------------------------------------------- 
    #  PCoA for Axis 1 and 2
    # --------------------------------------------------- 
    
    p <- plot_ordination(obj_dist, iMDS, color="description", axes = c(1, 2))
    
    df12 = p$data
    df12$donor_group <- factor(df12$donor_group, levels=c("PC", "PD", "HC"))
    
    p <- ggplot(df12, aes(Axis.1, Axis.2, fill = description, color=description))
    p <- p + geom_point(shape=21, size=3, alpha=0.7)
    ord <- p + 
      theme_bw() + 
      labs(fill="Donor Group") +
      xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1])*100, digits = 2), "%)")) +
      ylab(paste0("PCoA 2 (", round((iMDS$values$Relative_eig[2])*100, digits = 2), "%)")) +
      labs(fill="Donor Group") +
      scale_fill_manual(values = cols.pdpchc.dark) +
      scale_color_manual(values = cols.pdpchc.rim) +
      theme(plot.title = element_text(hjust = 0.5), 
            panel.grid = element_blank())
    
    # PCoA data to fit ridgeline Plots
    my.ggp.xrange <- ggplot_build(ord)$layout$panel_scales_x[[1]]$range$range # For PCoA1
    my.ggp.yrange2 <- ggplot_build(ord)$layout$panel_scales_y[[1]]$range$range # For PCoA2
    
    # --------------------------------------------------- 
    #  Boxplot - Axis 1
    # --------------------------------------------------- 
    
    r1 <- ggplot(df12, aes(x = Axis.1, y = description)) +
      geom_boxplot(aes(color = description, fill = description), alpha = 0.2) +
      xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1])*100, digits = 1), "%)")) +
      scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
      scale_x_continuous(limits = my.ggp.xrange) +
      coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
      scale_color_manual(values = cols.pdpchc.dark) +
      scale_fill_manual(values = cols.pdpchc) +
      theme_classic() +
      ggtitle(paste0(dist_label, " Distance PCoA")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.line.y = element_blank(),
            axis.line.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x=element_blank(),
            legend.position = "none")
    
    # --------------------------------------------------- 
    #  Boxplot - Axis 2
    # --------------------------------------------------- 
    
    r2 <- ggplot(df12, aes(x = Axis.2, y = description)) +
      geom_boxplot(aes(color = description, fill = description), alpha = 0.2) +
      scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
      scale_x_continuous(limits = my.ggp.yrange2) +
      coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
      scale_color_manual(values = cols.pdpchc.dark) +
      scale_fill_manual(values = cols.pdpchc) +
      theme_classic() +
      theme(axis.title.y=element_blank(),
            axis.line.y = element_blank(),
            axis.line.x = element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x = element_blank(),
            plot.title=element_blank(),
            legend.position = "none") +
      coord_flip()
    
    # --------------------------------------------------- 
    #  Beta Diversity Violin Plots
    # --------------------------------------------------- 
    p <- obj_dist
    s <- "donor_id"
    d <- "description"
    
    # Make Melted Distance Matrix 
    wu.m <- melt(as.matrix(iDist))
    # Exclude intra-sample distances 
    wu.m <- wu.m %>% filter(as.character(Var1) != as.character(Var2)) %>% 
      mutate_if(is.factor, as.character)
    # Pull metadata of interest
    sd <- meta(p) %>% dplyr::select(all_of(s), all_of(d)) %>% mutate_if(is.factor, as.character)
    # Add group name for Var1 Column
    colnames(sd) <- c("Var1", "Type1")
    wu.sd <- left_join(wu.m, sd, by = "Var1")
    # Add group name for Var2 Column
    colnames(sd) <- c("Var2", "Type2")
    wu.sd <- left_join(wu.sd, sd, by = "Var2")
    # Select only distances to Population control
    wu.sd <- filter(wu.sd, Type1 == "Population Control")
    
    # Specifying comparisons for analysis
    my_comparisons <-
      list(c("Household Control", "PD Patient"),
           c("Population Control", "PD Patient"))
    wu.sd$Type2 <-
      factor(wu.sd$Type2,
             levels = c("Population Control", "PD Patient", "Household Control"))
    
    if (cohort == "Merged"){
      beesize <- 0.2
      beescale <- 0.2
    } else {
      beesize <- 0.3
      beescale <- 0.5
    }
    
    v <- ggplot(wu.sd, aes(x = Type2, y = value)) + theme_minimal() + 
      geom_beeswarm(aes(color = Type2, fill = Type2), shape = 21, size= beesize, alpha = 0.5, cex = beescale) +
      geom_violin(aes(color = Type2), draw_quantiles = c(0.5), trim = T, alpha=0) +
      theme_classic() +
      ylab(paste(dist_label, "Distance")) +
      labs(fill="Group") +
      scale_color_manual(values = cols.pdpchc) +
      scale_fill_manual(values = cols.pdpchc) +
      scale_x_discrete(labels = c("Population Control" = "Population \nControl",
                                  "PD Patient" = "Parkinson's \nDisease", 
                                  "Household Control" = "Household \nControl")) +
      stat_compare_means(comparisons = my_comparisons, label = "p.signif", tip.length = 0.02, step.increase = 0) + 
      ggtitle(paste0(" Distance to Population Control\n", z[cnt], " Abundance")) +
      theme(plot.title = element_text(hjust = 0.5),
            panel.grid = element_blank(),
            axis.title.x=element_blank(),
            axis.text.x = element_text(angle=45, hjust=1),
            legend.position = "none")
    
    
    # --------------------------------------------------- 
    #  Assemble Plots
    # --------------------------------------------------- 
    cat(paste0("Assembling summary plots for: " , z[cnt], "\n"))
    ord.plot <- ord + theme(legend.position = "none")
    cow1 <- cowplot::plot_grid(r1, NULL, ord.plot, r2, nrow = 2, ncol = 2, rel_heights = c(1, 3.25), rel_widths = c(5, 1),
                               align = "vh", axis="tblr")
    cow2 <- cowplot::plot_grid(v, cow1, nrow = 1, rel_widths = c(1, 2.75), align = "h")
    print(cow2)
    
    ggsave(plot = cow2, 
           filename = paste0(
             "data/Community_Composition/Beta_Diversity_Analysis/",
             z[cnt], "_",
             cohort, "_",
             dist_label,
             ".svg"
           ),
           width = 8, height = 6)
    
    cnt <- cnt + 1
  }
}



dysbiosis_score <- function(obj, dist = "Aitchisons"){
  
  #' Function to calculate dysbiosis via the median distance to 
  #' population control samples 
  obj <- obj %>% 
    subset_samples(donor_id %ni% low_qc[[1]])
  
  if (dist == "Aitchisons"){
    obj_dist <- microbiome::transform(obj, "clr") 
    dist <- phyloseq::distance(obj_dist, method="euclidean")
  } else if (dist == "bray"){
    obj_dist <- microbiome::transform(obj, "compositional")
    dist <- phyloseq::distance(obj_dist, method="bray")
  }
  
  dist_prep <-
    dist %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column("donor_A") %>% 
    pivot_longer(!donor_A, names_to = "donor_B", values_to = "distance")
  
  dist_df <- 
    dist_prep %>% 
    group_col_from_ids(ids = dist_prep$donor_A) %>% 
    dplyr::rename(group_A = group) %>% 
    group_col_from_ids(ids = dist_prep$donor_B) %>%
    dplyr::rename(group_B = group) %>% 
    # Exclude same-sample distances
    filter(as.character(donor_A) != as.character(donor_B)) %>% 
    # Select distance to healthy population controls
    # filter(group_A == "PC")
    filter(group_A == "HC")
  
  #' Calculating Dysbiosis score:
  #' Determine median value of a given sample to all 
  #' samples within the PC group
  dysbiosis_score <-
    dist_df %>% 
    group_by(donor_B) %>% 
    dplyr::summarise(median = median(distance, na.rm = TRUE))
  dysbiosis_score <- 
    group_col_from_ids(dysbiosis_score, ids = dysbiosis_score$donor_B) %>% 
    dplyr::rename(donor_id = donor_B)
  return(dysbiosis_score)
}







# 
# # TROUBLESHOOTING QC
# # load_all_cohorts()
# 
# cohort = "Merged"
# dist = "Aitchisons"
# z <- c("Species", "Pathways", "Enzymes", "KOs", 
#        # "Eggnogs", "Pfams",
#        "Pathways.slim", "Enzymes.slim", "KOs.slim", "Eggnogs.slim", "Pfams.slim")
# i <- dat.EGGNOGs.slim
# cnt <- 8
# load_betadiv_colors() 
# 
# 
# load("files/low_quality_samples.RData")
# i <- i %>%
#   subset_samples(donor_id %ni% low_qc[[1]])
# 
# if (dist == "Aitchisons"){
#   obj_dist <- 
#     microbiome::transform(i, "clr")
#     # microbiome::transform(i, "compositional") %>%
#     # microbiome::transform("clr")
#   iDist <- phyloseq::distance(obj_dist, method="euclidean")
#   dist_label <- "Aitchison"
# } else if (dist == "bray"){
#   obj_dist <- microbiome::transform(i, "compositional") 
#   iDist <- phyloseq::distance(obj_dist, method="bray")
#   dist_label <- "Bray-Curtis"
# } else if (dist == "jaccard"){
#   obj_dist <- microbiome::transform(i, "compositional") 
#   iDist <- phyloseq::distance(obj_dist, method="jaccard")
#   dist_label <- "Jaccard"
# }
# 
# # Calculate ordination
# cat("Processing", dist_label, "Distance:", z[cnt], "\n")
# iMDS  <- phyloseq::ordinate(obj_dist, "MDS", distance=iDist)
# 
# # --------------------------------------------------- 
# #  PCoA for Axis 1 and 2
# # --------------------------------------------------- 
# 
# plot_ordination(obj_dist, iMDS, color="cohort", axes = c(1, 2))
# p <- plot_ordination(obj_dist, iMDS, color="description", axes = c(1, 2))
# 
# df12 = p$data
# df12$donor_group <- factor(df12$donor_group, levels=c("PC", "PD", "HC"))
# p <- ggplot(df12, aes(Axis.1, Axis.2, fill = description, color=description))
# p <- p + geom_point(shape=21, size=5, alpha=0.7)
# ord <- p + 
#   theme_bw() + 
#   labs(fill="Donor Group") +
#   xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1])*100, digits = 2), "%)")) +
#   ylab(paste0("PCoA 2 (", round((iMDS$values$Relative_eig[2])*100, digits = 2), "%)")) +
#   labs(fill="Donor Group") +
#   scale_fill_manual(values = cols.pdpchc.dark) +
#   scale_color_manual(values = cols.pdpchc.rim) +
#   theme(plot.title = element_text(hjust = 0.5), 
#         panel.grid = element_blank())
# 
# # PCoA data to fit ridgeline Plots
# my.ggp.xrange <- ggplot_build(ord)$layout$panel_scales_x[[1]]$range$range # For PCoA1
# my.ggp.yrange2 <- ggplot_build(ord)$layout$panel_scales_y[[1]]$range$range # For PCoA2
# 
# # --------------------------------------------------- 
# #  Boxplot - Axis 1
# # --------------------------------------------------- 
# 
# r1 <- ggplot(df12, aes(x = Axis.1, y = description)) +
#   geom_boxplot(aes(color = description, fill = description), alpha = 0.2) +
#   xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1])*100, digits = 1), "%)")) +
#   scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
#   scale_x_continuous(limits = my.ggp.xrange) +
#   coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
#   scale_color_manual(values = cols.pdpchc.dark) +
#   scale_fill_manual(values = cols.pdpchc) +
#   theme_classic() +
#   ggtitle(paste0(dist_label, " Distance PCoA")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.line.y = element_blank(),
#         axis.line.x = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title.x=element_blank(),
#         legend.position = "none")
# 
# # --------------------------------------------------- 
# #  Boxplot - Axis 2
# # --------------------------------------------------- 
# 
# r2 <- ggplot(df12, aes(x = Axis.2, y = description)) +
#   geom_boxplot(aes(color = description, fill = description), alpha = 0.2) +
#   scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
#   scale_x_continuous(limits = my.ggp.yrange2) +
#   coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
#   scale_color_manual(values = cols.pdpchc.dark) +
#   scale_fill_manual(values = cols.pdpchc) +
#   theme_classic() +
#   theme(axis.title.y=element_blank(),
#         axis.line.y = element_blank(),
#         axis.line.x = element_blank(),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x = element_blank(),
#         plot.title=element_blank(),
#         legend.position = "none") +
#   coord_flip()
# 
# # --------------------------------------------------- 
# #  Beta Diversity Violin Plots
# # --------------------------------------------------- 
# p <- obj_dist
# m <- "euclidean"
# s <- "donor_id"
# d <- "description"
# require("phyloseq");require("dplyr");require("reshape2");require("ggplot2")
# 
# # Make Melted Distance Matrix 
# wu.m <- melt(as.matrix(iDist))
# # Exclude intra-sample distances 
# wu.m <- wu.m %>% filter(as.character(Var1) != as.character(Var2)) %>% 
#   mutate_if(is.factor, as.character)
# # Pull metadata of interest
# sd <- meta(p) %>% dplyr::select(all_of(s), all_of(d)) %>% mutate_if(is.factor, as.character)
# # Add group name for Var1 Column
# colnames(sd) <- c("Var1", "Type1")
# wu.sd <- left_join(wu.m, sd, by = "Var1")
# # Add group name for Var2 Column
# colnames(sd) <- c("Var2", "Type2")
# wu.sd <- left_join(wu.sd, sd, by = "Var2")
# # Select only distances to Population control
# wu.sd <- filter(wu.sd, Type1 == "Population Control")
# 
# # Specifying comparisons for analysis
# my_comparisons <-
#   list(c("Household Control", "PD Patient"),
#        c("Population Control", "PD Patient"))
# wu.sd$Type2 <-
#   factor(wu.sd$Type2,
#          levels = c("Population Control", "PD Patient", "Household Control"))
# 
# v <- ggplot(wu.sd, aes(x = Type2, y = value)) + theme_minimal() + 
#   # geom_beeswarm(aes(color = Type2, fill = Type2), shape = 21, size=.3, alpha = .5, cex = .5) +
#   geom_violin(aes(color = Type2), draw_quantiles = c(0.5), trim = T, alpha=0) +
#   theme(axis.title.x=element_blank(),
#         legend.position = "none") +
#   ylab(paste(dist_label, "Distance")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   labs(fill="Group") +
#   scale_color_manual(values = cols.pdpchc) +
#   scale_fill_manual(values = cols.pdpchc) +
#   stat_compare_means(comparisons = my_comparisons, label = "p.signif", tip.length = 0.02, step.increase = 0) + 
#   ggtitle(paste0(dist_label," Distance to Population Control\n", z[cnt], " abundance"))
# 
# # --------------------------------------------------- 
# #  Assemble Plots
# # --------------------------------------------------- 
# cat(paste0("Assembling summary plots for: " , z[cnt], "\n"))
# ord.plot <- ord + theme(legend.position = "none")
# cow1 <- cowplot::plot_grid(r1, NULL, ord.plot, r2, nrow = 2, ncol = 2, rel_heights = c(1, 5), rel_widths = c(5, 1),  align = "vh")
# cow2 <- cowplot::plot_grid(v, cow1, nrow = 1, rel_widths = c(1, 2.75), align = "h")
# print(cow2)
# cow2
# # ggsave(plot = cow2, 
# #        filename = paste0(
# #          "data/Community_Composition/Beta_Diversity_Analysis/",
# #          z[cnt], "_",
# #          cohort, "_",
# #          dist_label,
# #          ".svg"
# #        ),
# #        width = 13, height = 9)

