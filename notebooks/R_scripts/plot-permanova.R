# Joe Boktor
# Caltech - Mazmanian Lab

source("src/_plot-functions.R")
base::load("data/Analyses/community_composition/PERMANOVA_UHGG_Genus.RData")
# permanova_df <- species_permanova_df # Delete me later
refDB <- "UHGG"
level <- "Genus"

#______________________________________________________________________________
#                    Visualize top PERMANOVA hits
#______________________________________________________________________________


permanova_df %>% 
  filter(FDR < 0.25) %>%
  filter(R2_perc < 100) %>%
  ggplot(aes(x = F.Model, y = R2_perc)) +
  geom_point(aes(size = -log10(FDR), fill = metadata_class), shape = 21, alpha = 0.8) +
  scale_fill_brewer(palette = "Set1") + 
  scale_x_log10() + scale_y_log10() +
  labs(size = expression('-log'[10]*'(FDR)'), y = expression(R^"2"~"(%)"),
       fill = "Metadata Category") +
  clean_theme() + legend_plus_fill(3)

permanova_df %>% 
  filter(FDR < 0.25) %>%
  filter(R2_perc < 100) %>%
  ggplot(aes(x = F.Model, y = -log10(FDR))) +
  geom_point(aes(size = R2_perc, fill = metadata_class), shape = 21, alpha = 0.8) +
  scale_fill_brewer(palette = "Set1") + 
  scale_x_log10() + #scale_y_log10() +
  # labs(size = expression('-log'[10]*'(FDR)'), y = expression(R^"2"~"(%)"),
  #      fill = "Metadata Category") +
  clean_theme() + legend_plus_fill(3)



perm1 <- permanova_df %>% 
  filter(FDR < 0.25) %>%
  ggplot(aes(x= R2_perc, y = -log10(FDR))) +
  geom_point(aes(size = n, fill = metadata_class), shape = 21, alpha = 0.8) +
  scale_fill_brewer(palette = "Set1") + 
  scale_x_log10() +
  labs(x = expression(R^"2"~"(%)"), y = expression('-log'[10]*'(FDR)'),
       fill = "Metadata Category") +
  clean_theme() + legend_plus_fill(3)
perm1

perm2 <- permanova_df %>% 
  filter(FDR < 0.25) %>%
  ggplot(aes(x = n, y = R2_perc)) +
  geom_point(aes(size = -log10(FDR), fill = metadata_class), shape = 21, alpha = 0.8) +
  scale_fill_brewer(palette = "Set1") + 
  scale_x_log10() +
  labs(size = expression('-log'[10]*'(FDR)'), y = expression(R^"2"~"(%)"),
       fill = "Metadata Category") +
  clean_theme() + legend_plus_fill(3)
perm2

perm3 <- permanova_df %>% 
  filter(FDR <= 0.05) %>%
  filter(R2_perc >= 1) %>%
  slice_max(n = 20, order_by = R2_perc, with_ties = F) %>% 
  ggplot(aes(x= R2_perc, y = reorder(metadata, R2))) +
  geom_segment(aes(x=0, xend=R2_perc, y=reorder(metadata, R2), yend=reorder(metadata, R2)),
               color="grey") +
  geom_point(aes(fill = metadata_class), shape = 21, size = 4) +
  scale_fill_brewer(palette = "Set1") + 
  labs(x = expression(R^"2"~"(%)"), y = "", fill = "Metadata Category") +
  clean_theme() + legend_plus_fill(3) +
  theme(legend.position = c(0.7, 0.4))
perm3


# ggsave(perm1, filename = paste0("figures/community_composition/permanova/PERMANOVA_", refDB, "_", level, "_bubble_plot_FDRbyR.png"),
#        width = 6.5, height = 4, dpi = 600)
# ggsave(perm2, filename = paste0("figures/community_composition/permanova/PERMANOVA_", refDB, "_", level, "_bubble_plot_RbyN.png"),
#        width = 6.5, height = 4, dpi = 600)
# ggsave(perm3, filename = paste0("figures/community_composition/permanova/PERMANOVA_", refDB, 
#                                 "_", level, "barplot_R2above1perc_sig.svg"),
#        width = 6, height = 3, dpi = 600)



