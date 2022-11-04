# Joe Boktor
# Caltech - Mazmanian Lab

source("src/_load_packages.R")
source("src/_plot-functions.R")
source("src/_misc_functions.R")

RefSeqPlusPF_reports_long <-
  readRDS("data/Kraken_reports/RefSeqPlusPF_reports_merged.rds") %>%
  slim_report_df() %>% 
  mutate(ref_DB = "RefSeqPlusPF")
UHGG_reports_long <-
  readRDS("data/Kraken_reports/UHGG_reports_merged.rds") %>%
  slim_report_df() %>% 
  mutate(ref_DB = "UHGG")
WoL_reports_long <-
  readRDS("data/Kraken_reports/WoL_reports_merged.rds") %>%
  slim_report_df() %>% 
  mutate(ref_DB = "WoL")
dat.species <- readRDS("data/Phyloseq_Objects/RefSeqPlusPF/Species_counts.rds")
meta_df <- meta(dat.species)



reports_all <- bind_rows(RefSeqPlusPF_reports_long,
                         WoL_reports_long,
                         UHGG_reports_long)


report_summary <- 
  reports_all %>% 
  filter(NCBI_taxon_ID %in% c(0, 1, 2, 4751, 10239)) %>% 
  ggplot(aes(x = ref_DB, y = count)) + 
  geom_point(aes(color = ref_DB), position = position_jitterdodge(), alpha = 0.6) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.4) +
  scale_y_log10() +
  facet_wrap(~NCBI_taxon_ID, nrow = 1, #scales = "free_y",
             labeller = labeller(
               NCBI_taxon_ID = c(
                 "0" = "Unclassified",
                 "1" = "Root",
                 "2" = "Bacteria",
                 "2157" = "Archaea",
                 "4751" = "Fungi",
                 "10239" = "Viruses"))) +
  scale_color_uchicago() + 
  clean_theme() +
  labs(x = NULL, y = "read counts", color = "Reference DB") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
report_summary

ggsave(report_summary,
       filename = "figures/mapping_stats/ref_DB_mapping_comparisons.png",
       width = 8.5, height = 3, dpi = 600)


report_summary_pd <- 
  reports_all %>% 
  filter(NCBI_taxon_ID %in% c(0, 1, 2, 4751, 10239)) %>% 
  left_join(meta_df, by = "participant_id") %>% 
  filter(case_control_other_latest != "Other") %>% 
  ggplot(aes(x = case_control_other_latest, y = count)) + 
  geom_point(aes(color = case_control_other_latest), position = position_jitterdodge(), alpha = 0.6) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.4) +
  scale_y_log10() +
  facet_grid(NCBI_taxon_ID~ref_DB, scales = "free",
             labeller = labeller(
               NCBI_taxon_ID = c(
                 "0" = "Unclassified",
                 "1" = "Root",
                 "2" = "Bacteria",
                 "4751" = "Fungi",
                 "10239" = "Viruses"))) +
  scale_color_aaas() + 
  stat_compare_means(method = "wilcox.test", label.y.npc = 0.95) +
  clean_theme() +
  labs(x = NULL, y = "read counts", color = "Group") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
report_summary_pd
ggsave(report_summary_pd,
       filename = "figures/mapping_stats/ref_DB_mapping_comparisons_PD.png",
       width = 7.5, height = 7.5, dpi = 600)


report_summary_sex <- 
  reports_all %>% 
  filter(NCBI_taxon_ID %in% c(0, 1, 2, 4751, 10239)) %>% 
  left_join(meta_df, by = "participant_id") %>% 
  ggplot(aes(x = sex, y = count)) + 
  geom_point(aes(color = sex), position = position_jitterdodge(), alpha = 0.6) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.4) +
  scale_y_log10() +
  facet_grid(NCBI_taxon_ID~ref_DB, scales = "free",
             labeller = labeller(
               NCBI_taxon_ID = c(
                 "0" = "Unclassified",
                 "1" = "Root",
                 "2" = "Bacteria",
                 "4751" = "Fungi",
                 "10239" = "Viruses"))) +
  scale_color_aaas() + 
  stat_compare_means(method = "wilcox.test", label.y.npc = 0.95) +
  clean_theme() +
  labs(x = NULL, y = "read counts", color = "Sex") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
report_summary_sex
ggsave(report_summary_sex,
       filename = "figures/mapping_stats/ref_DB_mapping_comparisons_sex.png",
       width = 7.5, height = 7.5, dpi = 600)

report_summary_study <- 
  reports_all %>% 
  filter(NCBI_taxon_ID %in% c(0, 1, 2, 4751, 10239)) %>% 
  left_join(meta_df, by = "participant_id") %>% 
  ggplot(aes(x = study, y = count)) + 
  geom_point(aes(color = study), position = position_jitterdodge(), alpha = 0.6) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.4) +
  scale_y_log10() +
  facet_grid(NCBI_taxon_ID~ref_DB, scales = "free",
             labeller = labeller(
               NCBI_taxon_ID = c(
                 "0" = "Unclassified",
                 "1" = "Root",
                 "2" = "Bacteria",
                 "4751" = "Fungi",
                 "10239" = "Viruses"))) +
  scale_color_aaas() + 
  stat_compare_means(method = "anova", label.y.npc = 0.95, label.x.npc = 0.1) +
  clean_theme() +
  labs(x = NULL, y = "read counts", color = "Study") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
report_summary_study
ggsave(report_summary_study,
       filename = "figures/mapping_stats/ref_DB_mapping_comparisons_study.png",
       width = 8, height = 8, dpi = 600)


report_summary_root_reads <- 
  reports_all %>% 
  filter(NCBI_taxon_ID == 1) %>% 
  left_join(meta_df, by = "participant_id") %>% 
  ggplot(aes(x = study, y = count)) + 
  geom_point(aes(color = study), position = position_jitterdodge(), alpha = 0.6) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.4) +
  scale_y_log10() +
  scale_color_aaas() + 
  stat_compare_means(method = "anova", label.y.npc = 0.95, label.x.npc = 0.1) +
  clean_theme() +
  labs(x = NULL, y = "read counts", color = "Study") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
report_summary_root_reads


ecdf_plot_study <- 
  reports_all %>% 
  filter(NCBI_taxon_ID == 1) %>% 
  left_join(meta_df, by = "participant_id") %>% 
  ggplot(aes(x=count, colour = study)) + 
  stat_ecdf(geom = "step", pad = FALSE, size = 1) + 
  scale_x_log10() +
  scale_color_aaas() + 
  clean_theme() +
  labs(x = "Mapped Reads", y = "eCDF", color = "Study") +
  facet_wrap(~ref_DB, nrow = 1) +
  theme(legend.position = c(0.95, 0.5))
ecdf_plot_study

ecdf_plot_study_fig <- 
  reports_all %>% 
  filter(NCBI_taxon_ID == 1) %>% 
  left_join(meta_df, by = "participant_id") %>% 
  filter(ref_DB == "RefSeqPlusPF") %>% 
  ggplot(aes(x=count, colour = study)) + 
  stat_ecdf(geom = "step", pad = FALSE, size = 1) + 
  scale_x_log10() +
  scale_color_aaas() + 
  clean_theme() +
  labs(x = "Mapped Reads", y = "eCDF", color = "Study") +
  # facet_wrap(~ref_DB, nrow = 1) +
  theme(legend.position = c(0.95, 0.5))
ecdf_plot_study_fig

ggsave(ecdf_plot_study_fig,
       width = 4, height = 4,
       filename = "figures/mapping_stats/ecdf_study_RefSeq_4_figure.svg")

save_me_cleanly(ggobj = ecdf_plot_study_fig, filename = paste0("figures/mapping_stats/ecdf_study_RefSeq_4_figure.svg"),
                plot_w = 3, plot_h = 3, leg_w = 8, leg_h = 8)


#_______________________________________________________________________________
#                       eCDF Plotting Helper Function
#_______________________________________________________________________________

plot_ecdf <- function(color_var){
  reports_all %>% 
    filter(NCBI_taxon_ID == 1) %>% 
    left_join(meta_df, by = "participant_id") %>% 
    ggplot(aes(x=count, colour = !!sym(color_var))) + 
    stat_ecdf(geom = "step", pad = FALSE, size = 1) + 
    scale_x_log10() +
    scale_color_aaas() + 
    clean_theme() +
    labs(x = "Mapped Reads", y = "eCDF", color = "") +
    facet_wrap(~ref_DB, nrow = 1) +
    theme(
      # legend.position = c(0.9, 0.4), 
          legend.background = element_blank(),
          legend.text  = element_text(size = 8))
}

plotMe <- c("study", "sex", "case_control_other_latest", "diagnosis_latest", "alcohol_prior_use")
for (plot_col in plotMe){
  p <- plot_ecdf(plot_col)
  print(p)
  save_me_cleanly(ggobj = p, filename = paste0("figures/mapping_stats/ecdf_", plot_col),
                  plot_w = 7, plot_h = 3, leg_w = 8, leg_h = 8)
}



