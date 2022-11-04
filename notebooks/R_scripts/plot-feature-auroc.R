
source("src/_load_packages.R")
source("src/_plot-functions.R")

base::load("data/Analyses/feature_AUROC/feature_AUROCs.RData")

aurocs_map <- aucs %>%
  dplyr::select(feature, data_level, ref_DB) %>%
  distinct()

aucs_shared_loose <- aucs %>%
  dplyr::mutate(auroc_center = auroc - 0.5) %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(mean_auroc = mean(auroc),
                   mean_auroc_cent = mean(auroc_center),
                   median = median(auroc),
                   n = n()) %>%
  left_join(aurocs_map, by = "feature") %>%
  dplyr::group_by(data_level, ref_DB) %>%
  filter(mean_auroc_cent > 0.1) %>%
  slice_max(order_by = abs(mean_auroc_cent), n = 50)


aucplot_errorbars_loose <- aucs %>%
  filter(feature %in% aucs_shared_loose$feature) %>%
  # filter(data_level == "Species") %>%
  ggplot(aes(x=fct_reorder(feature, -auroc), y = auroc, fill = ref_DB,
             ymin=ci_lower, ymax=ci_upper, group = feature)) +
  geom_pointrange(aes(group = feature),
                  position = position_jitterdodge(jitter.height = 0),
                  shape = 24, stroke = 0.2, colour="grey") +
  facet_wrap(ref_DB~data_level, scales = "free_y") +
  geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.8) +
  scale_fill_d3() +
  # theme(axis.text.y = element_blank()) +
  coord_flip() +
  clean_theme()
aucplot_errorbars_loose
# ggsave(aucplot_errorbars_loose,
#        filename = "figures/differential_abundance/feature_AUROCs_top_50_facets.svg",
#        width = 15, height = 8, dpi = 600)
