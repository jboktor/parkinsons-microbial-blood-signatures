# Joe Boktor
# Caltech - Mazmanian Lab

source("src/_load_packages.R")
source("src/_misc_functions.R")
library(broom.mixed)
library(lmerTest)

###' NOTE TO SELF!!!!: 
###' There some red flags in the longitudinal metadata with 
###' duplicate/missing sample_id info for the Biofind cohort

# load Data
res_quantiseq <- readRDS(file = "data/Immune_Cell_Proportions/quanTIseq.rds")
res_epic <- readRDS(file = "data/Immune_Cell_Proportions/EPIC.rds")
res_xcell <- readRDS(file = "data/Immune_Cell_Proportions/xCell.rds")
res_mcp_counter <- readRDS(file = "data/Immune_Cell_Proportions/MCP-counter.rds")
# load Metadata
longitudinal_metadata <- readRDS(file = "data/Metadata/longitudinal_metadata.rds")
rna_sample_inv <-
  read.csv(file = "input_files/2021_v2-5release_0510/rna_sample_inventory.csv",
           stringsAsFactors = F, header = TRUE)
rna_sample_inv %>% head()

# format metadata
longitudinal_metadata %<>%
  mutate(visit_name = gsub("#2", "", visit_name)) %>%
  left_join(rna_sample_inv, by = c("participant_id", "visit_month")) %>% 
  mutate_all(na_if,"")
longitudinal_metadata %>% dim()
# longitudinal_metadata$genetic_status_enrollment %>% table()
# longitudinal_metadata$genetic_status_wgs %>% table()

longitudinal_metadata %<>%
  mutate(
    age_at_baseline_factor =
      cut(
        longitudinal_metadata$age_at_baseline,
        breaks = quantile(longitudinal_metadata$age_at_baseline, na.rm = T),
        labels = c(1, 2, 3, 4),
        include.lowest = F
      )
  )

longitudinal_metadata %>% glimpse()
longitudinal_metadata$diagnosis_type %>% table()

# longitudinal_metadata$genetic_status_wgs %>% unique()

#_______________________________________________________________________________
#                  Loop through cell types with lme4 analysis
#_______________________________________________________________________________


lmer_modeling <- function(res_df) {
  
  list_out <-  vector("list")
  for (cell in unique(res_df$cell_type)) {
    print(cell)
    
    cell_data <- res_df %>%
      pivot_longer(!cell_type, names_to = "sample_id", values_to = "score") %>% 
      filter(cell_type == cell) %>% 
      right_join(longitudinal_metadata) %>% 
      filter(case_control_other_latest %in% c("Case", "Control"))

    model <- lmer(
      score ~ 
        case_control_other_latest +
        sex +
        race +
        education_level_years +
        age_at_baseline +
        (1 | study:participant_id),
      data = cell_data
    )
    model %>% summary() %>% print()
    list_out[[cell]] <- model
    # sjPlot::plot_model(tst_model, show.values = TRUE) %>% print()
  }
  return(list_out)
}

models_epic <- lmer_modeling(res_epic)
models_quantiseq <- lmer_modeling(res_quantiseq)
models_xcell <- lmer_modeling(res_xcell)
models_mcp_counter <- lmer_modeling(res_mcp_counter)

models_epic_df <- map_df(models_epic, tidy, .id = "name") %>% mutate(method = "EPIC")
models_quantiseq_df <- map_df(models_quantiseq, tidy, .id = "name")%>% mutate(method = "QuantiSeq")
models_xcell_df <- map_df(models_xcell, tidy, .id = "name") %>% mutate(method = "xCell")
models_mcp_counter_df <- map_df(models_mcp_counter, tidy, .id = "name") %>% mutate(method = "MCP-counter")

deconv_comparison <- 
  bind_rows(models_epic_df, models_quantiseq_df, models_xcell_df, models_mcp_counter_df)

# group_by(data_class, distance) %>% 
#   mutate(p_adj = p.adjust(p_value, method = 'BH')) %>%

tst <- models_quantiseq_df %>% 
  filter(effect == "fixed", 
         p.value < 0.05,
         term == "case_control_other_latestControl") 


pd_sig_df <- deconv_comparison %>% 
  filter(effect == "fixed", 
         p.value < 0.05,
         term == "case_control_other_latestControl") 

deconv_wide <- 
  deconv_comparison %>% 
  mutate(FDR = p.adjust(p.value, method = 'BH')) %>% 
  filter(term != "(Intercept)",
         effect == "fixed",
         method == "xCell", 
         term %in% c("age_at_baseline", "case_control_other_latestControl")) %>% 
  select(-c(effect, group)) %>% 
  pivot_wider(names_from = term, names_sep = "__", 
              values_from = c(estimate, std.error, statistic, df, p.value, FDR)) 
deconv_wide_pdsig <- deconv_wide %>% 
  filter(FDR__case_control_other_latestControl <= 1e-5) %>% 
  slice_min(n = 7, order_by = FDR__case_control_other_latestControl, with_ties = F)

age_by_disease_immune <- 
  deconv_wide %>% 
  ggplot(aes(x = -estimate__case_control_other_latestControl, 
             y = estimate__age_at_baseline,
             xmin = -estimate__case_control_other_latestControl -
               std.error__case_control_other_latestControl,
             xmax = -estimate__case_control_other_latestControl +
               std.error__case_control_other_latestControl,
             ymin = estimate__age_at_baseline -
               std.error__age_at_baseline ,
             ymax = estimate__age_at_baseline +
               std.error__age_at_baseline, 
             color = statistic__case_control_other_latestControl)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  geom_errorbar() +
  geom_errorbarh() +
  geom_point() +
  geom_label_repel(data = deconv_wide_pdsig, aes(label = name), 
                   size = 3, color = "black", 
                   nudge_x = -2e-3, nudge_y = 1e-4, 
                   segment.color	 =  "#e0e0e0", 
                   segment.curvature = 0.5, 
                   direction = "y") +
  clean_theme() +
  scale_color_viridis(option = "H", direction = -1) +
  labs(x = "Disease estimate", y = "age estimate", color = "Disease\nStatistic")
age_by_disease_immune

ggsave(age_by_disease_immune, 
       filename = "figures/immune_cell_prop/age_by_disease_immune.svg", 
       width = 5.5,
       height = 4.5)


deconv_comparison %>% 
  filter(term != "(Intercept)",
         effect == "fixed",
         method == "QuantiSeq", 
         term %in% c("age_at_baseline", "sexMale", "case_control_other_latestControl")) %>%
  ggplot(aes(x = estimate,
             y = name,
             fill = term)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_bar(size = 10^-1, stat = "identity", 
           width = 0.6, color = "black") +
  theme_bw() +
  labs(y = NULL, x = "Analyte Class Mean Estimate",
       fill = "Class") +
  facet_wrap( ~ term) +
  scale_fill_colorblind() +
  clean_theme()


pd_sig_df <- deconv_comparison %>% 
  filter(effect == "fixed", 
         p.value < 0.05,
         term == "case_control_other_latestControl") 
pd_sig_immune_prop <- 
  deconv_comparison %>% 
  filter(effect == "fixed", 
       term == "case_control_other_latestControl") %>% 
  ggplot(aes(x = -estimate, -log10(p.value))) +
  geom_point(aes(size = abs(statistic), fill = name), shape = 21, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05)) +
  # scale_color_aaas() +
  geom_text_repel(data = pd_sig_df,
                   aes(x = -estimate, y = -log10(p.value), label = name), size = 3) +
  facet_wrap(~method, scales = "free") +
  labs(x = "Estimate (PD - Controls)") +
  theme_bw()
pd_sig_immune_prop

ggsave(pd_sig_immune_prop, filename = "figures/immune_cell_prop/pd_sig_immune_prop.png",
       width = 17, height = 14)




sex_sig_df <- deconv_comparison %>% 
  filter(effect == "fixed", 
         p.value <= 0.05,
         term == "sexMale") 

sex_sig_immune_prop <- 
  deconv_comparison %>% 
  filter(effect == "fixed", 
         term == "sexMale") %>% 
  ggplot(aes(x = -estimate, -log10(p.value))) + # negative of estimate is used to switch direction 
  geom_point(aes(size = abs(statistic), fill = name), shape = 21, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05)) +
  # scale_color_aaas() +
  geom_text_repel(data = sex_sig_df,
                  aes(x = -estimate, y = -log10(p.value), label = name), size = 3) +
  facet_wrap(~method, scales = "free") +
  labs(x = "Estimate (Female - Male)") +
  theme_bw()
sex_sig_immune_prop

ggsave(sex_sig_immune_prop, 
       filename = "figures/immune_cell_prop/sex_sig_immune_prop.png",
       width = 17, height = 14)



age_sig_df <- deconv_comparison %>% 
  filter(term == "age_at_baseline",
         p.value <= 0.05)

age_sig_immune_prop <- 
  deconv_comparison %>% 
  filter(effect == "fixed", 
         term == "age_at_baseline") %>% 
  ggplot(aes(x = -estimate, -log10(p.value))) + # negative of estimate is used to switch direction 
  geom_point(aes(size = abs(statistic), color = name), alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05)) +
  # scale_color_aaas() +
  geom_text_repel(data = age_sig_df,
                  aes(x = -estimate, y = -log10(p.value), label = name), size = 3) +
  facet_wrap(~method, scales = "free") +
  labs(x = "Estimate (Age Regression)") +
  theme_bw()
age_sig_immune_prop

ggsave(sex_sig_immune_prop, filename = "figures/immune_cell_prop/age_sig_immune_prop.png",
       width = 17, height = 14)


# model_stats_plot <- 
#   deconv_comparison %>% 
#   filter(effect == "fixed", 
#          # name == "Neutrophil",
#          term != "(Intercept)") %>% #str()
#   mutate(p.value = as.numeric(p.value)) %>% 
#   ggplot(aes(x = estimate, y = term, color = method)) +
#   geom_point(alpha = 0.7) + 
#   geom_errorbar(aes(xmin = estimate-std.error, xmax = estimate+std.error), width = 0.2) +
#   facet_grid(cols = vars(method), 
#              rows = vars(name), 
#              scales = "free_x") +
#   geom_ysidecol(data = deconv_comparison, aes(x = -log10(p.value)) ) +
#   scale_color_aaas() +
#   scale_yfill_gradient(low ="#FFFFFF", high = "#0000FF") +
#   theme_bw()
# model_stats_plot
# ggsave(model_stats_plot, filename = "figures/immune_cell_prop/all_model_stats_plot.png",
#        width = 14, height = 80, limitsize = FALSE)


# neut_data <- res_xcell %>%
#   pivot_longer(!cell_type, names_to = "sample_id", values_to = "score") %>% 
#   filter(cell_type == "Neutrophil") %>% 
#   right_join(longitudinal_metadata) %>% 
#   filter(case_control_other_latest %in% c("Case", "Control"))
# 
# model <- lmer(
#   score ~ 
#     case_control_other_latest +
#     sex +
#     race +
#     education_level_years +
#     age_at_baseline_factor +
#     (1 | study) +
#     (1 | participant_id),
#   data = neut_data
# )



# plot_model(tst_model, show.values = TRUE, value.offset = .3)

tst_model
summary(tst_model)
plot(tst_model)
anova(tst_model)
ranef(tst_model)
tst <- tidy(tst_model)

# multicollinearity test
# Smallest possible value is one -  if value exceeds 5-10 then colinearity exists
library(car)
vif(tst_model)


tst_model
summary(tst_model)
plot(tst_model)

testvar <- "case_control_other_latest"
res_xcell_neutrophils <- 
  res_epic %>%
  pivot_longer(!cell_type, names_to = "sample_id", values_to = "score") %>% 
  left_join(longitudinal_metadata, by = "sample_id") %>% 
  filter(case_control_other_latest %in% c("Case", "Control"),
         cell_type == "Neutrophil") %>% 
  ggplot(aes(x=jitter(visit_month, factor = 10), y=score, color=!!sym(testvar))) +
  geom_point(alpha = 1, size = 0.5) +
  facet_wrap(~age_at_baseline_factor, scales="free_y", ncol=2) +
  geom_line(aes(group = participant_id), alpha = 0.1) +
  geom_smooth(method = "loess", alpha=0.6, fill = "white") + 
  geom_ysidedensity(aes(x=stat(density)), position = "stack") +
  geom_xsidedensity(aes(y=stat(density)), position = "stack") +
  scale_color_d3() +
  labs(x = "visit month", y = "cell-type score estimation", 
       title = "EPIC Bulk-RNASeq Deconvolution - Neutrophils", 
       color = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        # legend.position = c(0.8, 0.8),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))
res_xcell_neutrophils

longitudinal_metadata$age_at_baseline_factor



res_epic %>%
  pivot_longer(!cell_type, names_to = "sample_id", values_to = "score") %>% 
  left_join(longitudinal_metadata, by = "sample_id") %>% 
  filter(case_control_other_latest %in% c("Case", "Control") ) %>%  #,
         # visit_month == 0) %>% 
  ggplot(aes(x=study, y=score, color=case_control_other_latest)) +
  geom_boxplot() +
  facet_wrap(~cell_type, scales="free_y")
  




#_______________________________________________________________________________
# Plotting ----

res_epic_panels <- 
  res_epic %>%
  pivot_longer(!cell_type, names_to = "sample_id", values_to = "score") %>% 
  left_join(longitudinal_metadata, by = "sample_id") %>% 
  filter(case_control_other_latest %in% c("Case", "Control")) %>% 
  ggplot(aes(x=jitter(visit_month, factor = 10), y=score, color=case_control_other_latest)) +
  geom_point(size = 0.3) +
  # facet_grid(~cell_type, space = "free", scales = "free") +
  facet_wrap(~cell_type, scales="free_y") +
  geom_line(aes(color = case_control_other_latest, group = participant_id), alpha = 0.05) +
  geom_smooth(method = "loess", alpha=1, fill = "white") + 
  geom_ysidedensity(aes(x=stat(density)), position = "stack") +
  geom_xsidedensity(aes(y=stat(density)), position = "stack") +
  ggside(scales = "free") +
  labs(x = "visit month", y = "cell-type score estimation", 
       title = "EPIC Bulk-RNASeq Deconvolution") +
  scale_color_d3() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom",
        panel.grid = element_blank())
res_epic_panels

ggsave(res_epic_panels, filename = "figures/immune_cell_prop/EPIC_panels.pdf", 
       width = 12, height = 16)


res_mcp_counter_panels <- 
  res_mcp_counter %>%
  pivot_longer(!cell_type, names_to = "sample_id", values_to = "score") %>% 
  left_join(longitudinal_metadata, by = "sample_id") %>% 
  filter(case_control_other_latest %in% c("Case", "Control")) %>% 
  ggplot(aes(x=jitter(visit_month, factor = 10), y=score, color=case_control_other_latest)) +
  geom_point(size = 0.3) +
  facet_wrap(~cell_type, scales="free") +
  geom_line(aes(color = case_control_other_latest, group = participant_id), alpha = 0.05) +
  geom_smooth(method = "loess", alpha=1, fill = "white") + 
  geom_ysidedensity(aes(x=stat(density)), position = "stack") +
  geom_xsidedensity(aes(y=stat(density)), position = "stack") +
  ggside(scales = "free") +
  labs(x = "visit month", y = "cell-type score estimation", 
       title = "MCP-Counter Bulk-RNASeq Deconvolution") +
  scale_color_d3() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))
res_mcp_counter_panels

ggsave(res_mcp_counter_panels, filename = "figures/immune_cell_prop/MCP-counter_panels.pdf", 
       width = 14, height = 16)


res_quantiseq_panels <- 
  res_quantiseq %>%
  pivot_longer(!cell_type, names_to = "sample_id", values_to = "score") %>% 
  left_join(longitudinal_metadata, by = "sample_id") %>% 
  filter(case_control_other_latest %in% c("Case", "Control")) %>% 
  ggplot(aes(x=jitter(visit_month, factor = 10), y=score, color=case_control_other_latest)) +
  geom_point(size = 0.3) +
  facet_wrap(~cell_type, scales="free") +
  geom_line(aes(color = case_control_other_latest, group = participant_id), alpha = 0.05) +
  geom_smooth(method = "loess", alpha=1, fill = "white") + 
  geom_ysidedensity(aes(x=stat(density)), position = "stack") +
  geom_xsidedensity(aes(y=stat(density)), position = "stack") +
  ggside(scales = "free") +
  labs(x = "visit month", y = "cell-type score estimation", 
       title = "QuantiSeq Bulk-RNASeq Deconvolution") +
  scale_color_d3() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))
res_quantiseq_panels

ggsave(res_quantiseq_panels, filename = "figures/immune_cell_prop/Quantiseq_panels.pdf", 
       width = 14, height = 16)


res_xcell_panels <- 
  res_xcell %>%
  pivot_longer(!cell_type, names_to = "sample_id", values_to = "score") %>% 
  left_join(longitudinal_metadata, by = "sample_id") %>% 
  filter(case_control_other_latest %in% c("Case", "Control")) %>% 
  ggplot(aes(x=jitter(visit_month, factor = 10), y=score, color=case_control_other_latest)) +
  geom_point(size = 0.3) +
  facet_wrap(~cell_type, scales="free") +
  geom_line(aes(color = case_control_other_latest, group = participant_id), alpha = 0.05) +
  geom_smooth(method = "loess", alpha=1, fill = "white") + 
  geom_ysidedensity(aes(x=stat(density)), position = "stack") +
  geom_xsidedensity(aes(y=stat(density)), position = "stack") +
  ggside(scales = "free") +
  labs(x = "visit month", y = "cell-type score estimation", 
       title = "xCell Bulk-RNASeq Deconvolution") +
  scale_color_d3() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))
# res_xcell_panels

ggsave(res_xcell_panels, filename = "figures/immune_cell_prop/xCell_panels.pdf", 
       width = 18, height = 30)

#_______________________________________________________________________________

res_epic_neutrophils <- 
  res_epic %>%
  pivot_longer(!cell_type, names_to = "sample_id", values_to = "score") %>% 
  left_join(longitudinal_metadata, by = "sample_id") %>% 
  filter(case_control_other_latest %in% c("Case", "Control"),
         cell_type == "Neutrophil") %>% 
  ggplot(aes(x=jitter(visit_month, factor = 10), y=score, color=case_control_other_latest)) +
  geom_point(alpha = 1, size = 0.5) +
  # facet_wrap(~study, scales="free_y", ncol=3) +
  geom_line(aes(group = participant_id), alpha = 0.05) +
  geom_smooth(method = "loess", alpha=1, fill = "white") + 
  geom_ysidedensity(aes(x=stat(density)), position = "stack") +
  geom_xsidedensity(aes(y=stat(density)), position = "stack") +
  scale_color_d3() +
  labs(x = "visit month", y = "cell-type score estimation", 
       title = "EPIC Bulk-RNASeq Deconvolution - Neutrophils", 
       color = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = c(0.8, 0.8),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))
res_epic_neutrophils
ggsave(res_epic_neutrophils, filename = "figures/immune_cell_prop/Trashme_res_epic_neutrophils2.pdf", 
       width = 7, height = 7)

#_______________________________________________________________________________

# Fractional data visualization
longitudinal_metadata_essent <- 
  longitudinal_metadata %>% 
  dplyr::select(participant_id, sample_id, visit_month,
         case_control_other_latest, study, sex) %>%  #%>% get_dupes()
  drop_na(sample_id) %>% 
  distinct()

# dim(longitudinal_metadata_essent)
# dim(distinct(longitudinal_metadata_essent))
# dim(longitudinal_metadata)


# Res options (res_quantiseq or res_epic)
res_df <- res_epic
res_df_matrix <- res_df %>% 
  column_to_rownames("cell_type") %>% 
  as.matrix()
og_sample_order <- colnames(res_df_matrix)
og_cell_order <- colnames(t(res_df_matrix))

hc <- 
  res_df_matrix %>% 
  dist() %>% 
  hclust()
dd_col <- as.dendrogram(hc)
cell_order <- order.dendrogram(dd_col)
# dd_col.reorder <- reorder(dd_col, 10:1)
# plot(dd_col, main = "random dendrogram 'dd'")

hr <- 
  res_df_matrix %>% t() %>% 
  dist() %>%
  hclust()
dd_row <- as.dendrogram(hr)
sample_order <- order.dendrogram(dd_row)
# dd_row.reorder <- reorder(dd_row, 10:1)
# plot(dd_row, main = "random dendrogram 'dd'")

p1 <- res_df %>%
  pivot_longer(!cell_type, names_to = "sample_id", values_to = "fraction") %>% 
  left_join(longitudinal_metadata_essent, by = "sample_id") %>% 
  filter(case_control_other_latest %in% c("Case", "Control")) %>%
  mutate(sample_id = factor(sample_id, levels = og_sample_order[sample_order]),
         cell_type = factor(cell_type, levels = og_cell_order[cell_order])) %>% 
  # plot as stacked bar chart
  ggplot(aes(x=reorder(sample_id,-fraction, .fun='mean'), 
             y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  # coord_flip() +
  scale_fill_brewer(palette="Paired") +
  # scale_x_discrete(limits = rev(levels(res_quantiseq))) +
  facet_grid(~case_control_other_latest, scales = "free_x") +
  labs(x = NULL, y = "Relative Abundance") +
  theme_classic() +
  theme(axis.text.x = element_blank())


p2 <- res_df %>%
  pivot_longer(!cell_type, names_to = "sample_id", values_to = "fraction") %>% 
  left_join(longitudinal_metadata_essent, by = "sample_id") %>% 
  filter(case_control_other_latest %in% c("Case", "Control")) %>% 
  mutate(sample_id = factor(sample_id, levels = og_sample_order[sample_order]),
         cell_type = factor(cell_type, levels = og_cell_order[cell_order]),
         fraction = asin(sqrt(fraction))) %>% 
  ggplot(aes(y = fraction, x = cell_type, color = case_control_other_latest)) +
  scale_color_d3() +
  # scale_y_log10() +
  labs(x = NULL, y = "Normalized Abundance (AST)") +
  theme_bw() +
  geom_boxplot(notch = T, size = 0.3, outlier.alpha = 0.2) +
  facet_wrap(~cell_type, scales = "free", nrow = 1, strip.position = "right") +
  # geom_violin() +
  theme(axis.text.x = element_blank(), 
        axis.line.x.top = element_line(linetype = 1), 
        panel.grid = element_blank())
p2

c1 <- cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(2, 1), align = "h", axis = "lr")
ggsave(c1, filename = "figures/immune_cell_prop/res_epic_summary.pdf", 
       width = 14, height = 6)


