
source("src/_load_packages.R")
source("src/_misc_functions.R")
source("src/_plot-functions.R")
base::load("data/Phyloseq_Objects/Phyloseq_all_outliers_removed.RData")
base::load("data/Analyses/feature_AUROC/feature_AUROCs_sex.RData")
aucs_sex <- aucs # Temporary fix me later
base::load("data/Analyses/feature_AUROC/feature_AUROCs.RData")


# tstdat <- refDB_phyloseq_orm[[refDB]][[level]]
# tstab <- tstdat %>%
#   abundances() %>%
#   as.data.frame() %>% 
#   rownames_to_column()

# tstab$`Salmonella enterica`
#_______________________________________________________________________________
#                         MaAsLin2 Visualization
#_______________________________________________________________________________
refDB <- "RefSeqPlusPF"
level <- "Species"
fileinput <- (paste0("data/Analyses/differential_abundance/", refDB, 
                     "/Maaslin2_output_", level, "/all_results.tsv"))

aucs.masformat <- aucs %>% 
  filter(ref_DB == refDB) %>% 
  dplyr::mutate(auroc_center = auroc - 0.5) %>%
  mutate(feature = gsub("-", ".", feature)) %>% 
  mutate(feature = gsub(" ", ".", feature))
aucs_sex.masformat <- aucs_sex %>% 
  filter(ref_DB == refDB) %>% 
  dplyr::mutate(auroc_center = auroc - 0.5) %>%
  mutate(feature = gsub("-", ".", feature)) %>% 
  mutate(feature = gsub(" ", ".", feature))

maaslin2_df <- read_tsv(fileinput, col_names = T) %>% 
  left_join(aucs.masformat, by = "feature")
maaslin2_df_sex <- read_tsv(fileinput, col_names = T) %>% 
  left_join(aucs_sex.masformat, by = "feature")


plot_df <- maaslin2_df %>% 
  filter(metadata == "case_control_other_latest") %>% 
  mutate(hit = if_else(qval <= .05 & abs(auroc_center) > 0.1, "Yes", "No"))
plot_df_sex <- maaslin2_df_sex %>% 
  filter(metadata == "sex") %>% 
  mutate(hit = if_else(qval <= .05 & abs(auroc_center) > 0.1, "Yes", "No"))
  
# 
plot_df_sig <- plot_df %>%
  filter(hit == "Yes") %>%
  slice_min(order_by = qval, n = 10, with_ties = F)
plot_df_sex_sig <- plot_df_sex %>%
  filter(hit == "Yes") %>%
  slice_min(order_by = qval, n = 10, with_ties = F)

 # Create Output Dir ----
fig_output_dir_root <- paste0("figures/differential_abundance/", refDB)
fig_output_dir <- paste0(fig_output_dir_root, "/", level, "/")
dir.create(path = fig_output_dir_root, showWarnings = F)
dir.create(path = fig_output_dir, showWarnings = F)

#Plotting Figures ----
volc <- 
  plot_df %>% 
  ggplot(aes(x=auroc, y = -log10(qval + 1e-100))) +
  geom_point(aes(color = hit)) +
  scale_x_continuous(limits = c(0.2, 0.8)) +
  scale_color_manual(values = c("No" = "darkgrey", "Yes" = "red")) +
  labs(color = "Disease \nSignal") +
  geom_label_repel(
    data = plot_df_sig,
    aes(x = auroc, y = -log10(qval + 1e-100), label = feature), size = 3,
    seed = 42, force = 2) +
  clean_theme()
volc
# ggsave(volc, filename = paste0(fig_output_dir, "volcano-plot.png"), 
#        width = 5, height = 7, dpi = 600)

for(foi in plot_df_sig$feature){
  
  p <- boxplots_maaslin2(refDB = refDB, level = level, feature = foi)
  ggsave(p, filename = paste0(fig_output_dir, foi, ".png"),
         dpi = 400, width = 2, height = 3)
}



sex_volc <- 
  plot_df_sex %>% 
  ggplot(aes(x=auroc, y = -log10(qval + 1e-100))) +
  geom_point(aes(color = hit)) +
  # scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("No" = "darkgrey", "Yes" = "red")) +
  labs(color = "Sex \nSignal") +
  geom_label_repel(
    data = plot_df_sex_sig,
    aes(x = auroc, y = -log10(qval + 1e-100), label = feature), size = 3,
    seed = 42, force = 50, segment.size = 0.1,
    nudge_y = -5) +
  clean_theme()
sex_volc
# ggsave(sex_volc, filename = paste0(fig_output_dir, "volcano-plot_SEX.png"),
#        width = 7, height = 7, dpi = 600)


#_______________________________________________________________________________
#                         CORNCOB Visualization
#_______________________________________________________________________________

base::load("data/Analyses/differential_abundance/OLD/UHGG/Species_corncob_diagnosis_latest_Sex_DispControl.RData")
base::load("data/Analyses/differential_abundance/OLD/UHGG/Species_corncob_diagnosis_latest_DispControl.RData")

#_______________________________________________________________________________
summary(model_diagnosis_disp)
plot(model_diagnosis_disp)

summary(model_diagnosisSex_disp)
plot(model_diagnosisSex_disp)
#_______________________________________________________________________________

model_df <- tibble()
iterlength <- length(model_diagnosis_disp$all_models)
pb <- progress_bar$new(
  format = "  aggregating model data [:bar] :current/:total :percent eta: :eta",
  total = iterlength, clear = FALSE, width= 90)

for (Nmodel in 1:iterlength){
  est <- model_diagnosis_disp$all_models[[Nmodel]]$coefficients %>% 
    as.data.frame() %>% rownames_to_column(var = "parameter") %>% 
    mutate(model_n = Nmodel)
  model_df <- rbind(model_df, as.data.frame(est))
  pb$tick()
}

model_df_long <- model_df %>% 
  pivot_wider(names_from = parameter, names_sep = "__",
              values_from = c(Estimate, `Std. Error`, `t value`, `Pr(>|t|)`)) %>% 
  janitor::clean_names()

corcob_res_df <- data.frame(
  features = names(model_diagnosis_disp$p),
  p_val = model_diagnosis_disp$p,
  q_val = model_diagnosis_disp$p_fdr,
  model_n = 1:length(model_diagnosis_disp$p)
  )

model_stats <- corcob_res_df %>% 
  left_join(model_df_long, by = "model_n") %>% 
  as_tibble()
#_______________________________________________________________________________


model_stats %>% 
  ggplot(aes(x = estimate_mu_diagnosis_latest_no_pd_nor_other_neurological_disorder, 
             y = -log10(q_val), 
             color = model_n)) +
  clean_theme() +
  geom_point(alpha = 0.6)

model_stats %>% 
  ggplot(aes(x = estimate_phi_diagnosis_latest_no_pd_nor_other_neurological_disorder, 
             y = -log10(q_val), 
             color = model_n)) +
  clean_theme() +
  geom_point(alpha = 0.6)

model_stats %>% 
  filter(q_val <= 0.05) %>%
  ggplot(aes(x = estimate_phi_diagnosis_latest_no_pd_nor_other_neurological_disorder, 
             y = estimate_mu_diagnosis_latest_no_pd_nor_other_neurological_disorder,
             xmin = estimate_phi_diagnosis_latest_no_pd_nor_other_neurological_disorder - 
               std_error_phi_intercept, 
             xmax = estimate_phi_diagnosis_latest_no_pd_nor_other_neurological_disorder +
               std_error_phi_intercept,
             ymin = estimate_mu_diagnosis_latest_no_pd_nor_other_neurological_disorder - 
               std_error_mu_intercept, 
             ymax = estimate_mu_diagnosis_latest_no_pd_nor_other_neurological_disorder +
               std_error_mu_intercept,
             color = q_val)) +
  clean_theme() +
  geom_errorbar() +
  geom_errorbarh() +
  geom_point(alpha = 0.2) +
  scale_color_viridis(option = "I")


corcobstats_allsig <- model_stats %>% 
  filter(q_val <= .05) %>%
  # slice_min(q_val, n = 50, with_ties = F) %>%
  ggplot(aes(y = reorder(features, 
                         -q_val),
             x = estimate_mu_diagnosis_latest_no_pd_nor_other_neurological_disorder,
             xmin = estimate_mu_diagnosis_latest_no_pd_nor_other_neurological_disorder - 
               std_error_mu_intercept, 
             xmax = estimate_mu_diagnosis_latest_no_pd_nor_other_neurological_disorder +
               std_error_mu_intercept,
             color = q_val)) +
  geom_errorbar() +
  clean_theme() +
  labs(x = "Estimate Control vs PD", y = NULL) +
  scale_color_viridis(option = "magma") +
  geom_point() +
  theme(axis.text.y = element_text(size = 1))
corcobstats_allsig
# ggsave(corcobstats_allsig, filename = "figures/differential_abundance/allfdrsig_corncob.png",
#        width = 8, height = 9)

corcobstats <- model_stats %>% 
  filter(q_val <= 10^-8) %>%
  slice_min(q_val, n = 50, with_ties = F) %>%
  ggplot(aes(y = reorder(features, 
                         -q_val),
           x = estimate_mu_diagnosis_latest_no_pd_nor_other_neurological_disorder,
           xmin = estimate_mu_diagnosis_latest_no_pd_nor_other_neurological_disorder - 
             std_error_mu_intercept, 
           xmax = estimate_mu_diagnosis_latest_no_pd_nor_other_neurological_disorder +
             std_error_mu_intercept,
           color = q_val)) +
  geom_errorbar() +
  clean_theme() +
  labs(x = "Estimate Control vs PD") +
  scale_color_viridis(option = "magma") +
  geom_point()
corcobstats
# ggsave(corcobstats, filename = "figures/differential_abundance/top50corncob.png",
#        width = 8, height = 9)


