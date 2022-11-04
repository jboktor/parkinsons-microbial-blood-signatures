# Exploratory analysis

source("src/_load_packages.R")
source("src/_plot-functions.R")
# source("src/sample_stats.R")

df.plot.A <- mdsupdrs %>%
  left_join(df.sam, by = c("participant_id")) 

df.plot.plasma <- df.plot.A %>% 
  right_join(plasma_somalogic) 

df.plot.plasma %>% str()

df.plot.A %>%
  ggplot(aes(x=mds_updrs_part_i_summary_score, mds_updrs_part_iii_summary_score)) +
  geom_point(size = 0.3) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  # geom_smooth(method = "lm") +
  scale_fill_continuous(type = "viridis") +
  facet_wrap(~visit_month) +
  clean_theme()

tst <- df.plot.plasma %>%
  filter(test_name == "Interleukin-6")

tst %>%
  ggplot(aes(x=visit_month, y=test_value, group = visit_month)) +
  geom_boxplot() +
  geom_point(aes(fill=visit_month), position = position_jitterdodge())

df.plot %>%
  ggplot(aes(x=sex, mds_updrs_part_iii_summary_score)) +
  geom_boxplot() +
  # geom_point(aes(fill = cohort), position = position_jitterdodge(jitter.width = 0.5),
  #            shape=21, size=1.5, alpha = 0.9) +
  # geom_smooth(method = "loess") +
  stat_compare_means(aes(label = ..p.signif..)) +
  scale_fill_jco() +
  clean_theme()



ppmi_bio <-
  read.csv(file = "input_files/PPMI/Current_Biospecimen_Analysis_Results.csv",
           stringsAsFactors = F, header = TRUE)


ppmi_il6 <- 
  ppmi_bio %>% 
  janitor::clean_names() %>% 
  filter(testname == "IL-6") %>% 
  mutate(testvalue = as.numeric(testvalue)) %>% 
  mutate(grouping_col = paste(clinical_event, diagnosis, sep = "_")) %>% 
  ggplot(aes(x=clinical_event, y=log10(testvalue + 1), group = grouping_col)) +
  geom_boxplot() +
  clean_theme() +
  geom_point(aes(color=diagnosis), 
             position = position_jitterdodge(jitter.width = 0.2), alpha = 0.6) +
  scale_color_npg()



