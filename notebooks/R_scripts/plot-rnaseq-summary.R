# Joe Boktor
# Caltech - Mazmanian Lab

source("src/_load_packages.R")
source("src/_misc_functions.R")

# load Metadata
longitudinal_metadata <- readRDS(file = "data/Metadata/longitudinal_metadata.rds")
rna_sample_inv <-
  read.csv(file = "input_files/2021_v2-5release_0510/rna_sample_inventory.csv",
           stringsAsFactors = F, header = TRUE)

longitudinal_metadata %<>%
  mutate(visit_name = gsub("#2", "", visit_name)) %>%
  left_join(rna_sample_inv, by = c("participant_id", "visit_month"))
longitudinal_metadata %>% dim()

longitudinal_metadata_essent <-
  longitudinal_metadata %>%
  dplyr::select(participant_id, sample_id, visit_month, race, ethnicity,
                case_control_other_latest, age_at_baseline, study, sex) %>%
  drop_na(sample_id) %>%
  distinct()
longitudinal_metadata_essent %>% dim()

rna_sample_df <-
  rna_sample_inv %>%
  left_join(longitudinal_metadata_essent) %>%
  mutate(age = ((age_at_baseline*12)+visit_month)/12 )
rna_sample_df %>% dim()




rna_sample_df %>% 
  ggplot(aes(x = age, y = visit_month)) +
  geom_point(aes(color = case_control_other_latest), size = 1)





rna_sample_df %>% 
  drop_na(case_control_other_latest) %>% 
  # filter(age > 50, age < 55) %>%
  # filter(age > 50, age < 55) %>%
  mutate(sample_id = fct_reorder(sample_id, age, .fun='min', .desc = T )) %>%
  ggplot(aes(x = age, y = sample_id)) +
  geom_point(aes(color = case_control_other_latest), size = 1) +
  theme_classic() +
  # facet_grid(sex~race) +
  # scale_x_log10() +
  theme(
    # panel_grid = element_blank(),
        axis.text.y = element_blank())

# longitudinal_metadata$race

rna_sample_inv$visit_month %>% unique()
