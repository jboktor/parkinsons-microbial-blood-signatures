# ML Viz all

source("src/_load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
source("src/ml_models.R")
source("src/ml_plots.R")
load("files/low_quality_samples.RData")
library(fs)

file_paths <- fs::dir_ls("data/Machine_Learning_Analysis/model_stats")
file_paths

loso_files <- list()
s2s_files <- list()
for (i in seq_along(file_paths)){

  if (grepl("LOSO_", file_paths[[i]])){
    loso_files[[i]] <-
      read_csv(file = file_paths[[i]]) %>%
      mutate(filename = sub(".*model_stats/", "", file_paths[[i]] ))
  } else if (grepl("S2S_", file_paths[[i]])) {
    s2s_files[[i]] <-
      read_csv(file = file_paths[[i]]) %>%
      mutate(filename = sub(".*model_stats/", "", file_paths[[i]] ))
  }
}


df_s2s <- purrr::map_df(s2s_files, ~.x)

df_s2s <- df_s2s %>%
  mutate(model_perf = coalesce(mean, .estimate)) %>%
  filter(.metric == "roc_auc") %>% distinct() %>%
  pivot_wider(names_from = .metric, values_from = "model_perf") %>%
  mutate(rarefied = if_else(grepl("_Rarefied", filename), "Yes", "No" )) %>%
  mutate(filename = gsub("_Rarefied", "", filename)) %>%
  mutate(filename = gsub("all_features", "all", filename)) %>%
  mutate(filename = gsub(".csv", "", filename)) %>%
  separate(filename, c("analysis", "object", "model_type", "features"),
           sep = "_", remove = F)

df_s2s_slim <- df_s2s %>%
  select(train, test, filename, roc_auc, rarefied, std_err) %>%
  distinct() %>%
  dplyr::group_by(filename, rarefied) %>%
  pivot_wider(names_from = c("test"), values_from = "roc_auc")

df_s2s_slim %>%
  ungroup() %>%
  select(Shanghai, Bonn, TBC, Rush) %>%
  ggpairs()

# Calculate average test/prediction values on S2S analysis & distance from training values
all_cohorts <- c("Shanghai", "Bonn", "Rush", "TBC")
for (row in seq_along(df_s2s_slim$train)){
  group <- pull(df_s2s_slim[row, "train"])
  print(group)

  testvals <-
    df_s2s_slim[row,] %>%
    ungroup() %>%
    select(str_subset(all_cohorts, group, negate = T)) %>%
    as.double()

  test_avg <- testvals %>% mean()
  test_sd <- testvals %>% sd()
  dist_from_train <- pull(df_s2s_slim[row, group]) - test_avg
  print(test_avg)
  print(test_sd)
  print(dist_from_train)

  df_s2s_slim[row, "test_avg"] <- test_avg
  df_s2s_slim[row, "test_sd"] <- test_sd
  df_s2s_slim[row, "dist_from_train"] <- dist_from_train

}
tstlist <- c(1, 44, 5)
df_s2s_slim %>%
  separate(
    filename,
    c("analysis", "object", "model_type", "features"),
    sep = "_",
    remove = F
  ) %>%
  ggplot(aes(x = test_avg, y = dist_from_train)) +
  geom_rect(aes(xmin=-Inf, xmax=0.6, ymin=-Inf, ymax=Inf), fill = "#e9e9e9") +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0.2, ymax=Inf), fill = "#e9e9e9") +
  geom_point(aes(color = object, shape = rarefied)) +
  scale_color_npg() +
  my_clean_theme2()
