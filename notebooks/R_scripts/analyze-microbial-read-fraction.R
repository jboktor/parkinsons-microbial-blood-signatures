# Joe Boktor
# Caltech - Mazmanian Lab

source("src/_load_packages.R")
source("src/_misc_functions.R")

# Approx size of human genome
human_gen <- 2.9*10^9

wgs_stats <- read_tsv("input_files/2021_v2-5release_0510/releases_2021_v2-5release_0510_wgs_gatk_metrics_raw_wgs_metrics.tsv", 
                      col_names = T, show_col_types = FALSE)
wgs_stats %>% head()
average_coverage <- wgs_stats$MEAN_COVERAGE %>% mean()
average_coverage
seq_depth_clean <- average_coverage * human_gen
seq_depth_clean


base::load("data/Kraken_reports/WoL_reports_merged.RData")
WoL_reports <- kraken_df
WoL_reports %>% head()

WoL_reports_long <- slim_report_df(WoL_reports) 
WoL_map <- WoL_reports_long %>% 
  filter(feature %in% c("unclassified", "root"))
total_reads_input <- 
  WoL_map %>% 
  group_by(participant_id) %>% 
  summarise(count_sum = sum(count))
unmapped_reads <- total_reads_input %>% pull(count_sum) %>% mean
unmapped_reads

# On average, 7444922 reads an unclassified while on avareg

# unmapped_reads/(seq_depth_clean + unmapped_reads) * 100

# Stats for total reads that don't align to hg38
total_reads_input %>% 
  summarise(mean = mean(count_sum), median = median(count_sum))

# Stats for total reads that don't align to hg38
# stratified by root and unclass.
WoL_map %>% 
  group_by(feature) %>% 
  summarise(mean = mean(count), median = median(count))



