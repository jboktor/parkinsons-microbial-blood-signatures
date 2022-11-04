# Joe Boktor
# Caltech - Mazmanian Lab
# Dec 2021

source("src/_load_packages.R")
source("src/_plot-functions.R")
source("src/_misc_functions.R")

# Note: Compile flag stats before running this script with "process-flagstat.R"

qc_stats <-
  read.table(
    file = "input_files/2021_v2-5release_0510/releases_2021_v2-5release_0510_wgs_gatk_metrics_preBqsr_selfSM.tsv",
             stringsAsFactors = F, sep = "\t",
             header = TRUE)

flagstat_df <- readRDS("data/Sequencing_stats/unmapped_bam_flagstats.rds")
wgs_summary_stats <- readRDS("input_files/amp-pd-qc-data/wgs_duplicate_read_stats.rds")
wgs_summary_stats %<>% 
  dplyr::rename(participant_id = LIBRARY) %>% 
  mutate_at(vars(-participant_id), as.numeric)

RefSeqPlusPF_reports_long_trim <- 
  readRDS("data/Kraken_reports/RefSeqPlusPF_reports_merged.rds") %>% 
  slim_report_df() %>% 
  mutate(ref_DB = "RefSeqPlusPF") %>% 
  filter(NCBI_taxon_ID %in% c(0, 1, 2, 4751, 10239)) %>% 
  select(-c(NCBI_taxon_ID, taxon_rank_code)) %>% 
  pivot_wider(names_from = "feature", values_from = "count")

readqc_stats <- qc_stats %>% 
  select(participant_id, READS, FREEMIX, FREELK1, FREELK0) %>% 
  left_join(RefSeqPlusPF_reports_long_trim) %>% 
  left_join(flagstat_df) %>% 
  left_join(wgs_summary_stats)
readqc_stats %>% glimpse()

readqc_stats %>% 
  ggplot(aes(x=QC_passed_unmapped_reads, y= UNMAPPED_READS)) +
  geom_point() +
  geom_smooth(method = "lm") +
  # scale_color_viridis(option = "F") +
  clean_theme()

# tst <- readqc_stats[readqc_stats$participant_id %nin% wgs_summary_stats$participant_id, ]

#_______________________________________________________________________________

refseq_qc_plot <- 
  readqc_stats %>%
  select(participant_id, FREEMIX, READ_PAIRS_EXAMINED, QC_passed_unmapped_reads, root, Bacteria, Fungi, Viruses) %>%
  pivot_longer(-c(participant_id, FREEMIX), names_to = "metric") %>% 
  mutate(metric = factor(
    metric,
    levels = c(
      "READ_PAIRS_EXAMINED",
      "QC_passed_unmapped_reads",
      "root",
      "Bacteria",
      "Fungi",
      "Viruses"
    )
  )) %>%
  ggplot(aes(x = metric, y = value)) +
  geom_line(aes(group = fct_reorder(participant_id, FREEMIX), color = FREEMIX), alpha = 0.1) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.6, width = 0.3) +
  scale_color_viridis(option = "F", direction = -1) +
  scale_y_log10() +
  labs(y = "Reads/Counts", x = NULL) +
  scale_x_discrete(
    labels = c("READ_PAIRS_EXAMINED", "QC_passed_unmapped_reads" = "QC Passed Total\nUnmapped Reads", "root" = "Root")) +
  clean_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
refseq_qc_plot

ggsave(
  refseq_qc_plot,
  filename = "figures/quality_control/RefSeq_read_mapping.svg",
  width = 4, height = 5)



readqc_stats %>% 
  ggplot(aes(x=READS, y= root , color = FREEMIX)) +
  geom_point() +
  scale_color_viridis(option = "F") +
  clean_theme()
readqc_stats %>% 
  ggplot(aes(x=FREEMIX, y= unclassified , color = FREEMIX)) +
  geom_point() +
  scale_color_viridis(option = "F") +
  clean_theme()
readqc_stats %>% 
  ggplot(aes(x=FREELK1, y= root , color = FREEMIX)) +
  geom_point() +
  scale_color_viridis(option = "F") +
  clean_theme()


