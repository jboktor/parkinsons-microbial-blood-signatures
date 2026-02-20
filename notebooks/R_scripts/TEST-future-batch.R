
pdmbs_dir <- "/central/groups/MazmanianLab/joeB/PDMBS"
wgs_wkdir <- paste0(pdmbs_dir, "/workflow/WGS")
wkdir <- paste0(pdmbs_dir, "/parkinsons-microbial-blood-signatures")
source(paste0(wkdir, "/notebooks/R_scripts/_misc_functions.R"))
source(paste0(wkdir, "/notebooks/R_scripts/_load-core-pkgs.R"))
library(phyloseq)
library(microbiome)
library(ggsci)
library(aplot)

physeq_decon_dir <- glue(
  "{wkdir}/data/processed/",
  "phyloseq_objects/decontaminated"
)
dir.create(physeq_decon_dir, recursive = TRUE, showWarnings = FALSE)

# Platform unit naming convention
# {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}.
# ftp://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups

refDB <- "UHGG"
level <- "Species"
seq_method <- "WGS"

wgs_sequencing_meta <- readRDS(
  glue(
    "{wkdir}/data/interim/metadata/",
    "2022-12-31_sequencing_metadata.rds"
  )
)
wgs_sequencing_meta_formatted <- wgs_sequencing_meta %>%
  mutate(PLATFORM_UNIT = 
  if_else(grepl("LB-", SAMPLE), 
    gsub("_X", ".X", PLATFORM_UNIT) %>% gsub("_E", ".E", .),
    PLATFORM_UNIT)
  ) %>% 
  separate(PLATFORM_UNIT,
    c("FLOWCELL_BARCODE", "SAMPLE_BARCODE", "PLATFORM_UNIT_LANE"),
    sep = "\\.", remove = FALSE
  ) %>%
    mutate(SAMPLE = str_remove(SAMPLE, "SM:"))

meta_freq <-
  colnames(wgs_sequencing_meta_formatted) %>%
  purrr::map(
    ~ wgs_sequencing_meta_formatted %>%
      group_by(!!sym(.x)) %>%
      summarize(n = n())
  )

center_info <- wgs_sequencing_meta_formatted %>%
  select(SAMPLE, CENTER) %>%
  distinct()
flowcell_info <- wgs_sequencing_meta_formatted %>%
  get_dupes(FLOWCELL_BARCODE) %>%
  filter(dupe_count >= 25) %>%
  select(SAMPLE, FLOWCELL_BARCODE) %>%
  distinct() %>%
  mutate(detected = 1) %>%
  pivot_wider(
    names_from = FLOWCELL_BARCODE,
    values_from = detected, values_fill = 0
  )

dates <- wgs_sequencing_meta_formatted$DATE %>% unique

# Adding sequencing data metadata to phyloseq object
meta_wide <- wgs_sequencing_meta_formatted %>%
  get_dupes(DATE) %>%
  filter(dupe_count >= 25) %>%
  select(SAMPLE, DATE) %>%
  distinct() %>%
  mutate(detected = 1) %>%
  pivot_wider(names_from = DATE, values_from = detected, values_fill = 0) %>%
  dplyr::full_join(flowcell_info, by = "SAMPLE") %>%
  dplyr::full_join(center_info, by = "SAMPLE") %>%
  dplyr::rename(participant_id = SAMPLE) %>%
  clean_names() %>%
  mutate_if(is.numeric, as.logical) %>%
  mutate_if(is.logical, ~ replace_na(., FALSE))
date_timepoints <- meta_wide %>%
  select_if(is.logical) %>%
  select(contains("dt_")) %>% 
  colnames()
flowcell_barcodes <- meta_wide %>%
  select_if(is.logical) %>%
  select(contains("pu_")) %>% 
  colnames()


ps <- readRDS(
  glue(
    "{wkdir}/data/processed/phyloseq_objects/raw/",
    "2023-07-14_{seq_method}_{refDB}_{level}_phyloseq.rds"
  )
)
phyloseq::sample_data(ps) %>% dim
meta <- ps %>% meta()
meta %<>%
  dplyr::left_join(meta_wide, by = "participant_id") %>%
  mutate_if(is.logical, ~ replace_na(., FALSE))
rownames(meta) <- meta$participant_id
sample_data(ps) <- sample_data(meta)
phyloseq::sample_data(ps) %>% dim

# _________________________________________________________________
# Decontaminaton procedure ----
#' Iterate over cohort / flow cell / and sequencing dates
#' extract samples which have only been run in a single run
#' (check to see if this means we lose any cohorts or flowcells)
#' Then, calculate prevalence for each in-group / out-group sample set

#' Function take subsets of samples from a phyloseq object and return a
#' dataframe with prevalence stats on a per-feature basis for within and
#' without subset groups
# filter samples by cohort
# Sequencing date
# Sequencing center
# Flowcell barcode


# ______________________________________________________________________________
# Prevalence filter ----

calculate_batch_prevalence <- function(ps, subset_var) {
  require(magrittr)
  subset_data <- microbiome::meta(ps) %>%
    dplyr::pull(subset_var)
  batches <- subset_data %>% unique()
  prevalence_stats <- tibble::tibble()
  for (batch in batches) {
    print(batch)
    batch_samples <<- grepl(batch, subset_data)
    if (sum(batch_samples) < 25) next # skip if there are less than 3 samples
    ps_batch <- phyloseq::subset_samples(ps, batch_samples)
    print(ps_batch)
    batch_stats <- microbiome::prevalence(ps_batch) %>%
      stack() %>%
      dplyr::mutate(
        batch_subset = as.character(batch),
        batch_variable = subset_var,
        sample_n = phyloseq::nsamples(ps_batch)
      )
    prevalence_stats %<>% dplyr::bind_rows(batch_stats)
  }
  return(prevalence_stats)
}

batch_meta <- c(
  "study",
  "center",
  date_timepoints,
  flowcell_barcodes
)
ps_groups <- c(
  "Case" = phyloseq::subset_samples(
    ps,
    case_control_other_latest == "Case"
  ),
  "Control" = phyloseq::subset_samples(
    ps,
    case_control_other_latest == "Control"
  ),
  "Other" = phyloseq::subset_samples(
    ps,
    case_control_other_latest == "Other"
  )
)

# all possible combos of an x var and a y var
pair_combos <- expand_grid(
  ps_list = names(ps_groups),
  var_list = batch_meta
  ) %>%
  mutate(batch_id = glue("{ps_list}_{var_list}"))
ps_list <- pair_combos$ps_list
var_list <- pair_combos$var_list
names(ps_list) <- pair_combos$batch_id

tic()
future::plan("multisession", workers = 14)
prevalence_stats <-
  furrr::future_map2(
    ps_list,
    var_list,
    ~ calculate_batch_prevalence(ps_groups[[.x]], .y)
  ) %>%
  bind_rows(.id = "grouping_meta")
toc()


#' Select features with a max of at least 25% prevalence and
#' a minimum of < 2x for other batches
prevalence_stats_proc <- prevalence_stats %>%
  filter(batch_subset != "FALSE") %>%
  mutate(batch_subset = case_when(
    grepl("dt_|pu_", batch_variable) ~ batch_variable,
    TRUE ~ batch_subset
  ), batch_variable = case_when(
    grepl("dt_", batch_variable) ~ "date",
    grepl("pu_", batch_variable) ~ "flowcell_barcode",
    TRUE ~ batch_variable
  )) %>%
  mutate(case_control_other_latest = strex::str_before_first(
    grouping_meta, "_"
  ))

prevalence_stat_summary <-
  prevalence_stats_proc %>%
  # remove batches with fewer than 25 samples
  filter(sample_n > 25) %>%
  # double check that all batches to compare have at least a min and a max
  get_dupes(ind, batch_variable, case_control_other_latest) %>%
  filter(dupe_count > 1) %>% 
  # summarize the min and max prevalence for each microbe
  group_by(ind, batch_variable, case_control_other_latest) %>%
  dplyr::summarise(
    min = min(values), max = max(values),
    prevalence_delta = max - min
  )

  #' flag microbes that are present in at least 50% of samples in at least
  #' one batch and are at least 2x more abundant in that batch than in any other
blacklist_prev <- prevalence_stat_summary %>%
  filter(max > 0.5, min <= max / 2) %>%
  arrange(desc(prevalence_delta)) %>%
  pull(ind) %>%
  unique()

saveRDS(
  prevalence_stats_proc,
  glue(
    "{wkdir}/data/interim/decontamination/",
    "{Sys.Date()}_{seq_method}_{refDB}_{level}_prevalence_stats.rds"
  )
 )

saveRDS(
  blacklist,
  glue(
    "{wkdir}/data/interim/decontamination/",
    "{Sys.Date()}_{seq_method}_{refDB}_{level}_prevalence_blacklist.rds"
  )
 )

plot_df <- prevalence_stat_summary %>%
  mutate(blacklist = ind %in% blacklist) %>%
  arrange(blacklist, desc(prevalence_delta)) %>%
  mutate(ind = factor(ind))

p1 <- plot_df %>%
  ggplot(aes(x = ind, y = prevalence_delta, color = blacklist)) +
  geom_point(alpha = 0.7) +
  scale_color_d3() +
  facet_grid(batch_variable ~ case_control_other_latest) +
  theme_bw() +
  labs(x = "microbe", y = "max prevalence - min prevalence") +
    theme(axis.text.x = element_blank())

# How does blacklist detection compare across groups
p2 <- plot_df %>%
    ggplot(aes(
      x = case_control_other_latest,
      y = prevalence_delta,
      color = blacklist
    )) +
    geom_point(alpha = 0.3) +
    geom_line(aes(group = ind), alpha = 0.3) +
    facet_wrap(~batch_variable, ncol = 1) +
    scale_color_d3() +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(legend.position='none')

final_plot <- p1 %>% insert_right(p2, width = 0.5)
ggsave(
  glue(
    "{wkdir}/figures/decontamination/",
    "{Sys.Date()}_{refDB}_{seq_method}_{level}_prevalence_filter_EDA.png"
  ),
  final_plot, width = 12, height = 8
)


# ______________________________________________________________________________
# Read count filter ----

blacklist_readcounts <-
  filter_taxa(ps, function(x) max(x) < 100, TRUE) %>%
  taxa()
blacklist_readcounts

phyloseq::filter_taxa(ps, blacklist_readcounts)
all_taxa <- taxa(ps)
saveRDS(
  blacklist_readcounts,
  glue("{wkdir}/data/interim/decontamination/",
  "2023-01-09_{seq_method}_{refDB}_{level}_blacklist_readcounts.rds")
 )

ps_trim <- phyloseq::prune_taxa(
  all_taxa[all_taxa %nin% blacklist_readcounts], ps
)
saveRDS(
  ps_trim,
  glue(
    "{physeq_decon_dir}/",
    "{Sys.Date()}_{seq_method}_{refDB}_{level}_phyloseq.rds"
  )
)







#___________________
# Correlations to blacklist ----

corr_loop <- function(metadata, abundance, obj.name) {
  corr_output <- tibble()
  for (metavar in colnames(metadata)) {
    cat("Calculating correlations for: ", metavar[[1]], "\n")
    for (feature in colnames(abundance)) {
      # Calculate Spearman's Correlation
      spearman <-
        cor.test(
          x = metadata[[metavar]],
          y = abundance[[feature]],
          method = "spearman",
          na.action = na.exclude,
          alternative = "two.sided"
        )
      row2add <-
        cbind(
          "metadata" = metavar,
          "feature" = feature,
          "object_name" = obj.name,
          "rho" = spearman$estimate[[1]],
          "S" = spearman$statistic[[1]],
          "n" = length(na.omit(metadata[[metavar]])),
          "p" = spearman$p.value[[1]]
        )
      corr_output <- rbind(corr_output, row2add)
    }
  }
  # Remove NAs and add FDR (Benjamini Hochberg)
  statvars <- c("rho", "S", "n", "p")
  corr_output <-
    corr_output %>%
    na.omit() %>%
    mutate(
      across(all_of(statvars), as.character),
      across(all_of(statvars), as.numeric)
    ) %>%
    group_by(object_name, metadata) %>%
    mutate(q = p.adjust(p, method = "BH")) %>%
    ungroup()
  return(corr_output)
}


library(mRMRe)
# Correlation calculations
corr_abundance <-
  microbiome::transform(ps, "clr") %>% 
  microbiome::abundances() %>%
  t() %>% as.data.frame()

dd <- mRMR.data(corr_abundance[, blacklist])
blacklist_mRMRe <- mRMR.classic(
  data = dd,
  target_indices = c(1),
  # solution_count = 1,
  feature_count = 10
)
blacklist_mRMRe_features <-
  featureNames(blacklist_mRMRe)[solutions(blacklist_mRMRe)$`1`]

correlations_clr <-
  corr_loop(
    metadata = corr_abundance[, blacklist_mRMRe_features],
    abundance = corr_abundance,
    obj.name = "UHGG Species"
  )

saveRDS(
  correlations_clr,
  glue(
    "{wkdir}/data/interim/decontamination/",
    "{refDB}_{seq_method}_{level}_blacklist-correlations.rds"
  )
)

# correlations_clr %>%
#   filter(rho >= 0.6 & q <= 0.05 & feature %nin% blacklist_mRMRe_features) %>%
#   dim()

#   correlations_clr %>%
#   filter(rho >= 0.6 & q <= 0.05 ) %>%
#   pull(feature) %>% unique()

# correlations_clr %>%
#   filter(rho >= 0.5 & q <= 0.05 & feature %nin% blacklist) %>%
#   dim()

# correlations_clr %>%
#   filter(rho >= 0.6 & q <= 0.05) %>%
#   dim()

# ps_trim_corr_02 <- filter_taxa( ... )

corr_xy <- function(obj, corr_obj, feature_var, metadata_var) {
  #' Function creates a scatter plot of a given feature and a metadata column
  abund <- obj %>%
    abundances() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample_id")
  df.plot <- obj %>%
    meta() %>%
    left_join(abund, by = "sample_id")

  stat_col <-
    corr_obj %>%
    dplyr::filter(feature == feature_var) %>%
    dplyr::filter(metadata == metadata_var)
  stat_title <-
    paste0(
      "Spearman's Rho: ", round(stat_col$rho, digits = 3), "\n",
      "P-value: ", format(stat_col$p, digits = 3, scientific = T),
      "  FDR: ", format(stat_col$q, digits = 3, scientific = T)
    )
  cat("\nRho: ", stat_col$rho, ", ")
  cat("\nP-value: ", stat_col$p, ",  ")
  cat("\nQ-value: ", stat_col$q, "\n")

  df.plot %>%
    drop_na(metadata_var) %>%
    ggplot(aes(x = .data[[feature_var]], y = .data[[metadata_var]])) +
    geom_point(shape = 21, alpha = 1) +
    geom_smooth(method = lm, color = "darkgrey", linetype = "dotted", se = F) +
    theme_bw() +
    labs(x = feature_var, y = metadata_var, title = stat_title) +
    # scale_fill_manual(values = cols.pdpchc) +
    # scale_color_manual(values = cols.pdpchc.rim) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 12)
    )
}

ps_clr <- microbiome::transform(ps, "clr")
p_corr <- corr_xy(
  ps_clr, 
  corr_obj, 
  "S_s__Pseudomonas_E__massiliensis", 
  "S_s__Pseudomonas_E__fragi_B"
)
p_corr +
  labs(x = "Pseudomonas massiliensis", y = "Pseudomonas fragi B")
ggsave(
  glue("{wkdir}/figures/decontamination/{refDB}_{seq_method}_{level}_contamination-correlation-example.png"),
  width = 5,
  height = 5
)






























































#______________________________________________________________________________
library(ggforce)
# Load in Kraken Classification data ---
organisms_of_interest <- 
c("Bacteria",
  "Viruses",
  "Archaea",
  "Protozoa",
  "Fungi"
)
# read in wgs and kraken data
krak_rnaseq <- readRDS(
  glue("{wkdir}/data/interim/kraken_results/2023-01-04_RNASEQ_RefSeqPlusPF_kraken2_datatable.rds")
)
krak_wgs <- readRDS(
  glue("{wkdir}/data/interim/kraken_results/2023-01-04_WGS_RefSeqPlusPF_kraken2_datatable.rds")
)

sum_kraken_reads <- function(kraken_long, feature_list) {
  kraken_long %>%
    filter(feature %in% feature_list) %>%
    group_by(sample_id) %>%
    summarize(reads = sum(clade_counts))
}

classified_microbial_reads_wgs <- krak_wgs %>%
  sum_kraken_reads(organisms_of_interest) %>%
  mutate(
    seq_method = "WGS",
    analysis_stage = "Classified\n Microbial Reads"
  )
classified_microbial_reads_rnaseq <- krak_rnaseq %>%
  sum_kraken_reads(organisms_of_interest) %>%
  mutate(
    seq_method = "RNASEQ",
    analysis_stage = "Classified\n Microbial Reads"
  )


# Load in flagstat data ----
extract_penultimate_dir <- function(path) {
  path_split <- stringr::str_split(path, "/") %>% unlist()
  folder <- tail(path_split, 2)[1]
  return(unname(folder))
}
extract_penultimate_dir <- Vectorize(extract_penultimate_dir)

readqc_processed <- glue("{wkdir}/data/processed/readqc")
flagstat_dfs <-
  readRDS(glue(
    "{readqc_processed}/",
    "2023-01-06_flagstat-metrics_WGS-RNASEQ.rds"
  ))

flagstat_df_wgs <- flagstat_dfs[["WGS"]] %>%
  janitor::clean_names() %>% 
  mutate(flagstat_dir = extract_penultimate_dir(filepath)) %>%
  mutate(
    file = basename(filepath),
    sample_id = case_when(
      grepl("F3328", file) ~ gsub("_unmapped-F3328_flagstat.tsv", "", file),
      grepl("f4", file) ~ gsub("_unmapped-f4_flagstat.tsv", "", file),
      grepl("flagstat", file) ~ gsub("_flagstat.tsv", "", file),
      TRUE ~ "error"
    )
  )
flagstat_df_rnaseq <- flagstat_dfs[["RNASEQ"]] %>%
  janitor::clean_names() %>%
  mutate(flagstat_dir = extract_penultimate_dir(filepath)) %>%
  mutate(
    file = basename(filepath),
    sample_id = case_when(
      grepl("F3328", file) ~ gsub("_unmapped-F3328_flagstat.tsv", "", file),
      grepl("f4", file) ~ gsub("_unmapped-f4_flagstat.tsv", "", file),
      grepl("flagstat", file) ~ gsub("_flagstat.tsv", "", file),
      TRUE ~ "error"
    )
  )

# Load in BBDUK data ---
bbduk_summary_wgs <- readRDS(
  glue("{wkdir}/data/interim/readqc/2023-01-06_WGS_bbduk-stats-df.rds")) %>%
  select(sample_id, final_reads) %>%
  dplyr::rename(reads = final_reads) %>%
  mutate(analysis_stage = "04_bbduk-trim",
  reads = as.numeric(reads))

bbduk_summary_rnaseq <- readRDS(
  glue("{wkdir}/data/interim/readqc/2023-01-06_RNASEQ_bbduk-stats-df.rds")
) %>%
  select(sample_id, final_reads) %>%
  dplyr::rename(reads = final_reads) %>%
  mutate(analysis_stage = "04_bbduk-trim",
  reads = as.numeric(reads))

sample_metadata <- readRDS(
  glue("{wkdir}/data/interim/metadata/2023-01-06_phyloseq-metadata.rds")
)
sample_metadata %<>% purrr::map(
  ~ select(.,
  sample_id, participant_id, study, diagnosis_at_baseline,
  case_control_other_at_baseline, age_at_baseline,
  ethnicity, sex, race
  ))

# //TODO FIX ME HACK - resolve rnaseq metadata by joining with wgs static meta
core_meta <- sample_metadata[["WGS"]] %>% select(-sample_id)
sample_metadata[["RNASEQ"]] <-
  sample_metadata[["RNASEQ"]] %>%
  select(sample_id, participant_id) %>%
  left_join(core_meta)

reads_wgs <- flagstat_df_wgs %>%
  select(total_qc_passed_reads_qc_failed_reads, sample_id, flagstat_dir) %>%
  dplyr::rename(
    reads = total_qc_passed_reads_qc_failed_reads,
    analysis_stage = flagstat_dir
  ) %>%
  bind_rows(
    bbduk_summary_wgs,
    classified_microbial_reads_wgs
  ) %>% 
    left_join(sample_metadata[["WGS"]]) %>%
    mutate(seq_method = "WGS")

# sample_metadata[["WGS"]]$sample_id %>% unique()
# sample_metadata[["RNASEQ"]]$sample_id %>% unique()

reads_rnaseq <- flagstat_df_rnaseq %>%
  select(total_qc_passed_reads_qc_failed_reads, sample_id, flagstat_dir) %>%
  dplyr::rename(
    reads = total_qc_passed_reads_qc_failed_reads,
    analysis_stage = flagstat_dir
  ) %>%
  bind_rows(
    bbduk_summary_rnaseq,
    classified_microbial_reads_rnaseq
  ) %>% 
    left_join(sample_metadata[["RNASEQ"]], by = "sample_id") %>%
    mutate(seq_method = "RNASEQ")

mixed_reads_df <-
  bind_rows(
    reads_wgs,
    reads_rnaseq
  ) %>%
  mutate(
    analysis_stage = case_when(
      analysis_stage == "01_CRAM"
      ~ "All Reads",
      analysis_stage == "01_BAM"
      ~ "All Reads",
      analysis_stage == "02_BAM-unmapped"
      ~ "All Unmapped Reads",
      analysis_stage == "03_BAM-unmapped-3328-filtered"
      ~ "Unmapped Reads \n No Duplicates",
      analysis_stage == "04_bbduk-trim"
      ~ "Final Reads \n Quality Trimmed",
      TRUE ~ analysis_stage
    ),
    analysis_stage = factor(analysis_stage,
      levels =
        c(
          "All Reads",
          "All Unmapped Reads",
          "Unmapped Reads \n No Duplicates",
          "Final Reads \n Quality Trimmed",
          "Classified\n Microbial Reads"
        )
    ),
    seq_method = factor(seq_method, levels = c("WGS", "RNASEQ"))
  )

p_reads <- mixed_reads_df %>%
  ggplot(aes(
    x = analysis_stage,
    y = reads + 1
    )) +
  geom_jitter(aes(color = study), alpha = 0.5, size = 0.1) +
    geom_boxplot(alpha = 0.6, outlier.alpha = 0) +
    scale_color_nejm() +
    labs(x = "", y = "") +
    guides(
      color =
        guide_legend(override.aes = list(size = 4))
    ) +
    theme_bw() +
      facet_wrap(~seq_method) +
      scale_y_log10(breaks = c(1 %o% 10^(1:10))) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

# ggsave(p_reads,
#   filename =
#   glue("{wkdir}/figures/readqc/{Sys.Date()}_flagstat-bbdukqc-reads.png"),
#   dpi = 600, height = 6, width = 10
# )


#' Calculating the mean, median, and std dev for classified microbial reads
mixed_reads_df %>%
  group_by(seq_method, analysis_stage) %>%
  summarise(
    mean = mean(reads, na.rm = TRUE),
    median = median(reads, na.rm = TRUE),
    std_dev = sd(reads, na.rm = TRUE)
  )

# pivot data into wide format
pivot_bbduk_reads_df <- function(df) {
  df %>%
  pivot_wider(names_from = "analysis_stage", values_from = "reads") %>%
  janitor::clean_names() %>%
  mutate(filter_step_1 = all_reads - all_unmapped_reads,
  filter_step_2 =  all_unmapped_reads - unmapped_reads_no_duplicates,
  filter_step_3 = unmapped_reads_no_duplicates - final_reads_quality_trimmed,
  filter_step_4 = final_reads_quality_trimmed - classified_microbial_reads
  )
}

mixed_reads_rates <- mixed_reads_df %>% 
  pivot_bbduk_reads_df() %>% 
  mutate(perc_microbial = ((classified_microbial_reads * 2 * 100)/all_reads))


p_read_perc <- mixed_reads_rates %>%
  ggplot(aes(x = seq_method, y = perc_microbial)) +
  labs(x = "", y = "Percent of all reads classified as microbial") +
    facet_zoom(ylim = c(0, 0.03)) +
    theme_bw() +
    geom_jitter(aes(color = study), alpha = 1, size = 0.7) +
    scale_color_nejm() +
        guides(
      color =
        guide_legend(override.aes = list(size = 4))
    )

# ggsave(p_read_perc,
#   filename =
#   glue("{wkdir}/figures/readqc/{Sys.Date()}_final-read-percentages.png"),
#   dpi = 600, height = 6, width = 10
# )

rate_stat_summary <- mixed_reads_rates %>%
  group_by(seq_method) %>%
  summarise_at(
    vars(filter_step_1:filter_step_4),
    list(
      mean = mean,
      median = median,
      min = min,
      max = max
    ),  na.rm = TRUE
  )
View(rate_stat_summary)

mixed_reads_rates %>%
  ggplot(aes(x = ))


mixed_reads_rates %>%
  group_by(seq_method) %>%
  summarise(
    mean = mean(perc_microbial, na.rm = TRUE),
    median = median(perc_microbial, na.rm = TRUE),
    std_dev = sd(perc_microbial, na.rm = TRUE)
  )

#' Calculating the mean, median, and std dev for the fraction reads that are microbial
mixed_reads_rates %>%
  group_by(seq_method) %>%
  summarise(
    mean = mean(perc_microbial, na.rm = TRUE),
    median = median(perc_microbial, na.rm = TRUE),
    std_dev = sd(perc_microbial, na.rm = TRUE)
  )













# ggplot() +
#   scale_x_log10() +
#   geom_histogram(
#     data = mixed_reads_rates %>% filter(seq_method == "WGS"),
#     aes(((classified_microbial_reads * 2 * 100) / all_reads)),
#     bins = 500, alpha = 0.4, fill = "red"
#   ) +
#   geom_histogram(
#     data = mixed_reads_rates %>% filter(seq_method == "RNASEQ"),
#     aes(((classified_microbial_reads * 2 * 100) / all_reads)),
#     bins = 500, alpha = 0.4, fill = "blue"
#   )

# # Exploring microbal taxonomic assignment distributions
# ggplot() +
#   geom_histogram(
#     data = classified_microbial_reads_wgs,
#     aes(reads), bins = 500, fill = "red", alpha = 0.4
#   ) +
#   geom_histogram(
#     data = classified_microbial_reads_rnaseq,
#     aes(reads), bins = 500, fill = "blue", alpha = 0.4
#   ) +
#     scale_x_log10()


# ggsave(p_reads_rnaseq,
#   filename =
#   glue("{wkdir}/figures/readqc/{Sys.Date()}_RNASEQ_flagstat-bbdukqc-reads.png"),
#   dpi = 600, height = 8, width = 5
# )



#' Determine the fraction of final reads that are present
#' in original data and rates of removal for each step
#' THen, compare the fraction of unmapped reads (self-controlled) between WGS and RNA



pdmbs_dir <- "/central/groups/MazmanianLab/joeB/PDMBS"
wkdir <- paste0(pdmbs_dir, "/parkinsons-microbial-blood-signatures")
source(paste0(wkdir, "/notebooks/R_scripts/_misc_functions.R"))
source(paste0(wkdir, "/notebooks/R_scripts/_load-core-pkgs.R"))
library(janitor)
library(progress)

# library(readr)

run_status <- list()
slurm_out <- glue("{wkdir}/.cluster_runs/gcs_upload_HMP2")
stderr_paths <- list.files(slurm_out, full.names = TRUE) %>% keep(grepl(".err", .))
pb <- progress_bar$new(total = length(stderr_paths))

for (f in stderr_paths) {
  pb$tick()
  stderr <- read_tsv(f, col_names = FALSE,
    lazy = should_read_lazy(),
    num_threads = 4,
    show_col_types = FALSE, progress = FALSE)
  if (stderr$X1 %>% grepl("Operation completed over 1 objects", .) %>% any()) {
    run_status[[f]] <- "SUCCESS"
  } else {
    run_status[[f]] <- "ERROR"
  }
}

cat(
  "Workflow status:",
  "\n SUCCESS:",
  run_status %>% keep(grepl("SUCCESS", .)) %>% length(),
  "\n ERROR:",
  run_status %>% keep(grepl("ERROR", .)) %>% length(),
  "\n"
)











# final_assemblies_gc_stats$filename %>% unique()
# cdr3_gc_stats$filename %>% unique()






pdmbs_dir <- "/central/groups/MazmanianLab/joeB/PDMBS"
wgs_wkdir <- paste0(pdmbs_dir, "/workflow/WGS")
wkdir <- paste0(pdmbs_dir, "/parkinsons-microbial-blood-signatures")
source(paste0(wkdir, "/notebooks/R_scripts/_misc_functions.R"))
source(paste0(wkdir, "/notebooks/R_scripts/_load-core-pkgs.R"))




# nextflow run ./nf-core-ninjaindex/main.nf --genomes 's3://bucket/input/*.fna' --outdir 's3://bucket/output/' -profile aws

hcom2 <- readxl::read_xlsx(glue("{wkdir}/1-s2.0-S0092867422009904-mmc2 (2).xlsx"), sheet = 'hCom2')
gzfiles <- hcom2 %>% pull(Public_URI) %>% keep(grepl(".gz", .))
for (f in gzfiles) {
  shell_do(glue("wget -P /central/groups/MazmanianLab/shared/reference_genomes/hCom2/ {f}"))
}
gca_files <- read.table("gca_loc.txt", sep = "\t")
for (f in gca_files$V1) {
  shell_do(glue("wget {f} -P /central/groups/MazmanianLab/shared/reference_genomes/hCom2/"))
}






# assemblies_genbank <- read.delim("assembly_summary_genbank.txt",
#   stringsAsFactors = F, header = T, skip = 1
# )
# assemblies_refseq <- read.delim("assembly_summary_refseq.txt",
#   stringsAsFactors = F, header = T, skip = 1
# )
# assemblies_genbank %>% glimpse()
# assemblies_refseq %>% glimpse()
# assemblies <- bind_rows(assemblies_genbank, assemblies_refseq)


# library(SRAdb)
# if(!file.exists('/central/groups/MazmanianLab/joeB/Downloads/SRAmetadb.sqlite')) {
#   system("wget https://gbnci.cancer.gov/sra/SRAmetadb.sqlite.gz -P /central/groups/MazmanianLab/joeB/Downloads/")
#   system("gunzip /central/groups/MazmanianLab/joeB/Downloads/SRAmetadb.sqlite.gz")
#   sqlfile <<- 'SRAmetadb.sqlite'
# }
# sqlfile <<- 'SRAmetadb.sqlite'
# sra_con <- dbConnect(SQLite(), sqlfile )
# sra_tables <- dbListTables(sra_con)
# sra_tables

# dbListFields(sra_con,"study")
# conversion<-sraConvert(c('SRP001007','SRP000931'),sra_con=sra_con)






# gca_files <- list(
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/152/405/GCF_025152405.1_ASM2515240v1/GCF_025152405.1_ASM2515240v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/152/275/GCF_025152275.1_ASM2515227v1/GCF_025152275.1_ASM2515227v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/025/150/425/GCA_025150425.1_ASM2515042v1/GCA_025150425.1_ASM2515042v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/149/465/GCF_025149465.1_ASM2514946v1/GCF_025149465.1_ASM2514946v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/025/146/925/GCA_025146925.1_ASM2514692v1/GCA_025146925.1_ASM2514692v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/152/575/GCF_025152575.1_ASM2515257v1/GCF_025152575.1_ASM2515257v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/151/995/GCF_025151995.1_ASM2515199v1/GCF_025151995.1_ASM2515199v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/151/715/GCF_025151715.1_ASM2515171v1/GCF_025151715.1_ASM2515171v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/151/535/GCF_025151535.1_ASM2515153v1/GCF_025151535.1_ASM2515153v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/151/385/GCF_025151385.1_ASM2515138v1/GCF_025151385.1_ASM2515138v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/146/415/GCF_025146415.1_ASM2514641v1/GCF_025146415.1_ASM2514641v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/151/215/GCF_025151215.1_ASM2515121v1/GCF_025151215.1_ASM2515121v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/151/045/GCF_025151045.1_ASM2515104v1/GCF_025151045.1_ASM2515104v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/148/285/GCF_025148285.1_ASM2514828v1/GCF_025148285.1_ASM2514828v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/150/895/GCF_025150895.1_ASM2515089v1/GCF_025150895.1_ASM2515089v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/025/148/965/GCA_025148965.1_ASM2514896v1/GCA_025148965.1_ASM2514896v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/149/285/GCF_025149285.1_ASM2514928v1/GCF_025149285.1_ASM2514928v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/150/565/GCF_025150565.1_ASM2515056v1/GCF_025150565.1_ASM2515056v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/150/745/GCF_025150745.1_ASM2515074v1/GCF_025150745.1_ASM2515074v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/149/125/GCF_025149125.1_ASM2514912v1/GCF_025149125.1_ASM2514912v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/150/085/GCF_025150085.1_ASM2515008v1/GCF_025150085.1_ASM2515008v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/150/245/GCF_025150245.1_ASM2515024v1/GCF_025150245.1_ASM2515024v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/149/915/GCF_025149915.1_ASM2514991v1/GCF_025149915.1_ASM2514991v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/149/785/GCF_025149785.1_ASM2514978v1/GCF_025149785.1_ASM2514978v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/149/625/GCF_025149625.1_ASM2514962v1/GCF_025149625.1_ASM2514962v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/148/785/GCF_025148785.1_ASM2514878v1/GCF_025148785.1_ASM2514878v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/148/445/GCF_025148445.1_ASM2514844v1/GCF_025148445.1_ASM2514844v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/148/635/GCF_025148635.1_ASM2514863v1/GCF_025148635.1_ASM2514863v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/148/125/GCF_025148125.1_ASM2514812v1/GCF_025148125.1_ASM2514812v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/147/905/GCF_025147905.1_ASM2514790v1/GCF_025147905.1_ASM2514790v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/147/765/GCF_025147765.1_ASM2514776v1/GCF_025147765.1_ASM2514776v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/147/655/GCF_025147655.1_ASM2514765v1/GCF_025147655.1_ASM2514765v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/147/485/GCF_025147485.1_ASM2514748v1/GCF_025147485.1_ASM2514748v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/147/325/GCF_025147325.1_ASM2514732v1/GCF_025147325.1_ASM2514732v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/147/085/GCF_025147085.1_ASM2514708v1/GCF_025147085.1_ASM2514708v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/146/775/GCF_025146775.1_ASM2514677v1/GCF_025146775.1_ASM2514677v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/146/565/GCF_025146565.1_ASM2514656v1/GCF_025146565.1_ASM2514656v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/146/315/GCF_025146315.1_ASM2514631v1/GCF_025146315.1_ASM2514631v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/146/135/GCF_025146135.1_ASM2514613v1/GCF_025146135.1_ASM2514613v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/146/005/GCF_025146005.1_ASM2514600v1/GCF_025146005.1_ASM2514600v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/145/845/GCF_025145845.1_ASM2514584v1/GCF_025145845.1_ASM2514584v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/145/645/GCF_025145645.1_ASM2514564v1/GCF_025145645.1_ASM2514564v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/145/285/GCF_025145285.1_ASM2514528v1/GCF_025145285.1_ASM2514528v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/144/995/GCF_025144995.1_ASM2514499v1/GCF_025144995.1_ASM2514499v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/144/665/GCF_025144665.1_ASM2514466v1/GCF_025144665.1_ASM2514466v1_genomic.fna.gz",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/144/545/GCF_025144545.1_ASM2514454v1/GCF_025144545.1_ASM2514454v1_genomic.fna.gz"
# )









kraken_long <- list()
for (seq_method in c("RNASEQ", "WGS")) {
  for (refDB in c("UHGG", "RefSeqPlusPF")) {
    # for (taxa_rank_filter in names(taxa_subsets)) {
      taxa_rank_filter <- "All"
      kraken_long[[seq_method]][[refDB]][[taxa_rank_filter]] <- readRDS(
        glue(
          "{outputdir}/",
          "2023-01-04_{seq_method}_{refDB}_kraken2_datatable.rds"
        )
      )
    message(glue("{seq_method} {refDB} {taxa_rank_filter}"))
    kraken_long[[seq_method]][[refDB]][[taxa_rank_filter]] %>% dim %>% print
    # }
  }
}

kraken_long$RNASEQ$UHGG %>% names
kraken_long$RNASEQ$UHGG$Species %>% head()

# "All"     "Species" "Genus"   "Phylum"



tst <- readRDS(glue("{wkdir}/data/interim/gcs_locations/2022-12-23_URLs_unmapped-bam.rds"))



rnafqs <- list.files(
  glue("{pdmbs_dir}/workflow/RNASEQ/clean_reads"),
  full.names = TRUE
)

rnafqs_sizes <- rnafqs %>%
  purrr::set_names() %>%
  purrr::map(~ file.size(.))

View(rnafqs_sizes)

fq_counts <- rnafqs %>%
  basename() %>%
  strex::str_before_last("_") %>%
  table
  


run_list <- sampleIDs_rna %>%
  discard(
    ~ file.exists(
      file.path(data_out, glue("UHGG_reports/{.}_UHGG_report.tsv"))
    )
  )
  









# #' Since some samples have more than one sequencing runs we will need to
# #' reformat our data to have
# library(data.table)
# core_seq_meta <-
#   wgs_sequencing_meta_formatted %>%
#   select(SAMPLE, FLOWCELL_BARCODE) %>%
#   distinct() %>%
#   as.data.table()

# core_seq_meta_formatted <-
#   core_seq_meta[, lapply(.SD, paste0, collapse = ","), by = SAMPLE] %>%
#   mutate(SAMPLE = gsub("SM:", "", SAMPLE)) %>%
#   arrange(SAMPLE = sample_names(ps))

# # double check sample names align
# if (all(core_seq_meta_formatted$SAMPLE == sample_names(ps))){
#   sample_data(ps)$FLOWCELL_BARCODE <- core_seq_meta_formatted$FLOWCELL_BARCODE
# } else {
#   cat("\n ERROR in sample ordering!!! \n")
# }



#' visualize how various thresholds result in different numbers of
#' blacklist microbes

# blacklist_rarefaction <- tibble()
# for (n_threshold in seq(1, 125, 1)) {
#   prevalence_stat_summary <- prevalence_stats %>%
#     filter(sample_n > n_threshold) %>%
#     group_by(ind, batch_variable) %>%
#     dplyr::summarise(
#       min = min(values), max = max(values),
#       prevalence_delta = max - min
#     ) %>%
#     filter(max > 0.5, min <= max / 2)
#   blacklist <- prevalence_stat_summary %>%
#     pull(ind) %>%
#     unique()
#   blacklist_rarefaction %<>% bind_rows(
#     as_tibble_row(
#       list(
#         "group_size_threshold" = n_threshold,
#         "blacklist_microbes_n" = length(blacklist)
#       )
#     )
#   )
# }

# blacklist_rarefaction %>%
#   ggplot() +
#   geom_point(aes(x = group_size_threshold, y = blacklist_microbes_n))



# ______________________________________________________________________________
#

# blacklist <- readRDS(
#   glue("{wkdir}/data/interim/decontamination/",
#   "2023-01-09_{seq_method}_{refDB}_{level}_blacklist.rds")
#  )

# ps_trim_prev_01 <- microbiome::remove_taxa(ps, taxa = blacklist)

# prevalence_stats %>%
#   filter(batch_variable == "flowcell_barcode") %>%
#   filter(ind == "S_s__Pseudomonas_E__aeruginosa") %>%
#   ggplot(aes(fct_reorder(batch_subset, -values), values)) +
#   geom_col() +
#   # theme_bw() +
#   labs(y = "Pseudomonas aeruginosa Prevalence", x = "Flowcell barcode") +
#   theme(axis.text.x = element_blank())

# ggsave(
#   glue("{wkdir}/figures/decontamination/{refDB}_{seq_method}_{level}_S_s__Pseudomonas_E__aeruginosa.png"),
#   width = 7,
#   height = 4.5
# )

# Abundance vs Prevalance plot
# library(ggside)
# stats_prevalence <-
#   data.frame("prevalence" = microbiome::prevalence(ps)) %>%
#   rownames_to_column(var = "feature")
# stats_mean_clr <- ps %>%
#   microbiome::transform("clr") %>%
#   abundances() %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "feature") %>%
#   pivot_longer(!feature, names_to = "sample_id") %>%
#   group_by(feature) %>%
#   dplyr::summarise(clr_mean = mean(value), clr_sd = sd(value))
# stats_mean_log_counts <- ps %>%
#   microbiome::transform("log") %>%
#   abundances() %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "feature") %>%
#   pivot_longer(!feature, names_to = "sample_id") %>%
#   group_by(feature) %>%
#   dplyr::summarise(log_abund_mean = mean(value), log_abund_sd = sd(value))

# taxa_summary_stats <-
#   full_join(stats_prevalence, stats_mean_clr) %>%
#   full_join(stats_mean_log_counts) %>%
#   mutate(filtered = case_when(
#     feature %in% blacklist ~ "Prevalence Filtered",
#     feature %in% blacklist_readcounts ~ "Read Count Filtered",
#     TRUE ~ "No"
#     ))

# taxa_summary_stats %>%
#   ggplot(aes(x = prevalence, y = clr_mean+1)) +
#   geom_point(aes(fill = filtered, alpha = filtered),
#     size = 2, shape = 21, stroke = 0.2) +
#   geom_xsidedensity(aes(y = stat(density))) +
#   geom_ysidedensity(aes(x = stat(density))) +
#   scale_x_log10() +
#   scale_y_log10() +
#   scale_alpha_manual(values =
#   c("No" = 0.1, "Prevalence Filtered" = 0.1, "Read Count Filtered" = 1)) +
#   theme_bw()

# ggsave(
#   glue("{wkdir}/figures/decontamination/{refDB}_{seq_method}_{level}_abundance-vs-prevalence-filter-status.png"),
#   width = 7,
#   height = 5
# )

# bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | --sra-acc <acc> | b <bam>} -S [<sam>]




# bowtie2_cmd <- glue(
#   "bowtie2 -x {bowtiew_ind}",
#   " -b /central/groups/MazmanianLab/joeB/PDMBS/SY-PDZH104KR2.sam",
#   " -S /central/groups/MazmanianLab/joeB/PDMBS/REMAPPED_SY-PDZH104KR2.sam",
#   # " --align-paired-reads",
#   " --preserve-tags"
# )




# file_sizes <- fq_paths %>% file.size()/1000000
# file_sizes %>% hist(100)


# slurmstepd: error: Detected 1 oom-kill event(s) in StepId=36553876.batch. Some of your processes may have been killed by the cgroup out-of-memory handler.
# slurmstepd: error: poll(): Bad address







tst <- read.delim(
  glue("/central/groups/MazmanianLab/joeB/PDMBS/workflow/RNASEQ/stats_clean_reads_bowtie2_CHM13/PP-92834-SVM6T1_stderr.txt"),
  header = FALSE
)
suppressWarnings(has_error_message(
  glue("/central/groups/MazmanianLab/joeB/PDMBS/workflow/RNASEQ/stats_clean_reads_bowtie2_CHM13/PP-92834-SVM6T1_stderr.txt")
))

fsize <- sampleIDs %>%
  purrr::set_names() %>%
  purrr::map( ~ file.size(glue("/central/scratch/jbok/PDMBS/WGS/bowtie2_clean_reads/{.}_R1.fastq.gz")))

fsize %>% unlist() %>% sort()  %>%  head(200)



# # KrakenUniq test run command
# krakenuniq --db /central/groups/MazmanianLab/joeB/Downloads/RefDBs/KrakenUniq/MicrobialDB \
# --threads 8 \
# --paired \
# --preload \
# --output testrun_KrakenUniq_read-classification.tsv \
# --report-file testrun_KrakenUniq_report.tsv \
# --only-classified-output \
# --classified-out testrun_classified.fastq \
# /central/scratch/jbok/PDMBS/WGS/bowtie2_clean_reads/PP-40543_R1.fastq.gz \
# /central/scratch/jbok/PDMBS/WGS/bowtie2_clean_reads/PP-40543_R2.fastq.gz


# krakenuniq --db /central/groups/MazmanianLab/joeB/Downloads/RefDBs/KrakenUniq/MicrobialDB \
# --threads 8 \
# --paired \
# --preload \
# --output testrun2_KrakenUniq_read-classification.tsv \
# --report-file testrun2_KrakenUniq_report.tsv \
# /central/scratch/jbok/PDMBS/WGS/bowtie2_clean_reads/PP-40543_R1.fastq.gz \
# /central/scratch/jbok/PDMBS/WGS/bowtie2_clean_reads/PP-40543_R2.fastq.gz



# seq_method <- "WGS"






# kraken_uniq_results_df %>% View




# theme_clean <- function() {
#   font <- "Helvetica" # assign font family up front
#   theme_bw() %+replace% # replace elements we want to change
#     theme(
#       # grid elements
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       axis.ticks = element_blank(),
#       axis.title = element_text(
#         family = font,
#         size = 10
#       ),
#       axis.text = element_text(
#         family = font,
#         size = 9
#       ),
#       axis.text.x = element_text(
#         margin = margin(5, b = 10)
#       )
#     )
# }


blacklist_readcounts <- readRDS(
  glue(
    "{wkdir}/data/interim/decontamination/",
    "2023-10-19_WGS_KrakenUnique_blacklist_readcounts.rds"
  )
)
blacklist_prev <- readRDS(
  glue(
    "{wkdir}/data/interim/decontamination/",
    "2023-10-19_WGS_KrakenUnique_prevalence_blacklist.rds"
  )
) %>%
  as.character()

contaminant_lineages <- readRDS(
  glue(
    "{wkdir}/data/interim/decontamination/",
    "2024-02-05_reagent_blacklist_contaminant_lineages.rds"
  )
)
blacklist <- c(
  blacklist_prev,
  blacklist_readcounts,
  contaminant_lineages$taxid
) %>%
  unique()

kuq <- readRDS(
    glue(
        "{wkdir}/data/interim/kraken_results/",
        "2024-02-05_WGS_KrakenUniq_results_with_lineage.rds"
    )
)

kuq_decon <- readRDS(
    glue(
        "{wkdir}/data/interim/kraken_results/",
        "2024-02-05_WGS_KrakenUniq_results_with_lineage_decon.rds"
    )
)

# kraken_uniq_results_df %>% glimpse
# taxids <- kuq$taxID %>% unique()

sample_order <-
  kuq %>%
  filter(taxName == "unclassified") %>%
  arrange(desc(percent_reads)) %>%
  pull(sample_id)

kuq %<>%
  mutate(
    sample_id =
      factor(sample_id, levels = sample_order, ordered = TRUE)
  )

kuq_decon %<>%
  mutate(
    sample_id =
      factor(sample_id, levels = sample_order, ordered = TRUE)
  )
  
read_sum <- kuq %>%
  filter(taxName %in% c("root", "unclassified")) %>%
  group_by(sample_id) %>%
  summarize_at(vars(percent_reads, reads, kmers), sum)

p_readstat <- read_sum %>%
  ggplot(aes(y = reads, x = sample_id)) +
  geom_point(alpha = 0.2, size = 1) +
  theme_minimal() +
  scale_y_log10() +
    labs(y = "Quality \nReads", x = "Samples", fill = NULL) +
    theme(
      axis.text.x = element_blank(),
      panel.grid = element_blank()
    )


p_class_unclass <- kuq %>%
  filter(taxName %in% c("unclassified", "root")) %>%
  mutate(
    taxName =
      factor(taxName, levels = c("unclassified", "root"), ordered = TRUE)
  ) %>%
  ggplot(aes(
    x = sample_id,
    y = percent_reads, fill = taxName
  )) +
  geom_bar(stat = "identity") +
    scale_fill_manual(
      values =
        c("root" = "#1b1b1b", "unclassified" = "#7d7a7aea")
    ) +
    theme_minimal() +
    labs(y = "Reads (%)", x = NULL, fill = NULL) +
    theme(
      axis.text.x = element_blank(),
      panel.grid = element_blank()
    )


p_superking <- kuq %>%
  filter(rank == "superkingdom") %>%
  ggplot(aes(
    x = sample_id,
    y = percent_reads, fill = taxName
  )) +
  geom_bar(stat = "identity") +
    scale_fill_d3() +
    theme_minimal() +
    labs(y = "Reads (%)", x = NULL, fill = NULL) +
    theme(
      axis.text.x = element_blank(),
      panel.grid = element_blank()
    )

p_superking_decon <- kuq_decon %>%
  filter(rank == "superkingdom") %>%
  ggplot(aes(x = sample_id, y = percent_reads)) +
  geom_bar(stat = "identity", fill = "#FF7F0F") +
  theme_minimal() +
  labs(y = "Decontaminated\n Reads (%)", x = NULL, fill = NULL) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  )


clade_summary_stats <- kuq_decon %>%
  filter(grepl("; Bacteria;", lineage)) %>%
  filter(rank == "genus") %>%
  group_by(lineage) %>%
  summarize_at(vars(reads), sum)

top_clade_list <- clade_summary_stats %>%
  slice_max(order_by = reads, n = 10, with_ties = FALSE) %>%
  pull(lineage) %>%
  unique()

remaining_clade_list <- clade_summary_stats$lineage %>%
  setdiff(top_clade_list) %>%
  unique()


# loop through lineage ids, and create a df
#  summing each of the taxa within the clade level of interest

sid_clade_sum <- tibble()
taxid_list <- list()
tic()
for (clade in unique(top_clade_list)) {
  message("Processing: ", clade, "...")
  clade_data <- kuq_decon %>%
    filter(grepl(clade, lineage))
  taxid_list[[clade]] <- clade_data$taxid %>% unique()
  sid_clade_sum %<>% bind_rows(
    clade_data %>%
      group_by(sample_id) %>%
      dplyr::summarize(taxReads = sum(reads)) %>%
      mutate(clade_id = clade)
  )
}
toc()

all_other_clade_sum <- kuq_decon %>%
  filter(taxid %nin% unname(unlist(taxid_list))) %>%
  # filter clade levels of interest
  filter(rank %in% c(
    "order", "family", "genus",
    "species", "species group", "subspecies", "strain"
  )) %>%
  # filter out unidentified nucleotides
  filter(taxid  %nin% c("28384")) %>%
  group_by(sample_id) %>%
  dplyr::summarize(taxReads = sum(reads)) %>%
  mutate(clade_id = "Other")

clade_stats <-
  sid_clade_sum %>%
  mutate(
    clade_id =
      str_before_last(clade_id, ";") %>% str_after_last("; ")
  ) %>%
  bind_rows(all_other_clade_sum) %>%
  group_by(sample_id) %>%
  dplyr::mutate(
    clade_percent_reads = (taxReads / sum(taxReads)) * 100
  )

p_clade_relab <- clade_stats %>%
  ggplot(aes(
    x = sample_id,
    y = clade_percent_reads, fill = clade_id
  )) +
  geom_bar(stat = "identity") +
    scale_fill_d3(palette = "category20") +
    theme_minimal() +
    labs(y = "Bacteria Genera \n Relative Abundance", x = NULL, fill = NULL) +
    theme(
      axis.text.x = element_blank(),
      panel.grid = element_blank()
    )
# p_clade_relab

p_mapping_summary <- p_readstat %>%
  aplot::insert_top(p_class_unclass) %>%
  aplot::insert_top(p_superking) %>%
  aplot::insert_top(p_superking_decon) %>%
  aplot::insert_top(p_clade_relab, height = 3)

ggsave(
  glue(
    "{wkdir}/figures/readqc/",
    "{Sys.Date()}_KrakenUniq_mapping_summary.png"
    ),
  p_mapping_summary,
  width = 12, height = 9
)






# ______________________________________________________________________
# Domain level ECDFs for groups/ cohorts

kuq_decon <- readRDS(
  glue(
    "{wkdir}/data/interim/kraken_results/",
    "2024-02-25_WGS_KrakenUniq_results_with_lineage_decon.rds"
  )
)

kuq_dcon_labeled_clades <- kuq_decon %>%
  mutate(coi = case_when(
    grepl("; Bacteria;", lineage) | name == "Bacteria" ~ "Bacteria",
    grepl("; Archaea;", lineage) ~ "Archaea",
    grepl("; Viruses;", lineage) | name == "Viruses" ~ "Viruses",
    grepl("; Fungi;", lineage) | name == "Fungi" ~ "Fungi",
    TRUE ~ "OTHER"
  )) %>%
  filter(coi != "OTHER")

kuq_dcon_labeled_clades %>% glimpse

# read in sample summary stats across analysis levels
flagstat_df <- readRDS(
  glue(
    "{wkdir}/data/processed/readqc/",
    "2023-01-06_flagstat-metrics_WGS-RNASEQ.rds"
  )
)
wgs_flagstat_metrics <- flagstat_df[["WGS"]] %>%
  mutate(
    analysis_stage = filepath %>%
      str_before_last("/") %>% str_after_last("/")
  )  %>%
  mutate(
    sample_id = case_when(
      analysis_stage ==
        "01_CRAM" ~ basename(filepath) %>%
        str_remove("_flagstat.tsv"),
      analysis_stage ==
        "02_BAM-unmapped" ~ basename(filepath) %>%
        str_remove("_unmapped-f4_flagstat.tsv"),
      analysis_stage ==
        "03_BAM-unmapped-3328-filtered" ~ basename(filepath) %>%
        str_remove("_unmapped-F3328_flagstat.tsv"),
      TRUE ~ "ERROR"
    )
  ) %>%
    glimpse

wgs_flagstat_metrics %>%
  mutate(analysis_stage = glue("flagstat_{analysis_stage}")) %>%
  pivot_wider(
    names_from = analysis_stage,
    values_from = !sample_id,
      names_glue = "{analysis_stage}_{.value}",
  ) %>%
    janitor::clean_names() %>%
  glimpse()


wgs_flagstat_metrics %>% glimpse
wgs_flagstat_og <- wgs_flagstat_metrics %>%
  filter(analysis_stage == "01_CRAM")


# adding sample metadata to KrakenUniq report
ps_meta <- readRDS(
  glue("{wkdir}/data/interim/metadata/2023-07-14_phyloseq-metadata.rds")
)
wgs_meta <- ps_meta[["WGS"]] %>%
  dplyr::select(
    participant_id, sample_id, study,
    case_control_other_latest, sex,
    ethnicity, race, age_at_baseline
  )

coi_df <- kuq_dcon_labeled_clades %>%
  dplyr::group_by(sample_id, coi) %>%
  dplyr::summarise(clade_counts = sum(reads)) %>%
  left_join(wgs_meta, by = "sample_id") %>%
  left_join(wgs_flagstat_og, by = "sample_id") %>%
  ungroup() %>%
  mutate(clade_counts_by_humreads = clade_counts / mapped) %>%
  mutate(lognorm_cladecounts_by_hum = log10(clade_counts_by_humreads))

p_ecdf_study_clade <- coi_df %>%
  ggplot(aes(lognorm_cladecounts_by_hum,
    color = case_control_other_latest
  )) +
  stat_ecdf(geom = "point", alpha = 0.6, shape = 3) +
    theme_bw() +
    scale_color_d3() +
    labs(
      x = expression(log[10](italic("KrakenUniq read counts") /
        italic("BWA mapped human reads"))),
      y = "Empirical Cumulative \n Density Function",
      color = NULL
    ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  facet_grid(
    cols = vars(coi),
    rows = vars(study),
    scales = "free"
  )

ggsave(
  glue(
    "{wkdir}/figures/readqc/",
    "{Sys.Date()}_KrakenUniq_ECDFs_by-study-clade.png"
  ),
  p_ecdf_study_clade,
  width = 9, height = 9
)


p_ecdf_clade <- coi_df %>%
  ggplot(aes(lognorm_cladecounts_by_hum, color = case_control_other_latest)) +
  stat_ecdf(geom = "point", alpha = 0.4, shape = 3, size = 0.4) +
  theme_bw() +
  scale_color_d3() +
  labs(
    x = expression(log[10](italic("KrakenUniq read counts") /
      italic("BWA mapped human reads"))),
    y = "Empirical Cumulative \n Density Function",
    color = NULL
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  facet_grid(
    cols = vars(coi),
    scales = "free"
  ) +
  theme(legend.position = "top")

ggsave(
  glue(
    "{wkdir}/figures/readqc/",
    "{Sys.Date()}_KrakenUniq_ECDFs_by-clade.png"
  ),
  p_ecdf_clade,
  width = 7, height = 4
)



# length(blacklist) / (kuq$taxid %>% unique() %>% length()) * 100
# setdiff(kuq$taxid %>% unique(), blacklist) %>% length

kuq_dcon_labeled_clades %>%
  filter(coi != "Other") %>%
  pull(taxid) %>%
  unique() %>% 
  length

p_read_count_densitydist <- coi_df %>%
  ggplot(aes(clade_counts, color = case_control_other_latest)) +
  geom_density(linewidth = 1, alpha = 0.8) +
  scale_x_log10() +
  theme_bw() +
  scale_color_d3() +
  facet_grid(
    cols = vars(coi),
    scales = "free"
  ) +
  labs(x = "Read Counts", y = "Density", color = NULL) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave(
  glue(
    "{wkdir}/figures/readqc/",
    "{Sys.Date()}_KrakenUniq_density-dist-clade.png"
  ),
  p_read_count_densitydist,
  width = 7, height = 4
)



library(ggside)

pbase_ageby_reads <- coi_df %>%
  ggplot(aes(
    y = lognorm_cladecounts_by_hum,
    x = age_at_baseline, color = case_control_other_latest
  )) +
  geom_point(
    aes(fill = case_control_other_latest),
    shape = 21, alpha = 0.9, size = 0.9, color = "white") +
    labs(
      y = expression(log[10](italic("KrakenUniq read counts") /
        italic("BWA mapped human reads"))),
      x = "Age at Baseline",
      color = NULL
    ) +
    scale_fill_d3() +
    scale_color_d3() +
    geom_smooth(method = "lm") +
    geom_xsidedensity(aes(y = after_stat(density))) +
    geom_ysidedensity(aes(x = after_stat(density))) +
    ggside(collapse = "x") +
    theme_bw() +
    guides(color = NULL) +
    theme(legend.position = "right")

p_age_by_cladecount_norm  <-  pbase_ageby_reads +
    facet_grid(
      rows = vars(coi),
      # cols = vars(study),
      scales = "free_y"
    )
p_age_by_cladecount_norm_study  <-  pbase_ageby_reads +
    facet_grid(
      rows = vars(coi),
      cols = vars(study),
      scales = "free_y"
    )

ggsave(
  glue(
    "{wkdir}/figures/readqc/",
    "{Sys.Date()}_KrakenUniq_clade-counts-by-age.png"
  ),
  p_age_by_cladecount_norm,
  width = 6, height = 8
)
ggsave(
  glue(
    "{wkdir}/figures/readqc/",
    "{Sys.Date()}_KrakenUniq_clade-counts-by-age-cohort.png"
  ),
  p_age_by_cladecount_norm_study,
  width = 20, height = 8
)







library(lme4)
library(easystats)

# mamba install -y -c conda-forge r-datawizard
#  * datawizard  (0.7.1 -> 0.9.0)
#  * effectsize  (0.8.3 -> 0.8.6)
#  * insight     (0.19.2 -> 0.19.6)
#  * performance (0.10.4 -> 0.10.5)
#  * parameters  (0.21.1 -> 0.21.2)

	
# a specification for the model link function. This can be a name/expression, a literal character string, a length-one character vector, or an object of class "link-glm" (such as generated by make.link) provided it is not specified via one of the standard names given next. 
# The gaussian family accepts the links (as names) identity, log and inverse; the binomial family the links logit, probit, cauchit, (corresponding to logistic, normal and Cauchy CDFs respectively) log and cloglog (complementary log-log); the Gamma family the links inverse, identity and log; the poisson family the links log, identity, and sqrt; and the inverse.gaussian family the links 1/mu^2, inverse, identity and log.
# The quasi family accepts the links logit, probit, cloglog, identity, inverse, log, 1/mu^2 and sqrt, and the function power can be used to create a power link function.

model <- glmer(
  lognorm_cladecounts_by_hum ~ case_control_other_latest + age_at_baseline + (1 | study),
  family = binomial(link = "cauchit"),
  na.action = na.omit,
  data = coi_df %>% filter(coi == "Bacteria")
  )

report(model)
check_model(model)



# bphy_relab <- kuq %>%
#   filter(grepl("; Bacteria;", lineage)) %>%
#   filter(rank == "phylum") %>%
#   mutate(taxName = case_when(
#     taxName %in% top_10_bact_phy$taxName ~ taxName,
#     TRUE ~ "Other"
#   )) %>%
#   group_by(sample_id) %>%
#   dplyr::reframe(bact_percent_reads = (reads / sum(reads)) * 100, across())

# bphy_relab %>% glimpse

# bphy_relab_mat <- bphy_relab %>%
#   group_by(sample_id, taxName) %>%
#   dplyr::summarize(bact_percent_reads = sum(bact_percent_reads)) %>%
#   ungroup() %>%
#   pivot_wider(
#     id_cols = sample_id,
#     names_from = taxName,
#     values_from = bact_percent_reads,
#     values_fill = 0
#   ) %>%
#   column_to_rownames("sample_id") %>%
#   as.matrix()

# sid_order <- seriate_matrix_rows(bphy_relab_mat)
# phylum_order <- seriate_matrix_rows(t(bphy_relab_mat))

# saveRDS(
#   sid_order,
#   glue(
#     "{wkdir}/data/interim/kraken_results/",
#     "{Sys.Date()}_KrakenUniq-profiles-bacterial-phylum-relab-SampleID-ORDER-OLO.rds"
#   )
# )
# saveRDS(
#   phylum_order,
#   glue(
#     "{wkdir}/data/interim/kraken_results/",
#     "{Sys.Date()}_KrakenUniq-profiles-bacterial-phylum-relab-Phylum-ORDER-OLO.rds"
#   )
# )

saveRDS(
  sid_order,
  glue(
    "{wkdir}/data/interim/kraken_results/",
    "{Sys.Date()}_KrakenUniq-profiles-bacterial-phylum-relab-SampleID-ORDER-OLO.rds"
  )
)
saveRDS(
  phylum_order,
  glue(
    "{wkdir}/data/interim/kraken_results/",
    "{Sys.Date()}_KrakenUniq-profiles-bacterial-phylum-relab-Phylum-ORDER-OLO.rds"
  )
)


bphy_relab_mat %>% dim



meta_df



#________________________________________________________
# selecting pairs of samples

# Randomly select 5 PD and 5 control samples from each study
meta_df <- ps_analysis %>% microbiome::meta()
meta_df %>% glimpse



# impute data frame with 1/2 * minimum relab value
relab <- abundances(ps_analysis)
relab_values <- sort(relab) %>% unique
# minimum detected value 
impt <- relab_values[2]/2
relab <- relab + impt



kl_metric(
  relab[, "SY-PDZX943HWN"],
  relab[, "SY-PDZN008XD9"]
)


meta_df %>% 

ps_analysis

meta_df_sub <- meta_df %>%
  group_by(study, case_control_other_latest) %>%
  slice_sample(n=10)

basemat <- matrix(
  nrow = nrow(meta_df_sub),
  ncol = 1,
  rep(1, nrow(meta_df_sub))
) %>%
  set_rownames(meta_df_sub$participant_id)


dmat <- dist(basemat, method = "binary")
dmat

# group sample pairs of interest into a long-df
get_sample_pairs <- function(d, g, name) {
  dist_obj <- usedist::dist_groups(d, g)
  dist_obj %>%
    select(-Distance) %>%
    dplyr::rename_at(vars(Group1, Group2, Label), ~ glue("{.}_{name}"))
}

group_dist <- full_join(
  get_sample_pairs(
    dmat, meta_df_sub$case_control_other_latest, "case_control"
  ),
  get_sample_pairs(
    dmat, meta_df_sub$study, "study"
  )
)
group_dist %>% glimpse

# group_dist %>% glimpse
future::plan("multisession", workers = 8)
tic()
group_dist_df <- group_dist[1:5000, ] %>%
  mutate(kl = furrr::future_map2(
    .x = Item1,
    .y = Item2,
    ~ kl_metric(
      relab[, .x],
      relab[, .y]
    )
  )) %>%
  mutate(kl_divergence = purrr::map_dbl(kl, ~ sum(unlist(.))))
toc()
group_dist_df %>% glimpse

group_dist_df %>%
  ggplot(aes(x = Label_study, y = kl_divergence)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0) +
  geom_point(
    aes(color = Label_case_control),
    position = position_jitter(width = 0.2), alpha = 0.6) +
  labs(y = "KL Divergence", x = NULL) +
  theme_minimal()

# meta_df %>% filter(stud)

# group_dist_df <- tibble(
#   "Item1" = s,
#   "Item2" = s,
#   "Group1" = s,
#   "Group2" = s,
#   "Label" = s
# ) 




#________________________________________________________
# Testing out classifer method


ica_res <- readRDS(
  glue(
    "{analysis_res_dir}/REF_Tcell_PBS__TARGET_Tcell_WA1__ProstT5/",
    "ica_results_list.rds"
  )
)
ica_res %>% glimpse

gmm_obj <- readRDS(
  glue(
    "{analysis_res_dir}/REF_Tcell_PBS__TARGET_Tcell_WA1__ProstT5/",
    "GMMs.rds"
  )
)

ref_group <- "Tcell_PBS"
nrows_full <- nrow(ica_res[[ref_group]]$S)
rand_rows <- sample(1:nrows_full, nrows_full/20)
test_set_subset <- ica_res[[ref_group]]$S[rand_rows, ]
test_set_subset %>% dim

cat(glue(
    "{get_time()} Determining optimal",
    " number of GMM cluster..."
), "\n")



future::plan("multisession", workers = 24)
tic()
clust_search <- 1:100 %>%
  purrr::set_names() %>% 
  furrr::future_map(
  ~ mclust::Mclust(
    ica_res[[ref_group]]$S,
    G = .x,
    prior=priorControl(),
    modelNames = "VVV"
  )
)
toc()


bic_list <- clust_search %>%
  purrr::map(~ .x$BIC[[1]]) %>%
  bind_rows(.id = "cluster_n") %>%
  pivot_longer(
    cols = everything(),
    names_to = "cluster_size", 
    values_to = "BIC"
  )

bic_list %>%
  filter(BIC == max(BIC)) 
bic_list %>%
  mutate(cluster_size = as.numeric(cluster_size)) %>%
  ggplot(aes(x = cluster_size, y = BIC)) +
  geom_point() +
  geom_line() +
  theme_bw()

bic_list[1] 





gmms_open_fit <- mclust::Mclust(
    test_set_subset,
    G = NULL,
    modelNames = "VVV"
)
optimal_cluster_n <- gmms_open_fit$G




fit_gmms <- function(ica_list, ref_group, optimal_cluster_n = NULL) {
    if (is.null(optimal_cluster_n)) {
        nrows_full <- nrow(ica_list[[ref_group]]$S)
        rand_rows <- sample(1:nrows_full, nrows_full/10)
        test_set_subset <- ica_list[[ref_group]]$S[rand_rows, ]

        cat(glue(
            "{get_time()} Determining optimal",
            " number of GMM cluster..."
        ), "\n")

        gmms_open_fit <- mclust::Mclust(
            test_set_subset,
            G = NULL,
            modelNames = "VVV"
        )
        optimal_cluster_n <- gmms_open_fit$G
    }

    cat(glue(
        "{get_time()} Fitting GMMs with",
        " {optimal_cluster_n} clusters..."
    ), "\n")

    gmms <- names(ica_list) %>%
        purrr::set_names() %>%
        purrr::map(
            ~ Mclust(ica_list[[.]]$S,
                G = optimal_cluster_n, 
                prior=priorControl(),
                modelNames = "VVV"
            )
        )
    return(gmms)
}


ica_res[["Tcell_PBS"]] %>% glimpse
# independent components recovered
ica_res[["Tcell_PBS"]]$S %>% glimpse


ica_res[[2]]$S %>% glimpse



# library(mvtnorm)  # For dmvnorm, multivariate normal density function

ica_res_list <- list()
ica_res <- list.files(
  analysis_res_dir,
  recursive = TRUE,
  pattern = "ica_results_list.rds"
) %>%
  purrr::map(
    ~ readRDS(glue("{analysis_res_dir}/{.}")) %>%
      c(ica_res_list, .)
  )

ica_res_list %>% str
ica_res  %>% unlist %>% str
names(ica_res)



targets <- c("Tcell_WA1", "Tcell_B.1.351", "Tcell_B.1.617")
# Function to load and extract the "S" dataframe
load_S_dataframe <- function(target) {
  file_path <- glue(
    "{analysis_res_dir}/",
    "REF_Tcell_PBS__TARGET_{target}__ProstT5/ica_results_list.rds"
  )
  readRDS(file_path)[[target]]$S %>% 
    as.data.frame()
}

# ica_res %>% glimpse
ica_res <- purrr::map_df(targets, load_S_dataframe)
ica_res_cont <- readRDS(
  glue(
    "{analysis_res_dir}/REF_Tcell_PBS__TARGET_Tcell_WA1__ProstT5/",
    "ica_results_list.rds"
  ))[["Tcell_PBS"]]$S %>% as.data.frame()

ica_res_df <- ica_res %>% bind_rows(ica_res_cont) %>% glimpse()


# load in GMMs
load_gmms <- function(target) {
  file_path <- glue(
    "{analysis_res_dir}/",
    "REF_Tcell_PBS__TARGET_{target}__ProstT5/GMMs.rds"
  )
  readRDS(file_path)[[target]]
}

gmm_obj <- targets %>%
  purrr::set_names() %>%
  purrr::map(load_gmms)

gmm_obj[["Tcell_PBS"]] <- readRDS(
  glue(
    "{analysis_res_dir}/",
    "REF_Tcell_PBS__TARGET_Tcell_WA1__ProstT5/GMMs.rds"
  )
)[["Tcell_PBS"]]




sample_names_df <- data.frame(
  "feature" = rownames(ica_res_df)) %>%
  mutate(sid = strex::str_before_first(feature, "_")) %>%
    glimpse()

sample_names_df$sid %>% table

# for each ICA subset relavent samples (full set)
llr_metric <- list()
for (id in unique(sample_names_df$sid)) {
  message("Processing: ", id, "...")
  feat_list <- filter(sample_names_df, sid == id) %>% pull(feature)
  ica_set <- ica_res_df[feat_list, ]
  feat_group <- feat_list[1] %>% strex::str_after_last("__")

  llr_metric[[glue("{id}_{feat_group}")]] <- compute_gmm_llr(
    x_ref = ica_set,
    x_test = ica_set,
    gmm_ref = gmm_obj[["Tcell_PBS"]],
    gmm_test = gmm_obj[[feat_group]]
  )
}

llr_metrics_df <- data.frame(
  samples = names(llr_metric),
  llr = unlist(unname(llr_metric))
) %>%
  mutate(group_comp = str_after_first(samples, "_")) %>%
  glimpse()

llr_metrics_df %>%
  ggplot(aes(group_comp, y = llr)) +
  geom_point(position = position_jitter(width = 0.2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





compute_gmm_llr(
  x_ref = ica_set,
  x_test = ica_set,
  gmm_ref = gmm_obj[["Tcell_PBS"]],
  gmm_test = gmm_obj[["Tcell_WA1"]],
  randset_ref = sample(1:nrow(ica_res[["Tcell_PBS"]]$S), 10),
  randset_test = sample(1:nrow(ica_res[["Tcell_PBS"]]$S), 10)
)




# calibrate within sample set LLR robustness
# Define a function to wrap your existing code
compute_parallel_llr <- function(l, iter) {
  message("Processing: ", l, " samples. Iteration: ", iter)
  llr_value <- compute_gmm_llr(
    x_ref = ica_res[["Tcell_PBS"]]$S,
    x_test = ica_res[["Tcell_PBS"]]$S,
    gmm_ref = gmm_obj[["Tcell_PBS"]],
    gmm_test = gmm_obj[["Tcell_PBS"]],
    randset_ref = sample(1:nrow(ica_res[["Tcell_PBS"]]$S), l),
    randset_test = sample(1:nrow(ica_res[["Tcell_PBS"]]$S), l)
  )
  tibble(l = l, iter = iter, llr = llr_value)
}


# sample(1:nrow(ica_res[["Tcell_PBS"]]$S), 10)

future::plan("multisession", workers = 16)

# Create a dataframe of all combinations of l and iter
params_df <- expand.grid(
  l =
    c(5, 50, 100, 1000, 5000, 10000, 20000), iter = 1:1000
)

# Use future_map to run in parallel
res_df <- params_df %>%
  future_pmap_dfr(~ compute_parallel_llr(..1, ..2))

# Display or return the resulting dataframe
print(res_df)

p_llr_cont <- res_df %>%
  mutate(l = factor(l)) %>%
  ggplot(aes(x = l, y = llr)) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(alpha = 0.3) +
  labs(x = "Number of samples", y = "Log-likelihood ratio") +
  theme_bw()

ggsave(
  glue(
    "{wkdir}/figures/AIRR/benchmarks/",
    "LLR-control-subsamples.png"
  ),
  p_llr_cont,
  width = 7, height = 4
)

summmarize_iter <- c(1, 5, 10, 20, 50, 75, 100)

summary_boot_df <- tibble()
for (i2 in 1:100) {
  for (si in summmarize_iter) {
    set.seed(i2)
    message("Processing: ", si, " samples. Iteration: ", i2)
    subset_df <- res_df %>%
      group_by(l) %>%
      slice_sample(n = si) %>%
      dplyr::summarize(mean_llr = mean(llr)) %>%
      mutate(iter2 = i2, subset_size = si)
    summary_boot_df <- bind_rows(summary_boot_df, subset_df)
  }
}

p_llr_mean_summary <- summary_boot_df %>%
  mutate(subset_size = factor(subset_size)) %>%
  ggplot(aes(x = subset_size, y = mean_llr)) +
  geom_point(aes(color = l), position = position_jitter(width = 0.2)) +
  geom_boxplot(alpha = 0.3) +
  facet_wrap(~l) +
  # geom_boxplot(alpha = 0.3) +
  labs(x = "Number of samples", y = "Mean Log-likelihood ratio") +
  theme_bw()

ggsave(
  glue(
    "{wkdir}/figures/AIRR/benchmarks/",
    "LLR-control-_mean-subsamples.png"
  ),
  p_llr_mean_summary,
  width = 14, height = 12
)








