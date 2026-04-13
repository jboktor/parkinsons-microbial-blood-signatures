home_dir <- "/resnick/groups/MazmanianLab/jboktor"
pdmbs_dir <- paste0(home_dir, "/PDMBS")
ref_dir <- paste0(home_dir, "/Downloads/RefDBs")
wkdir <- paste0(pdmbs_dir, "/pdairr")
src_dir <- paste0(wkdir, "/notebooks")
source(paste0(src_dir, "/R_scripts/_load-core-pkgs.R"))
source(paste0(src_dir, "/R_scripts/_misc_functions.R"))
scratch_dir <- "/resnick/scratch/jbok/PDMBS"

library(ggsci)
library(aplot)

# loading metadata (this is a list with both RNASEQ and WGS metadata)
phymeta <- readRDS("data/interim/metadata/2023-07-14_phyloseq-metadata.rds")


# loading KrakenUnqiue Results
kraken_uniq_results_df <- readRDS(
  glue(
    "{wkdir}/data/interim/kraken_results/",
    "2023-10-15_WGS_KrakenUniq_results.rds"
  )
)
kraken_uniq_results_df %>% glimpse()

phy_abund_df <- kraken_uniq_results_df %>%
  dplyr::select(sample_id, taxReads, taxID) %>%
  pivot_wider(names_from = taxID, values_from = taxReads, values_fill = 0) %>%
  column_to_rownames("sample_id") %>%
  mutate_all(as.numeric)

# remove all zero columns
phy_abund_df <- phy_abund_df[, colSums(phy_abund_df) > 0]
phy_abund_df %>% dim

# formatting otu table
ps_counts_otu <- otu_table(t(phy_abund_df), taxa_are_rows = TRUE)
# formatting metadata
wgs_ps_meta <- phymeta$WGS %>%
  filter(participant_id %in% colnames(ps_counts_otu)) %>%
  as.data.frame()
rownames(wgs_ps_meta) <- wgs_ps_meta$participant_id
ps_meta <- phyloseq::sample_data(wgs_ps_meta)

ps <- phyloseq(ps_counts_otu, ps_meta)

saveRDS(
  ps,
  glue(
    "data/processed/phyloseq_objects/raw/",
    "{Sys.Date()}_WGS_KrakenUnique_phyloseq.rds"
  )
)
