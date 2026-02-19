# Script that 

source("notebooks/R_scripts/_misc_functions.R")
sample_info <- readRDS("data/interim/metadata/2022-11-25_static_metdata.rds")

library(magrittr)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tictoc)
library(glue)
library(foreach)
library(parallel)
library(doParallel)
library(future)
library(purrr)
library(future)
library(listenv)
library(batchtools)
library(future.batchtools)

# _______________________________________________________________________________
# Functions for reading in kraken2 tsv files ---
load_kraken_file <- function(reportdir, report) {
  sampleID <- str_split(report, "_")[[1]][1]
  df_report <- read_tsv(paste0(reportdir, report),
    col_names = FALSE,
    show_col_types = FALSE
  )

  if (nrow(df_report) > 0) {
    df_report %<>%
      as.data.table() %>%
      dplyr::rename(
        relative_abundance = X1,
        clade_counts = X2,
        taxon_counts = X3,
        taxon_rank_code = X4,
        NCBI_taxon_ID = X5,
        feature = X6
      ) %>%
      dplyr::mutate(sample_id = sampleID)
  }
}

aggregate_kraken_files <- function(input_files) {
  require(purrr)
  require(data.table)
  require(stringr)
  require(dplyr)
  require(readr)
  purrr::map(
    input_files,
    ~ load_kraken_file(report = ., reportdir = reportdir)
  ) %>% bind_rows()
}

#_______________________________________________________________________________
#              Merge Kraken Reports and create phyloseq objects
#_______________________________________________________________________________

# Notes: errors in taxa table creation function that need fixing
# Sequentially updating taxon list sometimes leads to incorrect/inconsistent
# taxon mapping for a species across varying samples ? Maybe related to taxon sub-classes (G1, P1, ..)?
# Also errors in progress bar refreshing properly

refDBlist <- c("RefSeqPlusPF", "UHGG")
outputdir <- glue("{wkdir}/data/interim/kraken_results")
refDB <- "RefSeqPlusPF"
kraken_df <- readRDS(
  glue("{outputdir}/2022-12-31_{refDB}_kraken2_datatable.rds")
)

kraken_df_wide <- kraken_df %>%
  select(feature, NCBI_taxon_ID, taxon_rank_code, sample_id, taxon_counts) %>% 
  pivot_wider(
    names_from = sample_id,
    values_from = taxon_counts,
    values_fill = 0
  )

#______________________________________________________________________________
  
  # save(kraken_df, file = paste0("data/Kraken_reports/", refDB, "_reports_merged.RData"))
  # save(tax_table_df, file = paste0("data/Kraken_reports/", refDB, "_tax_table_df_TEMP.RData"))
  
  #____________________________________________________________________________
  #                         Construct Phyloseq Objects
  #____________________________________________________________________________
  
  # Species
  species_counts <- kraken_df %>%
    filter(taxon_rank_code == "S" |
             grepl("s__", feature)) %>%
    dplyr::select(feature, contains("clade_counts")) %>%
    rename_at(vars(-c(feature)),
              function(x) gsub("clade_counts_", "", x)) %>%
    column_to_rownames(var = "feature") %>%
      replace(is.na(.), 0)
    
  physeq_meta <- sample_info %>%
    filter(participant_id %in% colnames(species_counts)) %>% sample_data()
    rownames(physeq_meta) <- physeq_meta$participant_id

  species_counts_otu <- otu_table(species_counts, taxa_are_rows=T)
  dat.species <- phyloseq(species_counts_otu, physeq_meta)
  dat.species
  save(dat.species, file = paste0(
    "data/Phyloseq_Objects/", refDB,
    "/Species_counts.RData"
  ))

  # Genus
  genus_counts <- kraken_df %>%
    filter(taxon_rank_code == "G" |
             grepl("g__", feature)) %>%
    dplyr::select(feature, contains("clade_counts")) %>%
    rename_at(vars(-c(feature)),
              function(x) gsub("clade_counts_", "", x)) %>%
    column_to_rownames(var = "feature") %>%
    replace(is.na(.), 0)

  genus_counts_otu <- otu_table(genus_counts, taxa_are_rows=T)
  dat.genus <- phyloseq(genus_counts_otu, physeq_meta)
  dat.genus
  save(dat.genus, file = paste0("data/Phyloseq_Objects/", refDB,
                                "/Genus_counts.RData"))

  # Phylum
  phylum_counts <- kraken_df %>%
    filter(taxon_rank_code == "P" |
             grepl("p__", feature)) %>%
    dplyr::select(feature, contains("clade_counts")) %>%
    rename_at(vars(-c(feature)),
              function(x) gsub("clade_counts_", "", x)) %>%
    column_to_rownames(var = "feature") %>%
    replace(is.na(.), 0)
  
  phylum_counts_otu <- otu_table(phylum_counts, taxa_are_rows=T)
  dat.phylum <- phyloseq(phylum_counts_otu, physeq_meta)
  dat.phylum
  save(dat.phylum, file = paste0("data/Phyloseq_Objects/", refDB,
                                 "/Phylum_counts.RData"))
  
  phyloseq_objs <- list(
    dat.species,
    dat.genus,
    dat.phylum
  )
  names(phyloseq_objs) <- c(
    "Species",
    "Genus",
    "Phylum"
  )
  save(phyloseq_objs, file = paste0("data/Phyloseq_Objects/", refDB,
                                    "/Phyloseq_", refDB, ".RData"))
  
# }


# # Create single R Object with all Phyloseq's in a list
# refDBlist <- c("RefSeqPlusPF", "UHGG", "WoL")
# base::load("data/Phyloseq_Objects/RefSeqPlusPF/Phyloseq_RefSeqPlusPF.RData")
# phyloseq_objs_RefSeqPlusPF <- phyloseq_objs
# base::load("data/Phyloseq_Objects/UHGG/Phyloseq_UHGG.RData")
# phyloseq_objs_UHGG <- phyloseq_objs
# base::load("data/Phyloseq_Objects/WoL/Phyloseq_WoL.RData")
# phyloseq_objs_WoL <- phyloseq_objs
# refDB_phyloseq <- list(
#   phyloseq_objs_RefSeqPlusPF,
#   phyloseq_objs_UHGG,
#   phyloseq_objs_WoL
# )
# names(refDB_phyloseq) <- c(
#   "RefSeqPlusPF",
#   "UHGG",
#   "WoL"
# )
# save(refDB_phyloseq, file = "data/Phyloseq_Objects/Phyloseq_all.RData")

#_______________________________________________________________________________


