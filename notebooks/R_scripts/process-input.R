# Joe Boktor
# Caltech - Mazmanian Lab
# Dec 2021

source("notebooks/R_scripts/_load_packages.R")
source("notebooks/R_scripts/_misc_functions.R")
sample_info <- readRDS("data/interim/metadata/2022-11-04_static_metdata.rds")

#_______________________________________________________________________________
#              Merge Kraken Reports and create phyloseq objects
# _______________________________________________________________________________

# Notes: errors in taxa table creation function that need fixing
# Sequentially updating taxon list sometimes leads to incorrect/inconsistent
# taxon mapping for a species across varying samples ? Maybe related to taxon sub-classes (G1, P1, ..)?
# Also errors in progress bar refreshing properly

# refDBlist <- c("RefSeqPlusPF", "UHGG", "WoL")
refDBlist <- c("WoL")

for (refDB in refDBlist){

  # reportdir <- paste0("input/OLD/", refDB, "_mapped 4/")
  reportdir <- paste0("data/input/WGX-Kraken/", refDB, "_mapped_paired/")
  filepaths <- list.files(path =reportdir)
  kraken_df <- tibble()
  tax_table_df <- tibble()

  # define progress bar
  iterlength <- length(filepaths)
  pb <- progress_bar$new(
    format = "  compiling kraken2 files [:bar] :current/:total :percent in :elapsed eta::eta",
    total = iterlength, clear = FALSE, width= 90, force = T )

  for (report in filepaths){
    # progress bar update
    pb$tick() 
    sampleID <- str_split(report, "__")[[1]][1]

    df_report <- read_tsv(
      paste0(reportdir, report),
      col_names = F, show_col_types = FALSE)
    
    if (nrow(df_report) > 0){
      df_report %<>%
        dplyr::rename(relative_abundance = X1,
                      clade_counts = X2,
                      taxon_counts = X3,
                      taxon_rank_code = X4,
                      NCBI_taxon_ID = X5,
                      feature = X6) %>%
        rename_at(vars(-c(feature, NCBI_taxon_ID, taxon_rank_code)),
                  function(x) paste(x, sampleID, sep = "_"))
      
      if (length(kraken_df) == 0){
        kraken_df <- df_report
      } else {
        kraken_df <- kraken_df %>%
          full_join(df_report, by = c("feature", "NCBI_taxon_ID", "taxon_rank_code"))
      }
      
    } else {
      cat("\nERROR: ", report, "is empty!!!\n")
    }
    
  }
  saveRDS(kraken_df, file = paste0("data/interim/kraken-reports/", refDB, "_reports_merged.rds"))
  saveRDS(tax_table_df, file = paste0("data/interim/kraken-reports/", refDB, "_tax_table_df.rds"))

  #_______________________________________________________________________________
  #                         Construct Phyloseq Objects
  #_______________________________________________________________________________

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
    filter(participant_id %in% colnames(species_counts)) %>% 
    sample_data()
  
  cat("Checking metadata length matches: ", 
      length(unique(physeq_meta$participant_id)) == nrow(physeq_meta))
  
  sample_names(physeq_meta) <- sample_data(physeq_meta)$participant_id

  species_counts_otu <- otu_table(species_counts, taxa_are_rows=T)
  dat.species <- phyloseq(species_counts_otu, physeq_meta)
  dat.species
  saveRDS(dat.species, file = paste0("data/Phyloseq_Objects/", refDB,
                                  "/Species_counts.rds"))

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
  saveRDS(dat.genus, file = paste0("data/Phyloseq_Objects/", refDB,
                                "/Genus_counts.rds"))


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
  saveRDS(dat.phylum, file = paste0("data/Phyloseq_Objects/", refDB,
                                 "/Phylum_counts.rds"))

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
  saveRDS(phyloseq_objs, file = paste0("data/Phyloseq_Objects/", refDB,
                                    "/Phyloseq_", refDB, ".rds"))

}


# Create single R Object with all Phyloseq's in a list
phyloseq_objs_RefSeqPlusPF <- readRDS("data/Phyloseq_Objects/RefSeqPlusPF/Phyloseq_RefSeqPlusPF.rds")
phyloseq_objs_UHGG <-  readRDS("data/Phyloseq_Objects/UHGG/Phyloseq_UHGG.rds")
phyloseq_objs_WoL <- readRDS("data/Phyloseq_Objects/WoL/Phyloseq_WoL.rds")
refDB_phyloseq <- list(
  phyloseq_objs_RefSeqPlusPF,
  phyloseq_objs_UHGG,
  phyloseq_objs_WoL)
names(refDB_phyloseq) <- c("RefSeqPlusPF", "UHGG", "WoL")
saveRDS(refDB_phyloseq, file = "data/Phyloseq_Objects/Phyloseq_all.rds")

#_______________________________________________________________________________


