#

source("notebooks/R_scripts/_load_packages.R")
source("notebooks/R_scripts/_misc_functions.R")
sample_info <- readRDS("data/interim/metadata/2022-11-04_static_metdata.rds")

#_______________________________________________________________________________
#              Merge Kraken Reports and create phyloseq objects
#_______________________________________________________________________________

# Notes: errors in taxa table creation function that need fixing
# Sequentially updating taxon list sometimes leads to incorrect/inconsistent
# taxon mapping for a species across varying samples ? Maybe related to taxon sub-classes (G1, P1, ..)?
# Also errors in progress bar refreshing properly


refDBlist <- c("RefSeqPlusPF", "UHGG", "WoL")

# for (refDB in refDBlist) {
#   # # TROUBLE
#   # refDB <- "WoL"
#   
#   reportdir <- paste0("input_files/OLD/", refDB, "_mapped 4/")
#   # reportdir <- paste0("input_files/", refDB, "_mapped/")
#   filepaths <- list.files(path = reportdir)
#   kraken_df <- tibble()
#   tax_table_df <- tibble()
#   
#   require(foreach)
#   require(doParallel)
#   start_time <- Sys.time()
#   
#   #setup parallel processing
#   start_time <- Sys.time()
#   cores = detectCores()
#   cl <- makeCluster(cores[1] - 2) # to prevent computer overload
#   registerDoParallel(cl)
#   
#   
#   readloop <-
#     foreach(
#       report = filepaths,
#       .combine = 'rbind',
#       .verbose = T,
#       .packages = c("readr", "dplyr", "stringr")
#     ) %dopar% {
#       sampleID <- str_split(report, "__")[[1]][1]
#       
#       df_report <- read_tsv(paste0(reportdir, report),
#                             col_names = F,
#                             show_col_types = FALSE)
#       
#       if (nrow(df_report) > 0) {
#         df_report %<>%
#           dplyr::rename(
#             relative_abundance = X1,
#             clade_counts = X2,
#             taxon_counts = X3,
#             taxon_rank_code = X4,
#             NCBI_taxon_ID = X5,
#             feature = X6
#           ) %>%
#           rename_at(vars(-c(
#             feature, NCBI_taxon_ID, taxon_rank_code
#           )),
#           function(x)
#             paste(x, sampleID, sep = "_"))
#         
#         if (length(kraken_df) == 0) {
#           kraken_df <- df_report
#         } else {
#           kraken_df <- kraken_df %>%
#             full_join(df_report,
#                       by = c("feature", "NCBI_taxon_ID", "taxon_rank_code"))
#         }
#         
#         kraken_df #%>% as.data.frame()
#       } else {
#         kraken_df #%>% as.data.frame()
#         
#       }
#       
#     }
# }

#______________________________________________________________________________

refDB <- "RefSeqPlusPF"
reportdir <- paste0("data/input/WGX-Kraken/", refDB, "_mapped_paired/")
filepaths <- list.files(path = reportdir)
kraken_df <- tibble()
tax_table_df <- tibble()

filepaths_trim <- filepaths[1:90]

require(foreach)
require(doParallel)
require(magrittr)
start_time <- Sys.time()

# setup parallel processing
## using local hardware
start_time <- Sys.time()
cores = detectCores()
cl <- makeCluster(cores[1] - 2) # to prevent computer overload
registerDoParallel(cl)


aggregate_kraken_files <- function(file_list) {
  aggregate_loop <-
    foreach(
      report = file_list,
      .combine = 'rbind',
      .verbose = T,
      .packages = c("readr", "dplyr", "stringr", "magrittr")
    ) %dopar% {
      sampleID <- str_split(report, "_")[[1]][1]
      df_report <- read_tsv(paste0(reportdir, report),
                            col_names = F,
                            show_col_types = FALSE)
      
      if (nrow(df_report) > 0) {
        df_report %<>%
          dplyr::rename(
            relative_abundance = X1,
            clade_counts = X2,
            taxon_counts = X3,
            taxon_rank_code = X4,
            NCBI_taxon_ID = X5,
            feature = X6
          ) %>%
          dplyr::mutate(participant_id = sampleID)
        kraken_df %<>% bind_rows(df_report)
      } else {
        next
      }
    }
  
  stopCluster(cl)
  end_time <- Sys.time()
  cat(
    "Kraken2 files aggregated in:",
    end_time - start_time,
    attr(end_time - start_time, "units"), "\n"
  )
  return(aggregate_loop)
}


load_kraken_file <- function(reportdir, report) {
  sampleID <- str_split(report, "_")[[1]][1]
  df_report <- read_tsv(paste0(reportdir, report),
                        col_names = F,
                        show_col_types = FALSE)
  
  if (nrow(df_report) > 0) {
    df_report %<>%
      dplyr::rename(
        relative_abundance = X1,
        clade_counts = X2,
        taxon_counts = X3,
        taxon_rank_code = X4,
        NCBI_taxon_ID = X5,
        feature = X6
      ) %>%
      dplyr::mutate(participant_id = sampleID)
    kraken_df %<>% bind_rows(df_report)
    return(kraken_df)
  }
}

  
# for loop implementation
  for (file in filepaths_trim[1:10]) {
    load_kraken_file(file, reportdir = reportdir) %>% head() %>% print()
  }

# purrr syntax
tst <- 
  purrr::map(filepaths_trim[1:10], ~ load_kraken_file(report = ., reportdir =  reportdir ) ) %>% bind_rows()
  

## using GCP cluster ----
library(future)
library(googleComputeEngineR)
library(tictoc)
library(doFuture)
library(furrr)

vm <- gce_vm_template("rstudio", predefined_type = "f1-micro",  name = "rstudiotest")
docker_run(vm, "rocker/rstudio", c("Rscript", "-e", "1+1"))
# gce_vm_delete(vm)



vms_nested <- gce_vm_cluster(
  vm_prefix = "kraken2", 
  cluster_size = 3,
  scheduling = list(preemptible = TRUE) ## Optional: Use cheaper, preemptible machines
  )

plan(list( 
  ## Topology 1: Use the cluster of remote VMs 
  tweak(cluster, workers = as.cluster(vms_nested)),  
  ## Topology 2: Use all CPUs on each VM
  tweak(multicore)  
))

# taken from: https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks/16275428#16275428
chunk_func <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))
# filepaths_trim_list <- chunk_func(filepaths_trim, length(vms_nested))

tictoc::tic()
kraken_df <-
  future_map(filepaths_trim[1:10],
             ~ load_kraken_file(report = ., reportdir = reportdir), ) %>% bind_rows()
tictoc::toc()

gce_vm_delete(vms_nested)



# library(googleComputeEngineR)
# library(future.apply)
# 
# ## Emulate a slow function that can be sped up in parallel
# slow_square <-
#   function(x = 1) {
#     x_sq <- x^2
#     Sys.sleep(5)
#     return(x_sq)
#   }
# 
# ## Set up our cluster: 3 VMs with 8 cores each
# vms_nested <-
#   gce_vm_cluster(
#     vm_prefix = "nested-cluster",
#     cluster_size = 3,
#     #docker_image = "rocker/r-parallel",  ## Default
#     predefined_type = "n1-highcpu-8",
#     scheduling = list(preemptible = TRUE) ## Optional: Use cheaper, preemptible machines
#   )
# 
# plan(list(
#   ## Topology 1: Use the cluster of remote VMs
#   tweak(cluster, workers = as.cluster(vms_nested)),
#   ## Topology 2: Use all CPUs on each VM
#   tweak(multiprocess)
# ))
# 
# chunk_func <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))
# 
# ## Input data (vector to be iterated over by our function)
# input_data <- 1:200
# 
# ## Run the function in nested parallel on our cluster and record timing
# tictoc::tic()
# ans <-
#   ## Parallelise over the three VMS
#   future_sapply(seq_along(vms_nested), function(v) {
# 
#     ## Split the input data into distinct chunks for each VM
#     input_chunk <- chunk_func(input_data, length(vms_nested))[[v]]
# 
#     ## Parallelise within each of the VMs
#     future_sapply(input_chunk, slow_square)
#   })
# tictoc::toc()
# 
# 
# ## Show that it worked
# as.vector(ans)
# gce_vm_delete(vms_nested)























  save(kraken_df, file = paste0("data/Kraken_reports/", refDB, "_reports_merged.RData"))
  save(tax_table_df, file = paste0("data/Kraken_reports/", refDB, "_tax_table_df_TEMP.RData"))
  
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
    filter(participant_id %in% colnames(species_counts)) %>% sample_data()
  rownames(physeq_meta) <- physeq_meta$participant_id
  
  species_counts_otu <- otu_table(species_counts, taxa_are_rows=T)
  dat.species <- phyloseq(species_counts_otu, physeq_meta)
  dat.species
  save(dat.species, file = paste0("data/Phyloseq_Objects/", refDB,
                                  "/Species_counts.RData"))
  
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
  
}


# Create single R Object with all Phyloseq's in a list
refDBlist <- c("RefSeqPlusPF", "UHGG", "WoL")
base::load("data/Phyloseq_Objects/RefSeqPlusPF/Phyloseq_RefSeqPlusPF.RData")
phyloseq_objs_RefSeqPlusPF <- phyloseq_objs
base::load("data/Phyloseq_Objects/UHGG/Phyloseq_UHGG.RData")
phyloseq_objs_UHGG <- phyloseq_objs
base::load("data/Phyloseq_Objects/WoL/Phyloseq_WoL.RData")
phyloseq_objs_WoL <- phyloseq_objs
refDB_phyloseq <- list(
  phyloseq_objs_RefSeqPlusPF,
  phyloseq_objs_UHGG,
  phyloseq_objs_WoL
)
names(refDB_phyloseq) <- c(
  "RefSeqPlusPF",
  "UHGG",
  "WoL"
)
save(refDB_phyloseq, file = "data/Phyloseq_Objects/Phyloseq_all.RData")

#_______________________________________________________________________________


