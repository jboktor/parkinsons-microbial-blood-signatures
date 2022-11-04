# Joe Boktor
# Caltech - Mazmanian Lab

source("src/_load_packages.R")
source("src/_misc_functions.R")
library(immunedeconv)
library(biomaRt)

#_______________________________________________________________________________
#  Data Pre-processing ----
#_______________________________________________________________________________

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

rnaseq_df <- read_tsv("input_files/2021_v2-5release_0510/releases_2021_v2-5release_0510_rnaseq_salmon_quantification_aggregated.genes.tsv", 
         col_names = T, show_col_types = FALSE)
rnaseq_df %>% head()

rnaseq_df_wide <- 
  rnaseq_df %>%
  dplyr::select(sample_id, Name, TPM) %>%
  pivot_wider(names_from = sample_id, values_from = TPM) %>% 
  mutate(Name = str_replace(Name, pattern = ".[0-9]+$", replacement = ""))
rnaseq_df_wide

gene_symbol_key <- 
  getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    values = rnaseq_df_wide$Name,
    mart = mart
  )
cat("Gene Key dimensions before filtering for NAs: ", gene_symbol_key %>% dim(), "\n")

gene_symbol_key %<>% 
  mutate_all(na_if,"") %>% 
  drop_na(hgnc_symbol)

cat("Gene Key dimensions after filtering: ", gene_symbol_key %>% dim(), "\n")

rnaseq_df_hugo <- 
  rnaseq_df_wide %>% 
  # Join hugo gene symbols to ensembl IDs
  dplyr::rename(ensembl_gene_id = Name) %>% 
  left_join(gene_symbol_key, by = "ensembl_gene_id") %>% 
  dplyr::select(-ensembl_gene_id) %>% 
  drop_na(hgnc_symbol) %>%
  # Sum TPM values for multiple Ensembl IDs that map to the same gene symbols
  group_by(hgnc_symbol) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames(var = "hgnc_symbol") %>% 
  ungroup()
rnaseq_df_hugo

saveRDS(rnaseq_df_hugo, file = "data/RNASeq/rnaseq_tpm.rds")
cat("TPM dataframe dimensions after summarizing: ", rnaseq_df_hugo %>% dim(), "\n")

#_______________________________________________________________________________
#  Immune cell-type proportion estimation ----
#_______________________________________________________________________________

gene_matrix <- rnaseq_df_hugo %>% as.matrix()

# deconvolution_methods

# important note: 
# All methods except CIBERSORT (relative mode) can support BETWEEN sample comparisons
# All methods except (xCell, TIMER, and MCP-counter) can be used for INTRA sample comparisons

# Deconvolution methods - inter and intra sample comparison
res_quantiseq <- deconvolute(gene_matrix, "quantiseq", tumor = F)
res_epic <- deconvolute(gene_matrix, "epic", tumor = F)
# Marker - gene based approaches (Inter sample comparisons only)
res_xcell <- deconvolute(gene_matrix, "xcell", tumor = F)
res_mcp_counter <- deconvolute(gene_matrix, "mcp_counter", tumor = F)

saveRDS(res_quantiseq, file = "data/Immune_Cell_Proportions/quanTIseq.rds")
saveRDS(res_epic, file = "data/Immune_Cell_Proportions/EPIC.rds")
saveRDS(res_xcell, file = "data/Immune_Cell_Proportions/xCell.rds")
saveRDS(res_mcp_counter, file = "data/Immune_Cell_Proportions/MCP-counter.rds")


 