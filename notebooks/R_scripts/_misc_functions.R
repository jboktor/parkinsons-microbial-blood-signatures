# Miscellanous functions

#______________________________________________________________________________
#                     Utility Functions
#______________________________________________________________________________


kraken2taxtable <- function(kraken_input){
  
  # # TROUBLESHOOTING
  # df_report <- read_tsv(
  #   "input_files/UHGG_mapped/PP-3003__report_UHGG.tsv",
  #   col_names = F, show_col_types = FALSE) %>%
  #   dplyr::rename(relative_abundance = X1,
  #                 clade_counts = X2,
  #                 taxon_counts = X3,
  #                 taxon_rank_code = X4,
  #                 NCBI_taxon_ID = X5,
  #                 feature = X6) %>%
  #   rename_at(vars(-c(feature, NCBI_taxon_ID, taxon_rank_code)),
  #             function(x) paste(x, "PP-12593", sep = "_"))
  # df_report
  # kraken_input <- df_report

  # Filter for classical [D, P, C, O, F, G, S] # removing kingdom
  taxon_levels <- c("D", "P", "C", "O", "F", "G", "S")
  # Create an empty vector for taxon classification 
  taxon_hierarchy <- vector(length = 7)
  names(taxon_hierarchy) = taxon_levels
  taxon_df <- tibble()
  
  if (any(taxon_levels %nin% unique(kraken_input$taxon_rank_code)) ) {
    cat("  double checking taxon rank annotation ...\n")
    kraken_input <- kraken_input %>%
      dplyr::select(taxon_rank_code, feature) %>%
      mutate(taxon_rank_code =
               if_else(str_detect(feature, "d__"),"D",
                       # if_else(str_detect(feature, "k__"),"K", # causes issues 
                       if_else(str_detect(feature, "p__"),"P",
                               if_else(str_detect(feature, "c__"),"C",
                                       if_else(str_detect(feature, "o__"),"O",
                                               if_else(str_detect(feature, "f__"),"F",
                                                       if_else(str_detect(feature, "g__"),"G",
                                                               if_else(str_detect(feature, "s__"), "S", taxon_rank_code
                                                               ))))))))
  }
  
  kraken_input_trim <- kraken_input %>%
    dplyr::select(taxon_rank_code, feature) %>%
    filter(taxon_rank_code %in% c("D", "P", "C", "O", "F", "G", "S"))
  
  for (feature in 1:nrow(kraken_input_trim)) {
    # continuously update taxon hierarchy while moving down list
    taxon_level <- kraken_input_trim[[feature, "taxon_rank_code"]]
    taxon_hierarchy[[taxon_level]] <- 
      gsub(paste0(tolower(taxon_level), "__"), "", kraken_input_trim[[feature, "feature"]])
    if (taxon_level == "S"){
      taxon_df <- rbind(taxon_df, unlist(taxon_hierarchy))
    }
  }
  colnames(taxon_df) <- c("Domain", "Phylum", "Class", "Order", 
                          "Family", "Genera", "Species")
  return(taxon_df)
}

#______________________________________________________________________________

# function to assist in data wrangling 
slim_report_df <- function(report){
  
  report_vars <- c("feature", "NCBI_taxon_ID", "taxon_rank_code")
  
  report %>% 
    dplyr::select(all_of(report_vars), 
                  contains("clade_counts")) %>% 
    rename_at(vars(-c(feature)), 
              function(x) gsub("clade_counts_", "", x)) %>% 
    pivot_longer(-all_of(report_vars), names_to = "participant_id", values_to = "count")
}


`%nin%` <- Negate( `%in%` )
