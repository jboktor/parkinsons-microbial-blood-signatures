# Joe Boktor
# Caltech - Mazmanian Lab
# Feb 2022

source("src/_load_packages.R")
source("src/_misc_functions.R")
sample_info <- readRDS("data/Metadata/static_metdata.rds")
devtools::source_url("https://raw.githubusercontent.com/jboktor/RIMA_pipeline/master/src/immune_repertoire/trust4_metric_functions.R")

#_______________________________________________________________________________
#                            Merge TRUST4 Reports 
# _______________________________________________________________________________

reportdir <- paste0("input_files/TRUST4_reports/")
filepaths <- list.files(path =reportdir)
airr_data <- c("CDR3nt", "CDR3aa", "V", "D", "J", "C")
airr_df <- tibble()
tcr <- tibble()
bcr.light <- tibble()
bcr.heavy <- tibble()


# define progress bar
iterlength <- length(filepaths)
pb <- progress_bar$new(
  format = "  compiling TRUST4 data [:bar] :current/:total :percent in :elapsed eta::eta",
  total = iterlength, clear = FALSE, width= 90, force = T )

for (report in filepaths){
  # progress bar update
  pb$tick() 
  sampleID <- gsub('_report.tsv', '', report)
  
  cdr3 <- read_tsv(
    paste0(reportdir, report),
    col_names = T, show_col_types = FALSE)
  
  if (nrow(cdr3) > 0){
    
    cdr3 %<>%
      dplyr::rename("count" = "#count") %>%
      filter(count > 0) %>% 
      mutate(sample = sampleID) #%>% 
      # mutate(
      #   is_complete = if_else(
      #     CDR3aa != "partial" 
      #     && CDR3aa != "out_of_frame" 
      #     && !grepl("^_", CDR3aa) 
      #     && !grepl("^\\?", CDR3aa), 
      #     "Y", "N"
      #   )
      # )
    
    cdr3$is_complete <- sapply(cdr3$CDR3aa, function(x) ifelse(x != "partial" && x != "out_of_frame" && !grepl("_",x) && !grepl("\\?", x),"Y","N"))
    
    
    # Partition TCRs and BCRs
    cdr3.bcr <- cdr3 %>% filter(grepl("^IG",V) | grepl("^IG",J) | grepl("^IG",C))
    cdr3.tcr <- cdr3 %>% filter(grepl("^TR",V) | grepl("^TR",J) | grepl("^TR",C))
    
    #add lib size and clinic traits
    cdr3.bcr %<>% mutate(lib.size = sum(count))
    cdr3.tcr %<>% mutate(lib.size = sum(count))
    
    # Partition BCR heavy chain and light chain dfs
    cdr3.bcr.heavy <- cdr3.bcr %>% filter(grepl("^IGH",V) | grepl("^IGH",J) | grepl("^IGH",C))
    cdr3.bcr.light <- cdr3.bcr %>% filter(grepl("^IG[K|L]",V) | grepl("^IG[K|L]",J) | grepl("^IG[K|L]",C))
    
    if (length(tcr) == 0){
      tcr <- cdr3.tcr
      bcr.heavy <- cdr3.bcr.heavy
      bcr.light <- cdr3.bcr.light
    } else {
      tcr %<>% bind_rows(cdr3.tcr)
      bcr.light %<>% bind_rows(cdr3.bcr.light)
      bcr.heavy %<>% bind_rows(cdr3.bcr.heavy)
    }
    
  } else {
    cat("\nERROR: ", report, "is empty!!!\n")
  }
  
}

tcr %>% head
tcr %>% glimpse
bcr.light %>% head
bcr.light %>% glimpse
bcr.heavy %>% head
bcr.heavy %>% glimpse

saveRDS(tcr, file = "data/AIRR/TRUST4_TCR.rds")
saveRDS(bcr.light, file ="data/AIRR/TRUST4_BCR_light.rds")
saveRDS(bcr.heavy, file ="data/AIRR/TRUST4_BCR_heavy.rds")

tcr$C %>% unique
bcr.light$C %>% unique
bcr.heavy$C %>% unique


