# Joe Boktor
# Caltech - Mazmanian Lab
# Dec 2021

source("src/_load_packages.R")
source("src/_plot-functions.R")
source("src/_misc_functions.R")


#_______________________________________________________________________________
# Aggregate Samtools Flagstat reports for Unmapped reads
#_______________________________________________________________________________

flagstat_df <- as_tibble()
flagstatdir <- paste0("input_files/readstats/")
filepaths <- list.files(path =flagstatdir)

for (flagstat in  filepaths){
  str_split(flagstat, "_flagstat")[[1]][1] %>% print()
  
  stats <-
    read.delim(
      file = paste0("input_files/readstats/", flagstat),
      stringsAsFactors = F, #sep = "\t",
      header = F)
  reads <- stats[1,] %>% str_split(" ") %>% unlist()
  row2add <-
    cbind(
      "participant_id" = str_split(flagstat, "_flagstat")[[1]][1],
      "QC_passed_unmapped_reads" = reads[1] %>% as.numeric(),
      "QC_failed_unmapped_reads" = reads[3] %>% as.numeric()
      
    )
  flagstat_df %<>% rbind(row2add) 
}
flagstat_df %<>%
  mutate_at(vars(QC_passed_unmapped_reads, QC_failed_unmapped_reads), 
            as.numeric)

saveRDS(flagstat_df, "data/Sequencing_stats/unmapped_bam_flagstats.rds")

#_______________________________________________________________________________

library(kableExtra)

flagstat_df %>%
  summarise_at(
    vars(QC_passed_unmapped_reads, QC_failed_unmapped_reads),
    list(
      Min = min,
      Mean = mean,
      Max = max,
      Sd = sd
    )
  ) %>% 
  kbl() %>% 
  kable_material(c("striped", "hover", "condensed"))

#_______________________________________________________________________________
