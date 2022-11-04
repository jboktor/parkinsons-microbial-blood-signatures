# Joe Boktor
# Caltech - Mazmanian Lab

#_______________________________________________________________________________
#####                   PERMANOVA Analysis Function                        ##### 
#_______________________________________________________________________________

phyloseq_permanova <- function(ps_object, metadata_list, nperm = 10){
  
  #' Function to analyze phyloseq object with euclidean distance PERMANOVA
  #' Input: a phyloseq object and list of metadata columns to test
  #' Output: data frame with analysis variables
  
  require(microbiome)
  require(vegan)
  require(dplyr)
  require(foreach)
  require(doParallel)
  permanova_df <- tibble()
  start_time <- Sys.time()
  
  # Transform count data with clr transformation
  clr_counts <- ps_object %>% 
    microbiome::transform("clr") %>% 
    microbiome::abundances() %>% t()
  
  # pull metadata of interest from ps object
  metadata_vars <- microbiome::meta(ps_object) %>% 
    dplyr::select(all_of(metadata_list))
  
  # remove metadata with less than 2 unique values
  columns2keep <- sapply(metadata_vars, function(x) length(unique(na.omit(x)))) > 2
  metadata_vars <- metadata_vars[, columns2keep]
  
  #setup parallel processing
  start_time <- Sys.time()
  cores = detectCores()
  cl <- makeCluster(cores[1] - 2) # to prevent computer overload
  registerDoParallel(cl)
  
  loop <- 
    foreach(i = 1:length(metadata_vars), 
            .combine = 'rbind', .verbose = T,
            .packages = c("vegan")) %dopar% {
              
              # filter NA values from metadata and abundance df
              a <- metadata_vars[,i]
              a.narm <- na.omit(a)
              if (any(is.na(a))) {
                clr_counts.narm <- clr_counts[-attr(a.narm, "na.action"), ]
              } else {
                clr_counts.narm <- clr_counts
              }
              
              # Calculate PERMANOVA
              meta_ano <- vegan::adonis(vegan::vegdist(
                clr_counts.narm, method = "euclidean") ~ a.narm, permutations = nperm)
              
              # update stats df
              data.frame(
                "metadata" = colnames(metadata_vars[i]),
                "Df" = meta_ano$aov.tab[1, ]$Df,
                "SumsOfSqs" = meta_ano$aov.tab[1, ]$SumsOfSqs,
                "MeanSqs" = meta_ano$aov.tab[1, ]$MeanSqs,
                "F.Model" = meta_ano$aov.tab[1, ]$F.Model,
                "R2" = meta_ano$aov.tab[1, ]$R2,
                "p_value" = meta_ano$aov.tab[1, ]$`Pr(>F)`,
                "distance" = "Aitchisons",
                "n" = length(a.narm),
                "permutations" = nperm
              )
            }
  stopCluster(cl)
  end_time <- Sys.time()
  cat("PERMANOVA calculated in : ",
      end_time - start_time, attr(end_time - start_time, "units"), "\n")
  
  return(loop)
}
#_______________________________________________________________________________
