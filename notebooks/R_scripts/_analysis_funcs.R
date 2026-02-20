# Joe Boktor
# Caltech - Mazmanian Lab

#_______________________________________________________________________________
#####                   PERMANOVA Analysis Function                        ##### 
#_______________________________________________________________________________

  #' Function to analyze phyloseq object with euclidean distance PERMANOVA
  #' Input: a phyloseq object and list of metadata columns to test
  #' Output: data frame with analysis variables

phyloseq_permanova <- function(ps_object, metadata_list, nperm = 10, threads = NULL){
  require(microbiome)
  require(vegan)
  require(dplyr)
  require(parallel)
  require(foreach)
  require(doParallel)
  start_time <- Sys.time()

  # Transform count data with clr transformation
  clr_counts <- ps_object %>%
    microbiome::transform("clr") %>%
    microbiome::abundances() %>%
    t()

  # pull metadata of interest from ps object
  metadata_vars <- microbiome::meta(ps_object) %>%
    dplyr::select(all_of(metadata_list))

  # remove metadata with less than 2 unique values
  columns2keep <- sapply(
    metadata_vars,
    function(x) length(unique(na.omit(x)))
  ) > 2
  metadata_vars <- metadata_vars[, columns2keep]

  #setup parallel processing
  if (is.null(threads)) {
    threads <- parallel::detectCores()[1] - 2
  }
  cl <- makeCluster(threads)
  registerDoParallel(cl)

  start_time <- Sys.time()
  loop <-
    foreach(
      i = 1:length(metadata_vars),
      .combine = "rbind", .verbose = T,
      .packages = c("vegan")
    ) %dopar% {
      # filter NA values from metadata and abundance df
      a <- metadata_vars[, i]
      a.narm <- na.omit(a)
      if (any(is.na(a))) {
        clr_counts.narm <- clr_counts[-attr(a.narm, "na.action"), ]
      } else {
        clr_counts.narm <- clr_counts
      }

      # Calculate PERMANOVA
      meta_ano <- vegan::adonis(vegan::vegdist(
        clr_counts.narm,
        method = "euclidean"
      ) ~ a.narm, permutations = nperm)

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



phyloseq_permanova_slurm <- function(ps_object, metadata_column, metadata, 
                                    output_path, nperm = 10, threads = NULL) {

  require(dplyr)
  require(magrittr)
  require(phyloseq)
  require(microbiome)
  require(vegan)
  tictoc::tic()

  # Extract taxa counts df
  taxa_counts <- microbiome::abundances(ps_object) %>% t()

  # pull metadata column of interest from ps object
  metadata_col <- microbiome::meta(ps_object) %>%
    dplyr::select(all_of(metadata_column))
  meta_vec <- pull(metadata_col, metadata_column)

  # skip if less than two unique metadata variables
  if (length(unique(meta_vec)) < 2) {
    return(NULL)
  }

  # drop NA values from metadata and abundance df
  meta_vec_narm <- na.omit(meta_vec)
  if (any(is.na(meta_vec))) {
    taxa_counts <- taxa_counts[!is.na(meta_vec), ]
  }

  # Calculate PERMANOVA
  permanova_res <- vegan::adonis2(
    distances::distances(taxa_counts) ~ meta_vec_narm,
    permutations = nperm,
    parallel = threads
  )

  tidy_res <- suppressWarnings(broom::tidy(permanova_res))[1, ]

  # populate a data frame with PERMANOVA results
  perm_stats <- data.frame(
    "metadata" = metadata_column,
    "Df" = tidy_res["df"][[1]],
    "SumsOfSqs" = tidy_res["SumOfSqs"][[1]],
    "F_stat" = tidy_res["statistic"][[1]],
    "R2" = tidy_res["R2"][[1]],
    "p_value" = tidy_res["p.value"][[1]],
    "n" = length(meta_vec_narm),
    "permutations" = nperm,
    "additional_metadata" = metadata
  )

  tictoc::toc()

  utils::write.table(perm_stats,
    file = output_path,
    row.names = FALSE,
    sep = "\t"
  )
  return(perm_stats)
  }

