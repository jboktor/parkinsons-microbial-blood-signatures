# Joe Boktor
# Caltech - Mazmanian Lab

pdmbs_dir <- "/central/groups/MazmanianLab/joeB/PDMBS"
wgs_wkdir <- paste0(pdmbs_dir, "/workflow/WGS")
wkdir <- paste0(pdmbs_dir, "/parkinsons-microbial-blood-signatures")
source(paste0(wkdir, "/notebooks/R_scripts/_load-core-pkgs.R"))
source(paste0(wkdir, "/notebooks/R_scripts/_misc_functions.R"))
source(paste0(wkdir, "/notebooks/R_scripts/_plot-functions.R"))
source(paste0(wkdir, "/notebooks/R_scripts/_analysis_funcs.R"))
library(phyloseq)
library(microbiome)

metadata_categories_df <- readRDS(
  glue(
    "{wkdir}/data/interim/metadata/",
    "2023-06-06_metadata_categories_dataframe.rds"
  )
)
metadata_categories <- readRDS(
    glue(
      "{wkdir}/data/interim/metadata/",
      "2023-06-06_metadata_categories.rds"
    )
  )

ps_trim <- readRDS(
    glue(
        "{wkdir}/data/processed/",
        "phyloseq_objects/decontaminated/",
        "2024-02-25_WGS_KrakenUnique_phyloseq.rds"
    )
)
# _______________________________________________________________________________
# PERMANOVA Analysis ----

# Example single environmental run
permanova_results_dir <- glue("{wkdir}/data/interim/permanova/distributed")
cohorts <- meta(ps_trim)$study %>% unique()

# split the phyloseq object by study into a list
ps_by_study <- cohorts %>%
  purrr::set_names() %>%
  purrr::map(function(x) {
    filtered_sample_data <- subset(sample_data(ps_trim), study == x)
    prune_samples(sample_names(filtered_sample_data), ps_trim)
  })


# essential_meta <- c(
#   metadata_categories$demographics$demographic_vars,
#   metadata_categories$dietary_behavioral$PDQ %>% keep(grepl("score", .)),
# metadata_categories$dietary_behavioral$caffeine_intake,
# metadata_categories$clincal_assessments  %>% 
#   unlist(recursive = TRUE, use.names = FALSE) %>% 
#   keep(grepl("score", .)),
# metadata_categories, 
# metadata_categories, 
# )


# ps_trim %>% meta() %>% glimpse

# tic()
# future::plan("multisession", workers = 24)
# permanova_analysis <- cohorts %>%
#   purrr::set_names() %>%
#   furrr::future_map2(
#     .y = unlist(metadata_categories$demographics, 
#       recursive = TRUE, use.names = FALSE
#       ),
#     seed=TRUE,
#     ~ phyloseq_permanova_slurm(
#       ps_object = ps_by_study[[.x]],
#       metadata_column = .y,
#       nperm = 1000,
#       threads = 4
#     )
#   ) %>%
#   bind_rows(.id = "study")
# toc()



ps_trim %>% meta() %>% glimpse

# tst <- phyloseq_permanova_slurm(
#       output_path = glue("{permanova_results_dir}/TEST3.tsv"),
#       ps_object = ps_by_study[["BioFIND"]],
#       metadata_column = "age_at_baseline",
#       nperm = 1000,
#       threads = 4,
#       metadata = "BioFIND"
#     )

ps_trim












dir.create(file.path("data/Analyses/community_composition/"), showWarnings = FALSE)
save(permanova_df, file = paste0("data/Analyses/community_composition/PERMANOVA_", 
                                         RefDB, "_", level, ".RData"))
openxlsx::write.xlsx(permanova_df,
                     file = paste0('data/Analyses/community_composition/PERMANOVA_', 
                                   RefDB, '_', level, '.xlsx'), overwrite = T)



library(listenv)
library(future)
library(batchtools)
library(future.batchtools)

# refDB <- "UHGG"
# level <- "Species"
seq_method <- "WGS"
decon_status <- "decontaminated"
refDBs <- c("UHGG", "RefSeqPlusPF")
taxa_levels <- c("Species", "Genus")

meta_df <- readRDS(glue(
  "{wkdir}/data/interim/metadata/",
  "2023-06-06_metadata_categories_dataframe.rds"
))

for (refDB in refDBs) {
  for (level in taxa_levels) {
    message(glue("Processing {refDB} {level}"))
    ps <- readRDS(
      glue(
        "{wkdir}/data/processed/phyloseq_objects/{decon_status}/",
        "2023-07-17_{seq_method}_{refDB}_{level}_phyloseq.rds"
      )) %>%
      prune_samples(sample_sums(.) > 1, .)
    env <- meta(ps)
    test_vars <- env %>%
      select(any_of(meta_df$metadata)) %>% colnames()

    future::plan(
      future.batchtools::batchtools_slurm,
      template = glue(
        "{wkdir}/batchtools_templates/",
        "batchtools.slurm.tmpl"
        ),
      resources = list(
        name = glue("{get_time()}_permanova"),
        memory = "5G",
        ncpus = 8,
        walltime = 3600
      )
    )
    n_jobs <- ceiling(length(test_vars) / 10)
    permanova_runs <- listenv()
    for (job in 1:n_jobs) {
      chunk <- chunk_func(test_vars, n_jobs)[[job]]
      permanova_runs[[job]] %<-% phyloseq_permanova(
        ps_object = ps,
        nperm = 9999,
        metadata_list = chunk
      )
    }
    permanova_res <- as.list(permanova_runs) %>%
      bind_rows()
    saveRDS(permanova_res, glue(
      "{wkdir}/data/processed/permanova/",
      "{Sys.Date()}_{seq_method}_{refDB}_{level}_permanova.rds"
    ))
  }
}




# chunk_func()
# permanova_analysis <- phyloseq_permanova(
#   ps_object = ps,
#   nperm = 99,
#   test_vars[1:3]
# )


# alpha_figures_dir <- glue("{wkdir}/figures/alpha-diversity")
# alpha_data_dir <- glue("{wkdir}/data/processed/alpha-diversity")
# dir.create(alpha_figures_dir, showWarnings = FALSE)
# dir.create(alpha_data_dir, showWarnings = FALSE, recursive = TRUE)

alpha_stats_df <- tibble()
ps_objs <- list()
for (refDB in refDBs) {
  for (level in taxa_levels) {
    ps <- readRDS(
      glue(
        "{wkdir}/data/processed/phyloseq_objects/{decon_status}/",
        "2023-07-17_{seq_method}_{refDB}_{level}_phyloseq.rds"
      )
    ) %>%
      prune_samples(sample_sums(.) > 1, .)
    
    ## Calculate Alpha Diversity Metrics and add cols to df
    env <- meta(ps)
    alpha_stats <- microbiome::alpha(abundances(ps), index = "observed") %>%
      rownames_to_column(var = "participant_id")
    stats_df <- env %>%
      left_join(alpha_stats) %>%
      dplyr::select(case_control_other_latest, colnames(alpha_stats)) %>%
      pivot_longer(!c(case_control_other_latest, participant_id),
        names_to = "alpha_metric"
      ) %>%
      mutate(DB = refDB, rank = level)
    alpha_stats_df %<>% bind_rows(stats_df)
  }
}