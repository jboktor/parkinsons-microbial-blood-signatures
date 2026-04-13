



# ______________________________________________________________________________

# Utility Functions
# Loading in data from invidual sample embedding results

# TCR full-length chain embeddings -------------
model_list <- c(
    "esm2_t12_35M",
    "esm2_t33_650M",
    "ProtT5-XL",
    "ProstT5",
    "Ankh",
    "Ankh-Large"
)

model_suffixes <- paste0(model_list, ".csv")
embedding_paths <- purrr::set_names(model_list) %>%
    purrr::map(
        ~ list.files(glue("{embedding_res}/TCR"),
            pattern = glue("_{.}_AVG.csv$"), full.names = TRUE
        )
    )
embedding_paths %>% purrr::map(length)

modeln <- "ProtT5-XL"
ref_group <- "Case"

df_latent <- embedding_paths[[modeln]] %>%
    keep(~ file.size(.) > 0) %>% 
    purrr::set_names() %>%
    purrr::map_dfr(~ data.table::fread(., header = TRUE), .id = "filepath") %>%
    mutate(sample_id = basename(filepath) %>%
        str_remove(glue("_{modeln}_AVG.csv"))) %>%
    left_join(rna_meta_essential, by = "sample_id") %>%
    mutate(Label = glue("{Label}_{sample_id}"))

df_latent %>% dim
df_latent %>% glimpse
df_latent %<>% select(where(~ sum(is.na(.)) < 1000 )) %>% glimpse

D <- df_latent %>%
    column_to_rownames("Label") %>%
    split_train_test(
        split_col = "case_control_other_at_baseline",
        remove_cols = c(
            "sample_id", "case_control_other_at_baseline", "filepath"
        )
    )

# print a preview of df
purrr::map(D, tibble)


# ______________________________________________________________________________

ica_list <- D %>% purrr::map(
    ~ as.data.frame(.) %>%
        fastICA(n.comp = 20)
)

# Visualizing ICA components
p_case <- plot_ica_heatmap(ica_list$Case) + ggtitle("Case")
p_cont <- plot_ica_heatmap(ica_list$Control) + ggtitle("Control")
ica_heatmap <- p_case + p_cont

# ica_list$Case$S %>% View

ica_df <- purrr::set_names(names(ica_list)) %>%
    purrr::map(~ ica_list[[.]]$S %>% as.data.frame()) %>%
    bind_rows(.id = "grouping") %>%
    dplyr::rename_all(  ~ gsub("V", "ICA_", .)) %>%
    rownames_to_column("Label") %>%
    dplyr::mutate(
        cid = strex::str_before_first(Label, "_"),
        sample_id = strex::str_after_first(Label, "_")
    ) %>% 
    glimpse()
ica_df %>% tibble


#______________________________________________________________________________
# Gaussian Mixture Models (GMMs)

# Determine optimal number of GMM cluster
gmm_cluster_fit <- Mclust(ica_list$Case$S, G = NULL)
optimal_cluster_n <- gmm_cluster_fit$G
message("Optimal number of clusters is: ", optimal_cluster_n)

# collecting data for BIC plot
bic_matrix <- gmm_cluster_fit$BIC
bic_data <- expand.grid(
  G = rownames(bic_matrix), 
  Model = colnames(bic_matrix)
)
bic_data$BIC <- as.vector(bic_matrix)
p_GMM_BIC <- bic_data %>%
    as.data.frame() %>%
    ggplot(aes(x = G, y = BIC, color = Model)) +
    geom_point(size = 2) +
    geom_line(aes(group = Model)) +
    labs(x = "Number of Clusters", y = "BIC") +
    theme_bw()
p_GMM_BIC

# Fit GMMs to reference and test sets
gmms <- ica_list %>% purrr::map(
    ~ Mclust(.x$S, G = optimal_cluster_n, modelNames = "VVV")
    )

# print a summary of the GMMs
purrr::map(gmms, summary)

# Computing Jeffrey's Divergence between each pair of mixture components of the two GMMs
jd_matrix <- align_gmm_distributions(gmms$Control, gmms$Case)
print(jd_matrix)

jd_pairs_long <- jd_matrix %>%
    as_tibble() %>%
    rownames_to_column("ref_components") %>%
    pivot_longer(
        !ref_components,
        names_to = "test_components",
        values_to = "jd"
    ) %>%
    mutate(
        ref_components = paste0("Mixture_", ref_components)) %>%
        mutate(test_components = gsub("V", "Mixture_", test_components)) %>%
    # Finding closest component from the TEST set to the REFERENCE set
    group_by(test_components) %>%
    dplyr::mutate(
        min_jd = min(jd),
        train_min = dplyr::case_when(
        jd == min_jd ~ "X",
        TRUE ~ ""
    ))

# List of test set GMM identities to reference set GMM identities
pair_map_df <- jd_pairs_long %>%
    filter(train_min == "X") %>%
    select(contains("components"))

gmm_class_map <- gmms %>%
    purrr::map_dfr(
        ~ .x$classification %>% as.data.frame(),
        .id = "grouping"
    ) %>%
    dplyr::rename("mixture" = ".") %>%
        mutate(
            mixture = paste0("Mixture_", mixture),
            # not actually the test group col,
            # temporarily using this to map to ref set
            test_components = mixture
        ) %>%
        rownames_to_column("Label") %>%
    dplyr::left_join(pair_map_df, by = "test_components") %>%
    dplyr::mutate(
        aligned_mixture = dplyr::case_when(
            grouping == ref_group ~ mixture,
            TRUE ~ ref_components
        )
    ) %>%
    dplyr::select(-c(test_components, ref_components)) %>%
    glimpse()



# function to select the closest matching component from the training set to the test
p_jd_heatmap <- jd_pairs_long %>%
    ggplot(aes(x = ref_components, y = test_components, fill = jd)) +
    geom_tile() +
    geom_text(aes(label = train_min), color = "white", size = 6) +
    theme_minimal() +
    theme_set(theme_bw()) +
    scale_fill_viridis_c(option = "mako") +
    labs(
        x = "Reference Components", y = "Test Components",
        fill = "Jeffrey's Divergence "
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top"
    )
p_jd_heatmap


fig_mapping_oveview <- (ica_heatmap | (p_GMM_BIC / p_jd_heatmap)) +
    plot_annotation(tag_levels = "A")

ggsave(
#     filename = glue("{wkdir}/figures/AIRR/gmm-alignment-oveview_{modeln}_{ref_group}_{Sys.Date()}.png"),
#     fig_mapping_oveview,
#     width = 11, height = 9
# )



## Dimensionality Reduction of Embeddings

# PCA
D_pca <- prcomp(bind_rows(D), rank. = 3)
pca_df <- D_pca[["x"]] %>%
    as.data.frame() %>%
    rownames_to_column("Label") %>% 
    glimpse()

# UMAP using default R settings
D_umap = umap(bind_rows(D), n_components = 3, method = "naive")

umap_df <- D_umap$layout %>%
    as.data.frame() %>%
    dplyr::rename_all( ~ gsub("V", "UMAP_", .)) %>%
    rownames_to_column("Label")

# Combining all receptor level data-fraes
gmm_meta_df <- gmm_class_map %>%
    left_join(ica_df) %>%
    left_join(trust4_reports_chains) %>%
    left_join(pca_df) %>%
    left_join(umap_df) %>%
    glimpse

gmm_meta_df %>% glimpse
dimred <- list(
    "pca_C" =  plotly_pca(gmm_meta_df, color_by = "C"),
    "pca_mixtures" = plotly_pca(gmm_meta_df, color_by = "mixture"),
    "pca_grouping" = plotly_pca(gmm_meta_df, color_by = "grouping"),
    "umap_C" = plotly_umap(gmm_meta_df, color_by = "C"),
    "umap_mixtures" = plotly_umap(gmm_meta_df, color_by = "mixture"),
    "umap_grouping" = plotly_umap(gmm_meta_df, color_by = "grouping")
)

# saveRDS(
#     dimred,
#     glue("{wkdir}/figures/AIRR/PCA-UMAP_{modeln}_{ref_group}_{Sys.Date()}.rds")
# )







# set.seed(27)
# x <- pp(matrix(runif(500),250,2))
# y <- pp(matrix(runif(500),250,2))
# wasserstein(x,x,p=1)
# wasserstein(y,y,p=1)
# wasserstein(x,y,p=1)
# wasserstein(x,y,p=3)


# requesting an interactive slurm session with a gpu for troubleshooting
# For an interactive job using 1 node and 1 tasks per node with a 4 hour run time on the short queue, this would look like:
# srun -N 1 -n 1 -t 02:00:00 --gpus-per-task=1 --mem-per-gpu=20G --pty bash

# gmm1$parameters$pro
# gmm1$parameters$mean
# gmm1$parameters$variance
# gmm1$z %>% glimpse
# gmm1$classification

# plot(gmm_cluster_fit, what = "BIC")
# plot(gmms$Control, what = "density")
# plot(gmms$Control, what = "uncertainty")
# plot(gmms$Control, what = "density", type = "persp")
# plot(gmms$Control, what = "boundaries", ngrid = 200)
# plot(gmms$Control, what = "pro")
# plot(gmms$Control, what = "mean")
# plot(gmms$Control, what = "density", type = "hdr", points.cex = 0.5)
# plot(gmms$Control, what = "density")

# gmm1 <- gmms$Control
# gmm2 <- gmms$Case
# Example of how to apply Jeffrey's Divergence to align distributions in two GMMs


# Assuming you have two GMMs fitted with mclust, named gmm1 and gmm2
# gmm1 <- Mclust(data1, G=3)
# gmm2 <- Mclust(data2, G=3)


# SCRAP ICA Analysis ----------------------
# reading in an ESM2 data-table

# # BCR light-chains
# bcr_light <- fread(
#     glue("{wkdir}/data/interim/airr/BCR_light_esm2_t30_150M.csv"),
#     header = TRUE
# )
# bcr_light %<>%
#     mutate(
#         particpant_id = str_after_nth(Label, "-", 2) %>%
#             str_before_last("_"),
#         case_control_other_latest = str_after_nth(Label, "-", 2) %>%
#             str_after_last("_")
#         )
# bcr_light %>% glimpse

# # random sample of 50 case and 50 controls
# set.seed(42)
# rand_samples <- bcr_light %>%
#     select(case_control_other_latest, particpant_id) %>%
#     distinct() %>%
#     filter(case_control_other_latest != "Other") %>%
#     group_by(case_control_other_latest) %>%
#     slice_sample(n = 25) %>%
#     pull(particpant_id)

# set.seed(42)
# bcr_light_subset <- bcr_light %>%
#     filter(particpant_id %in% rand_samples) %>%
#     group_by(particpant_id) %>%
#     slice_sample(n = 200) # 200 random receptors per sample
# bcr_light_subset %>% glimpse

# saveRDS(
#     bcr_light_subset,
#     glue("{wkdir}/data/interim/airr/bcr_light_subset.rds")
# )
# bcr_light_subset <- readRDS(
#     glue("{wkdir}/data/interim/airr/bcr_light_subset.rds")
# )
# bcr_light_subset %>% glimpse

# D <- bcr_light_subset %>%
#     column_to_rownames("Label") %>%
#     split_train_test(
#         split_col = "case_control_other_latest",
#         remove_cols = c("particpant_id", "case_control_other_latest")
#         )
# D[["Case"]] %>% tibble()
# D[["Control"]] %>% tibble()








qc_paths <- list.files(qc_stats_dir, full.names = TRUE)
readRDS(qc_paths[1])


contig_stats_df <- qc_paths %>%
    purrr::set_names(fs::path_ext_remove(basename(.))) %>%
    purrr::map(readRDS) %>%
    bind_rows(.id = "sid") %>%
    pivot_longer(!sid, names_to = "stat", values_to = "value") %>%
    mutate(receptor_type = case_when(
        grepl("tcr", stat) ~ "TCR",
        grepl("bcr", stat) ~ "BCR",
        TRUE ~ "ERROR"
    ), stat_rank = case_when(
        grepl("all", stat) ~ 1,
        grepl("qc1", stat) ~ 2,
        grepl("qc2", stat) ~ 3,
        grepl("qc3", stat) ~ 4
    )) %>% 
    glimpse

saveRDS(
    contig_stats_df,
    glue("{wkdir}/data/interim/airr/contig_qc_stats_df_{Sys.Date()}.rds"),
)

p_qc <- contig_stats_df %>%
    # filter(stat_rank > 1) %>%
    ggplot(aes(x = stat_rank, y = value)) +
    geom_point() +
    geom_line(aes(group = sid), alpha = 0.1) +
    geom_boxplot(aes(group = stat_rank), outlier.alpha = 0, alpha = 0.6) +
    facet_wrap(~receptor_type, scales = "free_y") +
    scale_y_log10() +
    theme_light() +
    theme(legend.position = "none")

ggsave(
    filename = glue("{wkdir}/figures/AIRR/contig_qc_stats_df_{Sys.Date()}.png"),
    p_qc,
    width = 11, height = 9

)






sid_pid_group_map <- rna_meta %>%
    select(sample_id, participant_id, case_control_other_at_baseline) %>%
    distinct() %>%
    glimpse()

df_latent %>% glimpse()
pids_50 <- df_latent %>%
    left_join(sid_pid_group_map) %>%
    pull(participant_id) %>%
    count() %>%
    filter(freq >= 50) %>%
    pull(x) %>%
    glimpse()

# computing all unique combinations of sample pairs
pid_pairs <- gtools::combinations(
  n = length(pids_50),
  r = 2, 
  v = pids_50,
  repeats.allowed = TRUE
  ) %>%
  as.data.frame() %>%
  dplyr::rename(
    sample_id1 = V1,
    sample_id2 = V2
  )

# 100 random subsamples with 50 rec per match
tic()
future::plan("multisession", workers = 12)
earth_movers_distance_df <- furrr::future_map2_dfr(
    sid_pairs$sample_id1, sid_pairs$sample_id2,
    ~ run_earth_movers_subsamples(pca_df_ot, s1 = .x, s2 = .y, n = 25)
)
toc()

earth_movers_distance_df <- readRDS(
    glue(
        "{wkdir}/data/interim/airr/",
        "earth_movers_distance_df_Ankh-Large_Case_2023-12-09.rds"
    )
)

all_pair_combos <- bind_rows(
    earth_movers_distance_df %>% select(sample_id1, sample_id2, min),
    earth_movers_distance_df %>%
        dplyr::select(sample_id1, sample_id2, min) %>%
        dplyr::rename(sample_id2 = sample_id1, sample_id1 = sample_id2)) %>%
    distinct() %>%
    mutate(min = as.numeric(min))
    
pair_matrix <- all_pair_combos %>%
    pivot_wider(names_from = sample_id2, values_from = min) %>%
    column_to_rownames("sample_id1")


sid1_order <- seriate_matrix_rows(as.matrix(pair_matrix))

p_earth_movers_heat <- all_pair_combos %>%
    # mutate_at(c("mean", "median", "sd", "min", "max"), as.numeric) %>%
    mutate(sample_id1 = factor(sample_id1, levels = sid1_order)) %>%
    mutate(sample_id2 = factor(sample_id2, levels = sid1_order)) %>%
    ggplot(aes(x = sample_id1, sample_id2, fill = min)) +
    geom_tile() +
    scale_fill_viridis_c(option = "mako") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_earth_movers_heat


ggsave(
    filename = glue("{wkdir}/figures/AIRR/earth-movers-heatmap_{Sys.Date()}.png"),
    p_earth_movers_heat,
    width = 10, height = 9
)




gmm_aligner_wrapper <- function(fpaths, modeln, out_name, wkdir) {
    require(glue)
    source(paste0(wkdir, "/notebooks/R_scripts/_misc_functions.R"))
    source(paste0(wkdir, "/notebooks/R_scripts/_airr-modeling_funcs.R"))
    future::plan("multisession", workers = 14)
    message("Running model: ", modeln)

    df_latent <- fpaths %>%
        purrr::keep(~ file.size(.) > 0) %>%
        purrr::set_names() %>%
        furrr::future_map(~ data.table::fread(., header = TRUE)) %>%
        dplyr::bind_rows(.id = "filepath") %>%
        dplyr::mutate(sample_id = basename(filepath) %>%
            stringr::str_remove(glue("_{modeln}_AVG.csv"))) %>%
        dplyr::left_join(rna_meta_essential, by = "sample_id") %>%
        dplyr::mutate(Label = glue("{Label}_{sample_id}"))

    message("Embeddings Dimensions: ", nrow(df_latent), " x ", ncol(df_latent))
    df_latent %<>%
        dplyr::select(where(~ sum(is.na(.)) < 1000))

    D <- df_latent %>%
        tibble::column_to_rownames("Label") %>%
        split_train_test(
            split_col = "case_control_other_at_baseline",
            remove_cols = c(
                "sample_id", "case_control_other_at_baseline", "filepath"
            )
        )
    rm(df_latent)
    purrr::map(D, tibble)

    trust4_reports_chains <- readRDS(
        glue("{wkdir}/data/interim/airr/",
        "2023-12-17_TRUST4-reports-chains.rds")
    )

    # print a preview of df
    gmm_res <- gmm_aligner(
        dlist = D,
        ref_group = "Case",
        test_group = "Control",
        t4_chains = trust4_reports_chains,
        threads = 14
    )

    saveRDS(
        gmm_res, glue("{wkdir}/data/interim/airr/{out_name}")
    )
}

possibly_align_gmms = possibly(
    .f = gmm_aligner_wrapper,
    otherwise = NULL, 
    quiet = FALSE
)

# # local run test
# possibly_align_gmms(
#     fpaths = embedding_paths[["esm2_t12_35M"]],
#     modeln = "esm2_t12_35M",
#     out_name = glue(
#         "GMMAligner_{rec_type}_{modeln}_Ref-{ref_group}",
#         "_Test-{test_group}_{Sys.Date()}.rds"
#     )
# )


# Initiate future.batchtools backend for parallel processing
future::plan(
  future.batchtools::batchtools_slurm,
  template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
  resources = list(
    name = "GMM-AligneR",
    memory = 200000, # 250GB
    ncpus = 12,
    walltime = 10800
  )
)

tic()
gmmAlignR_runs <- listenv()
for (rec_type in c("BCR", "TCR")) {
    ref_group <- "Case"
    test_group <- "Control"
    embedding_paths <- purrr::set_names(model_list) %>%
        purrr::map(
            ~ list.files(glue("{embedding_res}/{rec_type}/{.}"),
                pattern = glue("_{.}_AVG.csv$"), full.names = TRUE
            )
        )
    for (model in c("ProtT5-XL", "Ankh", "Ankh-Large")) {
        job <- glue("{rec_type}_{model}")
        file_out <- glue(
            "GMMAligner_{rec_type}_{model}",
            "_Ref-{ref_group}_Test-{test_group}_{Sys.Date()}.rds"
        )
        gmmAlignR_runs[[job]] %<-% possibly_align_gmms(
            fpaths = embedding_paths[[model]],
            modeln = model,
            wkdir = wkdir,
            out_name = file_out
        )
    }
}
toc()



gmmAlignR_runs_final <- as.list(gmmAlignR_runs)
names(gmmAlignR_runs_final) %>%
    purrr::map( ~ gmmAlignR_runs[[.]] %>% print)







# Earth Movers Distance Analysis ----------------------

meta_sid_pid <- rna_meta %>%
    dplyr::select(sample_id, participant_id)
meta_sid_pid %>% glimpse
modeln <- "esm2_t12_35M"


tst <-
    readRDS(
        glue(
            "{wkdir}/data/interim/airr/",
            "GMMAligner_esm2_t12_35M_Ref-Case_Test-Control_2023-12-19.rds"
        )
    )

tst$compiled %>% glimpse
pca_df_ot <- tst$compiled %>%
    dplyr::select(contains("PC"), cid, sample_id, Label) %>%
    left_join(meta_sid_pid) %>%
    glimpse

# Count the number of receptors per sample and make alist of all with > 20

quality_ids <- pca_df_ot %>%
    pull(participant_id) %>%
    count() %>%
    filter(freq >= 50) %>% # pull(freq) %>% hist(n = 100)
    pull(x) %>%
    sample(100)

# computing all unique combinations of sample pairs
quality_id_pairs <- gtools::combinations(
  n = length(quality_ids),
  r = 2, 
  v = quality_ids,
  repeats.allowed = TRUE
  ) %>% 
  as.data.frame() %>%
  dplyr::rename(
    id1 = V1,
    id2 = V2
  )


# 100 random subsamples with 50 rec per match
tic()
future::plan("multisession", workers = 14)
earth_movers_distance_df <- furrr::future_map2_dfr(
    quality_id_pairs$id1, quality_id_pairs$id2,
    ~ run_earth_movers_subsamples(pca_df_ot,
        s1 = .x, s2 = .y, n = 50,
        split_var = "participant_id",
        cols_to_filter = c(
            "Label",
            "participant_id",
            "sample_id",
            "cid"
        )
    )
)
toc()

earth_movers_distance_df %<>%
    mutate_at(c("mean", "median", "sd", "min", "max"), as.numeric)
earth_movers_distance_df %>% glimpse


# Visualizing results

meta_distance_df <- rna_meta %>%
    dplyr::select(sample_id, participant_id, case_control_other_at_baseline)
meta_distance_df %>% glimpse

pair_level <- "participant_id"
earth_movers_distance_df_annot <- earth_movers_distance_df %>%
  left_join(
    meta_distance_df %>% dplyr::rename_all(~ paste0("id1_", .)),
    by = c("id1" = glue("id1_{pair_level}")), 
    relationship = "many-to-many"
  ) %>%
  left_join(
    meta_distance_df %>% dplyr::rename_all(~ paste0("id2_", .)),
    by = c("id2" = glue("id2_{pair_level}")), 
    relationship = "many-to-many"
  ) %>%
  mutate(self = 
    case_when(
        # sample_id1 == sample_id2 ~ "Within Sample", 
        id1 == id2 ~ "Within Donor",
        TRUE ~ "Across Donors"
        )
    ) %>%
    mutate(
        group_comparison = case_when(
            id1_case_control_other_at_baseline == "Case" &
                id2_case_control_other_at_baseline == "Case" ~ "Within PD",
            id1_case_control_other_at_baseline == "Control" &
                id2_case_control_other_at_baseline == "Control" ~
                "Within Control",
            id1_case_control_other_at_baseline != id2_case_control_other_at_baseline ~
                "Across PD-Control",
            TRUE ~ "other"
        )
    )

p_emd_self <- earth_movers_distance_df_annot %>%
    ggplot(aes(y = median, x = self)) +
    geom_point(position = position_jitter(width = 0.4), alpha = 0.6) +
    geom_boxplot(alpha = 0.5, outlier.alpha = 0) +
    labs(x = NULL, y = "Median Earth Mover's Distance") +
    theme_set(theme_light()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_emd_group <- earth_movers_distance_df_annot %>%
    ggplot(aes(y = median, x = group_comparison)) +
    geom_point(position = position_jitter(width = 0.4), alpha = 0.6) +
    geom_boxplot(alpha = 0.5, outlier.alpha = 0) +
    labs(x = NULL, y = "Median Earth Mover's Distance") +
    theme_set(theme_light()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_emd <- p_emd_self + p_emd_group + p_emd_mediansd +
    plot_layout(widths = c(1, 1, 1.8)) +
    plot_annotation(tag_levels = "A")
p_emd


# saveRDS(
#     earth_movers_distance_df,
#     glue(
#         "{wkdir}/data/interim/airr/",
#         "earth_movers_distance_df_{modeln}_{ref_group}_{Sys.Date()}.rds"
#     )
# )






subset_test <- trust4_reports_clean_tcrs %>%
    sample_n(40) %>% 
    glimpse()

#-----------------------------------------------------------------------------

# stitchr analysis
# stitchr_wrapper <- function(V, J, CDR3nt, output_fmt, wkdir) {
#     require(glue)
#     source(glue("{wkdir}/notebooks/R_scripts/_misc_functions.R"))
#     shell_do(
#         glue("stitchr -v {V} -j {J} -cdr3 {CDR3nt} -m {output_fmt}"),
#         stdout_path = TRUE
#     )
# }



# Collecting import info on stop codons in the CDR3 region

thimble_stitchr_df_qc <- readRDS(
    glue("{wkdir}/data/interim/airr/thimble_stitchr_qc_seqs_2023-12-26.rds")
)

thimble_stitchr_df_qc %>% glimpse


# Function to find positions of '*' in a string
find_star_positions <- function(str) {
  # Split the string into individual characters
  chars <- unlist(strsplit(str, ""))
  # Find positions of '*'
  positions <- which(chars == "*")
  # Return NA if no '*' is found
  if (length(positions) == 0) {
    return(NA)
  } else {
    return(positions)
  }
}

# Apply the function to each element of the vector
tic()
star_positions <- thimble_stitchr_df_qc %>%
    filter(grepl("\\*", aa)) %>%
    pull(aa) %>%
    purrr::set_names() %>%
    purrr::map(., find_star_positions)
toc()

stop_codon_df <- as_tibble_row(
    list(
        "aa" = list(names(star_positions)),
        "stop_codon_indices" = list(unname(star_positions))
    )) %>%
    unnest(cols = c(aa, stop_codon_indices)) %>%
    mutate(earliest_stop_codon = purrr::map_dbl(stop_codon_indices, min)) %>%
    glimpse()


# stop_codon_df$earliest_stop_codon %>% as.numeric() %>% hist(n = 10)
# stop_codon_df$earliest_stop_codon %>% as.numeric() %>% unique

stitchr_stop_codon_indices <- thimble_stitchr_df_qc %>%
    left_join(stop_codon_df, by = "aa", relationship = "many-to-many") %>%
    dplyr::mutate(non_func = earliest_stop_codon < 100) %>%
        dplyr::select(
            filepath, unique_seq_id, warnings_errors,
            stop_codon_indices, earliest_stop_codon, non_func
        )

stitchr_stop_codon_indices %>% glimpse

saveRDS(
    stitchr_stop_codon_indices,
    glue(
        "{wkdir}/data/interim/airr/",
        "thimble_stitchr_qc_metrics_{Sys.Date()}.rds"
    )
)















err1 <- fread(
    glue(
        "{wkdir}/data/interim/airr/",
        "stiTChR_embeddings/stitchr_tcrs_chunk_1000_ProtT5-XL_AVG.csv"
        ), 
    header = TRUE
)

err1 %>% colnames
err1 %>% dim

err1 %>% glimpse


# Process for handeling data

# Read in all TRILL embeddings for a dataset and perform 1) PCA and 2) UMAP



# QC check on embeding results


# read in paired fasta file

# Match number of samples and number of dimensions


metadata <- readRDS(
    glue("{wkdir}/data/interim/metadata/2023-07-14_phyloseq-metadata.rds")
)
rna_meta <- metadata$RNASEQ
rna_meta_essential <- rna_meta %>%
    dplyr::select(sample_id, case_control_other_at_baseline)

trust4_reports_chains <- readRDS(
    glue("{wkdir}/data/interim/airr/",
    "2023-12-17_TRUST4-reports-chains.rds")
)
result_dims_df <- readRDS(
    glue("{wkdir}/data/interim/airr/embedding_dims_stiTChR_v2_2024-01-22.rds")
)
qc_embedding_paths <- result_dims_df %>%
    filter(correct_sample_n & correct_model_dim) %>%
    glimpse()

rec_type <- "TCR"
# modeln <- "Ankh-Large"
ref_group <- "Case"
test_group <- "Control"

for (modeln in model_list) {
    # Loading results
    df_latent <- qc_embedding_paths %>%
        filter(outer_list_name == modeln) %>%
        # sample_n(10) %>%
        pull(inner_list_name) %>%
        purrr::map(~ fread(
            glue("{wkdir}/data/interim/airr/stiTChR_embeddings_v2/{.}"),
            header = TRUE
        )) %>%
        bind_rows() %>% 
        left_join(rna_meta_essential, by = "sample_id") %>%
        glimpse()

    D <- df_latent %>%
        column_to_rownames("Label") %>%
        split_train_test(
            split_col = "case_control_other_at_baseline",
            remove_cols = c(
                "sample_id", "case_control_other_at_baseline"
            )
        )

    purrr::map(D, tibble)

    gmm_res <- gmm_aligner(
        dlist = D,
        ref_group = "Case",
        test_group = "Control",
        threads = 14
    )
}






metadata <- readRDS(
    glue("{wkdir}/data/interim/metadata/2023-07-14_phyloseq-metadata.rds")
)
rna_meta <- metadata$RNASEQ
rna_meta_essential <- rna_meta %>%
    dplyr::select(sample_id, case_control_other_at_baseline)

trust4_reports_chains <- readRDS(
    glue("{wkdir}/data/interim/airr/",
    "2023-12-17_TRUST4-reports-chains.rds")
)

# aggregate files and perform PCA and ICA

embd_res_dir <- glue("{wkdir}/data/interim/airr/stiTChR_embeddings")
embedding_paths <- model_list %>%
    purrr::set_names() %>%
    purrr::map(
        ~ list.files(embd_res_dir,
            pattern = glue("{.}_AVG.csv"), full.names = TRUE
        )
    )

embedding_paths %>% purrr::map(length)


ref_group = "Case"
test_group = "Control"
nthreads <- 28
set.seed(42)
rand_samps <- rna_meta_essential %>%
    filter(case_control_other_at_baseline != "Other") %>%
    group_by(case_control_other_at_baseline) %>%
    sample_n(500)

# 12gb  limit (1500*1024^2 = 1572864000) * 8
future::plan("sequential")
future::plan("multisession", workers = nthreads)
options(future.globals.maxSize = 12582912000)
for (modeln in model_list[2:6]){
    message("Running model: ", modeln)
    df_latent <- embedding_paths[[modeln]] %>%
        # sample(500) %>%
        furrr::future_map(~ data.table::fread(., header = TRUE)) %>%
        bind_rows() %>%
        mutate(sample_id = strex::str_before_nth(Label, "_", -2)) %>%
        # FILTER STEP
        filter(sample_id %in% rand_samps$sample_id) %>%
        left_join(rna_meta_essential, by = "sample_id") %>%
        column_to_rownames("Label")

    D <- df_latent %>%
        split_train_test(
            split_col = "case_control_other_at_baseline",
            remove_cols = c(
                "sample_id", "case_control_other_at_baseline"
            )
        )

    # print a preview of df
    purrr::map(D, tibble) %>% print()

    # Running GMM Aligner
    gmm_res <- gmm_aligner(
        dlist = D,
        ref_group = ref_group,
        test_group = test_group,
        threads = nthreads
    )

    saveRDS(
        gmm_res,
        glue("{wkdir}/data/interim/airr/GMMAligner/",
            "{modeln}_Ref-{ref_group}_Test-{test_group}",
            "_1ksubsample_{Sys.Date()}.rds"
        )
    )
    saveRDS(
        df_latent,
        glue("{wkdir}/data/interim/airr/stiTChR_embeddings_aggregated/",
            "{modeln}_1ksubsample_{Sys.Date()}.rds"
        )
    )
}



gmm_res <- readRDS(
        glue("{wkdir}/data/interim/airr/GMMAligner/",
            "ProstT5_Ref-Case_Test-Control_1ksubsample_2024-02-04.rds"
        )
    )

# selecting only receptors with ICA components within 4 sigma for
ica_vis <- gmm_res$compiled %>%
    dplyr::select(
        Label, grouping, mixture,
        aligned_mixture, sample_id, cid, C, contains("ICA_")
        #matches("^ICA_([1-8])$")
    ) %>% 
  filter(if_any(contains("ICA_"), ~ abs(. - mean(.)) >= 4 * sd(.))) %>%
 sample_n(50000) %>% 
 glimpse()

# ica_vis$Label %>% unique() %>% length
# ica_vis$Label %>% length

ica_vis_mat <- ica_vis %>%
    distinct() %>%
    column_to_rownames("Label") %>%
    select(contains("ICA_")) %>%
    as.matrix() %>%
    glimpse()

rec_order <- seriate_matrix_rows(ica_vis_mat)
component_order <- seriate_matrix_rows(t(ica_vis_mat))

ica_vis_df <- ica_vis %>%
    pivot_longer(
        cols = contains("ICA_"),
        names_to = "ICA_component",
        values_to = "ICA_value"
    ) %>%
    mutate(ICA_component = factor(ICA_component, levels = component_order)) %>%
    mutate(Label = factor(Label, levels = rec_order)) %>%
    glimpse

p_ica_vis <- ica_vis_df %>%
    ggplot(aes(x = ICA_component, y = Label)) +
    geom_tile(aes(color = ICA_value, fill = ICA_value)) +
    scale_fill_viridis_c(option = "cividis", limits = c(-10, 10)) +
    scale_color_viridis_c(option = "cividis", limits = c(-10, 10)) +
    theme_minimal() +
    labs(
        x = "ICA component", y = "Receptor",
        color = "weight", fill = "weight"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank()
    )

ggsave(
    filename = glue("{wkdir}/figures/AIRR/ICA_heatmap_{Sys.Date()}.png"),
    p_ica_vis,
    width = 3.5, height = 5
)


# Visualizing PCA Space

# gmm_res$compiled %>%
#     sample_n(1000) %>%
#     ggplot(aes(x = PC1, y = PC2, color = C)) +
#     geom_point() +
#     theme_light()

# figures/AIRR/gmm-param-summary_Ankh-Large_Case_2023-12-09.png

# gmm_res$compiled %>%
#     sample_n(50000) %>%
#     plotly_pca(color_by = "C")

gmm_res$compiled %>% glimpse



dimred <- list(
    "pca_C" = plotly_pca(sample_n(gmm_res$compiled, 50000),
        color_by = "C"
    ),
    "pca_mixtures" = plotly_pca(sample_n(gmm_res$compiled, 50000),
        color_by = "aligned_mixture"
    ),
    "pca_grouping" = plotly_pca(sample_n(gmm_res$compiled, 50000),
        color_by = "grouping"
    )
)
saveRDS(dimred,
    glue("{wkdir}/figures/AIRR/PCA_ProsT5-viztest.rds")
)



library(gg3D)
## An empty plot with 3 axes
qplot(x=0, y=0, z=0, geom="blank") + 
  theme_void() +
  axes_3D()

ggplot(
    sample_n(gmm_res$compiled, 50000),
    aes(x = PC1, y = PC2, z = PC3, color = aligned_mixture)) +
    scale_color_npg() +
    theme_void() +    
        axes_3D() +
        stat_3D()



#_______________________________________________________________________________
# Plotting GMM parameters

# Collecting stats on aligned distributions
pair_map_clean <- pair_map_df %>%
    mutate_at(
        vars(ref_components, test_components),
        ~ gsub("Mixture_", "", .) %>% as.numeric()
    ) %>%
    mutate(
        delta_mu = purrr::map2_dbl(.x = ref_components, .y = test_components,
                ~ gmm_mean_l2_norm(
                    gmm1 = gmms$Case, gmm2 = gmms$Control,
                    i = .x, j = .y
                )
            ), 
        delta_w = purrr::map2_dbl(
            .x = ref_components, .y = test_components,
            ~ gmm_abs_weight_diff(
                gmm1 = gmms$Case, gmm2 = gmms$Control,
                i = .x, j = .y
            )
        ),
        delta_sigma = purrr::map2_dbl(
            .x = ref_components, .y = test_components,
            ~ gmm_forstners_distance(
                gmm1 = gmms$Case, gmm2 = gmms$Control,
                i = .x, j = .y
            )
        )
    )
pair_map_clean


pair_map_clean %>%
    dplyr::mutate(test_components = as.character(test_components)) %>%
    ggplot(aes(x = test_components, y = delta_mu)) +
    geom_segment(aes(
        x = test_components, xend = test_components,
        y = 0, yend = delta_mu
    ), color = "grey") +
    geom_point( color="orange", size=4) +
    coord_flip() +
    labs( y = expression(Delta * mu[i])) +
    theme_set(theme_light())


gmm_meta_df %>% glimpse
C_freq_df <- gmm_meta_df %>%
    filter(grepl("TR", C)) %>%
    group_by(aligned_mixture, C) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(relative_freq = n / sum(n)) %>%
    mutate(aligned_mixture = gsub("Mixture_", "", aligned_mixture)) %>%
    dplyr::rename(test_components = aligned_mixture)

p_c_barplot <- C_freq_df %>%
    ggplot(aes(x = relative_freq, y = test_components)) +
    geom_bar(aes(fill = C), stat = "identity") +
    theme_light() +
    scale_fill_npg() +
    labs(y = "Aligned_Mixture", x = "Relative Abundance")

pair_map_clean_plot <- pair_map_clean %>%
    mutate(test_components = as.character(test_components))

# Example usage:
p_du <- plot_gmm_param(pair_map_clean_plot, "test_components", "delta_mu") +
    labs(y = expression(Delta * mu[i])) + labs(x = NULL) + theme(axis.text.y = element_blank())
p_dw <- plot_gmm_param(pair_map_clean_plot, "test_components", "delta_w") +
    labs(y = expression(Delta * omega[i])) + labs(x = NULL) + theme(axis.text.y = element_blank())
p_ds <- plot_gmm_param(pair_map_clean_plot, "test_components", "delta_sigma") +
    labs(y = expression(Delta * Sigma[i])) + labs(x = NULL) + theme(axis.text.y = element_blank())

gmm_param_summary <- p_c_barplot %>% 
    insert_right(p_du, width = 0.8) %>% 
    insert_right(p_dw, width = 0.8) %>%
    insert_right(p_ds, width = 0.8)

# ggsave(
#     filename = glue("{wkdir}/figures/AIRR/gmm-param-summary_{modeln}_{ref_group}_{Sys.Date()}.png"),
#     gmm_param_summary,
#     width = 9, height = 6
# )


# # Plotting and saving UMAPs
# p_ggumaps <- c("C", "mixture", "grouping") %>%
#     purrr::set_names() %>%
#     purrr::map(
#         ~ ggumap(gmm_res$compiled , color_by = .)
#     )
# names(p_ggumaps) %>%
#     purrr::walk(
#         ~ ggsave(
#             filename = glue(
#                 "{wkdir}/figures/AIRR/",
#                 "gmm_umap_{rec_type}_{modeln}_{Sys.Date()}_{.}.png"
#             ),
#             p_ggumaps[[.]], width = 11, height = 9
#         )
#     )

# # Converting files to plotly objects and saving
# p_ggumaps_plotly <- purrr::map(p_ggumaps,
#     ~ ggplotly(.x, tooltip = c("Label"))
# )

# names(p_ggumaps_plotly) %>%
#     purrr::walk(
#         ~ saveRDS(
#             p_ggumaps_plotly[[.]],
#             glue(
#                 "{wkdir}/figures/AIRR/",
#                 "gmm_umap_{modeln}_{Sys.Date()}_{.}.rds"
#             )
#         )
#     )


# saveRDS(
#     gmm_res,
#     glue("{wkdir}/data/interim/airr/",
#         "GMMAligner_{modeln}_Ref-{ref_group}",
#         "_Test-{test_group}_{Sys.Date()}.rds"
#     )
# )


# Perform ICA Analysis - Save data
# Compare PCA vs full embeddings input for ICA analysis


# Construct a loss function for ICA to determine the optimal number of components


# Build Gaussian Mixture Modelings on ICA and PCA embeddings and compare results
# Save Results


# Align Dtest to Dref using JDs

# Perform Parameter Analysis on aligned clusters










# # alternative tidyverse approach
# future::plan("multisession", workers = 100)
# tic()
# trust4_reports_clean_tcrs_stitched <- trust4_reports_clean_tcrs %>%
#     sample_n(10) %>%
#     mutate(
#         stitchr_aa = furrr::future_pmap(
#             list(V, J, CDR3nt),
#             ~ stitchr_wrapper(
#                 V = ..1,
#                 J = ..2,
#                 CDR3nt = ..3,
#                 output_fmt = "AA",
#                 wkdir = wkdir
#             )[1]
#         )
#     )
# toc()
# future::plan("sequential")



# trill testtune 1 finetune esm2_t12_35M {fasta} --epochs 100 --batch_size

# glue(
#         "cd {tmp_loc} &&",
#         " mamba run -n trill_v152",
#         " trill {run_name} {ngpu} embed {model} {fasta} --avg &&",
#         " mv * {output_dir}/"
#         )








# ______________________________________________________________________________






# Core analyses ---------------



library(RColorBrewer)
library(reshape)


ps <- readRDS(
  glue(
    "data/processed/phyloseq_objects/raw/",
    "2023-10-19_WGS_KrakenUnique_phyloseq.rds"
  )
)

ps_trim <- readRDS(
  glue(
    "data/processed/phyloseq_objects/decontaminated/",
    "2024-02-05_WGS_KrakenUnique_phyloseq.rds"
  )
)



# get a list of samples with fewer than 1000 counts total
ps_plot <- ps_trim %>% 
  microbiome::transform("compositional")

sample_sums <- ps_trim %>%
  microbiome::abundances() %>%
  as.data.frame() %>%
  colSums()



# # low_read_samples <- names(sample_sums[sample_sums < 1000])
# pseq <- analyses_ps %>% 
#   subset_samples(analysis_experiment.type == "amplicon") %>% 
#   subset_samples(analysis_accession %nin% low_read_samples) %>% 
#   microbiome::transform("compositional") %>%
#   core(detection = 1e-5, prevalence = 0.1)
# # get a list of samples with fewer than 1000 counts total 
# sample_sums_rnd2 <- microbiome::abundances(pseq) %>% as.data.frame() %>% colSums()
# low_read_samples_rnd2 <- names(sample_sums_rnd2[sample_sums_rnd2 == 0])
# ps_trim <- pseq %>% subset_samples(analysis_accession %nin% low_read_samples_rnd2)


prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(0.001), log10(1), length = 20), 3)

p_core_taxa <-
  ps_plot %>%
  plot_core(
    plot.type = "heatmap",
    colours = viridis::viridis_pal(option = "H")(8),
    prevalences = prevalences,
    detections = detections,
    min.prevalence =min(prevalence(ps_plot, sort = TRUE))
  ) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))", y= "Taxa") +
  #Adjusts axis text size and legend bar height
  theme(
    axis.text.y = element_blank(), # element_text(size = 8, face = "italic"),
    axis.text.x.bottom = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.key.height = unit(1, 'cm'),
    axis.ticks.y = element_blank()
  )

ggsave(
  glue("{wkdir}/figures/readqc/{Sys.Date()}_decon-Abundance-Prev-scatterplot.png"),
  p_core_taxa,
  width = 4, height = 4
)







# Prevalence / Abundance Plot
kuq_decon <- readRDS(
    glue(
        "{wkdir}/data/interim/kraken_results/",
        "2024-02-05_WGS_KrakenUniq_results_with_lineage_decon.rds"
    )
)
taxid_map <- kuq_decon %>%
  select(taxid, rank, name) %>%
  mutate(taxid = as.character(taxid)) %>%
  distinct()

prev_df <- data.frame("prevalence" = prevalence(ps_trim)) %>%
  rownames_to_column("taxa") %>%
  glimpse()

abund_prev_df <- ps_trim %>%
  transform("compositional") %>%
  microbiome::abundances() %>%
  as.data.frame() %>%
  rownames_to_column("taxa") %>%
  pivot_longer(-taxa, names_to = "sample_id", values_to = "counts") %>%
  group_by(taxa) %>%
  dplyr::summarize(mean_abundance = mean(counts)) %>%
  left_join(prev_df) %>%
  left_join(taxid_map, by = c("taxa" = "taxid")) %>%
  # filter(taxa %in% genus_ids) %>%
  glimpse()

p_abund_prev <- abund_prev_df %>%
  filter(rank %in% c("genus", "species", "strain")) %>%
  ggplot(aes(x = mean_abundance, y = prevalence)) +
  geom_point(alpha = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~rank, scales = "free") +
  labs(x = "Mean Relative Abundance", y = "Prevalence") +
  theme_minimal()

ggsave(
  glue("{wkdir}/figures/readqc/",
  "{Sys.Date()}_decon-Abundance-Prev-scatterplot.png"),
  p_abund_prev,
  width = 8, height = 3
)


plotly_pca <- function(df, color_by) {
    plot_ly(
        data = df,
        x = ~PC1,
        y = ~PC2,
        z = ~PC3,
        mode = "markers", marker = list(size = 3),
        text = ~ paste(
            "ID:", Label, "<br>", "C:", C, "<br>",
            "mixture:", mixture, "<br>", "grouping:", grouping
        )
    ) %>%
        add_markers(
            color = ~ .data[[color_by]],
            colors = pal_npg(palette = c("nrc"), alpha = 1)(9)
        ) %>%
        layout(scene = list(
            xaxis = list(title = "PC 1"),
            yaxis = list(title = "PC 2"),
            zaxis = list(title = "PC 3")
        ))
}


D <- ps_trim %>%
  microbiome::transform("clr") %>%
  abundances() %>%
  as.matrix() %>% 
  t()

D_pca <- prcomp(D, rank. = 4)
pca_df <- D_pca[["x"]] %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>% 
    glimpse()

meta_df <- ps_trim %>% 
  meta() %>%
  as.data.frame() %>%
  glimpse()

pca_plot_df <- pca_df %>%
  left_join(meta_df) 

p_pca1 <- pca_plot_df %>%
  ggplot(aes(x = PC1, y = PC2, color = case_control_other_latest)) +
  geom_point() +
  theme_light() +
  labs(
    x = glue("PC 1 [{(pcares$importance[2, 1] * 100)} %]"),
    y = glue("PC 2 [{(pcares$importance[2, 2] * 100)} %]"),
    color = ""
  ) +
  scale_color_npg() +
  theme(legend.position = "top")

p_pca2 <- pca_plot_df %>%
  ggplot(aes(x = PC1, y = PC2, color = study)) +
  geom_point() +
  theme_light() +
  labs(
    x = glue("PC 1 [{(pcares$importance[2, 1] * 100)} %]"),
    y = glue("PC 2 [{(pcares$importance[2, 2] * 100)} %]"),
    color = ""
  ) +
  scale_color_npg() +
  theme(legend.position = "top")

library(patchwork)
combined_pca <- p_pca1 + p_pca2
ggsave(
  glue("{wkdir}/figures/readqc/{Sys.Date()}_decon-PCA-plot.png"),
  combined_pca,
  width = 8, height = 4
)

# library(plotly)
# kraken_unique_pca3d <- plot_ly(
#     data = pca_plot_df,
#     x = ~PC1,
#     y = ~PC2,
#     z = ~PC3,
#     mode = "markers", marker = list(size = 3),
#     text = ~ paste(
#         "ID:", sample_id, "<br>", "study:", study, "<br>",
#         "group:", case_control_other_latest
#     )
# ) %>%
#     add_markers(
#         color = ~ case_control_other_latest,
#         colors = pal_npg(palette = c("nrc"), alpha = 1)(9)
#     ) %>%
#     layout(scene = list(
#         xaxis = list(title = "PC 1"),
#         yaxis = list(title = "PC 2"),
#         zaxis = list(title = "PC 3")
#     ))

# saveRDS(
#   kraken_unique_pca3d,
#   glue("{wkdir}/data/interim/")
# )





  
metadat <- meta(ps_trim)
ps_pd <- ps_trim %>%
  subset_samples(case_control_other_latest == "Case") %>%
  abundances()
ps_cont <- ps_trim %>%
  subset_samples(case_control_other_latest == "Control") %>%
  abundances()

require(foreach)
require(doParallel)
require(progress)

cores = detectCores()
cl <- makeCluster(cores[1] - 2)
registerDoParallel(cl)

aucs <- tibble()
roc2add <-
  foreach(
    feat = rownames(ps_pd),
    .combine = 'rbind',
    .packages = c('magrittr', 'pROC'),
    .options.snow = opts
  ) %dopar% {

    x <- ps_pd[feat, ] %>% t()
    y <- ps_cont[feat, ] %>% t()
    
    # AUROC
    rocdata <-
      c(roc(
        controls = y,
        cases = x,
        direction = '<',
        ci = TRUE,
        auc = TRUE
      )$ci)
    data.frame(
      "feature" = feat,
      "ci_lower" = rocdata[1],
      "auroc" = rocdata[2],
      "ci_upper" = rocdata[3]
    )
  }

stopCluster(cl)
aucs <- rbind(aucs, roc2add)


krakenq_res <- glue("{wkdir}/data/interim/KrakenUnique_analysis")
dir.create(krakenq_res, showWarnings = FALSE)

saveRDS(
  aucs,
  glue("{krakenq_res}/aucs_case_control_other_latest_{Sys.Date()}.rds")
)

aucs %>% glimpse

# aurocs_map <- aucs %>%
#   dplyr::select(feature, data_level, ref_DB) %>%
#   distinct()

aucs_shared_loose <- aucs %>%
  dplyr::mutate(auroc_center = auroc - 0.5) %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(mean_auroc = mean(auroc),
                   mean_auroc_cent = mean(auroc_center),
                   median = median(auroc),
                   n = n()) %>%
  left_join(aurocs_map, by = "feature") %>%
  dplyr::group_by(data_level, ref_DB) %>%
  filter(mean_auroc_cent > 0.1) %>%
  slice_max(order_by = abs(mean_auroc_cent), n = 50)


aucplot_full_dist <- aucs %>%
  ggplot(aes(
    x = fct_reorder(feature, -auroc), y = auroc,
    ymin = ci_lower, ymax = ci_upper,
    group = feature
  )) +
  geom_pointrange(aes(group = feature),
    shape = 24, stroke = 0.2, colour = "grey", fill = "red"
  ) +
  geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.8) +
    coord_flip() +
    theme_light() +
    labs(y = "AUROC", x = "Taxa") +
    theme(
      axis.text.y = element_blank(),
      panel.grid = element_blank()
    )

aucplot_tophits <- aucs %>%
  mutate(abs_diff = abs(auroc - 0.5)) %>%
  filter(abs_diff > 0.02) %>%
  left_join(taxid_map, by = c("feature" = "taxid")) %>%
  ggplot(aes(
    x = fct_reorder(name, -auroc), y = auroc,
    ymin = ci_lower, ymax = ci_upper,
    group = name
  )) +
  geom_pointrange(aes(group = name),
    shape = 24, stroke = 0.2, colour = "grey", fill = "red"
  ) +
  geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.8) +
    coord_flip() +
    theme_light() +
    labs(y = "AUROC", x = "Taxa") +
    theme(
      panel.grid = element_blank()
    )


ggsave(
  glue("{wkdir}/figures/readqc/{Sys.Date()}_decon-AUROC-plot.png"),
  aucplot_full_dist,
  width = 4, height = 5
)

ggsave(
  glue("{wkdir}/figures/readqc/{Sys.Date()}_decon-AUROC-plot-top-hits.png"),
  aucplot_tophits,
  width = 6, height = 4
)






aucs_studies <- tibble()
for (cohort in c("BioFIND", "HBS", "PDBP", "PPMI")) {
  cat(cohort, "\n")
  ps <- ps_trim %>%
    subset_samples(study %in% cohort)
  metadat <- meta(ps)
  ps_pd <- ps %>%
    subset_samples(case_control_other_latest == "Case") %>%
    abundances()
  ps_cont <- ps %>%
    subset_samples(case_control_other_latest == "Control") %>%
    abundances()

  cores = detectCores()
  cl <- makeCluster(cores[1] - 2)
  registerDoParallel(cl)

  roc2add <-
    foreach(
      feat = rownames(ps_pd),
      .combine = 'rbind',
      .packages = c('magrittr', 'pROC')
    ) %dopar% {

      x <- ps_pd[feat, ] %>% t()
      y <- ps_cont[feat, ] %>% t()

      # AUROC
      rocdata <-
        c(roc(
          controls = y,
          cases = x,
          direction = '<',
          ci = TRUE,
          auc = TRUE
        )$ci)
      data.frame(
        "feature" = feat,
        study = cohort,
        "ci_lower" = rocdata[1],
        "auroc" = rocdata[2],
        "ci_upper" = rocdata[3]
      )
    }

  stopCluster(cl)
  aucs_studies <- rbind(aucs_studies, roc2add)
}

aucs_studies %>% glimpse

saveRDS(
  aucs_studies,
  glue(
      "{krakenq_res}/",
      "aucs_COHORT-SPECIFIC_case_control_other_latest_{Sys.Date()}.rds"
  )

aucs <- aucs_studies

aurocs_map <- aucs %>%
  dplyr::select(feature, study) %>%
  distinct()

aucs_shared_loose <- aucs %>%
  dplyr::mutate(auroc_center = auroc - 0.5) %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(mean_auroc = mean(auroc),
                   mean_auroc_cent = mean(auroc_center),
                   median = median(auroc),
                   n = n()) %>%
  left_join(aurocs_map, by = "feature") %>%
  dplyr::group_by(study) %>%
  filter(mean_auroc_cent > 0.01) %>%
  slice_max(order_by = abs(mean_auroc_cent), n = 25)


aucplot_cohorts <- aucs %>%
    # filter(feature %in% aucs_shared_loose$feature) %>%
    # left_join(taxid_map, by = c("feature" = "taxid")) %>%
  ggplot(aes(x=fct_reorder(feature, -auroc), y = auroc, fill = study,
             ymin=ci_lower, ymax=ci_upper, group = feature)) +
  geom_pointrange(aes(group = feature),
                  position = position_jitterdodge(jitter.height = 0),
                  shape = 24, stroke = 0.2, colour="grey") +
  geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.8) +
      scale_fill_d3() +
      labs(y = "AUROC", x = "Taxa") +
      coord_flip() +
      theme_light() +
      theme(
        #   axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid = element_blank()
      )

ggsave(
  glue("{wkdir}/figures/readqc/",
  "{Sys.Date()}_decon-AUROC-plot-individual-cohorts.png"),
  aucplot_cohorts,
  width = 5, height = 6
)

aucplot_cohorts_tophits <- aucs %>%
    filter(feature %in% aucs_shared_loose$feature) %>%
    left_join(taxid_map, by = c("feature" = "taxid")) %>%
  ggplot(aes(x=fct_reorder(name, -auroc), y = auroc, fill = study,
             ymin=ci_lower, ymax=ci_upper, group = name)) +
  geom_pointrange(aes(group = name),
                  position = position_jitterdodge(jitter.height = 0),
                  shape = 24, stroke = 0.2, colour="grey") +
  geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.8) +
      scale_fill_d3() +
      labs(y = "AUROC", x = "Taxa") +
      coord_flip() +
      theme_light() +
      theme(
        #   axis.text.x = element_blank(),
          #   axis.text.y = element_blank(),
          panel.grid = element_blank()
      )

ggsave(
    glue(
        "{wkdir}/figures/readqc/",
        "{Sys.Date()}_decon-AUROC-plot-individual-cohorts-tophits.png"
    ),
    aucplot_cohorts_tophits,
    width = 9, height = 6
)


kuq %>%
    filter(grepl("Paenibacillus", name)) %>%
    pull(name) %>%
    unique()


kuq_decon %>%
    filter(grepl("Paenibacillus", name)) %>%
    pull(name) %>%
    unique()



ps_trim <- readRDS(
    glue(
        "{wkdir}/data/processed/",
        "phyloseq_objects/decontaminated/",
        "2024-02-25_WGS_KrakenUnique_phyloseq.rds"
    )
)

alpha_df <- ps_trim %>%
    microbiome::alpha("observed") %>%
    as.data.frame() %>%
    bind_cols(meta(ps_trim)) %>%
    glimpse



alpha_df %>%
    ggplot(aes(x = case_control_other_latest, y = observed)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
    scale_y_log10() +
    theme_light()

library(ggdist)
library(gghalves)
library(ggpubr)

grp_comparisons <- list(
    c("Control", "Case"),
    c("Case", "Other"),
    c("Control", "Other")
)

alpha_df$study %>% table
p_alpha_dist <- alpha_df %>%
    filter(study %in% c("BioFIND", "HBS", "PDBP", "PPMI")) %>%
  ggplot(aes(x = case_control_other_latest, y = observed)) +
  theme_bw() +
    scale_color_brewer(palette = "Greys") +
    stat_interval(
      .width = c(.25, 0.5, .75, 1),
      height = 5, show.legend = T
    ) +
    stat_halfeye(aes(fill = case_control_other_latest),
      .width = 0,
      scale = 0.5, size = 0.2, alpha = 0.4, height = 0.4,
      justification = -.25,
      point_alpha = 1, point_fill = "white",
      point_color = "white", shape = 23
    ) +
    geom_half_point(
      side = "l", range_scale = .5, alpha = 0.4,
      shape = 21, color = "black", size = 0.5
    ) +
    stat_compare_means(comparisons = grp_comparisons) +
    # stat_compare_means() +
  labs(x = NULL, y = "Observed Taxa") +
      guides(size = "none", fill = "none") +
      scale_fill_d3() +
      scale_y_log10() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.border = element_blank()
      )

ggsave(
  glue(
    "{wkdir}/figures/alpha-diversity/",
    "observed-taxa-distributions_{Sys.Date()}.png"),
    p_alpha_dist,
  width = 6, height = 6
)

# summary stats for observed taxa
alpha_df %>%
    filter(study %in% c("BioFIND", "HBS", "PDBP", "PPMI")) %>%
    group_by(case_control_other_latest) %>%
    dplyr::summarize(
        mean_observed = mean(observed),
        median_observed = median(observed)
    )

wgs_flagstat_meta <- readRDS(
  glue(
    "{wkdir}/data/processed/readqc/",
    "2024-02-25_WGS_flagstat-metrics_phyloseq-formatted.rds"
  )
)
all(wgs_flagstat_meta$sample_id == rownames(meta(ps_trim)))

# adding flagstat metrics to phyloseq object
phyloseq::sample_data(ps_trim)[colnames(wgs_flagstat_meta)] <-
    wgs_flagstat_meta

mapped_human_reads <- meta(ps_trim)$flagstat_01_cram_mapped

norm_abund <- ps_trim %>%
    abundances() %>% 
    as.data.frame()/mapped_human_reads

norm_abund %>% 
    colSums() %>%
    hist(100)
    summary



alpha_df %>% glimpse

alpha_df %>%
    ggplot(aes(observed)) +
    geom_density(aes(color = study)) +
    scale_x_log10() +
    theme_light()


alpha_df %>%
    ggplot(aes(x = age_at_baseline, y = observed)) +
    scale_y_log10() +
    geom_smooth(method = "lm") +
    geom_point() 

alpha_df %>%
    ggplot(aes(x = mds_updrs_part_ii_summary_score, y = observed)) +
    scale_y_log10() +
    geom_smooth(method = "lm") +
    geom_point() 

alpha_df %>%
    ggplot(aes(x = mds_updrs_part_iii_summary_score, y = observed)) +
    scale_y_log10() +
    geom_smooth(method = "lm") +
    geom_point() 

alpha_df %>%
    ggplot(aes(x = mds_updrs_part_iv_summary_score, y = observed)) +
    scale_y_log10() +
    geom_smooth(method = "lm") +
    geom_point() 






















# cd /resnick/groups/MazmanianLab/jboktor/Downloads/RefDBs/PalmDB

# kb count \
#     -t 8 \
#     -o /resnick/groups/MazmanianLab/jboktor/PDMBS/workflow/RNASEQ/test_palmdb_bwa_dlist_SINGLE \
#     --aa \
#     -m 4G \
#     -i index_GRCh38_PalmDB.idx \
#     -g palmdb_clustered_t2g.txt \
#     --parity single \
#     -x BULK \
#     /resnick/groups/MazmanianLab/jboktor/PDMBS/workflow/RNASEQ/clean_reads/PP-40737-BLM0T1_R1.fastq.gz








# mtx_data <- Seurat::ReadMtx(
#     mtx = glue("{results_path}/cells_x_genes.mtx"),
#     features = glue("{results_path}/cells_x_genes.genes.names.txt"),
#     cells = glue("{results_path}/cells_x_genes.barcodes.txt"),
#     feature.column = 1,
#     mtx.transpose = TRUE
# )


# palm_interim <- glue("{wkdir}/data/interim/PalmDB_mapping")
# # dir.create(palm_interim, showWarnings = FALSE)
# palmddb_res <- glue("{pdmbs_dir}/workflow/RNASEQ/palmdb_bwa_dlist")
# sample_id <- list.dirs(palmddb_res, recursive = FALSE) %>%
#     basename()


# progress_bar <- progress::progress_bar$new(
#     format = "  :spin [:bar] :percent eta: :eta",
#     total = length(sample_id)
# )
# palm_res <- tibble()
# for (id in sample_id) {
#     progress_bar$tick()
#     matrix <- readMM(
#         glue("{palmddb_res}/{id}/counts_unfiltered/cells_x_genes.mtx")
#     )
#     # matrix funct decreases index by 1
#     new_df <- data.frame(
#         direction = matrix@i + 1,
#         gene_id = matrix@j + 1,
#         count = matrix@x
#     ) %>%
#         dplyr::group_by(gene_id) %>%
#             dplyr::summarise(count = sum(count)) %>%
#             mutate(sample_id = id)
#     palm_res %<>% bind_rows(new_df)
# }

# count_df <- palm_res %>%
#     as.data.table() %>%
#     tidyr::pivot_wider(
#         names_from = gene_id,
#         values_from = count,
#         values_fill = 0
#     ) %>%
#     glimpse()

# saveRDS(
#     count_df,
#     glue("{palm_interim}/count_matrix_BWA-dlist-T2T-DNA_{Sys.Date()}.rds")
# )













# tst <- read_mtx(
tst <- Seurat::ReadMtx(
    mtx = "/resnick/groups/MazmanianLab/jboktor/PDMBS/workflow/RNASEQ/test_palmdb_bwa_dlist/counts_unfiltered/cells_x_genes.mtx",
    features = "/resnick/groups/MazmanianLab/jboktor/PDMBS/workflow/RNASEQ/test_palmdb_bwa_dlist/counts_unfiltered/cells_x_genes.genes.names.txt",
    cells = "/resnick/groups/MazmanianLab/jboktor/PDMBS/workflow/RNASEQ/test_palmdb_bwa_dlist/counts_unfiltered/cells_x_genes.barcodes.txt",
    feature.column = 1,
    mtx.transpose = TRUE
)


ps_taxtable <- palm_meta %>%
    mutate_all(
        ~ gsub("\\.", "", .) %>% tolower()
    ) %>%
    dplyr::select(-strandedness) %>%
    column_to_rownames("rep_ID") %>%
    as.matrix() %>%
    phyloseq::tax_table()

ps <- phyloseq(sotu, ps_meta, ps_taxtable)




# Use traqnsform function of microbiome to convert it to rel abun.
ps1.com.fam.rel <- microbiome::transform(ps, "compositional")

plot.composition.relAbun <- plot_composition(ps1.com.fam.rel,
    sample.sort = "case_control_other_latest"
    # x.label = "env_material"
)

plot.composition.relAbun <- plot.composition.relAbun + theme(legend.position = "bottom") 
plot.composition.relAbun <- plot.composition.relAbun + scale_fill_brewer("Family", palette = "Paired") + theme_bw() 
plot.composition.relAbun <- plot.composition.relAbun + theme(axis.text.x = element_text(angle = 90)) 
plot.composition.relAbun <- plot.composition.relAbun + ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))

print(plot.composition.relAbun)

ggsave(
    filename = glue("{wkdir}/figures/Viral/",
    "rel_abun_composition_{Sys.Date()}.png"),
    plot.composition.relAbun,
    width = 8, height = 6
)





library(infotheo)
library(mpmi)

metadata_df <- readRDS(
    glue("{wkdir}/data/interim/metadata/2023-07-14_phyloseq-metadata.rds")
)$RNASEQ %>% 
    glimpse()

?mmi.pw()
?cmi.pw()
?dmi.pw()


pairwise_mi <- function(x, y) {
  # Check if vectors are continuous (double) or discrete
  x_is_continuous <- is.double(x)
  y_is_continuous <- is.double(y)
  
  # Determine which function to use based on data types
  if (x_is_continuous && y_is_continuous) {
    result <- cmi.pw(x, y)
    func_name <- "cmi"
  } else if (x_is_continuous && !y_is_continuous) {
    result <- mmi.pw(x, y)
    func_name <- "mmi"
  } else if (!x_is_continuous && y_is_continuous) {
    # Swap x and y to ensure continuous is first for mmi
    result <- mmi.pw(y, x)
    func_name <- "mmi"
  } else {
    result <- dmi.pw(x, y)
    func_name <- "dmi"
  }
  # Add function name to result
  c(result, "func" = func_name)
}


# mutinformation(c(1, 2, 3), c(1, 2, 3))
# pairwise_mi(c(1, 2, 3), c(1, 2, 3))

# normalized_mi <- function(x, y) {
#     require(infotheo)
#     mutinformation(x, y) / sqrt(entropy(x) * entropy(y))
# }


immunarch_root <- glue("{wkdir}/data/interim/airr/immunarch")
immunarch_dir <- glue("{immunarch_root}/metrics")
immunarch_figs <- glue("{wkdir}/figures/AIRR/immunarch")


# # rarefaction vis
# imm_raref <- readRDS(
#     glue("{immunarch_dir}/rep-clonality-rare_2024-11-21.rds")
#     )

# p_raref <- immunarch::vis(imm_raref,
#     .by = "case_control_other_latest",
#     .meta = immdata$meta
# )

# ggsave(
#     p_raref,
#     filename = glue("{immunarch_dir}/rarefaction_2024-11-22.png"),
#     width = 6, height = 6
# )



# plot_df <- metadata_df %>% 
#     left_join(immunarch_res, by = c("sample_id" = "Sample")) %>%
#     glimpse()

# p_hill <- plot_df %>%
#     ggplot(aes(x = Q, y = Value, color = case_control_other_latest)) +
#     geom_point() +
#     geom_smooth(method = "loess") +
#     scale_color_d3() +
#     scale_y_log10() +
#     theme_bw() +
#     labs(x = "Diversity Estimation", y = "Q Value") +
#     theme(legend.position = "bottom")

# ggsave(
#     p_hill,
#     filename = glue("{immunarch_figs}/rarefaction_2024-11-22.png"),
#     width = 9, height = 6
# )













metric_paths <- list.files(immunarch_dir) %>%
    discard(~ grepl("overlap", .))

imm_res <- metric_paths %>% 
    purrr::set_names() %>% 
        purrr::map(
            ~ readRDS(glue("{immunarch_dir}/{.}"))
            )
imm_res %>% glimpse
imm_res %>% names


imm_res[["rep-diversity-div_2024-11-22.rds"]] %>% glimpse

imm_res_div_df <- imm_res[["rep-diversity-div_2024-11-22.rds"]] %>%
    as.data.frame() %>%
    dplyr::rename(sample_id = Sample) %>%
    left_join(metadata_df, by = c("sample_id")) %>%
    glimpse()

future::plan("multisession", workers = 16)
mmi_res_div <- colnames(metadata_df) %>% 
    purrr::set_names() %>%
    furrr::future_map(
        ~ possibly(pairwise_mi, otherwise = NA_real_)(
            imm_res_div_df$Value, 
            imm_res_div_df[[.]]
        )
)


div_screen_df <- enframe(mmi_res_div, name = "variable") %>%
 mutate(
    mi = map_dbl(value, 
        ~if(is.atomic(.x)) NA_real_ else .x$mi %||% NA_real_),
    bcmi = map_dbl(value, 
        ~if(is.atomic(.x)) NA_real_ else .x$bcmi %||% NA_real_),
    zvalue = map_dbl(value, 
        ~if(is.atomic(.x)) NA_real_ else .x$zvalue %||% NA_real_),
    func = map_chr(value, 
        ~if(is.atomic(.x)) NA_character_ else .x$func %||% NA_character_)
  ) %>%
  select(-value)

saveRDS(
    div_screen_df,
    glue("{immunarch_root}/div-screen_{Sys.Date()}.rds")
)


div_screen_df <- readRDS(
    glue("{immunarch_root}/div-screen_2024-11-25.rds")
)

meta_cats <- readRDS(
    glue("{wkdir}/data/interim/metadata/",
    "2023-06-06_metadata_categories_dataframe.rds")
) %>% 
    glimpse()

p_div_screen <- div_screen_df %>% 
    left_join(meta_cats, by = c("variable" = "metadata")) %>%
    mutate(label = ifelse(rank(-mi) <= 10, variable, "")) %>%
    ggplot(aes(x = mi, y = bcmi)) +
    geom_point(aes(color = metadata_class)) +
    ggrepel::geom_text_repel(aes(label = label), 
                             show.legend = FALSE,
                             max.overlaps = Inf) +
    theme_bw()

ggsave(
    p_div_screen,
    filename = glue("{immunarch_figs}/div-screen_{Sys.Date()}.png"),
    width = 9, height = 8
)


# Gene usage analysis
geneuseage_df <- imm_res[["gene-usage-HomoSapiens.TRBJ_2024-11-22.rds"]] %>% 
    filter(Names != ".") %>% 
    pivot_longer(
        !Names,
        names_to = "sample_id",
        values_to = "useage"
    )


metadata_df_minimal <- metadata_df %>% 
    dplyr::select(sample_id, participant_id, case_control_other_latest, study)

geneuseage_df_meta <- geneuseage_df %>% 
    left_join(metadata_df_minimal, by = c("sample_id")) %>%
    glimpse()


library(nlme)

unique(geneuseage_df$Names)

future::plan("multisession", workers = 16)
geneusage_lme_res <- unique(geneuseage_df$Names) %>% 
    purrr::set_names() %>%
    furrr::future_map(
        ~ possibly(
            ~ nlme::lme(
                useage ~ case_control_other_latest,
                random = ~ 1 | study/participant_id,
                data = geneuseage_df_meta %>% filter(Names == .x),
                na.action = na.omit
            ),
            otherwise = NULL
        )(.x)
    )

saveRDS(
    geneusage_lme_res,
    glue("{immunarch_root}/geneusage_lme_res_{Sys.Date()}.rds")
)



geneusage_lme_df <- names(geneusage_lme_res) %>% 
    purrr::set_names() %>% 
    purrr::map(
        ~ broom.mixed::tidy(geneusage_lme_res[[.]])
    ) %>% 
    bind_rows(.id = "gene") %>%
    glimpse()

geneusage_lme_df$term %>% unique

p_geneuseage_est_other <- geneusage_lme_df %>% 
    filter(term == "case_control_other_latestControl") %>%
    mutate(padj = p.adjust(p.value, method = "fdr")) %>%
    mutate(sig = padj <= 0.1) %>%
    mutate(rec_type = ifelse(grepl("IG", gene), "BCR", "TCR")) %>%
    ggplot(aes(x = estimate, y = fct_reorder(gene, estimate))) +
    geom_pointrange(aes(
        xmin = estimate - std.error,
        xmax = estimate + std.error,
        color = sig)) +
    facet_grid(rec_type ~ ., scales = "free", space = "free") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
        x = "Linear Mixed Model Estimate (Control vs Case)", 
        y = NULL) +
    theme_bw()

ggsave(
    p_geneuseage_est_other,
    filename = glue("{immunarch_figs}/geneuseage_est_other_{Sys.Date()}.png"),
    width = 9, height = 12
)

p_geneuseage_TRBJ1_1 <- geneuseage_df_meta %>% 
    filter(Names == "TRBJ1-1") %>%
    ggplot(aes(y = useage, x = case_control_other_latest)) +
    geom_point(
        position = position_jitter(width = 0.2, height = 0), 
        alpha = 0.3) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    labs(x = NULL, y = "TRBJ1-1 Gene Usage") +
    theme_bw()

ggsave(
    p_geneuseage_TRBJ1_1,
    filename = glue("{immunarch_figs}/geneuseage_TRBJ1_1_{Sys.Date()}.png"),
    width = 4, height = 6
)

# # ANALYZING BCR REPERTOIRE DATA
# immdata %>% glimpse

# immdata$data[[1]] %>% glimpse


# bcr_immdata <- repFilter(immdata, "by.clonotype",
#   list(V.name = include("IGH")),
#   .match = "substring"
# )

# # Filter for BCR data (immunoglobulin genes)
# # Try filtering with .strict parameter to handle mismatched samples
# bcr_immdata <- repFilter(immdata,
#   .method = "by.clonotype",
#   .query = list(V.name = include("TRB")),
#   .match = "startswith"
# )


# immdata$data %>% length
# immdata$meta %>% nrow

# # Examine the filtered data and metadata dimensions
# message("Data dimensions:")
# print(length(bcr_immdata$data))
# message("Metadata dimensions:")
# print(nrow(bcr_immdata$meta))
# bcr_immdata$data[[1]] %>% glimpse()






# data(bcrdata)

# bcrdata$data[[1]] %>% glimpse
# bcrdata$meta %>% nrow

# # take clusters that contain at least 1 sequence
# bcr_data <- bcrdata$data
# align_dt <- bcr_data %>%
#   seqCluster(seqDist(bcr_data, .col = 'CDR3.nt', .group_by_seqLength = TRUE),
#              .perc_similarity = 0.6) %>%
#   repGermline(.threads = 1) %>%
#   repAlignLineage(.min_lineage_sequences = 6, .align_threads = 2, .nofail = TRUE)




p1 <- vis(rep_eda$clones, .by = c("case_control_other_latest"), 
    .meta = immdata$meta
)


plot_df <- rep_eda$clones %>% 
    left_join(key_meta_df, by = c("Sample" = "sample_id")) %>%
    glimpse()

p_imarch_clonality <- plot_df %>% 
    ggplot(aes(x = case_control_other_latest, y = Clones)) +
    geom_point(position = 
        position_jitter(width = 0.2, height = 0), alpha = 0.3) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    theme_bw() +
    scale_y_log10() +
    labs(x = NULL, y = "Clonotypes")

ggsave(
    filename = glue("{wkdir}/figures/AIRR/",
    "TRUST4-clonal-metrics_{Sys.Date()}.png"),
    p_imarch_clonality,
    width = 4, height = 6
)


plot_df_volume <- rep_eda$volume %>% 
    left_join(key_meta_df, by = c("Sample" = "sample_id")) %>%
    glimpse()

p_imarch_volume <- plot_df_volume %>% 
    ggplot(aes(x = case_control_other_latest, y = Volume)) +
    geom_point(position = 
        position_jitter(width = 0.2, height = 0), alpha = 0.3) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    theme_bw() +
    scale_y_log10() +
    labs(x = NULL, y = "Clonal Volume")

ggsave(
    filename = glue("{wkdir}/figures/AIRR/",
    "TRUST4-clonal-volume_{Sys.Date()}.png"),
    p_imarch_volume,
    width = 4, height = 6
)




imm_gu <- geneUsage(immdata$data, "hs.trbv", .norm = T)
imm_gu_mds <- geneUsageAnalysis(imm_gu, .method = "mds", .verbose = F)

saveRDS(
    imm_gu,
    glue("{wkdir}/data/interim/airr/immunarch/gene-usage_{Sys.Date()}.rds")
)

saveRDS(
    imm_gu_mds,
    glue("{wkdir}/data/interim/airr/immunarch/gene-usage-mds_{Sys.Date()}.rds")
)


imm_gu <- readRDS(
    glue("{wkdir}/data/interim/airr/immunarch/gene-usage_2024-11-05.rds")
)
imm_gu %>% glimpse


imm_gu_mds <- readRDS(
    glue("{wkdir}/data/interim/airr/immunarch/gene-usage-mds_2024-11-05.rds")
)

# vis(geneUsageAnalysis(imm_gu, "cosine+hclust", .verbose = F))

imm_cl_mds <- geneUsageAnalysis(imm_gu, "js+mds+kmeans", .verbose = F)
p2 <- vis(imm_cl_mds, .plot = "clust")


imm_gu_mds %>% glimpse

imm_gu_mds$points[, 1]


# data.frame(
#     "sample_id" = rownames(imm_gu_mds$points),
#     "DimI" = imm_gu_mds$points[, 1],
#     "DimII" = imm_gu_mds$points[, 2]
# ) %>%

imm_cl_mds$x %>%
    rownames_to_column("sample_id") %>%
    left_join(key_meta_df, by = c("sample_id")) %>%
    ggplot(aes(x = DimI, y = DimII, color = case_control_other_latest)) +
    geom_point() +
    theme_bw()


# imm_gu_mds[is.na(imm_gu_mds)] <- 0


# imm_gu_js <- geneUsageAnalysis(imm_gu, .method = "js", .verbose = F)
# imm_gu_cor <- geneUsageAnalysis(imm_gu, .method = "cor", .verbose = F)

p1 <- vis(imm_gu_js, .title = "Gene usage JS-divergence", .leg.title = "JS", .text.size = 1.5)
p2 <- vis(imm_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size = 1.5)

# p1 + p2

p_full <- p1 + p2

ggsave(
    filename = glue("{wkdir}/figures/Viral/",
    "TRBV-gene-usage-JS-Cor_{Sys.Date()}.png"),
    p_full,
    limitsize = FALSE,
    width = 60, height = 60
)

imm_gu_js[is.na(imm_gu_js)] <- 0

vis(geneUsageAnalysis(imm_gu, "cosine+hclust", .verbose = F))


imm_cl_pca <- geneUsageAnalysis(imm_gu, "js+pca+kmeans", .verbose = F)
imm_cl_mds <- geneUsageAnalysis(imm_gu, "js+mds+kmeans", .verbose = F)
imm_cl_tsne <- geneUsageAnalysis(imm_gu, "js+tsne+kmeans", .perp = .01, .verbose = F)

p1 <- vis(imm_cl_pca, .plot = "clust")
p2 <- vis(imm_cl_mds, .plot = "clust")
p3 <- vis(imm_cl_tsne, .plot = "clust")
p1 + p2 + p3


c <- repClonality(immdata$data, "homeo") %>% vis()  # Visualise the relative abundance of clonotypes

repOverlap(immdata$data) %>% vis()  # Build the heatmap of public clonotypes shared between repertoires
geneUsage(immdata$data[[1]]) %>% vis()  # Visualise the V-gene distribution for the first repertoire
repDiversity(immdata$data) %>% vis(.by = "Status", .meta = immdata$meta)  # Visualise the Chao1 diversity of repertoires, grouped by the patient status

immdata %>% glimpse


#' # Example file path
file_path <- paste0(system.file(package = "immunarch"), "/extdata/db/vdjdb.example.txt")

# Load the database with human-only TRB-only receptors for all known antigens
db <- dbLoad(file_path, "vdjdb", "HomoSapiens", "TRB")

res <- dbAnnotate(immdata$data, db, "CDR3.aa", "cdr3")
res
















# Conclusion, looks like failure isn't correlated at all with size of the sample.


#------------------------------------------------------------------------------


# Extract wide format of clones shared within a sample


trust4_reports_meta_clean <- trust4_reports_meta %>%
    filter(C == "TRBC") %>%
    filter(V != "." & J != ".") %>%
    filter(CDR3aa != "out_of_frame") %>%
    filter(!grepl("_", CDR3aa)) %>%
    dplyr::rename(counts = `#count`) %>%
    glimpse()
# trust4_reports_meta_clean %>% glimpse()

saveRDS(
    trust4_reports_meta_clean,
    glue("{wkdir}/data/interim/airr/",
    "TRUST4-reports-meta-clean_{Sys.Date()}.rds")
)


# Get shared CDR3aa clones within participants
shared_clones <- trust4_reports_meta_clean %>%
  group_by(participant_id, CDR3aa) %>%
  summarise(
    n_samples = n_distinct(sample_id),
    .groups = "drop"
  ) %>%
  # Keep only clones found in multiple samples
  filter(n_samples > 1) %>%
  select(participant_id, CDR3aa)

saveRDS(
    shared_clones,
    glue("{wkdir}/data/interim/airr/",
    "TRUST4-shared-clones_{Sys.Date()}.rds")
)


# sid_pid_pairs <- key_meta_df %>% 
#     dplyr::select(sample_id, participant_id, raw_total_reads) %>% 
#     distinct() %>% 
#     glimpse()

#' Function to format clone stats for a given participant and clone
format_clone_stats <- function(pid, clone_id, wkdir) {
    require(glue)
    require(dplyr)
    require(tidyr)
    require(data.table)

    trust4_reports_meta_clean <- readRDS(
        glue("{wkdir}/data/interim/airr/",
        "TRUST4-reports-meta-clean_2024-12-27.rds")
    )
    shared_clones <- readRDS(
        glue("{wkdir}/data/interim/airr/",
        "TRUST4-shared-clones_2024-12-27.rds")
    )
    key_meta_df <- readRDS(
        glue("{wkdir}/data/interim/airr/key_meta_df_2024-12-27.rds")
    )

    # Filter and select in one step
    clone_stats <- trust4_reports_meta_clean[
        trust4_reports_meta_clean$participant_id == pid & 
        trust4_reports_meta_clean$CDR3aa == clone_id,
        c("sample_id", "counts", "frequency", "CDR3aa", "CDR3nt")
    ]
    
    # Filter pid_meta directly without pipe
    pid_meta <- key_meta_df[key_meta_df$participant_id == pid,]
    
    # Merge and mutate in one step
    res <- merge(pid_meta, clone_stats, by = "sample_id", all.x = TRUE) %>%
        mutate(
            counts = coalesce(counts, 0),
            frequency = coalesce(frequency, 0),
            CDR3aa = coalesce(CDR3aa, clone_id)
        )
    saveRDS(res,
        glue("{wkdir}/data/interim/airr/clone_formatting/",
        "{pid}__{clone_id}.rds")
    )
}

# format_clone_stats(
#     shared_clones$participant_id[1],
#     shared_clones$CDR3aa[1],
#     wkdir
# )

# Setup batchtools registry
cluster_run <- glue("{get_time()}_clone_formatting")
message("\n\nRUNNING: ", cluster_run, "\n")
breg <- makeRegistry(
    file.dir = glue(
        "{wkdir}/.cluster_runs/",
        cluster_run
    ),
    seed = 42
)
breg$cluster.functions <- batchtools::makeClusterFunctionsSlurm(
    template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
    scheduler.latency = 0.1,
    fs.latency = 1
)

# Create parameter list for each method
batchtools_params <- shared_clones %>% 
    mutate(wkdir = wkdir) %>% 
    dplyr::rename(pid = participant_id, clone_id = CDR3aa) %>% 
    mutate(expected_output = glue(
        "{wkdir}/data/interim/airr/clone_formatting/",
        "{pid}__{clone_id}.rds"),
        file_exists = file.exists(expected_output)
    ) %>% 
    filter(!file_exists) %>% 
    dplyr::select(pid, clone_id, wkdir) %>% 
    glimpse()

# Submit jobs
jobs <- batchMap(
    fun = format_clone_stats,
    args = batchtools_params,
    reg = breg
)

jobs[, chunk := chunk(job.id, chunk.size = 1)]
print(jobs[, .N, by = chunk])

submitJobs(jobs,
    resources = list(
        walltime = "10:00",
        memory = "10GB",
        ncpus = 1,
        max.concurrent.jobs = 9999
    )
)


pid_clone_paths <- list.files(
    glue("{wkdir}/data/interim/airr/clone_formatting"),
    full.names = TRUE
)


future::plan("multisession", workers = 30)
clone_df <- pid_clone_paths %>% 
    furrr::future_map(readRDS) %>% 
    bind_rows()

saveRDS(
    clone_df,
    glue("{wkdir}/data/interim/airr/",
    "TRUST4-cleaned-clones-per-pid_{Sys.Date()}.rds")
)



clone_df <- readRDS(
    glue("{wkdir}/data/interim/airr/",
    "TRUST4-cleaned-clones-per-pid_2024-12-27.rds")
) %>% 
    as.data.table() %>% 
    glimpse()



# Fit GLMM for each CDR3aa x participant_id combination
glmm_res <- clone_df %>%
  # Group by participant and CDR3aa sequence
  group_by(participant_id, CDR3aa) %>%
  # Fit model for each group
  group_modify(~{
    # Skip if less than 2 observations
    if(nrow(.x) < 2) return(NULL)
    
    # Fit GLMM with beta distribution for frequency response
    tryCatch({
      model <- glmmTMB::glmmTMB(
        frequency ~ visit_month + raw_total_reads + case_control_other_latest + (1|participant_id),
        data = .x,
        family = beta_family()
      )
      
      # Extract model summary
      summary <- broom.mixed::tidy(model)
      
      # Return coefficients and stats
      data.frame(
        n_timepoints = nrow(.x),
        visit_month_effect = coef(summary(model))$cond["visit_month","Estimate"],
        reads_effect = coef(summary(model))$cond["raw_total_reads","Estimate"],
        case_control_effect = coef(summary(model))$cond["case_control_other_latest","Estimate"]
      )
    }, error = function(e) NULL)
  }) %>%
  ungroup()



library(glmmTMB)
library(easystats)
# library(lme4)
# library(mgcv)



# simple_model <- glmmTMB::glmmTMB(
#     frequency ~  raw_total_reads,
#     data = clone_df
# )

model_nb <- glmmTMB(
  counts ~ visit_month + case_control_other_latest +
    (1 | participant_id) +           # random intercept for participant
    (1 | CDR3aa) +                   # optional random intercept for clone
    offset(log(raw_total_reads)),
  data = clone_df,
  family = "nbinom2"                 # Negative binomial
)

model_nb %>% summary()

saveRDS(
    model_nb,
    glue("{wkdir}/data/interim/airr/",
    "TRUST4-model-nb_2024-12-27.rds")
)



model_nb <- readRDS(
    glue("{wkdir}/data/interim/airr/",
    "TRUST4-model-nb_2024-12-27.rds")
)

model_nb_assumptions <- performance::check_model(model_nb)
# plot(parameters(model_nb))

key_meta_df$raw_total_reads %>% median



# pred_df <- expand.grid(
#   visit_month = seq(0, 24, by = 0.5), 
#   case_control_other_latest = c("Other"),
#   participant_id = NA,
#   CDR3aa = NA,
#   raw_total_reads = 489977824)  # log-scaled sequence of read depths
# )

pid_key_meta <- key_meta_df %>% 
    group_by(participant_id) %>% 
    dplyr::summarise(
        raw_total_reads = median(raw_total_reads),
        case_control_other_latest = first(case_control_other_latest)
    ) %>% 
    filter(participant_id %in% unique(clone_df$participant_id)) %>% 
    glimpse()


pred_df <- expand.grid(
  visit_month = seq(0, 24, by = 0.5),
  participant_id = pid_key_meta$participant_id,
  CDR3aa = NA
) %>% 
    left_join(pid_key_meta, by = "participant_id") %>% 
    glimpse()

# 2) Predict using the fitted model (population-level)
pred_df$pred_counts <- predict(
  model_nb,
  newdata = pred_df,
  type = "response",
  re.form = NA
)

pred_df %>% glimpse

# 3) Plot a heatmap
p_model_heatmap <- pred_df %>% 
    # mutate(raw_total_reads = as.character(raw_total_reads)) %>% 
    ggplot(aes(x = visit_month,
                    y = fct_reorder(participant_id, raw_total_reads),
                    fill = pred_counts)) +
  geom_tile(color = NA) +
  scale_fill_viridis_c() +
  labs(
    x = "Visit Month",
    # y = "Case/Control/Other",
    fill = "Predicted\nCounts"
  ) +
  theme_minimal()

ggsave(
    p_model_heatmap,
    filename = glue("{wkdir}/figures/AIRR/",
    "TRUST4-model-nb-heatmap_{Sys.Date()}.png"),
    width = 12, height = 16
)






model <- glmmTMB::glmmTMB(
    frequency ~ visit_month + raw_total_reads + case_control_other_latest + (1|participant_id),
    data = clone_df,
    family = beta_family()
)

clone_df %>% select(sample_id, participant_id, visit_month, raw_total_reads, case_control_other_latest, frequency, counts, CDR3aa) %>%
    head(20)


model <- glmmTMB::glmmTMB(
    frequency ~ visit_month + raw_total_reads + case_control_other_latest + (1|participant_id),
    data = clone_df,
    family = beta_family()
)

# Try fitting the model with error handling
  model <- mgcv::gam(
    frequency ~ visit_month + raw_total_reads + case_control_other_latest + s(participant_id, bs="re"),
    data = clone_df,
    family = mgcv::betar()
  )
  model <- mgcv::gam(
    frequency ~ visit_month + raw_total_reads + case_control_other_latest,
    data = clone_df,
    family = mgcv::betar()
  )

model %>% summary()


mod_lm = gam(frequency ~  participant_id + raw_total_reads + s(visit_month, bs="cr"),
    data = clone_df
)

summary(mod_lm)

# tic()
# clone_stats_summary <- purr   r::map2(
#   shared_clones$participant_id[1:10],
#   shared_clones$CDR3aa[1:10],
#   ~format_clone_stats(.x, .y),
#   .progress = TRUE
# ) %>% 
#     bind_rows()
# toc()


saveRDS(
    clone_stats_summary,
    glue(
        "{wkdir}/data/interim/airr/",
        "TRUST4-clone-stats-summary_{Sys.Date()}.rds")
)


clone_stats_summary %>% View

format_clone_stats(shared_clones$participant_id[3], shared_clones$CDR3aa[3])


clone_stats <- trust4_reports_meta_clean %>% 
    dplyr::filter(
        participant_id == shared_clones$participant_id[1], 
        CDR3aa == shared_clones$CDR3aa[1]) %>% 
        dplyr::select(sample_id, counts, frequency, CDR3aa, CDR3nt) %>% 
    glimpse()


pid_meta <- key_meta_df %>% filter(participant_id == shared_clones$participant_id[1])

pid_meta %>% 
    left_join(clone_stats) %>% 
    mutate(
        counts = replace_na(counts, 0),
        frequency = replace_na(frequency, 0),
        CDR3aa = replace_na(CDR3aa, CDR3aa[!is.na(CDR3aa)][1])
    ) %>% 
    View
    glimpse





# shared_clones$participant_id %>% table %>% hist(100)


# # Get counts of shared clones per participant-sample
# clone_counts <- trust4_reports_meta_clean %>%
#   # Keep only the shared clones identified above
#   semi_join(shared_clones, by = c("participant_id", "CDR3aa")) %>%
#   # Count occurrences per participant-sample
#   group_by(participant_id, sample_id) %>%
#   summarise(
#     shared_clone_count = n(),
#     .groups = "drop"
#   ) %>%
#   glimpse()





# emerson_files <- list.files(glue("{wkdir}/data/input/temp/emerson_2017_natgen"),
#     full.names = TRUE)

# tst <- read_tsv(emerson_files[1])
# tst %>% glimpse()
# View(tst)


# tst %>% filter(frame_type == "In", !is.na(amino_acid)) %>% glimpse





#------------------------------------------------------------------------------

# Reading in MultiQC data and comparing read stats to flagstat data 
# Also comparing HBS to other cohorts in read depth
library(tidyverse)
library(glue)
library(janitor)
library(strex)

home_dir <- "/resnick/groups/MazmanianLab/jboktor"
pdmbs_dir <- paste0(home_dir, "/PDMBS")
wkdir <- paste0(pdmbs_dir, "/pdairr")
qc_dir <- glue("{wkdir}/data/input/metadata/QC")
v4_1027_dir <- glue("{wkdir}/data/input/metadata/2023_v4release_1027")

sample_inv <- bind_rows(
    read_csv(glue("{v4_1027_dir}/rnaseq_WB-RWTS-VHBS_sample_inventory.csv")),
    read_csv(glue("{v4_1027_dir}/rnaseq_WB-RWTS_sample_inventory.csv"))) %>% 
    glimpse()

amppd_meta_core <- read_csv(
    glue("{v4_1027_dir}/amp_pd_participants.csv")) %>% 
    dplyr::select(participant_id, study) %>% 
    right_join(sample_inv) %>%
    glimpse()

main_cohort_stats <- read_tsv(glue(
    "{qc_dir}/rnaseq-WB-RWTS/multiqc_data/mqc_picard_aligned_reads_1.txt")) %>% 
    clean_names()

hbs_cohort_stats <- read_tsv(glue(
    "{qc_dir}/rnaseq-WB-RWTS-VHBS/multiqc_data/picard_alignment_summary_Aligned_Reads.txt")) %>% 
    clean_names()

picard_stats <- bind_rows(
    main_cohort_stats,
    hbs_cohort_stats) %>% 
    dplyr::rename(sample_id = sample) %>% 
    left_join(amppd_meta_core) %>%
    glimpse()

# ECDF plot of aligned reads by study
p_ecdf_aligned <- picard_stats %>%
    ggplot(aes(x = aligned_reads, color = study)) +
    stat_ecdf() +
    scale_x_log10() +
    labs(
        x = "Aligned Reads (log scale)",
        y = "Cumulative Probability",
        title = "Empirical CDF of Aligned Reads by Study",
        color = NULL
    ) +
    scale_color_jco() +
    theme_bw() +
    theme(legend.position = "bottom")

ggsave(
    filename = glue("{wkdir}/figures/readqc/",
    "picard_aligned_reads_ecdf_by_study_{Sys.Date()}.png"),
    p_ecdf_aligned,
    width = 6, height = 4
)

#------------------------------------------------------------------------------
# Saving sample table for HBS sample TRUST4 workflows
hbs_sample_sheet_trust4 <- read_csv(glue(
    "{wkdir}/data/input/amp-pd-sample-tables/",
    "releases_2023_v4release_1027_rnaseq-WB-RWTS-VHBS_rnaseq_WB-RWTS-VHBS_samples.csv"
    )) %>% 
    dplyr::relocate(sample_id, .before = 1)

# hbs_sample_sheet_trust4 %>% 
#     write_tsv(glue("{wkdir}/data/input/amp-pd-sample-tables/v4_hbs_samples.tsv"))


# After an inital run 1869 jobs succeeded, 279 failed (4 failures and 275 cost threshold exceeded)
# editing sample sheet to annotate failures for a re-run

hbs_sample_sheet_trust4 <- read_tsv(
    glue("{wkdir}/data/input/amp-pd-sample-tables/v4_hbs_samples.tsv")) %>% 
    glimpse()

first_round_results <- read_tsv(glue("{wkdir}/data/input/amp-pd-sample-tables/393ef746-ff5a-46b3-9010-fb742081d5cf.tsv")) %>%
  mutate(
    sample_id = str_extract(inputResolutions, 
        '"inputName":"TRUST4workflow\\.samplename","value":"([^"]+)"') %>%
                str_replace_all(
                    '"inputName":"TRUST4workflow\\.samplename","value":"', '') %>%
                str_replace_all('"', '')) %>%
  select(sample_id, status) %>%
  glimpse()

hbs_sample_sheet_trust4_round1 <- hbs_sample_sheet_trust4 %>% 
    left_join(first_round_results) %>% 
    glimpse()

hbs_sample_sheet_trust4_round1$status %>% table

hbs_sample_sheet_trust4_round1 %>% 
    write_tsv(glue("{wkdir}/data/input/amp-pd-sample-tables/v4_hbs_samples_round2.tsv"))


# JUST KDDING DONT NEED TO DO THIS - results were hiding in 2/3 attempts within a run
# # After the second round - one single sample was additionally run in a third round
# # upon checking available outputs - it appears some samples are missing and will need a 
# # fourth round

# hbs_t4_buckets <- paste0("GCS_TRUST4_HBS", 1:3)
# hbs_t4_buckets_paths <- hbs_t4_buckets %>% 
#   purrr::set_names() %>% 
#   purrr::map(
#     ~readRDS(
#       glue("{wkdir}/data/interim/gcs_locations/2025-06-27_URLs_TRUST4-HBS_{.}.rds")
#     ) %>% 
#       data.frame(path = .)
#   ) %>% 
#   bind_rows(.id = "bucket") %>% 
#   mutate(sample_id = basename(path) %>% str_before_last("_")) %>% 
#   mutate(filetype = path %>% str_after_last("_")) %>%
#   glimpse()

# hbs_t4_buckets_paths %>%
#   group_by(filetype, bucket) %>%
#   dplyr::summarise(n = n())

# successful_runs <- hbs_t4_buckets_paths %>% 
#     filter(filetype == "report.tsv") %>% 
#     pull(sample_id)


# hbs_sample_sheet_trust4_round3 <- hbs_sample_sheet_trust4 %>%
#     mutate(status =  case_when(
#         sample_id %in% successful_runs ~ "Succeeded",
#         TRUE ~ "Run Me"
#     )) %>% 
#     glimpse()

# hbs_sample_sheet_trust4_round3$status %>% table()

# hbs_sample_sheet_trust4_round3 %>% 
#     write_tsv(glue("{wkdir}/data/input/amp-pd-sample-tables/v4_hbs_samples_round4.tsv"))

#------------------------------------------------------------------------------

# Selecting HBS samples with the largest estimated T cell Reads
decon_meta <- readRDS(
    glue(
      "{wkdir}/data/interim/celltype_deconvolution/",
      "deconvoluted_cell_types_metadata_2025-06-23.rds")
)

decon_meta %>% glimpse
picard_stats %>% glimpse

hbs_samples_for_arcasHLA <- picard_stats %>% 
    filter(study == "HBS") %>%
    left_join(decon_meta) %>%
    mutate(tcell_reads = abis_t_cell * aligned_reads) %>% 
    group_by(participant_id) %>% 
    slice_max(order_by = tcell_reads, n = 1) %>% 
    pull(sample_id)

hbs_sample_sheet_arcasHLA <- hbs_sample_sheet_trust4 %>% 
    filter(sample_id %in% hbs_samples_for_arcasHLA)

hbs_sample_sheet_arcasHLA %>% 
    write_tsv(glue("{wkdir}/data/input/amp-pd-sample-tables/v4_hbs_samples_arcasHLA.tsv"))


#------------------------------------------------------------------------------


first_round_results_arcasHLA <- read_tsv(
    glue("{wkdir}/data/input/amp-pd-sample-tables/5bff65b9-e635-4aba-a29f-4ba9052880bc.tsv")
    ) %>%
  mutate(
    sample_id = str_after_first(
      inputResolutions, 
      "star/align-reads/"
    ) %>% str_before_first("\\/")
  ) %>%
  select(sample_id, status) %>%
#   filter(status != "Succeeded") %>%
  glimpse()


first_round_results_arcasHLA$status %>% table()

hbs_sample_sheet_arcasHLA %>% 
    left_join(first_round_results_arcasHLA) %>% 
    glimpse() %>%
    write_tsv(glue("{wkdir}/data/input/amp-pd-sample-tables/v4_hbs_samples_arcasHLA_round2.tsv"))



# gsutil ls -p amp-pd-gcp-joeb gs://fc-secure-e081ea07-d4aa-4d2f-bba9-c652349465e7/submissions/393ef746-ff5a-46b3-9010-fb742081d5cf


#------------------------------------------------------------------------------
# # Testing to see if all samples are downloaded 
v4_1027_dir <- glue("{wkdir}/data/input/metadata/2023_v4release_1027")

sample_inv <- bind_rows(
    read_csv(glue("{v4_1027_dir}/rnaseq_WB-RWTS-VHBS_sample_inventory.csv")),
    read_csv(glue("{v4_1027_dir}/rnaseq_WB-RWTS_sample_inventory.csv"))) %>% 
    glimpse()

hbs_samples_success <- 
    glue("{pdmbs_dir}/workflow/RNASEQ/AIRR/TRUST4_HBS/reports") %>% 
    list.files(full.names = TRUE, pattern = "report.tsv") %>% 
    basename() %>% 
    str_before_last("_")

tcrb_df <- fread(
    glue("{wkdir}/data/interim/airr/tcrb_soNNia_2025-05-31.csv")) %>%
    glimpse()

all_other_samples <- tcrb_df$sample_id %>% unique()

still_missing <- sample_inv %>% 
    filter(sample_id %nin% all_other_samples) %>% 
    filter(sample_id %nin% hbs_samples_success) %>% 
    glimpse()

case_control %>% glimpse()
case_control$diagnosis_at_baseline %>% table()
case_control$case_control_other_latest %>% table()
case_control$case_control_other_at_baseline %>% table()


case_control %>% 
    filter(case_control_other_at_baseline != case_control_other_latest) %>%
    View

longitudinal_data %>%
    left_join(sampdat) %>% 
    select(sample_id, participant_id, study, case_control_other_latest) %>% 
    filter(sample_id %in% still_missing$sample_id) %>% 
    pull(study) %>% table



# intersect(all_other_samples, hbs_samples_success)
sample_sheet_trust4_allothers <- read_csv(glue(
    "{wkdir}/data/input/amp-pd-sample-tables/",
    "releases_2023_v4release_1027_rnaseq-WB-RWTS_rnaseq_WB-RWTS_samples.csv"
    )) %>% 
    dplyr::relocate(sample_id, .before = 1) %>%
    glimpse()

sample_sheet_trust4_allothers %>% 
    filter(sample_id %in% still_missing$sample_id) %>% 
    write_tsv(glue("{wkdir}/data/input/amp-pd-sample-tables/v4_allothers_samples_2025-07-02.tsv"))


#------------------------------------------------------------------------------