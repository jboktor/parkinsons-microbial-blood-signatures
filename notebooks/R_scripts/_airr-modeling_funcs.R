#' Wrapper function to embed AIRR sequences using trill
trill_embed <- function(run_name, model, fasta, output_dir) {
    require(glue)
    wkdir <- glue(
        "/central/groups/MazmanianLab/joeB/PDMBS/",
        "parkinsons-microbial-blood-signatures"
    )
    source(glue("{wkdir}/notebooks/R_scripts/_misc_functions.R"))
    tmp_loc <- tempdir()
    dir.create(tmp_loc, recursive = TRUE)
    cmd <- glue(
        "cd {tmp_loc} &&",
        " mamba run -n trill_v152",
        " trill {run_name} 1 embed {model} {fasta} --avg &&",
        " mv * {output_dir}/"
        )
    shell_do(cmd)
    unlink(tmp_loc)
}

    # tmp_loc <- glue(
    #     "{wkdir}/tmp/",
    #     "{strex::str_after_last(tempdir(), '/')}"
    # )

trill_embed_parallel <- function(run_name, model, fasta, output_dir, ngpu) {
    require(glue)
    wkdir <- glue(
        "/central/groups/MazmanianLab/joeB/PDMBS/",
        "parkinsons-microbial-blood-signatures"
    )
    source(glue("{wkdir}/notebooks/R_scripts/_misc_functions.R"))
    tmp_loc <- glue("{wkdir}{tempdir()}")
    dir.create(tmp_loc, recursive = TRUE)
    cmd <- glue(
        "cd {tmp_loc} &&",
        " mamba run -n trill_v152",
        " trill {run_name} {ngpu} embed {model} {fasta} --avg &&",
        " mv * {output_dir}/"
        )
    shell_do(cmd)
    unlink(tmp_loc)
}

#' Function to split one large table by a column into a list of sub-dataframes
split_train_test <- function(m, split_col, remove_cols) {
    unique_vars <- m[[split_col]] %>% unique()
    m_list <- unique_vars %>%
        purrr::set_names() %>%
        purrr::map(~ dplyr::filter(m, !!sym(split_col) == .) %>%
            dplyr::select(-any_of(remove_cols)) %>%
            drop_na()
            )
    return(m_list)
}

# -----------------------------------------------------------------------------
# GMM-Aligner Functions

#' Function to run ICA on a list of datatables
map_ica <- function(df, n_ica_comp) {
    ica_list <- df %>%
        purrr::map(
            ~ as.data.frame(.) %>% fastICA(n.comp = n_ica_comp)
        )
    return(ica_list)
}

# TODO - cid column is not a generalizable analysis feature
format_ica_to_df <- function(ica_list) {
    purrr::set_names(names(ica_list)) %>%
        purrr::map(~ ica_list[[.]]$S %>% as.data.frame()) %>%
        dplyr::bind_rows(.id = "grouping") %>%
        dplyr::rename_all(~ gsub("V", "ICA_", .)) %>%
        tibble::rownames_to_column("Label") %>%
        dplyr::mutate(
            sample_id = strex::str_before_nth(Label, "_", -2),
            cid = strex::str_before_last(Label, "_") %>%
                strex::str_after_last("_")
        )
}

fit_gmms <- function(ica_list, ref_group, optimal_cluster_n = NULL) {
    if (is.null(optimal_cluster_n)) {
        nrows_full <- nrow(ica_list[[ref_group]]$S)
        rand_rows <- sample(1:nrows_full, nrows_full/10)
        test_set_subset <- ica_list[[ref_group]]$S[rand_rows, ]

        cat(glue(
            "{get_time()} Determining optimal",
            " number of GMM cluster..."
        ), "\n")

        gmms_open_fit <- mclust::Mclust(
            test_set_subset,
            G = NULL,
            modelNames = "VVV"
        )
        optimal_cluster_n <- gmms_open_fit$G
    }

    cat(glue(
        "{get_time()} Fitting GMMs with",
        " {optimal_cluster_n} clusters..."
    ), "\n")

    gmms <- names(ica_list) %>%
        purrr::set_names() %>%
        purrr::map(
            ~ Mclust(ica_list[[.]]$S,
                G = optimal_cluster_n, 
                prior=priorControl(),
                modelNames = "VVV"
            )
        )
    return(gmms)
}

#' Converts matrix of JD values to long format
make_jd_matrix_long <- function(jd_matrix) {
    jd_pairs_long <- tibble::as_tibble(jd_matrix) %>%
        tibble::rownames_to_column("ref_components") %>%
        tidyr::pivot_longer(
            !ref_components,
            names_to = "test_components",
            values_to = "jd"
        ) %>%
        dplyr::mutate(
            ref_components = paste0("Mixture_", ref_components)
        ) %>%
        mutate(test_components = gsub("V", "Mixture_", test_components)) %>%
        # Finding closest component from the TEST set to the REFERENCE set
        dplyr::group_by(test_components) %>%
        dplyr::mutate(
            min_jd = min(jd),
            train_min = dplyr::case_when(
                jd == min_jd ~ "X",
                TRUE ~ ""
            )
        )
    return(jd_pairs_long)
}

#' Function to align test set GMMs to reference set GMMs
align_gmms <- function(jd_matrix, gmms, ref_group) {
    message("Aligning GMMs -- Reference group is: ", ref_group)
    jd_pairs_long <- make_jd_matrix_long(jd_matrix)

    pair_map_df <- jd_pairs_long %>%
        filter(train_min == "X") %>%
        dplyr::select(tidyselect::contains("components"))

    gmm_class_map <- gmms %>%
        purrr::map_dfr(
            ~ .x$classification %>% as.data.frame(),
            .id = "grouping"
        ) %>%
        dplyr::rename("mixture" = ".") %>%
        dplyr::mutate(
            mixture = paste0("Mixture_", mixture),
            # not actually the test group col,
            # temporarily using this to map to ref set
            test_components = mixture
        ) %>%
        tibble::rownames_to_column("Label") %>%
        dplyr::left_join(pair_map_df, by = "test_components") %>%
        dplyr::mutate(
            aligned_mixture = dplyr::case_when(
                .data$grouping == ref_group ~ mixture,
                TRUE ~ ref_components
            )
        ) %>%
        dplyr::select(-c("test_components", "ref_components"))

    res <- list(
        "jd_pairs_long" = jd_pairs_long,
        "pair_map_df" = pair_map_df,
        "gmm_class_map" = gmm_class_map
    )
    return(res)
}

get_time <- function(){
  print(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"))
}


#' GMM aligner function
gmm_aligner <- function(dlist, ref_group, test_group,
                        optimal_cluster_n = 30,
                        threads = 6,
                        n_ica_comp = 30,
                        output_path,
                        verbose = TRUE,
                        overwrite = FALSE,
                        src = wkdir) {
  # Dependencies Check
  require(magrittr)
  require(glue)
  require(tibble)
  require(dplyr)
  require(purrr)
  source(glue("{src}/notebooks/R_scripts/_misc_functions.R"))

  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  if (verbose) {
    start_time <- Sys.time()
    cat(glue("{start_time} - Starting the analysis...\n\n"))
  }

  # Filtering dlist based on ref_group and test_group
  dlist_filtered <- dlist %>%
    purrr::keep(names(.) %in% c(ref_group, test_group))

  step_process <- function(data, step_name, step_function,
                           output_filename, overwrite = FALSE, ...) {
      output_filepath <- file.path(output_path, output_filename)
      if (!file.exists(output_filepath) | overwrite) {
          if (verbose) cat(glue("{Sys.time()} - {step_name}...\n\n"))
          result <- step_function(data, ...)
          saveRDS(result, output_filepath)
      } else {
          if (verbose) {
              cat(glue(
                  "{Sys.time()} - Loading precomputed {step_name}...\n\n"
              ))
          }
          result <- readRDS(output_filepath)
      }
      result
  }

  # PCA for test and train separately
  pca_results <- step_process(
    dlist_filtered, "Performing PCA",
    map_pca, "pca_results_list.rds",
    overwrite = overwrite
    )

  # ICA analysis on PCA results
  ica_results <- step_process(
    pca_results, "Performing ICA analysis",
    map_ica, "ica_results_list.rds",
    n_ica_comp = n_ica_comp,
    overwrite = overwrite
  )

  # Fitting Gaussian Mixture Models
  gmms <- step_process(
    ica_results, "GMM Analysis",
    fit_gmms, "GMMs.rds",
    overwrite = overwrite,
    ref_group = ref_group, 
    optimal_cluster_n = optimal_cluster_n
  )

  # Computing GMM component distances
  jd_matrix <- step_process(
    gmms[[ref_group]], "Computing Pairwise Jeffreys Divergence",
    pairwise_jeffreys_divergence, "jefferys_divergence_GMMs.rds",
    overwrite = overwrite,
    gmm2 = gmms[[test_group]]
  )

  # Aligning GMM components
  gmm_alignment <- step_process(
    jd_matrix, "Aligning GMM Components",
    align_gmms, "GMM-alignment.rds",
    overwrite = overwrite,
    gmms = gmms,
    ref_group = ref_group
  )

  if (verbose) {
    end_time <- Sys.time()
    cat(glue(
      "{end_time} - Analysis completed in",
      " {difftime(end_time, start_time, units = 'mins')} minutes.\n\n"
    ))
  }
}


# # main function
# gmm_aligner <- function(dlist, ref_group, test_group,
#                         optimal_cluster_n = 9,
#                         t4_chains = trust4_reports_chains,
#                         threads = 40,
#                         n_ica_comp = 30) {
#     require(magrittr)
#     require(glue)
#     wkdir <- glue(
#         "/central/groups/MazmanianLab/joeB/PDMBS/",
#         "parkinsons-microbial-blood-signatures"
#     )
#     source(paste0(wkdir, "/notebooks/R_scripts/_misc_functions.R"))

#     cat(glue("{get_time()}  Reference group is: {ref_group}"), "\n")
#     dlist %<>% purrr::keep(names(.) %in% c(ref_group, test_group))

#     cat(glue("{get_time()}  Beginning ICA analysis..."), "\n")
#     ica_list <- map_ica(dlist, n_ica_comp)

#     cat(glue("{get_time()} Formating ICA results..."), "\n")
#     ica_df <- format_ica_to_df(ica_list)

#     cat(glue("{get_time()}  Fitting GMMs..."), "\n")
#     gmms <- fit_gmms(ica_list, ref_group, optimal_cluster_n)

#     cat(glue("{get_time()} Computing Pairwise Jeffreys Divergence..."), "\n")
#     jd_matrix <- pairwise_jeffreys_divergence(
#         gmms[[ref_group]], gmms[[test_group]]
#     )

#     cat(glue("{get_time()} Optimizing TEST set matches to Ref..."), "\n")
#     gmm_alignment <- align_gmms(jd_matrix, gmms, ref_group)

#     # # # Cosine distance
#     # # cos_dist_matrix <- dist_cosine(as.matrix(dplyr::bind_rows(dlist)))
#     # cos_dist <- coop::cosine(
#     #     dplyr::bind_rows(dlist) %>%
#     #         sample_n(1000) %>%
#     #         as.matrix() %>%
#     #         t()
#     # )

#     # dist_data <- dplyr::bind_rows(dlist) %>%
#     #     sample_n(1000) %>%
#     #     as.matrix() %>%
#     #     t()

#     # euc_dist <- distances::distances(
#     #     dist_data,
#     #     id_variable = rownames(dist_data)
#     # )

#     cat(glue("{get_time()} Computing PCA..."), "\n")
#     d_pca <- prcomp(dplyr::bind_rows(dlist), rank. = 50)
#     pca_df <- d_pca[["x"]] %>%
#         as.data.frame() %>%
#         tibble::rownames_to_column("Label")

#     # cat(glue("{get_time()} Computing UMAP..."), "\n")
#     # d_umap <- uwot::umap(
#     #     dplyr::bind_rows(dlist),
#     #     n_components = 3,
#     #     n_threads = threads
#     # )

#     # umap_df <- d_umap %>%
#     #     as.data.frame() %>%
#     #     dplyr::rename_all(~ gsub("V", "UMAP_", .)) %>%
#     #     tibble::rownames_to_column("Label")

#     cat(glue("{get_time()} Compiling results..."), "\n")
#     # gmm_meta_df <- gmm_alignment$gmm_class_map %>%
#     #     dplyr::left_join(ica_df) %>%
#     #     dplyr::left_join(t4_chains) %>%
#     #     dplyr::left_join(pca_df) %>%
#     #     # dplyr::left_join(umap_df) %>%
#     #     pillar::glimpse()

#     # Prepare the list of results
#     res <- list(
#         "ica_list" = ica_list,
#         "gmms" = gmms,
#         "jd_matrix" = jd_matrix,
#         "jd_pairs_long" = gmm_alignment$jd_pairs_long,
#         "pair_map_df" = gmm_alignment$pair_map_df,
#         "gmm_class_map" = gmm_alignment$gmm_class_map,
#         "pca" = d_pca #,
#         # "umap" = d_umap,
#         # "compiled" = gmm_meta_df
#     )
#     # Return the list of results
#     return(res)
# }


#' Function to calculate Jeffrey's Divergence between two Gaussian distributions
jeffreys_divergence <- function(mu1, sigma1, mu2, sigma2) {
  det
  sigma1_inv <- solve(sigma1)
  sigma2_inv <- solve(sigma2)
  trace_term <- psych::tr(sigma1_inv %*% sigma2 + sigma2_inv %*% sigma1)
  mu_delta <- mu1 - mu2
  quadratic_term <- t(mu_delta) %*% (sigma1_inv + sigma2_inv) %*% mu_delta
  jd <- 0.5 * (trace_term + quadratic_term - length(mu1) * 2)
  return(jd)
}

#' Function to align the distributions of two GMMs 
#' by calculating the Jeffrey's Divergence between each pair of components
#' and returning a matrix of the divergences
pairwise_jeffreys_divergence <- function(gmm1, gmm2) {
    n_mixtures <- gmm1$G
    divergence_matrix <- matrix(, nrow = n_mixtures, ncol = n_mixtures)
    for (i in 1:n_mixtures) {
        for (j in 1:n_mixtures) {
            mu1 <- gmm1$parameters$mean[, i]
            mu2 <- gmm2$parameters$mean[, j]
            sigma1 <- gmm1$parameters$variance$sigma[, , i]
            sigma2 <- gmm2$parameters$variance$sigma[, , j]
            divergence_matrix[i, j] <-
                jeffreys_divergence(mu1, sigma1, mu2, sigma2)
        }
    }
    return(divergence_matrix)
}

map_distances <- function(dlist) {
  require(distances)
  require(magrittr)
  require(tibble)
  require(dplyr)
  require(purrr)

  dist_list <- purrr::map(
    dlist,
    ~ {
      dist_data <- as.matrix(.) %>% t()
      distances::distances(
        dist_data,
        id_variable = rownames(dist_data)
      )
    }
  )
  return(dist_list)
}

map_pca <- function(dist_list) {
  require(magrittr)
  require(dplyr)
  require(purrr)
  require(stats)

  pca_list <- purrr::map(
    dist_list,
    ~ prcomp(., rank. = 50)[["x"]]
  )
  return(pca_list)
}


#------------------------------------------------------------
#### Analysis of GMM Parameters 

# Calculating delta Mu (average change in state)
gmm_mean_l2_norm <- function(gmm1, gmm2, i, j) {
    delta_mu <- gmm1$parameters$mean[ ,i] - gmm2$parameters$mean[ ,j]
    norm(delta_mu, type = "2")
}

gmm_abs_weight_diff <- function(gmm1, gmm2, i, j) {
    abs(
        gmm1$parameters$pro[i] - gmm2$parameters$pro[j]
    )
}

gmm_forstners_distance <- function(gmm1, gmm2, i, j) {
    sigma_i <- gmm1$parameters$variance$sigma[,,i]
    sigma_j <- gmm2$parameters$variance$sigma[,,j]

  # Check if matrices are square and of the same dimension
  if (!is.matrix(sigma_i) || !is.matrix(sigma_j) || 
      nrow(sigma_i) != ncol(sigma_i) || nrow(sigma_j) != ncol(sigma_j) || 
      nrow(sigma_i) != nrow(sigma_j)) {
    stop("Both inputs must be square matrices of the same dimensions.")
  }
  
  # Compute the generalized eigenvalues for sigma_i and sigma_j
  gen_eigen <- function(A, B) {
    # Solve the generalized eigenvalue problem A * v = lambda * B * v
    ev <- eigen(solve(B, A))
    return(ev$values)
  }

  lambda <- gen_eigen(sigma_i, sigma_j)

 # Checking if any eigenvalues are negative and imputing if so
  min_eigenvalue <- min(lambda)
  if (min_eigenvalue <= 0) {
    # Add a small constant to all eigenvalues to make them positive
    shift_amount <- abs(min_eigenvalue) + .Machine$double.eps
    lambda <- lambda + shift_amount
  }
  # Compute the Forstner distance
  forstners_distance <- sqrt(sum(log(lambda)^2))
  return(forstners_distance)
}

#' Collecting stats on aligned distributions
calculate_gmm_param_analysis <- function(pair_map_df, gmm_ref, gmm_test) {
  pair_map_df %>%
    mutate_at(
      vars(ref_components, test_components),
      ~ gsub("Mixture_", "", .) %>% as.numeric()
    ) %>%
    mutate(
      delta_mu = purrr::map2_dbl(
        .x = ref_components, .y = test_components,
        ~ gmm_mean_l2_norm(
          gmm1 = gmm_ref, gmm2 = gmm_test,
          i = .x, j = .y
        )
      ),
      delta_w = purrr::map2_dbl(
        .x = ref_components, .y = test_components,
        ~ gmm_abs_weight_diff(
          gmm1 = gmm_ref, gmm2 = gmm_test,
          i = .x, j = .y
        )
      ),
      delta_sigma = purrr::map2_dbl(
        .x = ref_components, .y = test_components,
        ~ gmm_forstners_distance(
          gmm1 = gmm_ref, gmm2 = gmm_test,
          i = .x, j = .y
        )
      )
    )
}



# Plotting functions ----------------------------------------------------------

#' Function to plot heatmap of ICA components and receptor chains
plot_ica_heatmap <- function(fast_ica_obj) {
    require(ggplot2)
    S <- fast_ica_obj$S
    S_trans <- t(S)
    rownames(S_trans) <- paste0("V", 1:ncol(S))
    receptor_order <- seriate_matrix_rows(S, seriate_method = "MDS_angle")
    ica_order <- seriate_matrix_rows(S_trans, seriate_method = "MDS_angle")

    ind_comp <- as.data.frame(S) %>%
        tibble::rownames_to_column("Label") %>%
        as.data.table() %>%
        tidyr::pivot_longer(
            !Label,
            names_to = "component",
            values_to = "value") %>%
        dplyr::mutate(
            particpant_id = strex::str_after_nth(Label, "-", 2) %>%
                strex::str_before_last("_"),
            case_control_other_latest = str_after_nth(Label, "-", 2) %>%
                strex::str_after_last("_")) %>%
        dplyr::mutate(
            Label = factor(Label, levels = receptor_order, ordered = TRUE),
            component = factor(component, levels = ica_order, ordered = TRUE))

    p_ica_heatmap <- ind_comp %>%
        ggplot(aes(x = component, y = Label, fill = value)) +
        geom_tile(aes(color = value)) +
        theme_minimal() +
        theme_set(theme_bw()) +
        scale_fill_viridis_c(option = "mako") +
        scale_color_viridis_c(option = "mako") +
        guides(color = "none") +
        labs(x = "Components", y = "Receptor", fill = NULL) +
        theme(
            axis.text.x = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            legend.position = "top"
            )
    return(p_ica_heatmap)
}

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

plotly_umap <- function(df, color_by) {
    plot_ly(
        data = df,
        x = ~UMAP_1,
        y = ~UMAP_2,
        z = ~UMAP_3,
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
            xaxis = list(title = "UMAP 1"),
            yaxis = list(title = "UMAP 2"),
            zaxis = list(title = "UMAP 3")
        ))
}

plotly_umap2d <- function(df, color_by) {
    plot_ly(
        data = df,
        x = ~UMAP_1,
        y = ~UMAP_2,
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
            xaxis = list(title = "UMAP 1"),
            yaxis = list(title = "UMAP 2")
        ))
}

plot_gmm_param <- function(data, component_column, delta_mu_column) {
  data %>%
    # mutate(!!component_column := as.character(.[[component_column]])) %>%
    ggplot(aes_string(x = component_column, y = delta_mu_column)) +
    geom_segment(aes_string(
      x = component_column, xend = component_column,
      y = 0, yend = delta_mu_column
    ), color = "grey") +
    geom_point(color = "orange", size = 2) +
    coord_flip() +
    labs(x = "Aligned_Mixture") +
    theme_light()
}

#' ggplot function for 2D UMAP viz
ggumap <- function(df, color_by = "grouping") {
    ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = .data[[color_by]] )) +
        theme_light() +
        geom_point(alpha = 0.7, size = 0.5) +
        ggsci::scale_color_npg()
}

#' ggplot function for 2D UMAP viz
ggpca <- function(df, color_by = "grouping", x = "PC1", y = "PC2") {
    # require(Polychrome)
    # colpal = glasbey.colors(length(unique(df[[color_by]])))
    ggplot(df, aes(x = .data[[x]], y = .data[[y]], color = .data[[color_by]])) +
        theme_light() +
        geom_point(alpha = 0.7, size = 0.5) +
        guides(color = guide_legend(override.aes = list(size = 3)))
        # scale_color_manual(values = unname(colpal))
        # ggsci::scale_color_d3()
}


#### Optimal Transport Analysis -----------------------------------------------
# Splitting the data into a list of two dataframes S1 and S2
filter_and_split_samples <- function(df, s1, s2,
                                     cols_to_filter,
                                     split_var) {
    require(dplyr)
    df %>%
        filter(.data[[split_var]] %in% c(s1, s2)) %>%
        split_train_test(
            split_col = split_var,
            remove_cols = cols_to_filter
        )
}

# Function to randomly select n receptors from each dataframe in list
subsample_dfs <- function(split_dfs, n) {
    purrr::map(
        split_dfs,
        ~ dplyr::sample_n(.x, n) %>% tibble::tibble()
    )
}

# Function to calculate the Wasserstein distance between two dataframes
run_earth_movers_subsamples <- function(pca_df_ot, s1, s2, n,
                                        iter = 50,
                                        split_var = "sample_id",
                                        cols_to_filter = c(
                                            "Label",
                                            "sample_id",
                                            "cid"
                                        )) {
    sample_splits <- filter_and_split_samples(
        df = pca_df_ot,
        split_var = split_var,
        s1 = s1,
        s2 = s2,
        cols_to_filter = cols_to_filter
    )
    if (s1 == s2) {
        sample_splits <- c(sample_splits, sample_splits)
    }
    res <- c()
    for (i in 1:iter) {
        pp_pair <- subsample_dfs(sample_splits, n = n) %>%
            purrr::map(pp)
        res <- c(res, transport::wasserstein(
            pp_pair[[1]], pp_pair[[2]],
            p = 1, prob = TRUE
        ))
    }
    res_stats <- tibble::as_tibble_row(
        list(
            id1 = s1,
            id2 = s2,
            mean = mean(res),
            sd = sd(res),
            median = median(res),
            min = min(res),
            max = max(res),
            full_emd = list(res)
        ),
        .name_repair = "minimal"
    )
    return(res_stats)
}


#------------------------------------------------------------
# Thimble/Stitchr processing functions

#' Reading in Thimble outputs
load_thimble_results <- function(file, file_prefix = "trust4_") {
    filebase <- basename(file)
    c_reg <- filebase %>%
        strex::str_before_nth("_", 2) %>%
        gsub(file_prefix, "", .) %>%
        tolower()

    thimble_res <- data.table::fread(file, header = TRUE)
    formatted_df <- thimble_res %>%
        janitor::clean_names() %>%
        dplyr::select(tcr_name, contains(c_reg), warnings_errors) %>%
        dplyr::select(-matches("(_leader|_prime_seq|Link)")) %>%
        dplyr::mutate(constant = c_reg) %>%
        dplyr::rename_all(~ gsub(glue("{c_reg}_"), "", .)) %>%
        dplyr::rename_all(~ gsub(c_reg, "", .))
    return(formatted_df)
}


#------------------------------------------------------------
# Log-likelihood Ratio (LLR) functions

# Adjusted function to handle VVV covariance model
likelihood_gmm_mclust_vvv <- function(x, gmm) {
  # Mixture proportions
  pi <- gmm$parameters$pro
  # Means for each Gaussian component
  mu <- gmm$parameters$mean
  # number of mixture components
  k <- ncol(mu)
  # Number of data points
  n <- ncol(x)
  # Covariance matrices for each component
  sigma <- gmm$parameters$variance$sigma

  # Calculate the likelihood of the data points for each GMM component
  likelihood <- numeric(n)
  component_densities <- numeric(k)
  for (j in 1:k) {
    component_densities[j] <- mclust::dmvnorm(
      data = x,
      mean = mu[, j, drop = FALSE],
      sigma = sigma[, , j]
    ) * pi[j]
  }
  likelihood <- sum(component_densities)
  return(likelihood)
}

#' Calculate the log-likelihood ratio for a sample (vector)
#' between two pretrained Gaussian Mixture Models
#' @param x A numeric vector of the sample
#' @param gmm_obj A list of two mclust GMM objects
#' @return The log-likelihood ratio
compute_gmm_llr <- function(x_ref, x_test, gmm_ref, gmm_test) {
  sample_n <- nrow(x_ref)
  lref <- 1:sample_n %>%
    purrr::set_names(rownames(x_ref[., ])) %>%
    purrr::map_dfr(~ tibble(
      feature = names(.),
      LLR = likelihood_gmm_mclust_vvv(
        x = matrix(x_ref[.x, ], nrow = 1),
        gmm = gmm_ref
      )
    ), .id = "feature")

  ltest <- 1:sample_n %>%
    purrr::set_names(rownames(x_test[., ])) %>%
    purrr::map_dfr(~ tibble(
      feature = names(.),
      LLR = likelihood_gmm_mclust_vvv(
        x = matrix(x_test[.x, ], nrow = 1),
        gmm = gmm_test
      )
    ), .id = "feature")

  log_likelihoods <- log(lref$LLR / ltest$LLR)
  llr <- sum(log_likelihoods) / sample_n
  return(llr)
}

