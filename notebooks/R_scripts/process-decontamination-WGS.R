



for (refDB in refDBs) {
    for (level in taxa_levels) {
        message(glue("Processing {refDB} {level}"))
        ps <- readRDS(
            glue(
                "{wkdir}/data/processed/phyloseq_objects/raw/",
                "2023-07-14_{seq_method}_{refDB}_{level}_phyloseq.rds"
            )
        )
        # add sequencing metadata to phyloseq object
        meta <- ps %>% meta()
        meta %<>%
            dplyr::left_join(meta_wide, by = "participant_id") %>%
            mutate_if(is.logical, ~ replace_na(., FALSE))
        rownames(meta) <- meta$participant_id
        sample_data(ps) <- sample_data(meta)

        # Prevalence filter
        ps_groups <- c(
            "Case" = phyloseq::subset_samples(
                ps,
                case_control_other_latest == "Case"
            ),
            "Control" = phyloseq::subset_samples(
                ps,
                case_control_other_latest == "Control"
            ),
            "Other" = phyloseq::subset_samples(
                ps,
                case_control_other_latest == "Other"
            )
        )

        # all possible combos of an x var and a y var
        pair_combos <- expand_grid(
            ps_list = names(ps_groups),
            var_list = batch_meta
        ) %>%
            mutate(batch_id = glue("{ps_list}_{var_list}"))
        ps_list <- pair_combos$ps_list
        var_list <- pair_combos$var_list
        names(ps_list) <- pair_combos$batch_id

        tic()
        future::plan("multisession", workers = 14)
        prevalence_stats <-
        furrr::future_map2(
            ps_list,
            var_list,
            ~ calculate_batch_prevalence(ps_groups[[.x]], .y)
        ) %>%
        bind_rows(.id = "grouping_meta")
        toc()

        # Select features with a max of at least 25% prevalence and
        # a minimum of < 2x for other batches
        prevalence_stats_proc <- prevalence_stats %>%
            filter(batch_subset != "FALSE") %>%
            mutate(batch_subset = case_when(
                grepl("dt_|pu_", batch_variable) ~ batch_variable,
                TRUE ~ batch_subset
            ), batch_variable = case_when(
                grepl("dt_", batch_variable) ~ "date",
                grepl("pu_", batch_variable) ~ "flowcell_barcode",
                TRUE ~ batch_variable
            )) %>%
            mutate(case_control_other_latest = strex::str_before_first(
                grouping_meta, "_"
            ))

        prevalence_stat_summary <- prevalence_stats_proc %>%
            # remove batches with fewer than 25 samples
            filter(sample_n > 25) %>%
            # check that all batches have at least a min and a max
            get_dupes(ind, batch_variable, case_control_other_latest) %>%
            filter(dupe_count > 1) %>%
            # summarize the min and max prevalence for each microbe
            group_by(ind, batch_variable, case_control_other_latest) %>%
            dplyr::summarise(
                min = min(values), max = max(values),
                prevalence_delta = max - min
            )

        # flag microbes that are present in at least 50% of samples in at least
        # one batch and are at least 2x more abundant in
        # that batch than in any other
        blacklist_prev <- prevalence_stat_summary %>%
            filter(max > 0.5, min <= max / 2) %>%
            arrange(desc(prevalence_delta)) %>%
            pull(ind) %>%
            unique()

        saveRDS(
            prevalence_stats_proc,
            glue(
                "{wkdir}/data/interim/decontamination/",
                "{Sys.Date()}_{seq_method}_{refDB}_{level}_prevalence_stats.rds"
            )
        )
        saveRDS(
            blacklist_prev,
            glue(
                "{wkdir}/data/interim/decontamination/",
                "{Sys.Date()}_{seq_method}_{refDB}_{level}",
                "_prevalence_blacklist.rds"
            )
        )

        plot_df <- prevalence_stat_summary %>%
            mutate(blacklist = ind %in% blacklist_prev) %>%
            arrange(blacklist, desc(prevalence_delta)) %>%
            mutate(ind = factor(ind))

        p1 <- plot_df %>%
            ggplot(aes(x = ind, y = prevalence_delta, color = blacklist)) +
            geom_point(alpha = 0.7) +
            scale_color_d3() +
            facet_grid(batch_variable ~ case_control_other_latest) +
            theme_bw() +
            labs(x = "microbe", y = "max prevalence - min prevalence") +
                theme(axis.text.x = element_blank())

        # How does blacklist detection compare across groups
        p2 <- plot_df %>%
            ggplot(aes(
            x = case_control_other_latest,
            y = prevalence_delta,
            color = blacklist
            )) +
            geom_point(alpha = 0.3) +
            geom_line(aes(group = ind), alpha = 0.3) +
            facet_wrap(~batch_variable, ncol = 1) +
            scale_color_d3() +
            labs(x = NULL, y = NULL) +
            theme_bw() +
            theme(legend.position='none')

        final_plot <- p1 %>% insert_right(p2, width = 0.5)
        ggsave(
            glue(
                "{wkdir}/figures/decontamination/",
                "{Sys.Date()}_{refDB}_{seq_method}_{level}",
                "_prevalence_filter_EDA.png"
            ),
            final_plot, width = 12, height = 8
            )

        # Read count filter ----
        blacklist_readcounts <-
            filter_taxa(ps, function(x) max(x) < 100, TRUE) %>%
            taxa()
        saveRDS(
            blacklist_readcounts,
            glue(
                "{wkdir}/data/interim/decontamination/",
                "2023-01-09_{seq_method}_{refDB}_{level}",
                "_blacklist_readcounts.rds"
            )
        )
        final_blacklist <- c(blacklist_prev, blacklist_readcounts) %>%
            unique()
        all_taxa <- taxa(ps)
        ps_trim <- phyloseq::prune_taxa(
            all_taxa[all_taxa %nin% final_blacklist], ps
        )
        saveRDS(
            ps_trim,
            glue(
                "{physeq_decon_dir}/",
                "{Sys.Date()}_{seq_method}_{refDB}_{level}_phyloseq.rds"
            )
        )
    }
}
