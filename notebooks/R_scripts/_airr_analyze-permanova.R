
run_permanova <- function(dmat_path, meta_path, model_formula_string, outpath,
                          n_perms = 999, n_cores = 1) {
    # distance_matrix is used in the model formula
    distance_matrix <- readRDS(dmat_path)
    metadata <- readRDS(meta_path)
    perm_result <- vegan::adonis2(
        formula = as.formula(model_formula_string),
        data = metadata,
        permutations = n_perms,
        parallel = n_cores,
        na.action = "na.omit"
    )
    saveRDS(perm_result, outpath)
}
