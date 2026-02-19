# Adapted function that uses IMGT positions to score TCRs with HLA risk score beta-coefficients
# Able to score using all possible CDR3 phenotypes (stratified and unstratified)
score_tcrs_with_betacoefs <- function(x, phenotype_type, betacoef_list){
    require(dplyr)
    all_imgt_pos <- c(
        "P104","P105","P106","P107","P108","P109","P110","P111",
        "P111.1","P112.2","P112.1",
        "P112","P113","P114","P115","P116","P117","P118")
    align_imgt <- function(CDR3){
        AA <- unlist(strsplit(CDR3, ""));
        N_AA_fow <- nchar(CDR3) %/% 2 + nchar(CDR3) %% 2;
        N_AA_na <- 18 - nchar(CDR3);
        AAmod <- c(
            AA[1:(N_AA_fow)],
            rep("NA",N_AA_na),
            AA[(N_AA_fow + 1):nchar(CDR3)]
        );
        names(AAmod) <- all_imgt_pos
        return(t(AAmod))
    }
    
    len_cdr3 <- nchar(x);
    
    if( len_cdr3 >=12 & len_cdr3 <= 18 ){
        cdr3_aln <- align_imgt(x);
        df <- data.frame(ptarget = cdr3_aln[1,], pos = colnames(cdr3_aln)) %>% 
            mutate(tag = paste0(ptarget, ":", pos))
        
        if (grepl("_strat_", phenotype_type)){
            pos_beta <- subset(betacoef_list[[phenotype_type]], cdr3_len==len_cdr3)
        } else {
            pos_beta <- betacoef_list[[phenotype_type]]
        }
        pos_beta$tag <- paste0(pos_beta$cdr3_aa, ":", pos_beta$cdr3_pos)
        pos_beta <- pos_beta[,c("tag","estimate")]
        
        df <- merge( pos_beta, df, by="tag" );
        if ( nrow(df) == 0 ){
            score = 0
        } else {
            score = sum(df$estimate)
        };
    } else {
        score = 0
    }
    return( score )
}


score_tcrs_with_betacoefs_list <- function(x, phenotype_type, betacoef_list, filename) {
    require(dplyr)
    require(purrr)
    source("/central/groups/MazmanianLab/joeB/PDMBS/parkinsons-microbial-blood-signatures/notebooks/R_scripts/compute_CDR3-risk-score.R")
    
    res <- x %>% 
        purrr::set_names() %>% 
        purrr::map(
            ~ score_tcrs_with_betacoefs(x = ., 
                phenotype_type = phenotype_type,
                betacoef_list = betacoef_list)
        )
    saveRDS(res, filename)
    # return(res)
}
