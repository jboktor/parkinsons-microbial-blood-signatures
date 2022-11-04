# Alpha Diversity adaptable

##### Alpha Diversity Boxplots Script
######## Load Data & functions
source("src/_load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
source("src/community_composition_funcs.R")
source("src/stats_funcs.R")


alpha_diversity_summary <- function(cohort){

  #' Function to calculate a panel alpha diversity metrics and return plots
  #' Takes as input:
  #'
  #'


  load_alphadiv_colors()
  # x <- c(dat.species, dat.path, dat.ec, dat.KOs,
  #        dat.path.slim, dat.ec.slim, dat.KOs.slim, dat.EGGNOGs.slim, dat.PFAMs.slim)
  # z <- c("Species", "Pathways", "Enzymes", "KOs",
  #        "Pathways.slim", "Enzymes.slim", "KOs.slim", "Eggnogs.slim", "Pfams.slim")
  x <- c(dat.species)
  z <- c("Species")

  # ---------------------------------------------------
  #         Alpha Diversity Plotting Loop
  # --------------------------------------------------
  cnt <- 1
  for (i in x){

    cat("Processing input: ", z[cnt], "\n")
    i <- i %>%
      subset_samples(donor_id %ni% low_qc[[1]])
    print(i)
    cat("\n")
    dat_alpha <- i

    # Run Metadata pre-processing function
    env <- NULL
    env <- process_meta(dat_alpha, cohort = cohort)
    env$description <- factor(env$description, levels=c("PD Patient", "Population Control", "Household Control"))
    env$donor_group <- factor(env$donor_group, levels=c("PC", "PD", "HC"))
    ## Calculate Alpha Diversity Metrics and add cols to df
    env$Observed <- microbiome::alpha(abundances(dat_alpha), 'observed')$observed
    env$Shannon <- microbiome::alpha(abundances(dat_alpha), 'shannon')$diversity_shannon
    env$Evenness <- evenness(abundances(dat_alpha), 'simpson')$simpson
    cat("Alpha Div Calculated \n")

    # Create new paired column with only household pairs and NAs for rest
    env.pairs <- dplyr::filter(env, paired != "No")
    # Plot histograms to get a sense of data distribution
    par(mfrow = c(1, 3))
    hist(env$Observed, main="observed OTUs", xlab="", breaks=10)
    hist(env$Shannon, main="Shannon diversity", xlab="", breaks=10)
    hist(env$Evenness, main="Simpson's evenness", xlab="", breaks=10)


    # ---------------------------------------------------
    #                  Observed Species
    # ---------------------------------------------------

    ### STATS
    observed.PdPC <- lm.PdPc(metadf=env, metric="Observed", cohort = cohort)
    observed.PdHC <- lmm.PdHc(metadf=env.pairs, metric="Observed", cohort = cohort)
    ## Pull p-values
    PC.stats <- summ(observed.PdPC)
    observed.PdPC.pval <- PC.stats$coeftable["descriptionPopulation Control","p"]
    HC.stats <- summ(observed.PdHC)
    observed.PdHC.pval <- HC.stats$coeftable["descriptionHousehold Control","p"]

    p1 <- alpha_div_boxplots(df=env, x=env$donor_group, y=env$Observed,
                             df.pairs=env.pairs, df.pairs.x = env.pairs$donor_group, df.pairs.y=env.pairs$Observed,
                             pairs.column=env.pairs$paired,
                             cols=cols.pdpchc, cols.rim=cols.pdpchc.rim,
                             ylabel = paste0("Observed Counts: ", z[cnt]),
                             PDvPC.stat = observed.PdPC.pval, PDvHC.stat = observed.PdHC.pval)

    # ---------------------------------------------------
    #                  Shannon Diversity
    # ---------------------------------------------------

    ### STATS
    shannon.PdPC <- lm.PdPc(metadf=env, metric="Shannon", cohort = cohort)
    shannon.PdHC <- lmm.PdHc(metadf=env, metric="Shannon", cohort = cohort)
    ## Pull p-values
    PC.stats <- summ(shannon.PdPC)
    shannon.PdPC.pval <- PC.stats$coeftable["descriptionPopulation Control","p"]
    HC.stats <- summ(shannon.PdHC)
    shannon.PdHC.pval <- HC.stats$coeftable["descriptionHousehold Control","p"]

    p2 <- alpha_div_boxplots(df=env, x=env$donor_group, y=env$Shannon,
                             df.pairs=env.pairs, df.pairs.x = env.pairs$donor_group, df.pairs.y=env.pairs$Shannon,
                             pairs.column=env.pairs$paired,
                             cols=cols.pdpchc, cols.rim=cols.pdpchc.rim,
                             ylabel = paste0("Shannon's Diversity: ", z[cnt]),
                             PDvPC.stat = shannon.PdPC.pval, PDvHC.stat = shannon.PdHC.pval)

    # ---------------------------------------------------
    #                  Simpsons Eveness
    # ---------------------------------------------------

    ### STATS
    evenness.PdPC <- lm.PdPc(metadf=env, metric="Evenness", cohort = cohort)
    evenness.PdHC <- lmm.PdHc(metadf=env, metric="Evenness", cohort = cohort)
    ## Pull p-values
    PC.stats <- summ(evenness.PdPC)
    evenness.PdPC.pval <- PC.stats$coeftable["descriptionPopulation Control","p"]
    HC.stats <- summ(evenness.PdHC)
    evenness.PdHC.pval <- HC.stats$coeftable["descriptionHousehold Control","p"]

    p3 <- alpha_div_boxplots(df=env, x=env$donor_group, y=env$Evenness,
                             df.pairs=env.pairs, df.pairs.x = env.pairs$donor_group, df.pairs.y=env.pairs$Evenness,
                             pairs.column=env.pairs$paired,
                             cols=cols.pdpchc, cols.rim=cols.pdpchc.rim,
                             ylabel = paste0("Simpson's Evenness: ", z[cnt]),
                             PDvPC.stat = evenness.PdPC.pval, PDvHC.stat = evenness.PdHC.pval)

    cat("Alpha Div Plotted \n")


    ### MERGE PLOTS ###
    alpha_cow <- cowplot::plot_grid(p1, p2, p3, nrow = 1, align = "v")
    alpha_cow
    print(alpha_cow)
    ggsave(alpha_cow, filename =
             paste0("data/Community_Composition/Alpha_Diversity_Analysis/AlphaDiversity_BoxPlot_",
                    z[cnt], "_", cohort, "_Summary.svg"),
           height = 5, width = 7)

    cnt <- cnt + 1

  }
}
