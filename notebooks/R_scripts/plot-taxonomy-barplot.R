# Joe Boktor
# Caltech - Mazmanian Lab

source("src/_load_packages.R")
source("src/_plot-functions.R")
library(fantaxtic)
sample_info <- readRDS("data/Metadata/static_metdata.rds")

# refDB <- "WoL"
# ps_obj <- readRDS("data/Phyloseq_Objects/WoL/Species_counts.rds")
# taxa_names(ps_obj) %<>% str_remove("s__")
# 
# WoL_taxonomy <- 
#   read_tsv(file = "input_files/taxonomy/WoL/NCBI/ranks.tsv")
# WoL_taxonomy %<>% 
#   select(-genome) %>% 
#   filter(species %in% (ps_obj %>% taxa_names())) %>% 
#   distinct() %>% 
#   as.matrix()
# 
# rownames(WoL_taxonomy) <- WoL_taxonomy[,7]
# tax_table(ps_obj) <- WoL_taxonomy
# sample_data(ps_obj)$reads <- ps_obj %>% sample_sums() %>% as.numeric()



UHGG_taxonomy <- 
  read_tsv(file = "input_files/taxonomy/UHGG/genomes-nr_metadata.tsv")

ps_obj <- readRDS("data/Phyloseq_Objects/UHGG/Species_counts.rds")

UHGG_tax_df <- UHGG_taxonomy %>% 
  select(MGnify_accession, Lineage) %>% 
  separate(Lineage, c("domain", "phylum", "class", "order", "family", "genus", "species"), ";") %>% 
  distinct() %>% 
  mutate(genus = if_else(genus == "g__", paste0(genus, MGnify_accession), genus),
         species = if_else(species == "s__", paste0(species, MGnify_accession), species)) %>% 
  select(-MGnify_accession) %>% 
  distinct() %>% 
  as.matrix()

rownames(UHGG_tax_df) <- UHGG_tax_df[,7]
tax_table(ps_obj) <- UHGG_tax_df
sample_data(ps_obj)$reads <- ps_obj %>% sample_sums() %>% as.numeric()


dat.obj_r3 <- ps_obj %>% aggregate_taxa('order')

#______________________________________________________________________________

dat.top.N <- dat.obj_r3 %>% 
  core(detection = 0, prevalence = 1/1000) %>% 
  subset_samples(reads > 10000) %>% 
  get_top_taxa(n=30, relative = TRUE, discard_other = F, other_label = "Other")


phylum_cols <-
  c("#7fc97f",
   "#beaed4",
   "#fdc086",
   "#ffff99",
   "#386cb0",
   "#f0027f",
   "#e41a1c",
    "#377eb8",
   "#7fc97f",
   "#beaed4")

barplt2 <-
  dat.top.N %>%
  subset_samples(case_control_other_latest != "Other") %>%
  fantaxtic_bar(
    facet_type = "grid",
    color_by = "phylum",
    label_by = "order",
    other_label = "Other",
    facet_by = "case_control_other_latest",
    facet_cols = 2,
    order_alg = "hclust",
    palette = phylum_cols
  ) +
  labs(y = "Relative Abundance") +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())
barplt2


save_me_cleanly(barplt2,
                filename =
                  paste0("figures/community_composition/taxa_barplots/taxa_barplot_top10_classesbyphy_case-vs-control",
                                  refDB, "_", Sys.Date()),
                plot_w = 5, plot_h = 2, leg_w = 4, leg_h = 5,
                res = 1200, filetype = ".png")





# Plot all Samples
barplt1 <- 
  dat.top.N %>% 
  subset_samples(case_control_other_latest != "Other") %>% 
  fantaxtic_bar(
    color_by = "phylum",
    label_by = "order",
    other_label = "Other",
    facet_by = "case_control_other_latest",
    grid_by = "study",
    facet_cols = 2,
    order_alg = "hclust",
    base_color = "#5b9bd5",
  ) +
  labs(y = "Relative Abundance") +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())
barplt1

                                                                  save_me_cleanly(
  barplt1,
  filename = paste0(
    "figures/community_composition/taxa_barplots/taxa_barplot_top50_orders_",
    refDB,
    "_",
    Sys.Date()
  ),
  plot_w = 12,
  plot_h = 7,
  leg_w = 4,
  leg_h = 20,
  res = 1200
)


# dat.top.50 %>%
#   subset_samples(case_control_other_latest != "Other") %>%
#   fantaxtic_bar(
#     color_by = "Rank2",
#     label_by = "Rank6",
#     other_label = "Other",
#     facet_by = "upd2hy_hoehn_and_yahr_stage",
#     # grid_by = "study",
#     # facet_cols = 2,
#     order_alg = "hclust",
#     base_color = "#5b9bd5"
#   ) +
#   labs(y = "Relative Abundance") +
#   theme(axis.text.x = element_blank())



dat.obj.Viruses <- dat.obj_r6 %>% 
  subset_taxa(Rank6 != "g__") %>%  # Order level
  subset_taxa(Rank1 == "k__Viruses")
dat.obj.Viruses <- prune_samples(sample_sums(dat.obj.Viruses)>=5, dat.obj.Viruses)
dat.obj.Viruses <- dat.obj.Viruses %>% subset_samples(participant_id %in% sample_names(dat.quality))

# tst2 <- tax_table(dat.obj_r4) %>% as.data.frame()

# Create Metadata Column for study x Donor Group
sample_data(dat.obj.Viruses)$study_case_control_other_latest <-
  paste(sample_data(dat.obj.Viruses)$study, sample_data(dat.obj.Viruses)$case_control_other_latest)
# Abundance filter for top N features
dat.top.50.v <- dat.obj.Viruses %>% 
  get_top_taxa(n=500, relative = TRUE, discard_other = F, other_label = "Other")

# 
barplt_viruses <-
  dat.top.50.v %>%
  subset_samples(case_control_other_latest != "Other") %>%
  fantaxtic_bar(
    color_by = "Rank3",
    label_by = "Rank5",
    other_label = "Other",
    facet_by = "case_control_other_latest",
    # grid_by = "study",
    facet_cols = 2,
    order_alg = "hclust",
    base_color = "#5b9bd5"
  ) +
  labs(y = "Relative Abundance") +
  theme(axis.text.x = element_blank())
barplt_viruses

# save_me_cleanly(barplt_viruses,
#                 filename = paste0("figures/community_composition/taxa_barplots/taxa_barplot_top50_orders_case-vs-control_Viruses",
#                                   refDB, "_", Sys.Date()),
#                 plot_w = 12, plot_h = 3.5, leg_w = 4, leg_h = 20, res = 1200)




# biom_file <- paste0("input_files/biom_files/", refDB, "_kraken2_2021-09-04.biom")
# biom_otu_tax <- import_biom(biom_file)
# sample_names(biom_otu_tax) <- gsub(paste0("__report_", refDB), "", sample_names(biom_otu_tax))
