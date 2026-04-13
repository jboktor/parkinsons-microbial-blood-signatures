library(data.table)
library(tidyverse)
library(glue)

# This script is supposed to be run within a cloned git repo of pogorely/ALICE .. moved here for record keeping

wkdir <- "/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr"
load("VDJT.rda")
source("ALICE.R")

tcrb_df <- fread(
    glue("{wkdir}/data/interim/airr/tcrb_soNNia_2025-05-31.csv")) %>%
    mutate(sample_id_clean = gsub("-", "_", sample_id)) %>%
    glimpse()
# tcrb_df %>% glimpse()

sids <- tcrb_df$sample_id_clean %>% unique()
clonotype_df <- data.frame(table(tcrb_df$sample_id_clean)) %>% 
    dplyr::rename(sample_id_clean = Var1, count = Freq) %>%
    arrange(desc(count)) %>%
    glimpse()

# top_clonotype <- clonotype_df %>% slice_max(order_by = count, n = 10) %>%
#     pull(sample_id_clean)

tcrb_tst <- tcrb_df %>%
    dplyr::select(
    Read.count = "#count",
    Read.proportion = frequency,
    CDR3.nucleotide.sequence = CDR3nt,
    CDR3.amino.acid.sequence = CDR3aa, 
    bestVGene = v_gene,
    bestJGene = j_gene,
    sample_id_clean) %>%
    filter(bestVGene %in% segments$TRBV$V.alleles) %>%
    filter(bestJGene %in% segments$TRBJ$J.alleles) %>%
    # filter(sample_id_clean %in% top_clonotype) %>%
    glimpse()

tcrb_list <- split(tcrb_tst, tcrb_tst$sample_id_clean) %>%
    purrr::map(~ .x %>% dplyr::select(-sample_id_clean))

# tcrb_list %>% purrr::map(nrow)
# tcrb_list %>% glimpse

tcrb_alice <- ALICE_pipeline(
    DTlist=tcrb_list,
    folder="tcrb_alice_fullrun",
    cores=36,
    iter=100,
    nrec=1e6
    # iter=10,
    # nrec=5e5
    )
print(sapply(tcrb_alice, nrow))

output_dir <- glue(
    "/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr",
    "/data/interim/ALICE"
)
saveRDS(tcrb_alice,
    glue("{output_dir}/TEST_tcrb_iter-100_nrec-1e6_{Sys.Date()}.rds")
    )

#______________________________________________________

# # Tutorial data examples

# vgenes_ref <- segments$TRBV$V.alleles
# jgenes_ref <- segments$TRBJ$J.alleles

# # Process S1
# S1d15<-fread("sample/S1_d15_V9_J2_7.tsv")
# S1d0<-fread("sample/S1_d0_V9_J2_7.tsv")
# nrow(S1d0)
# nrow(S1d15)
# S1<-list(d0=S1d0,d15=S1d15)
# S1_chop <- list(
#     "d02" = S1d0 %>% slice_sample(n=40),
#     "d022" = S1d0 %>% slice_sample(n=400),
#     "d0" = S1d0,
#     "d152" = S1d15 %>% slice_sample(n=40),
#     "d1522" = S1d15 %>% slice_sample(n=400),
#     "d15" = S1d15
#     )

# names(S1_chop)

# S1_alice<-ALICE_pipeline(DTlist=S1,folder="S1_res",cores=1,iter=10,nrec=5e5)
# sapply(S1_alice,nrow)

# S1_alice_chop2 <- ALICE_pipeline(DTlist=S1_chop,folder="S1_res",cores=16,iter=10,nrec=5e5)
# sapply(S1_alice_chop2,nrow)

# S1_alice_olga<-ALICE_pipeline_OLGA(DTlist=S1,cores=16)
# sapply(S1_alice_olga,nrow)

# S1_alice_olga_choptest<-ALICE_pipeline_OLGA(DTlist=S1_chop,cores=16)
# sapply(S1_alice_olga_choptest,nrow)
