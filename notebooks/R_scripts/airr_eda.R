# Joe Boktor
# Caltech - Mazmanian Lab
# Feb 2022

source("src/_load_packages.R")
source("src/_misc_functions.R")
devtools::source_url("https://raw.githubusercontent.com/jboktor/RIMA_pipeline/master/src/immune_repertoire/trust4_metric_functions.R")

metadata_categories <- readRDS("data/Metadata/metadata_categories.rds")
static_metdata <- readRDS(file = "data/Metadata/static_metdata.rds")
rnaseq_inv <- read.csv(file = "input_files/2021_v2-5release_0510/rna_sample_inventory.csv", 
                stringsAsFactors = F, header = TRUE)
meta_RNASeq <- rnaseq_inv %>% left_join(static_metdata)
static_metdata %>% dim  #glimpse
rnaseq_inv %>% dim  #glimpse
meta_RNASeq %>% dim  #glimpse


tcr <- readRDS("data/AIRR/TRUST4_TCR.rds")
bcr.light <- readRDS("data/AIRR/TRUST4_BCR_light.rds")
bcr.heavy <- readRDS("data/AIRR/TRUST4_BCR_heavy.rds")

#_______________________________________________________________________________
#                         TRUST4 RIMA AIRR measure extraction
#_______________________________________________________________________________


#BCR clustering 
#Note that not every sample have BCR cluster
# sample_bcr_cluster <- BuildBCRlineage(sampleID = sample_id, Bdata = cdr3.bcr.heavy, start=3, end=10)

bcr_lineage <- tibble()
bcr_clonality <- tibble()
tcr_clonality <- tibble()

for (id in unique(bcr.heavy$sample) ){
  
  BCRLin <- BuildBCRlineage(sampleID = id,
                   Bdata = bcr.heavy,
                   start=3, end=10)
  BCRLin %>% print
  bcr_lineage %<>% bind_rows(BCRLin)
  
  #BCR clonality and entropy
  single_sample_bcr_clonality <- getClonality(id, bcr.heavy, start=3, end=10)
  bcr_clonality %<>% bind_rows(single_sample_bcr_clonality)

  #TCR clonality and entropy
  single_sample_tcr_clonality <- getClonalityTCR(id, tcr)
  tcr_clonality %<>% bind_rows(single_sample_tcr_clonality)

}

bcr_lineage %>% glimpse()
bcr_clonality %>% glimpse()
tcr_clonality %>% glimpse()




#________________________________________________________________________________
# TROUBLE SHOOTING

 BuildBCRlineage(sampleID = id,
                          Bdata = bcr.heavy,
                          start=1, end=10)
# bcr_lineage %<>% bind_rows(BCRLin)

 
 
 
 
 
 
 BuildBCRlineage2 <- 
   function(sampleID, Bdata = BCRdata, start=3, end=10) {
     
   # # Given sample ID, start and end position of complete CDR3,  return all the lineages in the sample
   # sampleID <- "BF-1253-SVM0_5T1"
   # Bdata = bcr.heavy
   # start=3
   # end=10
   # 
   
   tmp.dd.ss = subset(Bdata, sample == sampleID)
   tmp.dd.ss = tmp.dd.ss[!duplicated(tmp.dd.ss[,"CDR3nt"]),]
   if (is.null(dim(tmp.dd.ss)))
     return(NA)
   tmp.comp.vv <- which(tmp.dd.ss[, "is_complete"] == "Y")
   comp.CDR3.ss = data.frame(CDR3aa = tmp.dd.ss[tmp.comp.vv, "CDR3aa"])
   if (length(comp.CDR3.ss) == 0)
     return(NA)
   tmp.tt = table(substr(comp.CDR3.ss$CDR3aa, start, end))
   tmp.tt = sort(tmp.tt, decreasing = T)
   tmp.tt <- tmp.tt[which(nchar(names(tmp.tt))==(end-start+1))]
   tmp.motifs = MergeMotifs(names(tmp.tt))
   count = 0
   BCRlineage = c() ## a list of BCR lineage trees
   kept.motifs = c()
   for (mm in tmp.motifs) {
     
     # # TROUBLE
     # mm <- "ADPPRSGY"
     
     
     mm.list = CreateMotifList(mm)
     tmp.vv.ss = c()
     for (tmp.mm in mm.list) {
       tmp.vv.ss = c(tmp.vv.ss, grep(tmp.mm, tmp.dd.ss$CDR3aa))
     }
     tmp.vv.ss = unique(tmp.vv.ss)
     if (length(tmp.vv.ss) < 2)
       # print("motif search stop")
       next
     SEQs = unique(as.character(tmp.dd.ss[tmp.vv.ss, "CDR3nt"]))
     #SEQs = SEQs$CDR3nt
     tmp.dd0 = tmp.dd.ss[tmp.vv.ss, ]
     setDF(tmp.dd0)   ###format as dataframe
     rownames(tmp.dd0) = tmp.dd0$CDR3nt   ####same cdr3dna, same cdr3aa, different Ig gene and totaldna
     tmp.dd0 = tmp.dd0[SEQs, ]
     if (length(SEQs) < 3)
       next
     MSAalign = msa(DNAStringSet(SEQs), 'ClustalW')
     seqs = as.character(attributes(MSAalign)$unmasked)
     seqs0 = gsub('-', '', seqs)
     tmp.dd0 = tmp.dd0[match(seqs0, SEQs),]
     tmp = MergeSeqs(seqs, tmp.dd0)
     seqs = tmp$SS
     tmp.dd0 = tmp$DD
     if (is.null(dim(tmp.dd0)))
       next
     nn = nrow(tmp.dd0)
     if (nn <= 3)
       next
     sDist = matrix(0, nn, nn)
     for (ii in 1:nn) {
       for (jj in ii:nn) {
         if (jj == ii)
           next
         tmp.dist = SeqDist(seqs[ii], seqs[jj])
         sDist[ii, jj] = sDist[ii, jj] + tmp.dist
       }
     }
     kept.motifs = c(kept.motifs, mm)
     rownames(tmp.dd0) = NULL
     lineage.obj = list(distMat = sDist,
                        Sequences = seqs,
                        data = tmp.dd0)
     BCRlineage = c(BCRlineage, list(lineage.obj))
     count = count + 1
   }
   names(BCRlineage) = kept.motifs
   return(BCRlineage)
 }
 

 #_________________________________________________________________________________
 #_________________________________________________________________________________
 
 
 
 
 
 
 
 
 
 
 


#TCR CPK
tcr_cpk <- aggregate(CDR3aa ~ sample+lib.size, tcr, function(x) length(unique(x))) %>%
  mutate(CPK = signif(CDR3aa/(lib.size/1000),4))
tcr_cpk

bcr.heavy_cpk <- aggregate(CDR3aa ~ sample+lib.size, bcr.heavy, function(x) length(unique(x))) %>%
  mutate(CPK = signif(CDR3aa/(lib.size/1000),4))
bcr.heavy_cpk

bcr.light_cpk <- aggregate(CDR3aa ~ sample+lib.size, bcr.light, function(x) length(unique(x))) %>%
  mutate(CPK = signif(CDR3aa/(lib.size/1000),4))
bcr.light_cpk


bcr_clonality %>% 
  mutate(clonality = as.numeric(clonality)) %>% 
  dplyr::rename(sample_id = sample) %>% 
  left_join(meta_RNASeq) %>% 
  filter(case_control_other_latest != "Other") %>% 
  ggplot(aes(x = case_control_other_latest, y = clonality)) +
  geom_boxplot() +
  geom_point(aes(fill = case_control_other_latest), position = position_jitterdodge()) +
  stat_compare_means() +
  theme_light()

tcr_clonality %>% 
  mutate(clonality = as.numeric(clonality)) %>% 
  dplyr::rename(sample_id = sample) %>% 
  left_join(meta_RNASeq) %>% 
  filter(case_control_other_latest != "Other") %>%
  ggplot(aes(x = case_control_other_latest, y = clonality)) +
  geom_boxplot() +
  geom_point(aes(fill = case_control_other_latest), position = position_jitterdodge()) +
  stat_compare_means() +
  theme_light()

bcr.heavy_cpk %>% 
  dplyr::rename(sample_id = sample) %>% 
  left_join(meta_RNASeq) %>% 
  filter(case_control_other_latest != "Other") %>%
  ggplot(aes(x = case_control_other_latest, y = CPK)) +
  geom_boxplot() +
  geom_point(aes(fill = case_control_other_latest), position = position_jitterdodge()) +
  stat_compare_means() +
  theme_light()
bcr.light_cpk %>% 
  dplyr::rename(sample_id = sample) %>% 
  left_join(meta_RNASeq) %>% 
  filter(case_control_other_latest != "Other") %>%
  ggplot(aes(x = case_control_other_latest, y = CPK)) +
  geom_boxplot() +
  geom_point(aes(fill = case_control_other_latest), position = position_jitterdodge()) +
  stat_compare_means() +
  theme_light()
tcr_cpk %>% 
  dplyr::rename(sample_id = sample) %>% 
  left_join(meta_RNASeq) %>% 
  filter(case_control_other_latest != "Other") %>%
  ggplot(aes(x = case_control_other_latest, y = CPK)) +
  geom_boxplot() +
  geom_point(aes(fill = case_control_other_latest), position = position_jitterdodge()) +
  stat_compare_means() +
  theme_light()



# Ig isotype frequency
ig_freq <- bcr.heavy %>% 
  group_by(sample) %>%
  mutate(est_clonal_exp_norm = frequency/sum(frequency)) %>%   #as.numeric(sample.clones[filename,2])
  dplyr::filter(C != "*" & C !=".") %>%
  group_by(sample, C) %>% 
  dplyr::summarise(Num.Ig = sum(est_clonal_exp_norm))
ig_freq



# Sample and Ig Isotype clustering
ig_freq_matrix <- ig_freq %>% 
  pivot_wider(values_from = "Num.Ig", names_from = "sample") %>% 
  column_to_rownames("C") %>% 
  replace(is.na(.), 0) %>% 
  as.matrix()

og_sample_order <- colnames(ig_freq_matrix)
og_Ig_order <- colnames(t(ig_freq_matrix))

hc <- 
  ig_freq_matrix %>% 
  dist() %>% 
  hclust()
dd_col <- as.dendrogram(hc)
Ig_order <- order.dendrogram(dd_col)

hr <- 
  ig_freq_matrix %>% t() %>% 
  dist() %>%
  hclust()
dd_row <- as.dendrogram(hr)
sample_order <- order.dendrogram(dd_row)



ig_freq %>% 
  dplyr::rename(sample_id = sample) %>% 
  left_join(meta_RNASeq) %>% 
  filter(case_control_other_latest %in% c("Case", "Control")) %>%
  mutate(sample_id = factor(sample_id, levels = og_sample_order[sample_order])#,
         # C = factor(C, levels = og_Ig_order[Ig_order])
         ) %>% 
  ggplot(aes(x = sample_id, y = Num.Ig, fill = C)) +
  geom_bar(stat='identity', width = 1) +
  facet_grid(~case_control_other_latest, scales = "free_x", space = "free") +
  scale_fill_brewer(palette = "Paired") +
  labs(x=NULL, y = "IGH Isotype Frequency") +
  # geom_ysidedensity(aes(x = stat(density), color = C)) +
  # scale_ysidex_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.3))) +
  theme_classic() +
  theme(axis.text.x = element_blank())







# file <- read.table("TRUST4/SRR8281233.sorted.bam.stat.txt",sep = "\t",row.names = 1)
# mapped.reads <- file["reads mapped:","V2"]
# individual.stats <- cbind.data.frame(sample = sample, map.reads = mapped.reads)
# 
# #---------fraction of BCR reads------------------
# ##extract library size
# lib.size <- cdr3.bcr.heavy %>% group_by(sample) %>%
#   dplyr::summarise(lib = mean(lib.size))
# 
# ##combine stats and library size
# bcr.lib.reads <- merge(individual.stats,lib.size,by = "sample") %>% 
#   mutate(Infil = signif(as.numeric(lib)/as.numeric(map.reads),4))
# 
# #------------fraction of TCR reads-----------------
# ##extract library size
# lib.size <- cdr3.tcr %>% group_by(sample) %>%
#   dplyr::summarise(lib = mean(lib.size)) 
# ##combine stats and library size
# tcr.lib.reads <- merge(individual.stats,lib.size,by = "sample") %>% 
#   mutate(Infil = signif(as.numeric(lib)/as.numeric(map.reads),4))
# 
# combined <- rbind(bcr.lib.reads, tcr.lib.reads) %>% mutate(Type = c("BCR", "TCR"))
# combined


# # Somatic hypermutation rate of BCR
# SHM.ratio <- getSHMratio(sample_bcr_cluster)





