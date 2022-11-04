
# install.packages("tidyverse",  dependences = TRUE, 
#                  lib = "/central/groups/MazmanianLab/joeB/R-Libraries/")

library(glue)
library(purrr)
library(tibble)
library(magrittr)
library(dplyr)
library(stringr)

# Download and save bbduk stdout QC file
pdmbs_dir = "/central/groups/MazmanianLab/joeB/PDMBS/"
setwd(glue("{pdmbs_dir}parkinsons-microbial-blood-signatures") )

# create list of all file paths
clean_fqs <- glue("{pdmbs_dir}workflow/WGS/clean_fastqs")
readqc_stdout <- list.files(clean_fqs, full.names = T) %>% 
  keep(grepl("bbduk_stdout", .))


readqc_stats <- tibble()

for (f in readqc_stdout){
  stdout <- read.delim(f, header = F, sep = "\t")[13:31,]
  stats_list <- stdout %>% str_split(., pattern = " ")
  
  sample_id <- f %>%
    str_split("/") %>%
    unlist() %>% 
    tail(1) %>%
    str_replace("_bbduk_stdout.txt", "")
  
  sample_stdout <- try(
    readqc_stats %<>% bind_rows(as_tibble_row(
      c(
        "sample_id" = sample_id,
        "input_reads" = stats_list[[2]][1],
        "input_bases" = stats_list[[4]][1],
        "qtrimmed_reads" = stats_list[[6]][1],
        "qtrimmed_bases" = stats_list[[7]][1],
        "low_phred_Q10_reads" = stats_list[[9]][1],
        "low_phred_Q10_bases" = stats_list[[10]][1],
        "low_entropy_reads" = stats_list[[12]][1],
        "low_entropy_bases" = stats_list[[13]][1],
        "removed_reads" = stats_list[[15]][1],
        "removed_bases" = stats_list[[16]][1],
        "final_reads" = stats_list[[18]][1],
        "final_bases" = stats_list[[19]][1]
      )
    ))
  )
  
  if (inherits(sample_stdout, "try-error")) {
    cat("\nERROR: ", sample_id, f, "\n")
    next
  }
}


#_____________________________________________________________________________________
# Visual exploration of data 
library(ggplot2)

df_readqc_stats <-
  readqc_stats %>%
  column_to_rownames(var = "sample_id") %>%s
  mutate_all(as.numeric) %>%
  mutate(
    percent_reads_removed = 100 * ((input_reads - final_reads)/input_reads) ,
    percent_bases_removed = 100 * ((input_bases - final_bases)/input_bases)
  )

df_readqc_stats %>% glimpse()

df_readqc_stats %>% 
  ggplot(aes(input_reads, low_entropy_reads)) +
  geom_point(shape = 21) +
  # scale_x_log10() +
  theme_bw()

df_readqc_stats %>% 
  ggplot(aes(input_bases, percent_bases_removed)) +
  geom_point(shape = 21) +
  scale_x_log10() +
  theme_bw()

df_readqc_stats %>% 
  ggplot(aes(input_reads, percent_reads_removed)) +
  geom_point(shape = 21) +
  scale_x_log10() +
  theme_bw()

df_readqc_stats %>% 
  ggplot(aes(final_reads+ 1)) +
  geom_histogram(bins = 300) +
  geom_vline(xintercept = 1e6) +
  scale_x_log10()



#_________________________________________________________
# install.packages("rslurm")
library(rslurm)

test_func <- function(par_mu, par_sd) {
  samp <- rnorm(10^6, par_mu, par_sd)
  c(s_mu = mean(samp), s_sd = sd(samp))
}
pars <- data.frame(par_mu = 1:100,
                   par_sd = seq(0.1, 1, length.out = 10
                                )
                   )
head(pars, 3)
map2(pars$par_mu, pars$par_sd, test_func) %>% 
  bind_rows()

library(rslurm)
sjob <- slurm_apply(test_func, pars, jobname = 'test_apply',
                    nodes = 2, cpus_per_node = 2, submit = TRUE)



slist.files('_rslurm_test_apply', 'results')

res <- get_slurm_out(sjob, outtype = 'table', wait = FALSE)
head(res, 3)
