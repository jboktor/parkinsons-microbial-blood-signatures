library(dplyr)
library(magrittr)
library(future)
library(batchtools)
library(future.batchtools)
library(readr)
library(listenv)
library(glue)
library(logr)
library(logger)

pdmbs_dir <<- "/central/groups/MazmanianLab/joeB/PDMBS"
wkdir <<- glue("{pdmbs_dir}/parkinsons-microbial-blood-signatures")
source(glue("{wkdir}/notebooks/R_scripts/_misc_functions.R"))
setwd(wkdir)

kraken2_alignment <- function(sampleID, nworkers = 2) {
  library(logger)
  log_threshold(TRACE)
  refdbs_kraken <- list(
    "RefSeqPlusPF" =
      "/central/groups/MazmanianLab/joeB/Downloads/refseq_pluspf_v4/",
    "UHGG" =
      "/central/groups/MazmanianLab/joeB/Downloads/uhgg_kraken2-db/",
    "WoL" =
      "/central/groups/MazmanianLab/joeB/WebOfLife/databases/kraken2/"
  )
  reads_dir <-
  "/central/groups/MazmanianLab/joeB/PDMBS/workflow/WGS/clean_fastqs/"
  temp_reads_dir <-
    "/central/scratch/jbok/krakenScratch/"
  output_dir <-
    "/central/groups/MazmanianLab/joeB/PDMBS/workflow/WGS/results/kraken2/"
  quality_reads <-
    glue("{temp_reads_dir}{sampleID}.fq.gz")
  
  log_info(glue("Clean Reads Dir: {reads_dir}"))
  log_info(glue("Temp Reads Dir: {temp_reads_dir}"))
  log_info(glue("Output Dir: {output_dir}"))
  log_info(glue("Quality Reads: {quality_reads}"))

  for (refDB in names(refdbs_kraken)) {
    classified_seqs_out <-
      glue("{output_dir}{refDB}_mapped/{sampleID}__classified_{refDB}.fastq")
    report_out <-
      glue("{output_dir}{refDB}_mapped/{sampleID}__report_{refDB}.tsv")
    log_info(glue(
      "_____________________________________________________________\n\n",
      "Executing: {sampleID} with {refDB} \n\n",
    ))
    log_info(glue("Classified Reads: {classified_seqs_out}"))
    log_info(glue("Kraken2 Report: {report_out}"))
    # concatneate f/r/s reads
    log_info("Concatenating F/R/S reads")
    system2(
      command = "cat",
      args = c(
        glue("{reads_dir}{sampleID}_1.fq.gz"),
        glue("{reads_dir}{sampleID}_2.fq.gz"),
        glue("{reads_dir}{sampleID}_single.fq.gz"),
        ">", quality_reads
      )
    )
    # execute kraken2
    log_info("Sumbitting Kraken2 Alignment")
    system2(
      command = "kraken2",
      args = c(
        "--db", refdbs_kraken[[refDB]],
        "--threads", nworkers,
        "--classified-out", classified_seqs_out,
        "--report", report_out,
        "--gzip-compressed",
        quality_reads
      ), stdout = FALSE
    )
    # compress classified reads
    system2(command = "gzip", args = classified_seqs_out)
  }
  # delete concatenated reads
    log_debug("Removing concatenated reads")
    system2(command = "rm", args = quality_reads)
}

#______________________________________________________________________________

future::plan(
  future.batchtools::batchtools_slurm,
  template = glue("{wkdir}/batchtools_templates/batchtools-kraken2.slurm.tmpl"),
  resources = list(
    name = "Kraken2",
    memory = "30G",
    ncpus = 2,
    chunks.as.arrayjobs = FALSE
  )
)

sampleIDs <- readRDS(glue("{wkdir}/data/interim/metadata/sample-id_WGS.rds"))
kraken_jobs <- list()
logfile <- log_open(
  glue("{wkdir}/{Sys.Date()}_kraken2_slurm-job-submission.log")
)

for (sample in sampleIDs) {
  log_print(glue("Executing Alignment: {sample}"))
  job <- future(kraken2_alignment(
    sampleID = sample,
    nworkers = 2
  ), packages = c("glue", "logger"))
  # collect job info
  kraken_jobs[[sample]] <- capture.output(job$config$reg)
}

saveRDS(
    kraken_jobs,
    glue(
        "{wkdir}/data/interim/kraken_results/",
        "{Sys.Date()}_kraken2_batchtools-registry.rds"
    )
)


#' Go through job folders and collect log files
#' (will populate when jobs start) run this step when all jobs are finished
while (count_slurm_jobs(params = c('-u', 'jboktor', "--name=Kraken2")) > 0) {
  log_print(glue(
  "Slurm Jobs Running: ",
  "{count_slurm_jobs(params = c('-u', 'jboktor', '--states=RUNNING'))}",
  "  Pending: ",
  "{count_slurm_jobs(params = c('-u', 'jboktor', '--states=PENDING'))}"
  ))
  Sys.sleep(300)
}
log_close()
writeLines(readLines(logfile))

kraken_log_loc <- list()
for (job in names(kraken_jobs)) {
  job_path <- kraken_jobs[[job]][3] %>% gsub("  File dir : ", "", .)
  log_file_name <- unlist(list.files(file.path(job_path, "logs/")))[1]
  kraken_log_loc[[job]] <- file.path(job_path, "logs", log_file_name)
}
saveRDS(
    kraken_log_loc,
    glue(
        "{wkdir}/data/interim/kraken_results/",
        "{Sys.Date()}_kraken2_logfile-locations.rds"
    )
)
