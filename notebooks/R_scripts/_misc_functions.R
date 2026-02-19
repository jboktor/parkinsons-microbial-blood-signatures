# Miscellanous functions
#' This script contains function that are used
#' across multiple scripts notebooks 

# ______________________________________________________________________________
#                     Utility Functions
# ______________________________________________________________________________

# inverse of logical statement
`%nin%` <- Negate(`%in%`)


get_time <- function(){
  print(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"))
}

#' Function for executing shell command in R with STDOUT and STEDRR control
#' specifying TRUE for stdout_path redirects the output of the command 
#' into an R object
shell_do <- function(command_string, stdout_path = "", stderr_path = "") {
  full_command <- glue::glue("bash -c 'source ~/.bashrc && {command_string}'")
  inputs <- unlist(stringr::str_split(full_command, " "))
  system2(
    command = inputs[1],
    args = inputs[-1],
    stdout = stdout_path,
    stderr = stderr_path
  )
}

breg_setup <- function(cluster_name, log_dir, template_path, seed = 42) {
  require(batchtools)
  require(glue)
  
  # Create a unique subdirectory for this registry
  cluster_run <- glue("{get_time()}_{cluster_name}")
  reg_dir <- file.path(log_dir, cluster_run)
  
  # Create the registry
  breg <- batchtools::makeRegistry(
    file.dir = reg_dir,
    seed = seed
  )
  
  # Set up SLURM cluster functions
  breg$cluster.functions <- batchtools::makeClusterFunctionsSlurm(
    template = template_path
  )
  
  return(breg)
}

easy_batch <- function(cluster_name, fun, params, template_path, log_dir,
                       walltime = "10:00:00", memory = "50GB",
                       ncpus = 4, max_jobs = 9999, seed = 42, ngpus = FALSE) {
  require(batchtools)
  require(glue)
  
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  cluster_run <- glue("{get_time()}_{cluster_name}")
  
  # Use breg_setup to create the registry
  breg <- breg_setup(
    cluster_name = cluster_name,
    log_dir = log_dir,
    template_path = template_path,
    seed = seed
  )
  
  jobs <- batchtools::batchMap(fun = fun, args = params, reg = breg)
  
  if (!ngpus) {
  batchtools::submitJobs(jobs,
    resources = list(
      walltime = walltime,
      memory = memory,
      ncpus = ncpus,
      max.concurrent.jobs = max_jobs
    )
  )
  } else {
    batchtools::submitJobs(jobs,
      resources = list(
        walltime = walltime,
        memory_gpu = memory,
        ncpus = ncpus,
        ngpus = ngpus,
        max.concurrent.jobs = max_jobs
      )
    )
  }
  return(breg)
}

# Function to run any shell command using slurm and future.batchtools
slurm_shell_do <- function(cmd,
                           jobname = glue("slurm-shell-{get_time()}"),
                           working_dir = wkdir,
                           template_path = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
                           memory = "1G",
                           ncpus = 1,
                           walltime = 3600) {
  require(magrittr)
  require(future)
  require(future.batchtools)
  # Initiate future.batchtools backend for parallel processing
  future::plan(
    future.batchtools::batchtools_slurm,
    template = template_path,
    resources = list(
      name = jobname,
      memory = memory,
      ncpus = ncpus,
      walltime = walltime
    )
  )
  job %<-% shell_do(cmd)
}

# function to count ongoing slurm jobs
count_slurm_jobs <- function(params = c("-u", "jboktor")) {
  queue <- system2(
    command = "squeue",
    args = params,
    stdout = TRUE
  )
  length(queue) - 1
}

#' sleep until a condition is true
wait_until <- function(conditional, interval = 2) {
  # Keep looping until the condition is met
  while (!conditional()) {
    Sys.sleep(interval)
  }
}

check_slurm_overload <- function(njobs = 9999) {
  wait_until(function() {
    count_slurm_jobs() < njobs
  })
}

#' Function to rapidly extract the number of files matching a search string
#' while integrating a list of sample names.
#' Input:
#' search_pattern   - a string with glue syntax {.}, not yet glued
#'                      to indicate placement of sampleIDs
#' name_list        - list of names
#' nworkders        - number of threads to use for
#'                    parallel processing
#' Output:
#' A named list of samples with the number of files
#' matching the search critera for each sample
search_file_n <- function(search_pattern, name_list, nworkers = 2, ...) {
  plan(multisession, workers = nworkers)
  file_n <- name_list %>%
    purrr::set_names() %>%
    purrr::map(~glue(search_pattern)) %>%
    furrr::future_map(
      ~ system2(
        command = "ls",
        args = c(., "|", "wc", "-l"),
        stdout = TRUE
      ) %>% as.numeric(),
      .progress = TRUE
    )
  return(file_n)
  }

# chunking function from https://stackoverflow.com/a/16275428
chunk_func <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))

has_error_message <- function(stderr_file) {
  lines <- readLines(stderr_file)
  for (line in lines) {
    if (grepl("error|ERR|quota exceeded", line)) {
      return(TRUE)
    }
  }
  return(FALSE)
}

#' Splits lists of fasta STRINGS and headers into chunked fasta files outputs,
#' with a given number of chunks
write_fasta_chunks <- function(seqs, names, outdir, filename, nchunks) {
  # future::plan("multisession", workers = parallel::detectCores()/2)
  seq_chunks <- chunk_func(seqs, nchunks)
  name_chunks <- chunk_func(names, nchunks)
  purrr::walk(
    1:nchunks,
    ~ seqinr::write.fasta(
      sequences = as.list(seq_chunks[[.]]),
      names = name_chunks[[.]],
      file.out = glue("{outdir}/{filename}_chunk_{.}.fasta"),
      as.string = TRUE
    )
  )
}


#______________________________________________________________________________


slim_report_df <- function(report) {
  report_vars <- c("feature", "NCBI_taxon_ID", "taxon_rank_code")

  report %>%
    dplyr::select(
      all_of(report_vars),
      contains("clade_counts")
    ) %>%
    rename_at(
      vars(-c(feature)),
      function(x) gsub("clade_counts_", "", x)
    ) %>%
    pivot_longer(-all_of(report_vars), names_to = "participant_id", values_to = "count")
}


# Function to download files via wget and submit to slurm
wget_download_slurm <- function(jobname,
                                download_link,
                                slurm_out,
                                output_dir,
                                threads = 1,
                                walltime = "72:00:00",
                                mem_per_cpu = "1G") {
  shell_do(
    glue(
      "sbatch",
      " --job-name={jobname}",
      " --ntasks={threads}",
      " --output={slurm_out}/{jobname}.out",
      " --error={slurm_out}/{jobname}.err",
      " --time={walltime}",
      " --mem-per-cpu={mem_per_cpu}",
      " /central/home/jboktor/slurm_wget.sh",
      " -u {download_link}",
      " -o {output_dir}"
    )
  )
}

get_time <- function(){
  print(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"))
}

shell_do_bowtie2 <- function(bowtiew_ind,
                             inpath_f,
                             inpath_r,
                             outpath_sam,
                             outpath_fqs,
                             outpath_metrics,
                             outpath_stderr,
                             threads,
                             wkdir) {
  require(glue)
  source(glue("{wkdir}/notebooks/R_scripts/_misc_functions.R"))
  bowtie2_cmd <- glue(
    "bowtie2 -x {bowtiew_ind}",
    " -1 {inpath_f}",
    " -2 {inpath_r}",
    " -S {outpath_sam}",
    " --un-conc-gz {outpath_fqs}",
    " --met-file {outpath_metrics}",
    " --time",
    " --threads {threads}",
    " 2> {outpath_stderr}"
  )
  shell_do(bowtie2_cmd)
  # delete large temp sam file
  unlink(outpath_sam)
}

#' R wrapper for KrakenUniq classification of reads
shell_do_krakenunique <- function(id,
                                  inpath_f,
                                  inpath_r,
                                  outputdir,
                                  threads,
                                  wkdir) {
  require(glue)
  source(glue("{wkdir}/notebooks/R_scripts/_misc_functions.R"))
  message("Processing: ", id, "\n")
  ku_cmd <- glue(
    "krakenuniq",
    " --db /central/groups/MazmanianLab/joeB/Downloads/RefDBs/KrakenUniq/MicrobialDB",
    " --threads {threads}",
    " --paired",
    " --preload",
    " --output {outputdir}/{id}_KrakenUniq_read-classification.tsv",
    " --report-file {outputdir}/{id}_KrakenUniq_report.tsv",
    " {inpath_f} {inpath_r}",
    " > /dev/null"
  )
  shell_do(ku_cmd)
}

#  [1] "ARSA"          "BBURCG"        "BBWRCG"        "GW"           
#  [5] "GW_average"    "GW_complete"   "GW_single"     "GW_ward"      
#  [9] "HC"            "HC_average"    "HC_complete"   "HC_single"    
# [13] "HC_ward"       "Identity"      "MDS"           "MDS_angle"    
# [17] "MDS_metric"    "MDS_nonmetric" "OLO"           "OLO_average"  
# [21] "OLO_complete"  "OLO_single"    "OLO_ward"      "QAP_2SUM"     
# [25] "QAP_BAR"       "QAP_Inertia"   "QAP_LS"        "R2E"          
# [29] "Random"        "SA"            "Spectral"      "Spectral_norm"
# [33] "SPIN_NH"       "SPIN_STS"      "TSP"           "VAT"          

seriate_matrix_rows <- function(mat,
                                seriate_method = "HC_average",
                                dist_method = "euclidean",
                                nthreads = 8) {
  order <- mat %>%
    parallelDist::parDist(method = dist_method, threads = nthreads) %>%
    # stats::dist(method = dist_method) %>%
    seriation::seriate(method = seriate_method) %>%
    seriation::get_order()
  ranked_order <- rownames(mat)[order]
  return(ranked_order)
}

# seriate_matrix_rows <- function(
#     mat, seriate_method = "OLO", dist_method = "euclidean") {
#   order <- mat %>%
#     stats::dist(method = dist_method) %>%
#     seriation::seriate(method = seriate_method) %>%
#     seriation::get_order()
#   ranked_order <- rownames(mat)[order]
#   return(ranked_order)
# }

# Clustering samples for distance matrix and cluster membership matrices
seriate_dist_obj <- function(dist_obj, seriate_method = "OLO") {
  order <- seriation::seriate(dist_obj, method = seriate_method) %>%
    seriation::get_order()
  ranked_order <- labels(dist_obj)[order]
  return(ranked_order)
}
