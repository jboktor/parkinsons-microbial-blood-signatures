library(data.table)
library(glue)
library(tictoc)

n_workers <- 8 # Var to set
data.table::setDTthreads(n_workers)

# ---- Inputs ----
alignment_file <- "/resnick/groups/MazmanianLab/clarayu/full_run/fastqs/tsv_outputs/ERP000108/ERR011117.tsv"
mmseq2_file <- "/central/scratch/clarayu/mmseq2_clustering/clusters/cluster50_cluster.tsv"
mmseq2_level <- "50%"
read_ID <- "ERR011117"

# ---- Logging function ----
log_time <- function(msg) {
  message(glue("\n[{format(Sys.time(), '%Y-%m-%d %H:%M:%S')}] {msg}\n"))
}

# ---- Read data ----
log_time("Reading alignment data...")
alignments <- fread(alignment_file)
log_time(glue("Loaded {format(nrow(alignments), big.mark=',')} alignment records"))

log_time("Reading mmseq2 clustering data...")
mmseq2_map <- fread(mmseq2_file, col.names = c("cluster_rep", "cluster_member"))
log_time(glue("Loaded {format(nrow(mmseq2_map), big.mark=',')} cluster mappings"))

# Set keys for fast joins
log_time("Setting up data.table keys for optimized joins...")
setkey(mmseq2_map, cluster_member)
setkey(alignments, qname)
log_time("Keys set successfully")

log_time("Using fully vectorized approach...")
tic("Vectorized processing")

# Step 1: Join primary alignments with mmseq2 clusters
log_time("Step 1: Processing primary alignments...")
primary_alignments <- alignments[is_primary == TRUE & !is.na(primary_rname)]
log_time(glue("Found {format(nrow(primary_alignments), big.mark=',')} primary alignments"))

primary_with_clusters <- primary_alignments[mmseq2_map, on = c("primary_rname" = "cluster_member"), nomatch = 0]
log_time(glue("Joined primary alignments with clusters: {format(nrow(primary_with_clusters), big.mark=',')} records"))

# Step 2: Join secondary alignments with mmseq2 clusters
log_time("Step 2: Processing secondary alignments...")
secondary_alignments <- alignments[is_primary == FALSE & !is.na(rname)]
log_time(glue("Found {format(nrow(secondary_alignments), big.mark=',')} secondary alignments"))

secondary_with_clusters <- secondary_alignments[mmseq2_map, on = c("rname" = "cluster_member"), nomatch = 0]
log_time(glue("Joined secondary alignments with clusters: {format(nrow(secondary_with_clusters), big.mark=',')} records"))

# Step 3: Create secondary cluster lists using a simpler approach
log_time("Step 3: Aggregating secondary clusters by qname...")
secondary_summary <- secondary_with_clusters[, .N, by = .(qname, cluster_rep)]
log_time(glue("Created secondary summary with {format(nrow(secondary_summary), big.mark=',')} cluster counts"))

log_time("Converting secondary clusters to list format...")
secondary_cluster_lists <- tapply(
  seq_len(nrow(secondary_summary)), 
  secondary_summary$qname, 
  function(idx) {
    subset <- secondary_summary[idx, ]
    as.list(setNames(subset$N, subset$cluster_rep))
  }
)
log_time("Converted secondary clusters to list format")

# Step 4: Create final result with vectorized operations
log_time("Step 4: Creating final result table...")
results_dt <- primary_with_clusters[, .(
  read_ID = read_ID,
  primary_read_cluster = cluster_rep,
  MMSeq2_seqID = mmseq2_level
), by = qname]
log_time(glue("Created base result table with {format(nrow(results_dt), big.mark=',')} reads"))

log_time("Adding secondary cluster information...")
# Create a lookup table for vectorized access
secondary_lookup_dt <- data.table(
  qname = names(secondary_cluster_lists),
  secondary_clusters = secondary_cluster_lists
)
setkey(secondary_lookup_dt, qname)

# Use data.table join for maximum speed
results_dt <- results_dt[secondary_lookup_dt, on = "qname", secondary_read_clusters := i.secondary_clusters]

# Handle missing values efficiently
results_dt[is.na(secondary_read_clusters), secondary_read_clusters := list(list())]
log_time("Secondary clusters added successfully")

toc()

log_time(glue("Found {format(nrow(results_dt), big.mark=',')} valid primary reads."))

# ---- Save output ----
log_time("Saving results to RDS file...")
saveRDS(results_dt, glue("read_cluster_summary_{read_ID}_vectorized.rds"))
log_time("✅ Saved result successfully!")
# log_time("=== SUMMARY STATISTICS ===")
log_time(glue("Total alignment records: {format(nrow(alignments), big.mark=',')}"))
log_time(glue("Primary alignments: {format(nrow(primary_alignments), big.mark=',')}"))
log_time(glue("Secondary alignments: {format(nrow(secondary_alignments), big.mark=',')}"))
log_time(glue("Final processed reads: {format(nrow(results_dt), big.mark=',')}"))
log_time(glue("Unique clusters in mmseq2: {format(length(unique(mmseq2_map$cluster_rep)), big.mark=',')}"))
log_time("=== END SUMMARY ===")
