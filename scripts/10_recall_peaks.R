#!/usr/bin/env Rscript
# =============================================================================
# 10_recall_peaks.R  |  Pipeline Step 10: Cell-Type-Aware Peak Recall
# =============================================================================
# Recalls MACS2 peaks per sample × cell type group on the annotated merged
# object. Signac's CallPeaks() handles fragment subsetting per group internally
# using cell barcodes — no manual splitting needed.
#
# The new cell-type-specific peaks are unioned with the existing peak set,
# filtered, and recounted across all cells, replacing the peaks assay.
#
# Usage:
#   Rscript scripts/10_recall_peaks.R configs/samplesheet.csv \
#       --project_prefix myproject \
#       --RDS_file_in output/RDS-files/myproject-09-annotate-obj.RDS \
#       --annotation_file configs/annotations.csv \
#       --macs2_path /path/to/macs2
#
# Input:  output of 09_annotate.R  (--RDS_file_in)
# Output: output/RDS-files/{project_prefix}-10-recall-peaks-obj.RDS
# =============================================================================

library(argparser,    quietly = TRUE)
library(yaml,         quietly = TRUE)
library(Seurat,       quietly = TRUE)
library(Signac,       quietly = TRUE)
library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Arguments
# -----------------------------------------------------------------------------
p <- arg_parser("Step 10: Cell-Type-Aware Peak Recall")
p <- add_argument(p, "samplesheet",
                  help = "Path to samplesheet CSV", type = "character")
p <- add_argument(p, "--pipeline_config",
                  help    = "Path to pipeline_config.yml",
                  default = "configs/pipeline_config.yml", type = "character")
p <- add_argument(p, "--project_prefix",
                  help    = "Output file prefix (no .RDS extension)",
                  default = "multiome", type = "character")
p <- add_argument(p, "--RDS_file_in",
                  help = "Input RDS — output of 09_annotate.R", type = "character")
p <- add_argument(p, "--annotation_file",
                  help    = "Path to annotations.csv",
                  default = "configs/annotations.csv", type = "character")
p <- add_argument(p, "--macs2_path",
                  help = "Path to MACS2 executable", type = "character")
argv <- parse_args(p)

if (is.na(argv$RDS_file_in) || argv$RDS_file_in == "NA") {
    stop("--RDS_file_in is required.")
}
if (!file.exists(argv$RDS_file_in)) {
    stop("RDS file not found: ", argv$RDS_file_in)
}
if (is.na(argv$macs2_path) || argv$macs2_path == "NA") {
    stop("--macs2_path is required.")
}
if (!file.exists(argv$annotation_file)) {
    stop("Annotation file not found: ", argv$annotation_file)
}

# -----------------------------------------------------------------------------
# 2. Load config and helper functions
# -----------------------------------------------------------------------------
pipeline_config <- yaml::read_yaml(argv$pipeline_config)

source("scripts/functions.R")
source("scripts/genome_utils.R")

# -----------------------------------------------------------------------------
# 3. Load genome resources and input object
# -----------------------------------------------------------------------------
genome <- load_genome(pipeline_config)

message("Starting Step 10: recall_peaks")
message("Loading: ", argv$RDS_file_in)
seu_obj <- readRDS(argv$RDS_file_in)

# -----------------------------------------------------------------------------
# 4. Determine annotation column (second column of annotations.csv)
# -----------------------------------------------------------------------------
anno_df   <- read.csv(argv$annotation_file, stringsAsFactors = FALSE, check.names = FALSE)
annot_col <- colnames(anno_df)[2]
message("Using annotation column: ", annot_col)

# -----------------------------------------------------------------------------
# 5. Create sample × cell type grouping column
# -----------------------------------------------------------------------------
seu_obj$Sample_Celltype <- paste0(
    seu_obj$orig.ident, "_",
    seu_obj[[annot_col, drop = TRUE]]
)
n_groups <- length(unique(seu_obj$Sample_Celltype))
message("Sample × cell type groups (", n_groups, "): ",
        paste(sort(unique(seu_obj$Sample_Celltype)), collapse = ", "))

# -----------------------------------------------------------------------------
# 6. Recall peaks per sample × cell type group
#    Signac's CallPeaks() subsets fragments per group by barcode internally —
#    no manual splitting needed. New peaks are unioned with existing peaks,
#    filtered, and recounted across all cells.
# -----------------------------------------------------------------------------
n_peaks_before <- nrow(seu_obj[["peaks"]])
message("Peaks before recall: ", n_peaks_before)
message("Recalling peaks across ", n_groups, " groups (sample × cell type)...")
seu_obj <- CallMyPeaks(
    seu_obj,
    grouping.var  = "Sample_Celltype",
    my.macs2.path = argv$macs2_path,
    my.annotation = genome$annotation,
    my.blacklist  = genome$blacklist
)
n_peaks_after <- nrow(seu_obj[["peaks"]])
message("Peak recall complete.")
message("  Peaks before: ", n_peaks_before)
message("  Peaks after:  ", n_peaks_after)
message("  New peaks added: ", n_peaks_after - n_peaks_before)

# -----------------------------------------------------------------------------
# 7. Save
# -----------------------------------------------------------------------------
out_rds <- file.path("output/RDS-files",
                     paste0(argv$project_prefix, "-10-recall-peaks-obj.RDS"))
saveRDS(seu_obj, file = out_rds)
message("Saved: ", out_rds)
message("\nStep 10 complete.")
