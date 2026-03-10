#!/usr/bin/env Rscript
# =============================================================================
# 06_merge.R  |  Pipeline Step 6: Consensus Peaks + Merge
# =============================================================================
# Builds a consensus peak set across all filtered samples, recounts each sample
# against the unified peak universe, then merges all objects into one.
#
# Why consensus peaks: each sample has its own MACS2-called peak set from step 3.
# For multi-sample ATAC analysis, all samples must share the same peak feature
# space. This step reduces all per-sample peaks to a union set, recounts, then
# merges so downstream LSI/SVD uses a consistent feature matrix.
#
# NOTE: JoinLayers() is NOT called here. In Seurat v5, RNA data is stored as
# separate per-sample layers after merge. SCTransform in step 7 runs on each
# layer independently — JoinLayers() is called after SCTransform.
#
# Usage:
#   Rscript scripts/06_merge.R configs/samplesheet.csv \
#       --project_prefix myproject \
#       --RDS_file_in output/RDS-files/myproject-05-filter-obj-list.RDS
#
# Input:  output of 05_filter.R  ({prefix}-05-filter-obj-list.RDS)
# Output: output/RDS-files/{prefix}-06-merge-obj.RDS
# =============================================================================

library(argparser,     quietly = TRUE)
library(yaml,          quietly = TRUE)
library(Seurat,        quietly = TRUE)
library(Signac,        quietly = TRUE)
library(GenomicRanges, quietly = TRUE)

# 1. Parse arguments
p <- arg_parser("Step 6: Build consensus peaks and merge filtered objects")
p <- add_argument(p, "samplesheet",
                  help = "Path to samplesheet CSV", type = "character")
p <- add_argument(p, "--pipeline_config",
                  help    = "Path to pipeline_config.yml",
                  default = "configs/pipeline_config.yml", type = "character")
p <- add_argument(p, "--project_prefix",
                  help    = "Output file prefix",
                  default = "multiome", type = "character")
p <- add_argument(p, "--RDS_file_in",
                  help = "Input RDS file — output of 05_filter.R", type = "character")
argv <- parse_args(p)

if (is.na(argv$RDS_file_in) || argv$RDS_file_in == "NA") {
    stop("--RDS_file_in is required. Provide the path to output of 05_filter.R.")
}

# 2. Load config and helper functions
pipeline_config <- yaml::read_yaml(argv$pipeline_config)
samplesheet     <- read.csv(argv$samplesheet)

source("scripts/functions.R")
source("scripts/genome_utils.R")

# 3. Load genome (annotation used when rebuilding ChromatinAssay)
genome <- load_genome(pipeline_config)

# 4. Load filtered object list from step 5
message("Starting Step 6: merge")
message("Loading: ", argv$RDS_file_in)
seu_obj_list <- readRDS(argv$RDS_file_in)

# 5. Build consensus peak set and recount per sample
#    Each object is stripped to RNA + new consensus peaks assay
message("Building consensus peaks and recounting per sample...")
seu_obj_list <- create_consensus_peaks(seu_obj_list, genome)
message("Consensus peaks complete.")

# 6. Merge all objects into one
#    Seurat v5: RNA data stored as separate layers per sample after merge
#    DO NOT call JoinLayers() here — SCTransform in step 7 runs per layer
message("Merging ", length(seu_obj_list), " objects...")
seu_obj <- merge(seu_obj_list[[1]], seu_obj_list[-1])
message("  Merge complete. Total cells: ", ncol(seu_obj))

# 7. Save merged object (single object — not a list)
out_path <- paste0("output/RDS-files/", argv$project_prefix, "-06-merge-obj.RDS")
saveRDS(seu_obj, file = out_path)
message("Step 6 complete. Saved: ", out_path)
