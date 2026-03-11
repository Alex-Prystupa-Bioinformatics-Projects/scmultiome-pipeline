#!/usr/bin/env Rscript
# =============================================================================
# 06_merge.R  |  Pipeline Step 6: Consensus Peaks + SCTransform + Merge
# =============================================================================
# Builds a consensus peak set across all filtered samples, recounts each sample
# against the unified peak universe, runs SCTransform per sample, selects
# integration features, then merges all objects into one.
#
# Why consensus peaks: each sample has its own MACS2-called peak set from step 3.
# For multi-sample ATAC analysis, all samples must share the same peak feature
# space. This step reduces all per-sample peaks to a union set, recounts, then
# merges so downstream LSI/SVD uses a consistent feature matrix.
#
# Why SCTransform here (before merge): running SCTransform per sample gives
# per-sample variable features. SelectIntegrationFeatures() then picks the top
# nfeatures genes ranked by consistency across all samples — a principled
# consensus. This avoids PCA being driven by genes variable in only one sample.
# JoinLayers() is called in step 7 after PCA variable features are set.
#
# Usage:
#   Rscript scripts/06_merge.R configs/samplesheet.csv \
#       --project_prefix myproject \
#       --RDS_file_in output/RDS-files/myproject-05-filter-obj-list.RDS \
#       --nfeatures 3000
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
p <- arg_parser("Step 6: Build consensus peaks, SCTransform per sample, and merge")
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
p <- add_argument(p, "--nfeatures",
                  help    = "Number of integration features for SelectIntegrationFeatures()",
                  default = 3000L, type = "integer")
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

# 6. SCTransform per sample
#    return.only.var.genes = FALSE so all genes are present after merge
#    Per-sample models stored in SCTModel.list for downstream use
message("Running SCTransform per sample...")
seu_obj_list <- lapply(seu_obj_list, function(seu_obj) {
    DefaultAssay(seu_obj) <- "RNA"
    SCTransform(seu_obj, return.only.var.genes = FALSE, verbose = FALSE)
})

# 7. Select integration features across all samples
#    Ranks genes by how consistently variable they are across samples —
#    a principled weighted intersection rather than sample-biased union
message("Selecting integration features (nfeatures = ", argv$nfeatures, ")...")
var_features <- SelectIntegrationFeatures(
    object.list = seu_obj_list,
    nfeatures    = argv$nfeatures
)
message("  Integration features selected: ", length(var_features))

# 8. Merge all objects into one
#    merge.data = TRUE preserves per-sample SCT normalized data across layers
message("Merging ", length(seu_obj_list), " objects...")
seu_obj <- merge(seu_obj_list[[1]], seu_obj_list[-1], merge.data = TRUE)
message("  Merge complete. Total cells: ", ncol(seu_obj))

# Set consensus variable features on the merged SCT assay
VariableFeatures(seu_obj) <- var_features

# 9. Save merged object (single object — not a list)
out_path <- paste0("output/RDS-files/", argv$project_prefix, "-06-merge-obj.RDS")
saveRDS(seu_obj, file = out_path)
message("Step 6 complete. Saved: ", out_path)
