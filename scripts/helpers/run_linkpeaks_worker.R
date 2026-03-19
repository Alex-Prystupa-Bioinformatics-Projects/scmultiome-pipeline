#!/usr/bin/env Rscript
# =============================================================================
# run_linkpeaks_worker.R  |  LinkPeaks Worker (per-group)
# =============================================================================
# Called by 11_linkpeaks.R — one job per group defined by --grouping_col.
# Runs Signac LinkPeaks on a single subset Seurat object and saves the result.
#
# Usage (called by generated LSF job scripts):
#   Rscript scripts/helpers/run_linkpeaks_worker.R \
#       <group_name> <group_rds> <pipeline_config> <output_dir>
#
# Input:  Single-group Seurat object RDS (output/RDS/linkpeaks/split/{group}.RDS)
# Output: {output_dir}/{group_name}-linkpeaks.RDS
# =============================================================================

library(Seurat,  quietly = TRUE)
library(Signac, quietly = TRUE)

# -----------------------------------------------------------------------------
# 1. Args
# -----------------------------------------------------------------------------
argv        <- commandArgs(trailingOnly = TRUE)
group_name  <- argv[1]
group_rds   <- argv[2]
config_path <- argv[3]
output_dir  <- argv[4]

if (any(is.na(c(group_name, group_rds, config_path, output_dir)))) {
    stop("Usage: run_linkpeaks_worker.R <group_name> <group_rds> <pipeline_config> <output_dir>")
}
if (!file.exists(group_rds)) {
    stop("Group RDS not found: ", group_rds)
}

message("=== LinkPeaks Worker ===")
message("Group:      ", group_name)
message("Input RDS:  ", group_rds)
message("Output dir: ", output_dir)

# -----------------------------------------------------------------------------
# 2. Load genome resources and helper functions
# -----------------------------------------------------------------------------
source("scripts/functions.R")
source("scripts/genome_utils.R")

pipeline_config <- yaml::read_yaml(config_path)
genome          <- load_genome(pipeline_config)

# -----------------------------------------------------------------------------
# 3. Load group Seurat object
# -----------------------------------------------------------------------------
message("Loading object for group: ", group_name)
seu_obj <- readRDS(group_rds)
message("  Cells: ", ncol(seu_obj), " | Peaks: ", nrow(seu_obj[["peaks"]]))

# -----------------------------------------------------------------------------
# 4. Run LinkPeaks
# -----------------------------------------------------------------------------
linkpeaks_genes <- pipeline_config$linkpeaks_genes
if (is.null(linkpeaks_genes)) linkpeaks_genes <- "all"

genes_use <- NULL
if (linkpeaks_genes == "variable_features") {
    genes_use <- seu_obj@assays$SCT@var.features
    message("Running LinkMyPeaks on SCT variable features only (", length(genes_use), " genes)...")
} else {
    message("Running LinkMyPeaks on all genes...")
}

seu_obj <- LinkMyPeaks(
    seu_obj,
    peak.genome     = genome$peak.genome,
    genes           = genes_use,
    distance.to.use = 250001
)
n_links <- length(Links(seu_obj[["peaks"]]))
message("LinkPeaks complete. Links identified: ", n_links)

# -----------------------------------------------------------------------------
# 5. Save
# -----------------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
safe_name <- gsub("[^A-Za-z0-9._-]", "_", group_name)
out_rds   <- file.path(output_dir, paste0(safe_name, "-linkpeaks.RDS"))
saveRDS(seu_obj, file = out_rds)
message("Saved: ", out_rds)
message("\nWorker complete for group: ", group_name)
