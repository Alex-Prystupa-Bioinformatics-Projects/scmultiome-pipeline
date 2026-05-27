#!/usr/bin/env Rscript
# =============================================================================
# merge_linkpeaks.R  |  LinkPeaks Merge Step
# =============================================================================
# Runs after all per-group linkpeaks child jobs finish successfully.
# Collects per-group Seurat objects from output/RDS/linkpeaks/, extracts
# peak-gene links from each, combines them into a single GRanges, and
# assigns the combined links back onto the unsplit step-10 object.
#
# Usage (called by scripts/helpers/job_merge.sh):
#   Rscript scripts/helpers/merge_linkpeaks.R \
#       output/RDS-files/{prefix}-10-recall-peaks-obj.RDS \
#       output/RDS/linkpeaks \
#       configs/pipeline_config.yml
#
# Input:  output/RDS-files/{prefix}-10-recall-peaks-obj.RDS  (step-10 object)
#         output/RDS/linkpeaks/*-linkpeaks.RDS               (per-group worker outputs)
# Output: output/RDS-files/{prefix}-11-linkpeaks-obj.RDS
# =============================================================================

library(Seurat,         quietly = TRUE)
library(Signac,         quietly = TRUE)
library(dplyr,          quietly = TRUE)
library(GenomicRanges,  quietly = TRUE)

# -----------------------------------------------------------------------------
# 1. Args
# -----------------------------------------------------------------------------
argv          <- commandArgs(trailingOnly = TRUE)
step10_rds    <- argv[1]
linkpeaks_dir <- argv[2]
config_path   <- argv[3]

project_prefix <- sub("-10-recall-peaks-obj\\.RDS$", "", basename(step10_rds))

message("=== LinkPeaks Merge ===")
message("Step-10 RDS:    ", step10_rds)
message("Project prefix: ", project_prefix)

# -----------------------------------------------------------------------------
# 3. Load unsplit merged object from step 10
# -----------------------------------------------------------------------------
message("Loading step-10 object...")
seu_obj <- readRDS(step10_rds)
message("  Cells: ", ncol(seu_obj))

# -----------------------------------------------------------------------------
# 4. Collect all per-group linkpeaks RDS files (excludes split/ subdirectory)
# -----------------------------------------------------------------------------
rds_files <- list.files(linkpeaks_dir, pattern = "-linkpeaks\\.RDS$", full.names = TRUE, recursive = FALSE)
message("Per-group RDS files found: ", length(rds_files))

# -----------------------------------------------------------------------------
# 5. Extract links from each group, combine into one data frame
# -----------------------------------------------------------------------------
message("Extracting and combining links across all groups...")

links_df <- dplyr::bind_rows(lapply(rds_files, function(f) {
    group_name <- sub("-linkpeaks\\.RDS$", "", basename(f))
    obj        <- readRDS(f)
    gr         <- obj[["peaks"]]@links

    if (is.null(gr) || length(gr) == 0) {
        message("  ", group_name, ": no links detected — skipping")
        return(NULL)
    }

    df          <- as.data.frame(gr)
    df$Celltype <- group_name
    message("  ", group_name, ": ", nrow(df), " links")
    df
}))

message("Total links combined: ", nrow(links_df))

# -----------------------------------------------------------------------------
# 6. Convert combined links df back to GRanges
# -----------------------------------------------------------------------------
links_gr <- GenomicRanges::makeGRangesFromDataFrame(
    links_df,
    keep.extra.columns = TRUE,
    seqnames.field     = "seqnames",
    start.field        = "start",
    end.field          = "end"
)

# -----------------------------------------------------------------------------
# 7. Assign combined links onto the unsplit step-10 object and save
# -----------------------------------------------------------------------------
seu_obj[["peaks"]]@links <- links_gr

out_rds <- file.path("output/RDS-files", paste0(project_prefix, "-11-linkpeaks-obj.RDS"))
saveRDS(seu_obj, file = out_rds)
message("Saved: ", out_rds)
message("=== merge_linkpeaks complete ===")
