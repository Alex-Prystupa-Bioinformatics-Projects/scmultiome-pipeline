#!/usr/bin/env Rscript
# =============================================================================
# 09_annotate.R  |  Pipeline Step 9: Cell Type Annotation
# =============================================================================
# Reads configs/annotations.csv, left-joins annotation columns onto Seurat
# metadata by cluster ID, removes cells labelled "Discard" in any annotation
# column, and saves the annotated object.
#
# annotations.csv format:
#   wsnn_res.0.4,Annotation_Broad,Annotation_Narrow,...
#   - First column name   = wsnn_res.* column to join on
#   - Remaining columns   = annotation layers (any name, any number)
#   - Any cell labelled "Discard" in any column will be removed
#
# Usage:
#   Rscript scripts/09_annotate.R configs/samplesheet.csv \
#       --project_prefix myproject \
#       --RDS_file_in output/RDS-files/myproject-07-normalize-reduce-obj.RDS \
#       --annotation_file configs/annotations.csv
# =============================================================================

library(argparser, quietly = TRUE)
library(Seurat,    quietly = TRUE)
library(dplyr,     quietly = TRUE)
library(tibble,    quietly = TRUE)


# -----------------------------------------------------------------------------
# 1. Arguments
# -----------------------------------------------------------------------------
p <- arg_parser("Step 09: Cell Type Annotation")
p <- add_argument(p, "samplesheet",       help = "Path to samplesheet CSV")
p <- add_argument(p, "--project_prefix",  help = "Output file prefix", default = "multiome")
p <- add_argument(p, "--RDS_file_in",     help = "Input RDS (step 7 or 8 output)")
p <- add_argument(p, "--annotation_file", help = "Path to filled-in annotations.csv",
                  default = "configs/annotations.csv")
argv <- parse_args(p)

if (is.na(argv$RDS_file_in))            stop("--RDS_file_in is required.")
if (!file.exists(argv$RDS_file_in))     stop("RDS file not found: ",        argv$RDS_file_in)
if (!file.exists(argv$annotation_file)) stop("Annotation file not found: ", argv$annotation_file)

# -----------------------------------------------------------------------------
# 2. Load object and annotations
# -----------------------------------------------------------------------------
message("Loading RDS: ", argv$RDS_file_in)
seu_obj <- readRDS(argv$RDS_file_in)

anno_df    <- read.csv(argv$annotation_file, stringsAsFactors = FALSE, check.names = FALSE)
res_col    <- colnames(anno_df)[1]
annot_cols <- colnames(anno_df)[-1]
message("Joining on: ", res_col)
message("Annotation columns: ", paste(annot_cols, collapse = ", "))

anno_df[[res_col]] <- as.factor(anno_df[[res_col]])

# -----------------------------------------------------------------------------
# 3. Left-join annotations onto metadata
# -----------------------------------------------------------------------------
seu_obj@meta.data <- seu_obj@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(anno_df, by = res_col) %>%
    tibble::column_to_rownames("barcodes")

# -----------------------------------------------------------------------------
# 4. Remove cells labelled "Discard" in any annotation column
# -----------------------------------------------------------------------------
n_before <- ncol(seu_obj)

keep <- Reduce("&", lapply(annot_cols, function(col) {
    seu_obj@meta.data[[col]] != "Discard"
}))

seu_obj <- seu_obj[, keep]
message("Removed ", n_before - ncol(seu_obj), " Discard cells; ",
        ncol(seu_obj), " cells remaining.")

# -----------------------------------------------------------------------------
# 5. Save
# -----------------------------------------------------------------------------
out_rds <- file.path("output/RDS-files",
                     paste0(argv$project_prefix, "-09-annotate-obj.RDS"))
saveRDS(seu_obj, file = out_rds)
message("Saved: ", out_rds)
message("\nStep 9 complete.")
