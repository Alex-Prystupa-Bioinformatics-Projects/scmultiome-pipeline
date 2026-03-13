#!/usr/bin/env Rscript
# =============================================================================
# 08_cluster_markers.R  |  Pipeline Step 8: Cluster Marker Detection
# =============================================================================
# Loads the step 7 object and runs cluster marker detection at every wsnn_res.*
# resolution found in the object. Uses custom_all_markers_function() (Wilcoxon,
# positive markers only) with Libra-formatted output columns.
#
# For each resolution, three output files are written:
#   Raw-Markers      — all markers pre-filtering              (.csv)
#   Top-Markers      — top N markers per cluster (pivot)      (.csv)
#   Filtered-Markers — per-cluster data frames, one tab each  (.xlsx)
#
# Usage:
#   Rscript scripts/08_cluster_markers.R configs/samplesheet.csv \
#       --project_prefix myproject \
#       --RDS_file_in output/RDS-files/myproject-07-normalize-reduce-obj.RDS
#
# Input:  output of 07_normalize_reduce.R  ({prefix}-07-normalize-reduce-obj.RDS)
# Output: output/markers/{prefix}-08-markers/{res_label}/
# =============================================================================

library(argparser,  quietly = TRUE)
library(yaml,       quietly = TRUE)
library(Seurat,     quietly = TRUE)
library(Signac,     quietly = TRUE)
library(dplyr,      quietly = TRUE)
library(tibble,     quietly = TRUE)
library(tidyr,      quietly = TRUE)
library(Matrix,     quietly = TRUE)
library(openxlsx,   quietly = TRUE)
library(ggplot2,    quietly = TRUE)

# 1. Parse arguments
p <- arg_parser("Step 8: Cluster marker detection")
p <- add_argument(p, "samplesheet",
                  help = "Path to samplesheet CSV", type = "character")
p <- add_argument(p, "--pipeline_config",
                  help    = "Path to pipeline_config.yml",
                  default = "configs/pipeline_config.yml", type = "character")
p <- add_argument(p, "--project_prefix",
                  help    = "Output file prefix",
                  default = "multiome", type = "character")
p <- add_argument(p, "--RDS_file_in",
                  help = "Input RDS file — output of 07_normalize_reduce.R", type = "character")
p <- add_argument(p, "--p_adj_max",
                  help    = "Bonferroni-adjusted p-value cutoff for filtering",
                  default = 0.05, type = "double")
p <- add_argument(p, "--lfc_min",
                  help    = "Minimum avg_log2FC for filtering",
                  default = 1.0, type = "double")
p <- add_argument(p, "--pct_min",
                  help    = "Minimum fraction of cells in cluster expressing the gene",
                  default = 0.25, type = "double")
p <- add_argument(p, "--top_n",
                  help    = "Top N markers per cluster for Top-Markers pivot table",
                  default = 50L, type = "integer")
argv <- parse_args(p)

if (is.na(argv$RDS_file_in) || argv$RDS_file_in == "NA") {
    stop("--RDS_file_in is required. Provide the path to output of 07_normalize_reduce.R.")
}

# 2. Load config and source helper functions
pipeline_config <- yaml::read_yaml(argv$pipeline_config)
source("scripts/functions.R")
source("scripts/genome_utils.R")

# 3. Load step 7 object
message("Starting Step 8: cluster markers")
message("Loading: ", argv$RDS_file_in)
seu_obj <- readRDS(argv$RDS_file_in)

# 4. Preprocess for DE
#    NormalizeData on RNA (step 7 uses SCTransform; log-normalised counts are
#    needed for the expression display columns in custom_all_markers_function)
message("Preprocessing: NormalizeData + converting character metadata to factors...")
DefaultAssay(seu_obj) <- "RNA"
seu_obj <- NormalizeData(seu_obj, verbose = FALSE)
seu_obj@meta.data <- seu_obj@meta.data %>%
    dplyr::mutate(dplyr::across(where(is.character), as.factor))

# 5. Detect all wsnn_res.* clustering columns present in the object
res_cols <- grep("^wsnn_res\\.", colnames(seu_obj@meta.data), value = TRUE)
if (length(res_cols) == 0) {
    stop("No wsnn_res.* columns found in object metadata. Ensure step 7 completed successfully.")
}
message("Found ", length(res_cols), " resolution(s): ", paste(res_cols, collapse = ", "))

# 6. Create output directories
markers_dir   <- file.path("output/markers", paste0(argv$project_prefix, "-08-markers"))
annotation_dir <- file.path("output/plots", paste0(argv$project_prefix, "-annotation"))
dir.create(markers_dir,                            recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(annotation_dir, "dotplot"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(annotation_dir, "heatmap"),   recursive = TRUE, showWarnings = FALSE)
message("Markers directory: ", markers_dir)
message("Plots directory:   ", annotation_dir)

# Scale RNA using variable features determined by SCT in step 6
message("Scaling RNA data for heatmap plotting...")
DefaultAssay(seu_obj) <- "RNA"
seu_obj <- ScaleData(seu_obj, features = VariableFeatures(seu_obj), assay = "RNA", verbose = FALSE)

# 7. Loop over each resolution — run markers and save outputs
for (res_col in res_cols) {

    # Derive a clean label e.g. wsnn_res.0.2 -> res0.2
    res_label <- sub("^wsnn_res\\.", "res", res_col)
    message("\n--- Resolution: ", res_col, " (", res_label, ") ---")

    # Set proper numeric cluster ordering once — used by both markers and plots
    cluster_levels <- as.character(sort(as.numeric(unique(as.character(seu_obj@meta.data[[res_col]])))))
    seu_obj@meta.data[[res_col]] <- factor(seu_obj@meta.data[[res_col]], levels = cluster_levels)

    res_dir <- file.path(markers_dir, res_label)
    dir.create(res_dir, showWarnings = FALSE)

    # Run marker detection at this resolution
    message("  Running custom_all_markers_function...")
    markers <- custom_all_markers_function(
        seu_obj         = seu_obj,
        ident_col       = res_col,
        assay           = "RNA",
        p_adj_max       = argv$p_adj_max,
        lfc_min         = argv$lfc_min,
        pct_min         = argv$pct_min,
        top_n           = argv$top_n,
        background_name = "background"
    )

    # Save Raw-Markers CSV
    raw_path <- file.path(res_dir,
        paste0(argv$project_prefix, "-", res_label, "-raw-markers.csv"))
    write.csv(markers[["Raw-Markers"]], file = raw_path, row.names = FALSE)
    message("  Saved: ", raw_path)

    # Save Top-Markers CSV
    top_path <- file.path(res_dir,
        paste0(argv$project_prefix, "-", res_label, "-top-markers.csv"))
    write.csv(markers[["Top-Markers"]], file = top_path, row.names = FALSE)
    message("  Saved: ", top_path)

    # Save Filtered-Markers Excel workbook — one tab per cluster
    xlsx_path <- file.path(res_dir,
        paste0(argv$project_prefix, "-", res_label, "-filtered-markers.xlsx"))
    wb <- openxlsx::createWorkbook()
    filtered_list <- markers[["Filtered-Markers"]]
    for (ct in names(filtered_list)) {
        openxlsx::addWorksheet(wb, sheetName = ct)
        openxlsx::writeData(wb, sheet = ct, x = filtered_list[[ct]])
    }
    openxlsx::saveWorkbook(wb, file = xlsx_path, overwrite = TRUE)
    message("  Saved: ", xlsx_path,
            " (", length(filtered_list), " cluster tab(s))")

    # Extract top genes from Top-Markers (already sorted + filtered)
    top_genes_dot <- markers[["Top-Markers"]] %>%
        head(4) %>%
        dplyr::select(-rank) %>%
        unlist() %>% na.omit() %>% unique() %>% as.character()

    top_genes_heat <- markers[["Top-Markers"]] %>%
        head(8) %>%
        dplyr::select(-rank) %>%
        unlist() %>% na.omit() %>% unique() %>% as.character()

    n_clusters <- length(cluster_levels)
    Idents(seu_obj) <- seu_obj@meta.data[[res_col]]

    # DotPlot — top 4 markers per cluster
    if (length(top_genes_dot) > 0) {
        DefaultAssay(seu_obj) <- "RNA"
        p <- DotPlot(seu_obj, features = top_genes_dot) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            ggtitle(paste0("Top 4 Markers per Cluster (", res_label, ")"))
        pdf(file.path(annotation_dir, "dotplot",
                      paste0(argv$project_prefix, "-", res_label, "-dotplot.pdf")),
            width  = max(8, length(top_genes_dot) * 0.5),
            height = max(6, n_clusters * 0.6))
        print(p)
        dev.off()
        message("  Saved dotplot: ", res_label)
    }

    # DoHeatmap — top 8 markers per cluster
    if (length(top_genes_heat) > 0) {
        DefaultAssay(seu_obj) <- "RNA"
        p <- DoHeatmap(seu_obj, features = top_genes_heat) +
            ggtitle(paste0("Top 8 Markers per Cluster (", res_label, ")"))
        pdf(file.path(annotation_dir, "heatmap",
                      paste0(argv$project_prefix, "-", res_label, "-heatmap.pdf")),
            width  = max(18, n_clusters * 1.0),
            height = max(6, length(top_genes_heat) * 0.4))
        print(p)
        dev.off()
        message("  Saved heatmap: ", res_label)
    }
}

# 8. Summary
message("\nStep 8 complete.")
message("  Resolutions processed: ", paste(res_cols, collapse = ", "))
message("  Outputs written to:    ", markers_dir)
