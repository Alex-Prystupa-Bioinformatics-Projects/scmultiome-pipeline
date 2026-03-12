#!/usr/bin/env Rscript
# =============================================================================
# 07_normalize_reduce.R  |  Pipeline Step 7: Dimensionality Reduction + Cluster
# =============================================================================
# Takes the merged object from step 6 (which already has SCTransform run
# per sample and consensus variable features set via SelectIntegrationFeatures)
# and runs the full dimensionality reduction and clustering pipeline.
#
# RNA:  JoinLayers → PCA (using pre-set consensus variable features) →
#       optional Harmony → UMAP
# ATAC: TF-IDF → FindTopFeatures → SVD/LSI → optional Harmony → UMAP
# WNN:  FindMultiModalNeighbors → WNN UMAP → multi-resolution FindClusters
#
# Harmony: optional batch correction on PCA + LSI. Pass one or more metadata
# column names as a comma-separated string via --harmony_vars.
#
# Usage:
#   Rscript scripts/07_normalize_reduce.R configs/samplesheet.csv \
#       --project_prefix myproject \
#       --RDS_file_in output/RDS-files/myproject-06-merge-obj.RDS \
#       --harmony_vars Institution
#
# Input:  output of 06_merge.R  ({prefix}-06-merge-obj.RDS)
# Output: output/RDS-files/{prefix}-07-normalize-reduce-obj.RDS
# =============================================================================

library(argparser,  quietly = TRUE)
library(yaml,       quietly = TRUE)
library(Seurat,     quietly = TRUE)
library(Signac,     quietly = TRUE)
library(harmony,    quietly = TRUE)
library(ggplot2,    quietly = TRUE)
library(patchwork,  quietly = TRUE)

# 1. Parse arguments
p <- arg_parser("Step 7: Dimensionality reduction and clustering")
p <- add_argument(p, "samplesheet",
                  help = "Path to samplesheet CSV", type = "character")
p <- add_argument(p, "--pipeline_config",
                  help    = "Path to pipeline_config.yml",
                  default = "configs/pipeline_config.yml", type = "character")
p <- add_argument(p, "--project_prefix",
                  help    = "Output file prefix",
                  default = "multiome", type = "character")
p <- add_argument(p, "--RDS_file_in",
                  help = "Input RDS file — output of 06_merge.R", type = "character")
p <- add_argument(p, "--harmony_vars",
                  help    = "Comma-separated metadata columns for Harmony batch correction (e.g. 'Institution' or 'Institution,Sequencing_Date'). Use 'none' to skip.",
                  default = "none", type = "character")
p <- add_argument(p, "--rna_pcs",
                  help = "Number of PCs for RNA dimensionality reduction", default = 30, type = "integer")
p <- add_argument(p, "--atac_pcs",
                  help = "Number of PCs for ATAC dimensionality reduction", default = 30, type = "integer")
p <- add_argument(p, "--resolution",
                  help    = "Comma-separated clustering resolutions",
                  default = "0.2,0.4,0.6,0.8,1.0", type = "character")
argv <- parse_args(p)

if (is.na(argv$RDS_file_in) || argv$RDS_file_in == "NA") {
    stop("--RDS_file_in is required. Provide the path to output of 06_merge.R.")
}

# 2. Parse multi-value args
run_harmony    <- argv$harmony_vars != "none"
harmony_vars   <- if (run_harmony) strsplit(argv$harmony_vars, ",")[[1]] else NULL
resolutions    <- as.numeric(strsplit(argv$resolution, ",")[[1]])
rna_pcs        <- argv$rna_pcs
atac_pcs       <- argv$atac_pcs

# 3. Load config and helper functions
pipeline_config <- yaml::read_yaml(argv$pipeline_config)
samplesheet     <- read.csv(argv$samplesheet)

source("scripts/functions.R")
source("scripts/genome_utils.R")

# 4. Load merged object from step 6
#    SCTransform has already been run per sample and consensus variable features
#    are pre-set via SelectIntegrationFeatures() in step 6
message("Starting Step 7: reduce + cluster")
message("Loading: ", argv$RDS_file_in)
seu_obj <- readRDS(argv$RDS_file_in)
message("  Consensus variable features loaded: ", length(VariableFeatures(seu_obj)))

# ── RNA ───────────────────────────────────────────────────────────────────────

# 5. JoinLayers on RNA before PCA
message("Joining layers...")
seu_obj <- JoinLayers(seu_obj, assay = "RNA")

# 6. PCA using consensus variable features set in step 6
message("Running PCA...")
seu_obj <- RunPCA(seu_obj, npcs = rna_pcs, verbose = FALSE)

# 7. Optional Harmony on PCA (RNA)
rna_reduction <- "pca"
if (run_harmony) {
    message("Running Harmony on RNA (vars: ", paste(harmony_vars, collapse = ", "), ")...")
    seu_obj <- RunHarmony(
        seu_obj,
        group.by.vars  = harmony_vars,
        reduction      = "pca",
        assay.use      = "SCT",
        dims.use       = 1:rna_pcs,
        reduction.save = "rna.harmony",
        project.dim    = FALSE,
        max_iter       = 25,
        verbose        = FALSE
    )
    rna_reduction <- "rna.harmony"
}

# 8. RNA UMAP
message("Running RNA UMAP...")
seu_obj <- RunUMAP(seu_obj,
                   dims           = 1:rna_pcs,
                   reduction      = rna_reduction,
                   reduction.name = "umap.rna",
                   reduction.key  = "rnaUMAP_",
                   verbose        = FALSE)

# ── ATAC ──────────────────────────────────────────────────────────────────────

# 9. TF-IDF normalization + top features + LSI/SVD
message("Running TF-IDF + SVD (LSI)...")
DefaultAssay(seu_obj) <- "peaks"
seu_obj <- RunTFIDF(seu_obj, verbose = FALSE)
seu_obj <- FindTopFeatures(seu_obj, min.cutoff = "q50", verbose = FALSE)
seu_obj <- RunSVD(seu_obj, n = atac_pcs, verbose = FALSE)

# 10. Optional Harmony on LSI (ATAC)
#     Exclude first LSI dimension (correlated with sequencing depth)
atac_reduction <- "lsi"
atac_dims      <- 2:atac_pcs
if (run_harmony) {
    message("Running Harmony on ATAC...")
    seu_obj <- RunHarmony(
        seu_obj,
        group.by.vars  = harmony_vars,
        reduction      = "lsi",
        assay.use      = "peaks",
        dims.use       = 2:atac_pcs,
        reduction.save = "atac.harmony",
        project.dim    = FALSE,
        max_iter       = 25,
        verbose        = FALSE
    )
    atac_reduction <- "atac.harmony"
    atac_dims      <- 1:(atac_pcs - 1)
}

# 11. ATAC UMAP
message("Running ATAC UMAP...")
seu_obj <- RunUMAP(seu_obj,
                   dims           = atac_dims,
                   reduction      = atac_reduction,
                   reduction.name = "umap.atac",
                   reduction.key  = "atacUMAP_",
                   verbose        = FALSE)

# ── WNN ───────────────────────────────────────────────────────────────────────

# 12. WNN graph: combines RNA + ATAC modalities
message("Constructing WNN graph...")
seu_obj <- FindMultiModalNeighbors(
    seu_obj,
    reduction.list = list(rna_reduction, atac_reduction),
    dims.list      = list(1:rna_pcs, atac_dims),
    verbose        = FALSE
)

# 13. WNN UMAP (joint embedding)
message("Running WNN UMAP...")
seu_obj <- RunUMAP(seu_obj,
                   nn.name        = "weighted.nn",
                   reduction.name = "wnn.umap",
                   reduction.key  = "wnnUMAP_",
                   verbose        = FALSE)

# 14. Multi-resolution clustering on WNN graph
message("Running FindClusters at resolutions: ", paste(resolutions, collapse = ", "))
for (res in resolutions) {
    seu_obj <- FindClusters(seu_obj,
                            graph.name = "wsnn",
                            algorithm  = 3,
                            resolution = res,
                            verbose    = FALSE)
}

# ── PLOTS ─────────────────────────────────────────────────────────────────────

# 15. Create output directories for UMAP plots
plot_dir <- paste0("output/plots/", argv$project_prefix, "-umaps")
dir.create(file.path(plot_dir, "orig.ident"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(plot_dir, "metadata"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(plot_dir, "clustering"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(plot_dir, "qc-metrics"), recursive = TRUE, showWarnings = FALSE)
message("Plot output directory: ", plot_dir)

# 16. Combined 3-panel orig.ident UMAP (RNA + ATAC + WNN)
message("Generating combined orig.ident UMAP panel...")
p_rna  <- DimPlot(seu_obj, reduction = "umap.rna",  group.by = "orig.ident", label = FALSE) + ggtitle("RNA UMAP")
p_atac <- DimPlot(seu_obj, reduction = "umap.atac", group.by = "orig.ident", label = FALSE) + ggtitle("ATAC UMAP")
p_wnn  <- DimPlot(seu_obj, reduction = "wnn.umap",  group.by = "orig.ident", label = FALSE) + ggtitle("WNN UMAP")
pdf(file.path(plot_dir, paste0(argv$project_prefix, "-all-reduction-umaps.pdf")), width = 24, height = 8)
print(p_rna | p_atac | p_wnn)
dev.off()

# 17. Individual orig.ident UMAPs (RNA, ATAC, WNN)
message("Generating individual orig.ident UMAPs...")
modalities <- list(
    list(name = "rna",  reduction = "umap.rna",  title = "RNA UMAP"),
    list(name = "atac", reduction = "umap.atac", title = "ATAC UMAP"),
    list(name = "wnn",  reduction = "wnn.umap",  title = "WNN UMAP")
)
for (mod in modalities) {
    p <- DimPlot(seu_obj, reduction = mod$reduction, group.by = "orig.ident", label = FALSE) +
        ggtitle(mod$title)
    pdf(file.path(plot_dir, "orig.ident",
                  paste0(argv$project_prefix, "-umap-", mod$name, "-orig.ident.pdf")),
        width = 10, height = 8)
    print(p)
    dev.off()
}

# 18. Metadata UMAPs (WNN UMAP, one per samplesheet metadata column)
#     ss_cols derived from samplesheet already loaded in step 3
ss_cols <- setdiff(colnames(samplesheet), c("SampleID", "path"))
message("Generating metadata UMAPs for: ", paste(ss_cols, collapse = ", "))
for (col in ss_cols) {
    p <- DimPlot(seu_obj, reduction = "wnn.umap", group.by = col, label = FALSE) +
        ggtitle(col)
    pdf(file.path(plot_dir, "metadata",
                  paste0(argv$project_prefix, "-umap-wnn-", col, ".pdf")),
        width = 10, height = 8)
    print(p)
    dev.off()
}

# 19. Clustering UMAPs (WNN UMAP, one per resolution)
cluster_cols <- paste0("wsnn_res.", resolutions)
message("Generating clustering UMAPs at resolutions: ", paste(resolutions, collapse = ", "))
for (i in seq_along(resolutions)) {
    p <- DimPlot(seu_obj, reduction = "wnn.umap", group.by = cluster_cols[i], label = TRUE) +
        ggtitle(paste0("WNN Clusters (res ", resolutions[i], ")"))
    pdf(file.path(plot_dir, "clustering",
                  paste0(argv$project_prefix, "-umap-wnn-res", resolutions[i], ".pdf")),
        width = 10, height = 8)
    print(p)
    dev.off()
}

# 20. QC metrics FeaturePlot on WNN UMAP — 3x2 grid, one plot
message("Generating QC metrics FeaturePlot...")
qc_features <- c("nCount_RNA", "nFeature_RNA", "percent.mt",
                  "nCount_ATAC", "TSS.enrichment", "nucleosome_signal")
p_qc <- FeaturePlot(seu_obj, features = qc_features, reduction = "wnn.umap", ncol = 3)
pdf(file.path(plot_dir, "qc-metrics",
              paste0(argv$project_prefix, "-umap-wnn-qc-metrics.pdf")),
    width = 18, height = 10)
print(p_qc)
dev.off()
message("All UMAP plots saved to: ", plot_dir)

# 21. Save output
out_path <- paste0("output/RDS-files/", argv$project_prefix, "-07-normalize-reduce-obj.RDS")
saveRDS(seu_obj, file = out_path)
message("Step 7 complete. Saved: ", out_path)
