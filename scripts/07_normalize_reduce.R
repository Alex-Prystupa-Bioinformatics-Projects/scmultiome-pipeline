#!/usr/bin/env Rscript
# =============================================================================
# 07_normalize_reduce.R  |  Pipeline Step 7: Normalize, Reduce, Cluster
# =============================================================================
# Takes the merged object from step 6 and runs the full normalization,
# dimensionality reduction, and clustering pipeline.
#
# RNA:  SCTransform (per layer in Seurat v5) → consensus variable features →
#       JoinLayers → PCA → optional Harmony → UMAP
# ATAC: TF-IDF → FindTopFeatures → SVD/LSI → optional Harmony → UMAP
# WNN:  FindMultiModalNeighbors → WNN UMAP → multi-resolution FindClusters
#
# Consensus variable features: after SCTransform runs per layer (per sample),
# the intersection of per-layer variable genes is used for PCA to avoid
# PCA being driven by sample-specific variable genes.
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

library(argparser, quietly = TRUE)
library(yaml,      quietly = TRUE)
library(Seurat,    quietly = TRUE)
library(Signac,    quietly = TRUE)
library(harmony,   quietly = TRUE)

# 1. Parse arguments
p <- arg_parser("Step 7: Normalize, dimensionality reduction, and clustering")
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
message("Starting Step 7: normalize + reduce")
message("Loading: ", argv$RDS_file_in)
seu_obj <- readRDS(argv$RDS_file_in)

# ── RNA ───────────────────────────────────────────────────────────────────────

# 5. SCTransform on merged object
#    In Seurat v5, runs per layer (per sample) automatically
message("Running SCTransform (per layer)...")
DefaultAssay(seu_obj) <- "RNA"
seu_obj <- SCTransform(seu_obj, verbose = FALSE)

# 6. Compute consensus variable features across layers
#    Intersect per-layer variable genes so PCA uses only cross-sample variable features
message("Computing consensus variable features...")
sct_layers        <- Layers(seu_obj[["SCT"]])
consensus_var_features <- Reduce(intersect, lapply(sct_layers, function(l) {
    VariableFeatures(seu_obj[["SCT"]], layer = l)
}))
message("  Consensus variable features: ", length(consensus_var_features))
VariableFeatures(seu_obj[["SCT"]]) <- consensus_var_features

# 7. JoinLayers after SCTransform
message("Joining layers...")
seu_obj <- JoinLayers(seu_obj)

# 8. PCA using consensus variable features
message("Running PCA...")
seu_obj <- RunPCA(seu_obj, npcs = rna_pcs, verbose = FALSE)

# 9. Optional Harmony on PCA (RNA)
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

# 10. RNA UMAP
message("Running RNA UMAP...")
seu_obj <- RunUMAP(seu_obj,
                   dims           = 1:rna_pcs,
                   reduction      = rna_reduction,
                   reduction.name = "umap.rna",
                   reduction.key  = "rnaUMAP_",
                   verbose        = FALSE)

# ── ATAC ──────────────────────────────────────────────────────────────────────

# 11. TF-IDF normalization + top features + LSI/SVD
message("Running TF-IDF + SVD (LSI)...")
DefaultAssay(seu_obj) <- "peaks"
seu_obj <- RunTFIDF(seu_obj, verbose = FALSE)
seu_obj <- FindTopFeatures(seu_obj, min.cutoff = "q50", verbose = FALSE)
seu_obj <- RunSVD(seu_obj, n = atac_pcs, verbose = FALSE)

# 12. Optional Harmony on LSI (ATAC)
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

# 13. ATAC UMAP
message("Running ATAC UMAP...")
seu_obj <- RunUMAP(seu_obj,
                   dims           = atac_dims,
                   reduction      = atac_reduction,
                   reduction.name = "umap.atac",
                   reduction.key  = "atacUMAP_",
                   verbose        = FALSE)

# ── WNN ───────────────────────────────────────────────────────────────────────

# 14. WNN graph: combines RNA + ATAC modalities
message("Constructing WNN graph...")
seu_obj <- FindMultiModalNeighbors(
    seu_obj,
    reduction.list = list(rna_reduction, atac_reduction),
    dims.list      = list(1:rna_pcs, atac_dims),
    verbose        = FALSE
)

# 15. WNN UMAP (joint embedding)
message("Running WNN UMAP...")
seu_obj <- RunUMAP(seu_obj,
                   nn.name        = "weighted.nn",
                   reduction.name = "wnn.umap",
                   reduction.key  = "wnnUMAP_",
                   verbose        = FALSE)

# 16. Multi-resolution clustering on WNN graph
message("Running FindClusters at resolutions: ", paste(resolutions, collapse = ", "))
for (res in resolutions) {
    seu_obj <- FindClusters(seu_obj,
                            graph.name = "wsnn",
                            algorithm  = 3,
                            resolution = res,
                            verbose    = FALSE)
}

# 17. Save output
out_path <- paste0("output/RDS-files/", argv$project_prefix, "-07-normalize-reduce-obj.RDS")
saveRDS(seu_obj, file = out_path)
message("Step 7 complete. Saved: ", out_path)
