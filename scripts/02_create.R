#!/usr/bin/env Rscript
# =============================================================================
# 02_create.R  |  Pipeline Step 2: Create Seurat Objects
# =============================================================================
# Creates a multiome Seurat object (RNA + ATAC assays) for each sample in the
# samplesheet. Adds metadata columns from the samplesheet and computes:
#
#   1. Basic QC metrics: percent.mt, NucleosomeSignal, TSSEnrichment
#      - These are fragment-based metrics computed on the ATAC assay
#      - Values are independent of the peak set used later in step 3
#
#   2. Doublet detection via scDblFinder (RNA only)
#      - Run per sample on raw RNA counts before any normalization
#      - Adds scDblFinder.score (0-1 doublet probability) and
#        scDblFinder.class ("singlet" / "doublet") to cell metadata
#      - Doublet filtering is optional and controlled in 05_filter.R
#        via the remove_doublets flag in configs/qc_config.yml
#
# Usage:
#   Rscript scripts/02_create.R configs/samplesheet.csv \
#       --project_prefix myproject
#
# Input:  data/raw/{SampleID}/  (output of 01_init.R)
# Output: output/RDS-files/{project_prefix}-create-obj-list.RDS
# =============================================================================

library(argparser,           quietly = TRUE)
library(yaml,                quietly = TRUE)
library(future,              quietly = TRUE)
library(future.apply,        quietly = TRUE)
library(Seurat,              quietly = TRUE)
library(Signac,              quietly = TRUE)
library(scDblFinder,         quietly = TRUE)
library(SingleCellExperiment,quietly = TRUE)

# 1. Parse arguments
p <- arg_parser("Step 2: Create Seurat objects, compute QC metrics, detect doublets")
p <- add_argument(p, "samplesheet",
                  help = "Path to samplesheet CSV", type = "character")
p <- add_argument(p, "--pipeline_config",
                  help    = "Path to pipeline_config.yml",
                  default = "configs/pipeline_config.yml", type = "character")
p <- add_argument(p, "--project_prefix",
                  help    = "Output file prefix (no .RDS extension)",
                  default = "multiome", type = "character")
argv <- parse_args(p)

# 2. Load config, samplesheet, and helper functions
pipeline_config <- yaml::read_yaml(argv$pipeline_config)
samplesheet     <- read.csv(argv$samplesheet)
samples         <- samplesheet$SampleID

source("scripts/functions.R")
source("scripts/genome_utils.R")

# 3. Set up parallelization
options(future.globals.maxSize = Inf)
options(future.rng.onMisuse   = "ignore")
plan("multicore", workers = future::availableCores())

# 4. Load genome resources (annotation, blacklist, mt pattern) based on species
genome <- load_genome(pipeline_config)

dir.create("output/RDS-files", recursive = TRUE, showWarnings = FALSE)

message("Starting Step 2: create | Samples: ", paste(samples, collapse = ", "))

# 5. Create a Seurat object for each sample
#    - Creates RNA and ATAC assays from the 10X cellranger output
#    - Appends any extra metadata columns defined in the samplesheet
#    - Computes percent.mt, NucleosomeSignal, TSSEnrichment on the ATAC assay
#    - Runs scDblFinder on raw RNA counts to flag putative doublets
seu_obj_list <- lapply(samples, function(sample) {

    message("  Creating object for: ", sample)
    seu_obj <- CreateMultiomeSeurat(
        data.dir      = paste0("data/raw/", sample),
        my.annotation = genome$annotation
    )
    seu_obj@project.name <- sample
    seu_obj$orig.ident   <- sample

    # Add any extra metadata columns from the samplesheet (columns after 'path')
    path.col <- grep("^path$", colnames(samplesheet))
    if (path.col < ncol(samplesheet)) {
        for (md.col in colnames(samplesheet)[(path.col + 1):ncol(samplesheet)]) {
            seu_obj[[md.col]] <- samplesheet[samplesheet$SampleID == sample, md.col]
        }
    }

    # QC metrics — fragment-based, independent of peak set
    seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = genome$mt_pattern)
    seu_obj <- NucleosomeSignal(seu_obj, assay = "ATAC")
    seu_obj <- TSSEnrichment(seu_obj,    assay = "ATAC")

    # Doublet detection on raw RNA counts
    # scDblFinder requires a SingleCellExperiment with raw counts
    # Results are pulled back onto the Seurat object as metadata
    sce     <- as.SingleCellExperiment(seu_obj)
    sce     <- scDblFinder(sce)
    seu_obj$scDblFinder.score <- sce$scDblFinder.score
    seu_obj$scDblFinder.class <- sce$scDblFinder.class

    n_doublets <- sum(seu_obj$scDblFinder.class == "doublet")
    message("  Done: ", sample,
            " | Cells: ",    ncol(seu_obj),
            " | Doublets: ", n_doublets,
            " (", round(n_doublets / ncol(seu_obj) * 100, 1), "%)")
    return(seu_obj)
})
names(seu_obj_list) <- samples

# 6. Save output
out_path <- paste0("output/RDS-files/", argv$project_prefix, "-02-create-obj-list.RDS")
saveRDS(seu_obj_list, file = out_path)
message("Step 2 complete. Saved: ", out_path)
