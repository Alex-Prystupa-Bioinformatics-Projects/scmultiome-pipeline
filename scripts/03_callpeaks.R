#!/usr/bin/env Rscript
# =============================================================================
# 03_callpeaks.R  |  Pipeline Step 3: Call Peaks with MACS2
# =============================================================================
# Runs MACS2 peak calling on each sample object and creates a new "peaks" assay.
# This replaces the 10X cellranger peak set with a more sensitive, custom one.
#
# Optionally groups peak calling by a metadata variable (e.g. cell type
# annotation) to generate cell-type-aware peaks.
#
# Usage:
#   Rscript scripts/03_callpeaks.R configs/samplesheet.csv \
#       --project_prefix myproject \
#       --RDS_file_in output/RDS-files/myproject-create-obj-list.RDS \
#       --macs_path /path/to/macs2
#
# Optional:
#   --grouping_var  metadata column to group peak calling by (e.g. cell.lineage)
#
# Input:  output of 02_create.R  (--RDS_file_in)
# Output: output/RDS-files/{project_prefix}-03-callpeaks-obj-list.RDS
# =============================================================================

library(argparser,    quietly = TRUE)
library(yaml,         quietly = TRUE)
library(future,       quietly = TRUE)
library(future.apply, quietly = TRUE)
library(Seurat,       quietly = TRUE)
library(Signac,       quietly = TRUE)

# 1. Parse arguments
p <- arg_parser("Step 3: Call peaks with MACS2")
p <- add_argument(p, "samplesheet",
                  help = "Path to samplesheet CSV", type = "character")
p <- add_argument(p, "--pipeline_config",
                  help    = "Path to pipeline_config.yml",
                  default = "configs/pipeline_config.yml", type = "character")
p <- add_argument(p, "--project_prefix",
                  help    = "Output file prefix (no .RDS extension)",
                  default = "multiome", type = "character")
p <- add_argument(p, "--RDS_file_in",
                  help = "Input RDS file — output of 02_create.R", type = "character")
p <- add_argument(p, "--macs_path",
                  help = "Path to MACS2 executable", type = "character")
p <- add_argument(p, "--grouping_var",
                  help    = "Metadata column to group peak calling by (optional)",
                  default = "NA", type = "character")
argv <- parse_args(p)

if (is.na(argv$RDS_file_in) || argv$RDS_file_in == "NA") {
    stop("--RDS_file_in is required. Provide the path to output of 02_create.R.")
}
if (is.na(argv$macs_path) || argv$macs_path == "NA") {
    stop("--macs_path is required. Provide the path to the MACS2 executable.")
}

# 2. Load config, samplesheet, and helper functions
pipeline_config <- yaml::read_yaml(argv$pipeline_config)
samplesheet     <- read.csv(argv$samplesheet)

source("scripts/functions.R")
source("scripts/genome_utils.R")

# 3. Set up parallelization
options(future.globals.maxSize = Inf)
options(future.rng.onMisuse   = "ignore")
plan("multicore", workers = future::availableCores())

# 4. Load genome resources based on species in pipeline_config
genome <- load_genome(pipeline_config)

# 5. Load input objects from previous step
message("Starting Step 3: callpeaks")
message("Loading: ", argv$RDS_file_in)
seu_obj_list <- readRDS(argv$RDS_file_in)

grouping_var <- if (argv$grouping_var == "NA") NULL else argv$grouping_var

dir.create("data/raw/macs-peaks", recursive = TRUE, showWarnings = FALSE)
dir.create("output/RDS-files",    recursive = TRUE, showWarnings = FALSE)

# 6. Run MACS2 peak calling for each sample
#    - Calls peaks using MACS2 (optionally grouped by a metadata variable)
#    - Creates a new "peaks" assay on each object with the custom peak set
#    - The original "ATAC" assay (10X peaks) remains on the object
seu_obj_list <- lapply(seu_obj_list, function(seu_obj) {
    message("  Calling peaks for: ", seu_obj@project.name)

    seu_obj <- CallMyPeaks(
        seu_obj,
        grouping.var  = grouping_var,
        my.macs2.path = argv$macs_path,
        my.annotation = genome$annotation,
        my.blacklist  = genome$blacklist
    )

    message("  Done: ", seu_obj@project.name, " | Peaks: ", nrow(seu_obj[["peaks"]]))
    return(seu_obj)
})

# 7. Save output
out_path <- paste0("output/RDS-files/", argv$project_prefix, "-03-callpeaks-obj-list.RDS")
saveRDS(seu_obj_list, file = out_path)
message("Step 3 complete. Saved: ", out_path)
