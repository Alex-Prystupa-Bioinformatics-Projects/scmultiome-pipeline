#!/usr/bin/env Rscript
# =============================================================================
# 05_filter.R  |  Pipeline Step 5: Apply QC Filters
# =============================================================================
# Applies per-sample QC thresholds from configs/qc_config.yml to each Seurat
# object. Optionally removes predicted doublets if remove_doublets is set to
# TRUE in qc_config.yml.
#
# Run this after reviewing QC plots from 04_qc.R and editing qc_config.yml
# with your final thresholds.
#
# Usage:
#   Rscript scripts/05_filter.R configs/samplesheet.csv \
#       --project_prefix myproject \
#       --RDS_file_in output/RDS-files/myproject-03-callpeaks-obj-list.RDS
#
# Input:  output of 03_callpeaks.R  ({prefix}-03-callpeaks-obj-list.RDS)
#         configs/qc_config.yml     (written by 04_qc.R, edited by user)
# Output: output/RDS-files/{prefix}-05-filter-obj-list.RDS
# =============================================================================

library(argparser, quietly = TRUE)
library(yaml,      quietly = TRUE)
library(Seurat,    quietly = TRUE)
library(Signac,    quietly = TRUE)

# 1. Parse arguments
p <- arg_parser("Step 5: Apply QC filters from qc_config.yml")
p <- add_argument(p, "samplesheet",
                  help = "Path to samplesheet CSV", type = "character")
p <- add_argument(p, "--pipeline_config",
                  help    = "Path to pipeline_config.yml",
                  default = "configs/pipeline_config.yml", type = "character")
p <- add_argument(p, "--project_prefix",
                  help    = "Output file prefix",
                  default = "multiome", type = "character")
p <- add_argument(p, "--RDS_file_in",
                  help = "Input RDS file — output of 03_callpeaks.R", type = "character")
p <- add_argument(p, "--qc_config",
                  help    = "Path to qc_config.yml",
                  default = "configs/qc_config.yml", type = "character")
argv <- parse_args(p)

if (is.na(argv$RDS_file_in) || argv$RDS_file_in == "NA") {
    stop("--RDS_file_in is required. Provide the path to output of 03_callpeaks.R.")
}

# 2. Load config and helper functions
pipeline_config <- yaml::read_yaml(argv$pipeline_config)
samplesheet     <- read.csv(argv$samplesheet)

source("scripts/utils.R")

# 3. Load input objects from previous step
message("Starting Step 5: filter")
message("Loading: ", argv$RDS_file_in)
seu_obj_list <- readRDS(argv$RDS_file_in)

# 4. Load QC thresholds from qc_config.yml
message("Loading QC config: ", argv$qc_config)
qc_configs <- yaml::read_yaml(argv$qc_config)

# 5. Apply per-sample QC filters
#    filter_seu_list_by_qc() applies cell-level thresholds and optionally
#    removes predicted doublets (remove_doublets flag in qc_config.yml)
message("Applying QC filters...")
seu_obj_list <- filter_seu_list_by_qc(seu_obj_list, qc_configs)

# 6. Save filtered object list
out_path <- paste0("output/RDS-files/", argv$project_prefix, "-05-filter-obj-list.RDS")
saveRDS(seu_obj_list, file = out_path)
message("Step 5 complete. Saved: ", out_path)
