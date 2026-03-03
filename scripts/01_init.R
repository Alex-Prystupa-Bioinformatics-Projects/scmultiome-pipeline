#!/usr/bin/env Rscript
# 01_init.R
# Pipeline Step 1: Initialize directory structure from samplesheet.
#
# Reads the samplesheet and copies raw cellranger-arc output files
# (atac_fragments + GEX filtered matrix) into data/raw/{SampleID}/.
#
# Usage:
#   Rscript scripts/01_init.R configs/samplesheet.csv
#
# Outputs:
#   data/raw/{SampleID}/ populated for each sample in the samplesheet

library(argparser, quietly = TRUE)
library(yaml,      quietly = TRUE)

p <- arg_parser("Step 1: Initialize raw data directories from samplesheet")
p <- add_argument(p, "samplesheet",     help = "Path to samplesheet CSV",              type = "character")
p <- add_argument(p, "--pipeline_config", help = "Path to pipeline_config.yml",
                  default = "configs/pipeline_config.yml", type = "character")

argv <- parse_args(p)

pipeline_config <- yaml::read_yaml(argv$pipeline_config)
samplesheet     <- read.csv(argv$samplesheet)
samples         <- samplesheet$SampleID

source("scripts/functions.R")

message("Starting Step 1: init")
message("Samples to initialize: ", paste(samples, collapse = ", "))

lapply(samples, function(sample) {
    data.dir <- samplesheet$path[samplesheet$SampleID == sample]
    new.dir  <- paste0("data/raw/", sample)

    dir.create(new.dir, recursive = TRUE, showWarnings = FALSE)
    message("Initializing sample: ", sample, " from ", data.dir)

    CombineDirectories(data.dir, new.dir, sample)
})

message("Step 1 complete: all sample directories initialized under data/raw/")
