#!/bin/bash
# =============================================================================
# setup_dirs.sh  |  Initialize Project Directory Structure
# =============================================================================
# Run once when setting up a new project. Creates all required directories
# and initializes config files with the correct headers.
#
# Usage:
#   bash run/runmultiome init
# =============================================================================

# 1. Create all required directories
mkdir -p data/raw \
         data/raw/macs-peaks \
         data/scenicplus/cisTarget_dbs \
         scripts/outs \
         output/RDS-files \
         output/plots \
         output/tables \
         configs

# 2. Ensure all pipeline scripts are executable
chmod +x routes/*.sh run/runmultiome

# 3. Initialize samplesheet with required column headers
#    Add one row per sample: SampleID, path, plus any metadata columns
echo "SampleID,path" > configs/samplesheet.csv

# 4. Remind user to configure pipeline_config.yml before running
echo "Edit configs/pipeline_config.yml to set your species and genome before running."

echo "Project directory structure initialized."
echo ""
echo "Before running the pipeline, restore the R environment:"
echo "  module load R/4.4.1 && Rscript -e \"renv::restore()\""
echo ""
echo "Next steps:"
echo "  1. Fill in configs/samplesheet.csv with your samples"
echo "  2. Edit configs/pipeline_config.yml with your species/genome settings"
echo "  3. Run: bash run/runmultiome preprocess"
