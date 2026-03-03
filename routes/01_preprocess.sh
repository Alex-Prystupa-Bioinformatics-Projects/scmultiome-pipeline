#!/bin/bash
# =============================================================================
# 01_preprocess.sh  |  Phase 1: Full Preprocessing Pipeline (Steps 1-4)
# =============================================================================
# Runs the complete preprocessing pipeline in sequence:
#   Step 1 (init)      — copies raw cellranger data into data/raw/
#   Step 2 (create)    — creates Seurat objects, computes QC metrics, detects doublets
#   Step 3 (callpeaks) — calls peaks with MACS2, creates peaks assay
#   Step 4 (qc)        — generates QC plots, doublet summary, and qc_config.yml
#
# After this script completes:
#   1. Review output/plots/{project_prefix}-qc-plots.pdf
#   2. Review output/tables/{project_prefix}-doublet-summary.csv
#   3. Edit configs/qc_config.yml with your per-sample thresholds
#   4. Run routes/02_filter.sh
# =============================================================================

set -e

# 1. Load appropriate R module based on scheduler
if [ "$SCHEDULER" == "slurm" ]; then
    module load R/4.4.1
    module load macs2
elif [ "$SCHEDULER" == "lsf" ]; then
    module load R/4.4.1
else
    echo "No scheduler set, assuming modules already loaded."
fi

# 2. Set renv library path explicitly to avoid platform mismatch errors.
#
#    Background: renv organizes its package library by OS + R version, e.g.:
#      renv/library/linux-rocky-9.7/R-4.4/x86_64-pc-linux-gnu
#
#    On Minerva, login nodes and compute nodes run slightly different OS minor
#    versions (e.g. rocky-9.7 vs rocky-9.6). renv detects the platform at
#    runtime, so if the compute node returns a different string than what was
#    used when renv::restore() was run on the login node, renv cannot find its
#    library and tries to re-bootstrap from the internet — which fails because
#    compute nodes have no outbound internet access.
#
#    Fix: dynamically resolve the platform prefix on the compute node itself
#    (using --no-init-file to bypass .Rprofile) and set RENV_PATHS_LIBRARY to
#    the exact path. This guarantees renv always finds the correct library
#    regardless of which OS minor version the compute node is running.
export RENV_PATHS_LIBRARY="$(pwd)/renv/library/$(Rscript --no-init-file -e 'cat(renv:::renv_platform_prefix())')"

# 3. Step 1: Initialize raw data directories from samplesheet
echo "[Step 1/4] Initializing raw data directories..."
Rscript scripts/01_init.R \
    configs/samplesheet.csv \
    --pipeline_config configs/pipeline_config.yml

# 4. Step 2: Create Seurat objects, compute QC metrics, detect doublets
echo "[Step 2/4] Creating Seurat objects..."
Rscript scripts/02_create.R \
    configs/samplesheet.csv \
    --pipeline_config configs/pipeline_config.yml \
    --project_prefix  $project_prefix

# 5. Step 3: Call peaks with MACS2
echo "[Step 3/4] Calling peaks with MACS2..."
Rscript scripts/03_callpeaks.R \
    configs/samplesheet.csv \
    --pipeline_config configs/pipeline_config.yml \
    --project_prefix  $project_prefix \
    --RDS_file_in     output/RDS-files/${project_prefix}-02-create-obj-list.RDS \
    --macs_path       $my_macs_path

# 6. Step 4: Generate QC plots and write qc_config.yml
echo "[Step 4/4] Generating QC plots..."
Rscript scripts/04_qc.R \
    configs/samplesheet.csv \
    --pipeline_config configs/pipeline_config.yml \
    --project_prefix  $project_prefix \
    --RDS_file_in     output/RDS-files/${project_prefix}-03-callpeaks-obj-list.RDS

echo ""
echo "Preprocessing complete. Next steps:"
echo "  1. Review output/plots/${project_prefix}-qc-plots.pdf"
echo "  2. Review output/tables/${project_prefix}-doublet-summary.csv"
echo "  3. Edit configs/qc_config.yml with your thresholds"
echo "  4. Run routes/02_filter.sh"
