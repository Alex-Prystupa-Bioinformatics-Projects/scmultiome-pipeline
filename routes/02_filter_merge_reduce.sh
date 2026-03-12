#!/bin/bash
# =============================================================================
# 02_filter_merge_reduce.sh  |  Phase 2: Filter, Merge, and Dimensionality Reduction
# =============================================================================
# Runs the filter, merge, and normalization pipeline in sequence:
#   Step 5 (filter)            — applies QC thresholds from qc_config.yml
#   Step 6 (merge)             — consensus peaks + SCTransform per sample +
#                                SelectIntegrationFeatures + merge
#   Step 7 (normalize_reduce)  — JoinLayers, PCA, LSI, Harmony, WNN, UMAP, clustering
#
# Prerequisites:
#   1. Run routes/01_preprocess.sh
#   2. Review output/plots/{project_prefix}-qc/
#   3. Edit configs/qc_config.yml with final per-sample thresholds
#
# Configure batch correction below:
#   harmony_vars — comma-separated metadata column(s) to batch correct on
#                  e.g. "Institution" or "Institution,Sequencing_Date"
#                  Set to "none" to skip Harmony
# =============================================================================

set -e

# 1. Load appropriate R module based on scheduler
if [ "$SCHEDULER" == "slurm" ]; then
    module load R/4.4.1
elif [ "$SCHEDULER" == "lsf" ]; then
    module load R/4.4.1
else
    echo "No scheduler set, assuming modules already loaded."
fi

# 2. Set renv library path to avoid platform mismatch on compute nodes
#    See routes/01_preprocess.sh for full explanation of this fix.
export RENV_PATHS_LIBRARY="$(pwd)/renv/library"
export RENV_PATHS_PREFIX="$(cut -d'/' -f1 .renv_platform)"

# 3. Steps 5-6: Filter and merge (skipped if merged RDS already exists)
MERGE_RDS="output/RDS-files/${project_prefix}-06-merge-obj.RDS"

if [ -f "$MERGE_RDS" ]; then
    echo "[Steps 5-6] Skipping — merged RDS already exists: $MERGE_RDS"
    echo "            Delete it to force a full re-filter/re-merge."
else
    echo "[Step 5/3] Applying QC filters..."
    Rscript scripts/05_filter.R \
        configs/samplesheet.csv \
        --pipeline_config configs/pipeline_config.yml \
        --project_prefix  $project_prefix \
        --qc_config       configs/qc_config.yml \
        --RDS_file_in     output/RDS-files/${project_prefix}-03-callpeaks-obj-list.RDS

    echo "[Step 6/3] Building consensus peaks and merging objects..."
    Rscript scripts/06_merge.R \
        configs/samplesheet.csv \
        --pipeline_config configs/pipeline_config.yml \
        --project_prefix  $project_prefix \
        --RDS_file_in     output/RDS-files/${project_prefix}-05-filter-obj-list.RDS
fi

# 5. Step 7: Normalize, dimensionality reduction, clustering
echo "[Step 7/3] Normalizing, reducing, and clustering..."
Rscript scripts/07_normalize_reduce.R \
    configs/samplesheet.csv \
    --pipeline_config configs/pipeline_config.yml \
    --project_prefix  $project_prefix \
    --harmony_vars    $(Rscript --no-init-file -e "cat(yaml::read_yaml('configs/pipeline_config.yml')[['harmony_vars']])") \
    --RDS_file_in     output/RDS-files/${project_prefix}-06-merge-obj.RDS

# 6. Step 8: Cluster marker detection (runs at all wsnn_res.* resolutions)
echo "[Step 8/4] Running cluster marker detection..."
Rscript scripts/08_cluster_markers.R \
    configs/samplesheet.csv \
    --pipeline_config configs/pipeline_config.yml \
    --project_prefix  $project_prefix \
    --RDS_file_in     output/RDS-files/${project_prefix}-07-normalize-reduce-obj.RDS

echo ""
echo "filter_merge_reduce complete."
echo "  Final object: output/RDS-files/${project_prefix}-07-normalize-reduce-obj.RDS"
echo "  Markers:      output/markers/${project_prefix}-08-markers/"
echo "  Next steps:"
echo "    1. Inspect WNN UMAP and clustering"
echo "    2. Review cluster markers and annotate cell types"
