#!/bin/bash
# =============================================================================
# 04_linkpeaks.sh  |  Phase 4: Per-Group Parallel LinkPeaks
# =============================================================================
# Runs step 11: splits the annotated merged object by the user-chosen
# grouping column and submits one LSF job per group to run LinkPeaks in
# parallel via scripts/helpers/run_linkpeaks_worker.R.
#
# The grouping_col is flexible — set my_linkpeaks_grouping_col in runmultiome:
#   orig.ident        → per sample
#   Annotation_Broad  → per cell type
#   Sample_Celltype   → per sample × cell type
#
# Prerequisites:
#   1. Run routes/03_annotate_recall_peaks.sh
#   2. Set my_linkpeaks_grouping_col in run/runmultiome
#
# Output: output/RDS/linkpeaks/{group}-linkpeaks.RDS  (one per group, via child jobs)
#         output/RDS-files/{prefix}-11-linkpeaks-obj.RDS  (merged, via merge job)
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
export RENV_PATHS_LIBRARY="$(pwd)/renv/library"
export RENV_PATHS_PREFIX="$(cut -d'/' -f1 .renv_platform)"

# 3. Step 11: Submit per-group LinkPeaks jobs
echo "[Step 11] Splitting object and submitting per-group LinkPeaks jobs..."
Rscript scripts/11_linkpeaks.R \
    configs/samplesheet.csv \
    --project_prefix  $project_prefix \
    --pipeline_config configs/pipeline_config.yml \
    --RDS_file_in     output/RDS-files/${project_prefix}-10-recall-peaks-obj.RDS

echo ""
echo "linkpeaks orchestration complete."
echo "  Per-group jobs submitted — monitor with: bjobs -J 'linkpeaks_*'"
echo "  Logs:    scripts/outs/linkpeaks/"
echo "  Output:  output/RDS/linkpeaks/ (per-group), output/RDS-files/ (merged)"
