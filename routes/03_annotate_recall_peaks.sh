#!/bin/bash
# =============================================================================
# 03_annotate_recall_peaks.sh  |  Phase 3: Annotate Cell Types
# =============================================================================
# Runs step 9 (annotate) — applies user-defined cell type annotations from
# configs/annotations.csv onto the Seurat object and removes Discard cells.
#
# Note: recall_peaks logic will be added to this route in a future step.
#
# Prerequisites:
#   1. Run routes/02_filter_merge_reduce.sh
#   2. Review output/reports/{project_prefix}-clustering-report.pptx and
#      output/markers/{project_prefix}-08-markers/ to identify cluster identities
#   3. Fill in configs/annotations.csv:
#        - Rename first column to your chosen resolution (e.g. wsnn_res.0.4)
#        - Fill in Annotation_Broad, Annotation_Narrow (or your own column names)
#        - Label any clusters to exclude as "Discard"
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

# 3. Step 9: Annotate cell types
echo "[Step 9] Applying cell type annotations..."
Rscript scripts/09_annotate.R \
    configs/samplesheet.csv \
    --project_prefix  $project_prefix \
    --RDS_file_in     output/RDS-files/${project_prefix}-07-normalize-reduce-obj.RDS \
    --annotation_file configs/annotations.csv

echo ""
echo "annotate_recall_peaks complete."
echo "  Annotated object: output/RDS-files/${project_prefix}-09-annotate-obj.RDS"
