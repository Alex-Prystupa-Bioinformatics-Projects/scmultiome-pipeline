#!/bin/bash
# =============================================================================
# 03_annotate_recall_peaks.sh  |  Phase 3: Annotate + Recall Peaks
# =============================================================================
# Runs step 9 (annotate) then step 10 (cell-type-aware peak recall).
#
# Step 9: applies user-defined cell type annotations from configs/annotations.csv
#         onto the Seurat object and removes Discard cells.
# Step 10: splits the annotated object by sample × cell type, recalls MACS2
#          peaks per group, builds a consensus peak set, and re-merges.
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

# 4. Step 10: Cell-type-aware peak recall
echo "[Step 10] Recalling peaks per sample × cell type..."
Rscript scripts/10_recall_peaks.R \
    configs/samplesheet.csv \
    --project_prefix  $project_prefix \
    --RDS_file_in     output/RDS-files/${project_prefix}-09-annotate-obj.RDS \
    --annotation_file configs/annotations.csv \
    --macs2_path      $my_macs_path

echo ""
echo "annotate_recall_peaks complete."
echo "  Annotated object:      output/RDS-files/${project_prefix}-09-annotate-obj.RDS"
echo "  Recall-peaks object:   output/RDS-files/${project_prefix}-10-recall-peaks-obj.RDS"
