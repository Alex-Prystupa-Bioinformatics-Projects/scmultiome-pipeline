#!/bin/bash
# =============================================================================
# job_merge.sh  |  LinkPeaks Merge Job (static, submitted by 11_linkpeaks.R)
# =============================================================================
# Submitted after all per-group linkpeaks child jobs. Waits for done(linkpeaks_*)
# then runs merge_linkpeaks.R to combine per-group links back onto the full object.
#
# Dependency: all LSF jobs named linkpeaks_* must finish successfully.
# Output:     output/RDS-files/{prefix}-11-linkpeaks-obj.RDS
# =============================================================================
#BSUB -J linkpeaks_merge
#BSUB -q premium
#BSUB -P acc_naiklab
#BSUB -w "done(linkpeaks_*)"
#BSUB -o scripts/outs/linkpeaks/merge-%J.out
#BSUB -e scripts/outs/linkpeaks/merge-%J.err
#BSUB -n 8
#BSUB -R "rusage[mem=64000]"
#BSUB -R "span[hosts=1]"
#BSUB -W 6:00

module load R/4.4.1

export RENV_PATHS_LIBRARY="$(pwd)/renv/library"
export RENV_PATHS_PREFIX="$(cut -d'/' -f1 .renv_platform)"

Rscript scripts/helpers/merge_linkpeaks.R \
    "$RDS_FILE_IN" \
    output/RDS/linkpeaks \
    configs/pipeline_config.yml
