#!/usr/bin/env Rscript
# =============================================================================
# 11_linkpeaks.R  |  Pipeline Step 11: Per-Group Parallel LinkPeaks
# =============================================================================
# Splits the annotated merged object by a metadata column defined in
# pipeline_config.yml (linkpeaks_grouping_col) and submits one LSF job per
# group to run LinkPeaks in parallel via scripts/helpers/run_linkpeaks_worker.R.
# After all child jobs are submitted, a merge job is submitted with LSF
# dependency (done) on all child job IDs — runs merge_linkpeaks.R to
# re-assemble results back onto the unsplit step-10 object.
#
# Set linkpeaks_grouping_col in configs/pipeline_config.yml:
#   none                → run on full object (one job)
#   orig.ident          → one job per sample
#   Annotation_Broad    → one job per cell type
#   Sample_Celltype     → one job per sample × cell type combination
#
# Usage:
#   Rscript scripts/11_linkpeaks.R configs/samplesheet.csv \
#       --project_prefix  myproject \
#       --RDS_file_in     output/RDS-files/myproject-10-recall-peaks-obj.RDS
#
# Input:  output of 10_recall_peaks.R (--RDS_file_in)
# Output: output/RDS/linkpeaks/{group}-linkpeaks.RDS  (one per group, via LSF jobs)
#         output/RDS-files/{prefix}-11-linkpeaks-obj.RDS  (merged, via merge job)
# =============================================================================

library(argparser, quietly = TRUE)
library(yaml,      quietly = TRUE)
library(Seurat,    quietly = TRUE)

# Hardcoded LSF job settings
LSF_QUEUE <- "premium"
LSF_PROJ  <- "acc_naiklab"
LSF_CORES <- 16
LSF_MEM   <- 64000
LSF_TIME  <- "36:00"
DISTANCE  <- 250001

# -----------------------------------------------------------------------------
# 1. Arguments
# -----------------------------------------------------------------------------
p <- arg_parser("Step 11: Per-Group Parallel LinkPeaks")
p <- add_argument(p, "samplesheet",
                  help = "Path to samplesheet CSV", type = "character")
p <- add_argument(p, "--pipeline_config",
                  help    = "Path to pipeline_config.yml",
                  default = "configs/pipeline_config.yml", type = "character")
p <- add_argument(p, "--project_prefix",
                  help    = "Output file prefix",
                  default = "multiome", type = "character")
p <- add_argument(p, "--RDS_file_in",
                  help = "Input RDS — output of 10_recall_peaks.R", type = "character")
argv <- parse_args(p)

if (is.na(argv$RDS_file_in) || argv$RDS_file_in == "NA") {
    stop("--RDS_file_in is required.")
}
if (!file.exists(argv$RDS_file_in)) {
    stop("RDS file not found: ", argv$RDS_file_in)
}

# -----------------------------------------------------------------------------
# 2. Load config and read grouping_col
# -----------------------------------------------------------------------------
pipeline_config <- yaml::read_yaml(argv$pipeline_config)

grouping_col <- pipeline_config$linkpeaks_grouping_col
if (is.null(grouping_col)) grouping_col <- "none"

message("Starting Step 11: linkpeaks")
message("linkpeaks_grouping_col: ", grouping_col)

# -----------------------------------------------------------------------------
# 3. Load merged object
# -----------------------------------------------------------------------------
message("Loading: ", argv$RDS_file_in)
seu_obj <- readRDS(argv$RDS_file_in)
message("Object loaded — ", ncol(seu_obj), " cells")

# -----------------------------------------------------------------------------
# 4. Build named list of groups to process
# -----------------------------------------------------------------------------
if (grouping_col == "none") {
    # Run on full object — one job
    message("grouping_col = none: submitting one job for full merged object")
    seu_obj_list <- list("full" = seu_obj)
} else {
    # Validate column exists
    if (!grouping_col %in% colnames(seu_obj@meta.data)) {
        stop("linkpeaks_grouping_col '", grouping_col, "' not found in object metadata.\n",
             "Available columns: ", paste(colnames(seu_obj@meta.data), collapse = ", "))
    }
    groups   <- sort(unique(seu_obj[[grouping_col, drop = TRUE]]))
    n_groups <- length(groups)
    message("Groups (", n_groups, "): ", paste(groups, collapse = ", "))
    seu_obj_list <- SplitObject(seu_obj, split.by = grouping_col)
}

# -----------------------------------------------------------------------------
# 5. Create output directories and save per-group RDS files
# -----------------------------------------------------------------------------
split_dir <- "output/RDS/linkpeaks/split"
out_dir   <- "output/RDS/linkpeaks"
jobs_dir  <- "scripts/jobs/linkpeaks"
logs_dir  <- "scripts/outs/linkpeaks"

dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(jobs_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir,  recursive = TRUE, showWarnings = FALSE)

message("Saving per-group RDS files to: ", split_dir)
for (grp in names(seu_obj_list)) {
    safe_grp <- gsub("[^A-Za-z0-9._-]", "_", grp)
    grp_rds  <- file.path(split_dir, paste0(safe_grp, ".RDS"))
    saveRDS(seu_obj_list[[grp]], file = grp_rds)
    message("  Saved: ", grp_rds, " (", ncol(seu_obj_list[[grp]]), " cells)")
}

# -----------------------------------------------------------------------------
# 6. Read renv platform prefix for compute nodes
# -----------------------------------------------------------------------------
renv_platform <- ""
if (file.exists(".renv_platform")) {
    platform_raw  <- readLines(".renv_platform", n = 1)
    renv_platform <- strsplit(platform_raw, "/")[[1]][1]
}

# -----------------------------------------------------------------------------
# 7. Submit one LSF job per group
# -----------------------------------------------------------------------------
message("\nSubmitting LSF jobs...")

renv_lines <- paste0(
    "export RENV_PATHS_LIBRARY=\"$(pwd)/renv/library\"\n",
    if (nchar(renv_platform) > 0)
        paste0("export RENV_PATHS_PREFIX=\"", renv_platform, "\"\n")
    else ""
)

for (grp in names(seu_obj_list)) {
    safe_grp    <- gsub("[^A-Za-z0-9._-]", "_", grp)
    grp_rds     <- file.path(split_dir, paste0(safe_grp, ".RDS"))
    job_name    <- paste0("linkpeaks_", safe_grp)
    script_path <- file.path(jobs_dir, paste0("job_", safe_grp, ".sh"))
    out_log     <- file.path(logs_dir, paste0(safe_grp, "-%J.out"))
    err_log     <- file.path(logs_dir, paste0(safe_grp, "-%J.err"))

    job_script <- paste0(
        "#!/bin/bash\n",
        "#BSUB -J ", job_name, "\n",
        "#BSUB -q ", LSF_QUEUE, "\n",
        "#BSUB -P ", LSF_PROJ, "\n",
        "#BSUB -o ", out_log, "\n",
        "#BSUB -e ", err_log, "\n",
        "#BSUB -n ", LSF_CORES, "\n",
        "#BSUB -R 'rusage[mem=", LSF_MEM, "]'\n",
        "#BSUB -R span[hosts=1]\n",
        "#BSUB -W ", LSF_TIME, "\n",
        "\n",
        "module load R/4.4.1\n",
        "\n",
        renv_lines,
        "\n",
        "Rscript scripts/helpers/run_linkpeaks_worker.R \\\n",
        "    \"", grp, "\" \\\n",
        "    ", grp_rds, " \\\n",
        "    ", argv$pipeline_config, " \\\n",
        "    ", out_dir, "\n"
    )

    writeLines(job_script, script_path)
    system(paste0("bsub < ", script_path))
    message("  Submitted: ", job_name)
}

# -----------------------------------------------------------------------------
# 8. Submit merge job — pre-written script, depends on done(linkpeaks_*)
# -----------------------------------------------------------------------------
system(paste0('bsub -env "RDS_FILE_IN=', argv$RDS_file_in, '" < scripts/helpers/job_merge.sh'))
message("  Submitted: linkpeaks_merge (waits on all linkpeaks_* jobs)")

# -----------------------------------------------------------------------------
# 9. Summary
# -----------------------------------------------------------------------------
message("\n=== Step 11 submission complete ===")
message("Groups submitted: ", length(seu_obj_list))
message("Monitor jobs:     bjobs -J 'linkpeaks_*'")
message("Logs:             ", logs_dir, "/")
message("Per-group outputs (on child job completion):")
for (grp in names(seu_obj_list)) {
    safe_grp <- gsub("[^A-Za-z0-9._-]", "_", grp)
    message("  ", file.path(out_dir, paste0(safe_grp, "-linkpeaks.RDS")))
}
message("Final merged output (on merge job completion):")
message("  output/RDS-files/", argv$project_prefix, "-11-linkpeaks-obj.RDS")
