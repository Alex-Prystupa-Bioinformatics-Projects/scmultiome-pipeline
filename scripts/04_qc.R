#!/usr/bin/env Rscript
# =============================================================================
# 04_qc.R  |  Pipeline Step 4: QC Plots and Threshold Configuration
# =============================================================================
# Generates QC diagnostic plots for each sample and writes a qc_config.yml
# with default thresholds for the user to review and edit before filtering.
#
# No normalization or clustering is performed here — objects are raw counts
# with QC metrics attached (from steps 2 and 3).
#
# Plots generated:
#   - Violin plots per sample: nCount_RNA, nFeature_RNA, nCount_peaks,
#     nFeature_peaks, percent.mt, nucleosome_signal, TSS.enrichment
#   - Density scatter per sample: nCount_peaks vs TSS.enrichment
#   - Doublet score violin: scDblFinder.score per sample
#
# After reviewing the plots, edit configs/qc_config.yml to set your
# per-sample thresholds, then run 05_filter.R.
#
# Usage:
#   Rscript scripts/04_qc.R configs/samplesheet.csv \
#       --project_prefix myproject \
#       --RDS_file_in output/RDS-files/myproject-callpeaks-obj-list.RDS
#
# Input:  output of 03_callpeaks.R  ({project_prefix}-03-callpeaks-obj-list.RDS)
# Output: output/plots/{project_prefix}-qc-plots.pdf
#         output/tables/{project_prefix}-doublet-summary.csv
#         configs/qc_config.yml  (default thresholds — edit before filtering)
# =============================================================================

library(argparser, quietly = TRUE)
library(yaml,      quietly = TRUE)
library(Seurat,    quietly = TRUE)
library(Signac,    quietly = TRUE)
library(dplyr,     quietly = TRUE)
library(ggplot2,    quietly = TRUE)

# 1. Parse arguments
p <- arg_parser("Step 4: Generate QC plots and write qc_config.yml")
p <- add_argument(p, "samplesheet",
                  help = "Path to samplesheet CSV", type = "character")
p <- add_argument(p, "--pipeline_config",
                  help    = "Path to pipeline_config.yml",
                  default = "configs/pipeline_config.yml", type = "character")
p <- add_argument(p, "--project_prefix",
                  help    = "Output file prefix (no .RDS extension)",
                  default = "multiome", type = "character")
p <- add_argument(p, "--RDS_file_in",
                  help = "Input RDS file — output of 03_callpeaks.R", type = "character")
argv <- parse_args(p)

if (is.na(argv$RDS_file_in) || argv$RDS_file_in == "NA") {
    stop("--RDS_file_in is required. Provide the path to output of 03_callpeaks.R.")
}

# 2. Load config and helper functions
pipeline_config <- yaml::read_yaml(argv$pipeline_config)
samplesheet     <- read.csv(argv$samplesheet)

source("scripts/utils.R")

# 3. Load input objects from previous step
message("Starting Step 4: qc")
message("Loading: ", argv$RDS_file_in)
seu_obj_list <- readRDS(argv$RDS_file_in)

# 4. Merge all objects into one for multi-sample violin plots
#    Allows VlnPlot to group by orig.ident across all samples at once
seu_obj_merged <- merge(seu_obj_list[[1]], seu_obj_list[-1])

dir.create("output/plots",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

pdf_path <- paste0("output/plots/", argv$project_prefix, "-qc-plots.pdf")
pdf(file = pdf_path, height = 8, width = 12)

# 5. Violin plots of QC metrics grouped by sample (orig.ident)
#    - Counts:    nCount_RNA, nCount_peaks
#    - Features:  nFeature_RNA, nFeature_peaks
#    - Mito:      percent.mt
#    - Chromatin: nucleosome_signal, TSS.enrichment
qc_feature_groups <- list(
    counts    = c("nCount_RNA",   "nCount_peaks"),
    features  = c("nFeature_RNA", "nFeature_peaks"),
    mito      = "percent.mt",
    chromatin = c("nucleosome_signal", "TSS.enrichment")
)

violin_plots <- lapply(names(qc_feature_groups), function(group_name) {
    VlnPlot(seu_obj_merged,
            features = qc_feature_groups[[group_name]],
            group.by = "orig.ident",
            pt.size  = 0,
            log      = group_name %in% c("counts", "features")) +
        NoLegend() +
        ggtitle(group_name)
})
print(violin_plots)

# 6. Density scatter plots per sample: nCount_peaks vs TSS.enrichment
#    - Helps identify low-quality cells with poor chromatin accessibility
#    - Quantile lines help spot where natural cutoffs fall
density_plots <- lapply(seu_obj_list, function(seu_obj) {
    DensityScatter(seu_obj, x = "nCount_peaks", y = "TSS.enrichment",
                   log_x = TRUE, quantiles = TRUE) +
        ggtitle(seu_obj@project.name)
})
print(density_plots)

# 7. Doublet score violin plot across all samples
#    - scDblFinder.score: 0 = confident singlet, 1 = confident doublet
print(
    VlnPlot(seu_obj_merged,
            features = "scDblFinder.score",
            group.by = "orig.ident",
            pt.size  = 0) +
        NoLegend() +
        ggtitle("Doublet Score by Sample (scDblFinder)")
)

dev.off()
message("  QC plots saved: ", pdf_path)

# 8. Save doublet summary CSV: total cells, doublet count, and % per sample
doublet_summary <- bind_rows(lapply(seu_obj_list, function(seu_obj) {
    n_total    <- ncol(seu_obj)
    n_doublets <- sum(seu_obj$scDblFinder.class == "doublet")
    data.frame(
        sample          = seu_obj@project.name,
        total_cells     = n_total,
        n_doublets      = n_doublets,
        pct_doublets    = round(n_doublets / n_total * 100, 2)
    )
}))

doublet_csv_path <- paste0("output/tables/", argv$project_prefix, "-doublet-summary.csv")
write.csv(doublet_summary, file = doublet_csv_path, row.names = FALSE)
message("  Doublet summary saved: ", doublet_csv_path)

# 9. Write qc_config.yml with default thresholds for each sample
#    Review the plots above, then edit configs/qc_config.yml to set
#    your actual per-sample cutoffs before running 05_filter.R
#    Note: remove_doublets defaults to FALSE — set to TRUE to filter doublets
create_qc_yaml_file(seu_obj_list)
message("  Default thresholds written to: configs/qc_config.yml")
message("  --> Review QC plots and edit configs/qc_config.yml before running 05_filter.R")

message("Step 4 complete.")
